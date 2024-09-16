import itertools
from argparse import Namespace
from pathlib import Path

import numpy as np
import polars as pl
import scipy
from xlsxwriter import Workbook

from cat.constants import DELIMITER
from cat.dataset import DatasetDiff


def generate_tables(diff: DatasetDiff, sigma_th: float) -> dict[str, str]:
    """
    Takes the raw results from the CAT routine function and turns them into nice table per cluster.

    Parameters
    ----------
    diff : DatasetDiff
        Contains mean and std distance matrix (N x N)
    sigma_th: float
        Cutoff filter

    Returns
    -------
    dict[str, str]
        A dictionary containing tables for each cluster.
        The dataframe can be access in the following way:

            tables[dataset_name_from][dataset_name_to]["example_cluster"]

        In the example above, we get the distances from "example_cluster" in the dataset with name: "dataset_name_from", compared
        to all clusters in the dataset with the name: dataset_name_to.
    """
    ds_comparisons = list(itertools.product((diff.ds1_name, diff.ds2_name), repeat=2))
    dist_mean = diff.mean.to_pandas().set_index("cat_cluster")
    idxs = dist_mean.columns

    tables = {name.split(DELIMITER)[0]: {} for name in idxs}
    for ds1_name, ds2_name in ds_comparisons:
        tables[ds1_name][ds2_name] = {}

        subset_cols = idxs[idxs.str.startswith(ds2_name)]

        for cluster in subset_cols:
            # remove self-loop
            subset_rows = (pl.col("cat_cluster").str.starts_with(ds1_name)) & (
                pl.col("cat_cluster") != cluster
            )
            mean_per_cluster = diff.mean.filter(subset_rows).select(
                [cluster, "cat_cluster"]
            )
            std_per_cluster = diff.std.filter(subset_rows).select(cluster)

            tables[ds1_name][ds2_name][cluster] = (
                pl.DataFrame(
                    {
                        "dist_mean": mean_per_cluster.get_column(cluster),
                        "dist_std": std_per_cluster,
                        "cat_cluster": mean_per_cluster.get_column("cat_cluster"),
                    }
                )
                .sort(by="dist_mean")
                .with_columns(
                    diff_to_closest=pl.col.dist_mean - pl.col.dist_mean.get(0)
                )
                .with_columns(
                    diff_uncertainty=np.sqrt(
                        pl.col.dist_std**2 + pl.col.dist_std.get(0) ** 2
                    )
                )
                .with_columns(
                    diff_sigma_away=pl.col.diff_to_closest
                    / pl.col.diff_uncertainty.get(0)
                )
                .with_columns(
                    diff_sigma_away_p=pl.col.diff_sigma_away.map_elements(
                        lambda x: scipy.stats.norm.sf(x), return_dtype=pl.Float32
                    )
                )
                .with_columns(significant=pl.col.diff_sigma_away < sigma_th)
            )

    return tables


def to_excel(
    dashboard: pl.DataFrame, tables: dict[str, str], output: str, distance: str
) -> None:
    """Generate CAT results to Excel format

    Parameters
    ----------
    dashboard : pl.DataFrame
        :py:class:`polars.DataFrame` Dataframe containing specified parameters
    tables : dict[str, str]
        Pairwise comparisons
    output : str
        Output folder
    distance : str
        Used distance metric
    """
    for ds_from in tables:
        for ds_to in tables[ds_from]:
            filename = f"{output}/{ds_from}_{ds_to}_{distance}.xlsx"

            with Workbook(filename) as wb:
                dashboard.write_excel(workbook=wb, worksheet="Dashboard")
                for cluster in tables[ds_from][ds_to]:
                    sheet_name = cluster.replace(" ", "_").replace(":", ".")
                    tables[ds_from][ds_to][cluster].write_excel(
                        workbook=wb, worksheet=sheet_name
                    )


def save_tables(args: Namespace, tables: dict[str, str]):
    """Save results into specified format

    Parameters
    ----------
    args
        :py:class:`argparse.Namespace` Cli arguments
    tables : dict[str, str]
        Tables from generate_tables
    """
    Path(args.output).mkdir(parents=True, exist_ok=True)

    dashboard = pl.DataFrame(
        {
            "dataset1": [args.ds1, args.ds1_name, args.ds1_cluster, args.ds1_genes],
            "dataset2": [args.ds2, args.ds2_name, args.ds2_cluster, args.ds2_genes],
            "features": args.features,
            "distance": args.distance,
            "sigma": args.sigma,
            "iterations": args.n_iter,
        }
    )

    match args.format:
        case "excel":
            to_excel(dashboard, tables, args.output, args.distance)
        case "html":
            ...
        case "pdf":
            ...
        case _:
            raise ValueError(f"Output format `{args.format}` not supported!")
