import itertools
import json
from argparse import Namespace
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import polars as pl
import scipy
from rich.progress import track
from xlsxwriter import Workbook

from cat.constants import CLUSTER_FIELD, DELIMITER
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

            tables[dataset_name_from][dataset_name_to][CLUSTER_FIELD]

        In the example above, we get the distances from "example_cluster" in the dataset with name: "dataset_name_from", compared
        to all clusters in the dataset with the name: dataset_name_to.
    """
    ds_comparisons = list(itertools.product((diff.ds1_name, diff.ds2_name), repeat=2))
    dist_mean = diff.mean.to_pandas().set_index(CLUSTER_FIELD)
    idxs = dist_mean.columns

    tables = {name.split(DELIMITER)[0]: {} for name in idxs}
    for ds1_name, ds2_name in ds_comparisons:
        tables[ds1_name][ds2_name] = {}

        subset_cols = idxs[idxs.str.startswith(ds1_name)]

        for cluster in subset_cols:
            # remove self-loop
            subset_rows = (pl.col(CLUSTER_FIELD).str.starts_with(ds2_name)) & (
                pl.col(CLUSTER_FIELD) != cluster
            )

            tables[ds1_name][ds2_name][cluster] = (
                pl.DataFrame(
                    {
                        CLUSTER_FIELD: diff.mean.filter(subset_rows).get_column(
                            CLUSTER_FIELD
                        ),
                        "dist_mean": diff.mean.filter(subset_rows).get_column(cluster),
                        "dist_std": diff.std.filter(subset_rows).get_column(cluster),
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


def to_excel(args: Namespace, tables: dict[str, str]) -> None:
    """Generate CAT results to Excel format

    Parameters
    ----------
    args
        :py:class:`argparse.Namespace` Cli arguments
    tables : dict[str, str]
        Pairwise comparisons
    """

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

    for ds_from in tables:
        for ds_to in tables[ds_from]:
            filename = f"{args.output}/{ds_from}_{ds_to}_{args.distance}.xlsx"

            with Workbook(filename, {"nan_inf_to_errors": True}) as wb:
                dashboard.write_excel(workbook=wb, worksheet="Dashboard")
                for cluster in tables[ds_from][ds_to]:
                    sheet_name = cluster.replace(" ", "_").replace(":", ".")
                    tables[ds_from][ds_to][cluster].write_excel(
                        workbook=wb, worksheet=sheet_name
                    )


def to_dict(
    args: Namespace, tables: dict[str, str], significant_only: bool = False
) -> dict[str, str | dict]:
    """Gerate CAT results to Python dictionaries

    Parameters
    ----------
    args : Namespace
        :py:class:`argparse.Namespace` Cli arguments
    tables : dict[str, str]
        Pairwise comparisons
    significant_only : bool, optional
        Subset only significant interactions, by default False

    Returns
    -------
    dict[str, str | dict]
        CAT comparisons
    """

    comparisons = {}
    for ds_from in tables:
        comparisons[ds_from] = {}
        for ds_to in tables[ds_from]:
            comparisons[ds_from][ds_to] = {
                cluster: (
                    tables[ds_from][ds_to][cluster]
                    .filter(pl.col.significant)
                    .to_dicts()
                    if significant_only
                    else tables[ds_from][ds_to][cluster].to_dicts()
                )
                for cluster in tables[ds_from][ds_to]
            }

    return {
        "dashboard": {
            "dataset1": {"name": args.ds1_name, "file": args.ds1},
            "dataset2": {"name": args.ds2_name, "file": args.ds2},
            "features": args.features,
            "sigma": args.sigma,
            "interations": args.n_iter,
            "distance": args.distance,
        },
        "comparisons": comparisons,
    }


def to_json(args: Namespace, tables: dict[str, str]) -> None:
    """Gerate CAT results to JSON format

    Parameters
    ----------
    args : Namespace
        :py:class:`argparse.Namespace` Cli arguments
    tables : dict[str, str]
        Pairwise comparisons
    """

    filename = f"{args.output}/cat_comparisons.json"
    with open(filename, "w") as f:
        data = to_dict(args, tables)
        json.dump(data, f)


def to_html(args: Namespace, tables: dict[str, str]) -> None:
    """Gerate CAT results to HTML format

    Parameters
    ----------
    args : Namespace
        :py:class:`argparse.Namespace` Cli arguments
    tables : dict[str, str]
        Pairwise comparisons
    """

    data = to_dict(args, tables, significant_only=True)

    comparisons = data["comparisons"]
    for ds_from in track(comparisons, description="Saving plots..."):
        for ds_to in comparisons[ds_from]:
            sources, targets, values = [], [], []

            # create all labels
            labels = set(comparisons[ds_from][ds_to].keys())
            for _, item_to in comparisons[ds_from][ds_to].items():
                missing_labels = (
                    set([x.get(CLUSTER_FIELD, "N/A") for x in item_to]) - labels
                )
                labels = labels.union(missing_labels)
            labels = list(labels)

            # Sankey plots
            for item_from, item_to in comparisons[ds_from][ds_to].items():
                sources += [labels.index(item_from) for _ in item_to]
                targets += [labels.index(x.get(CLUSTER_FIELD, "N/A")) for x in item_to]
                values += [x.get("dist_mean", "N/A") for x in item_to]

            go.Figure(
                data=[
                    go.Sankey(
                        arrangement="snap",
                        node=dict(
                            pad=15,
                            thickness=20,
                            line=dict(color="black", width=0.5),
                            label=labels,
                        ),
                        link=dict(
                            source=sources,
                            target=targets,
                            value=values,
                            hovertemplate="From: %{source.label}<br />"
                            + "To: %{target.label}<br />Value: %{value}",
                        ),
                    )
                ]
            ).write_html(f"{args.output}/{ds_from}_{ds_to}_{args.distance}.html")


def save_tables(args: Namespace, tables: dict[str, str]):
    """Save results into multiple formats

    Parameters
    ----------
    args : Namespace
        :py:class:`argparse.Namespace` Cli arguments
    tables : dict[str, str]
        Pairwise comparisons
    """
    Path(args.output).mkdir(parents=True, exist_ok=True)

    to_html(args, tables)
    to_json(args, tables)
    to_excel(args, tables)
