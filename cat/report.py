import itertools
from pathlib import Path
from typing import Any, Dict

import numpy as np
import pandas as pd
import scipy


def highlight_max(
    data: pd.DataFrame,
    is_significant: bool,
    color1: str = "yellow",
    color2: str = "green",
):
    """
    highlight the maximum in a Series or DataFrame
    """

    attr1 = "color: {}".format(color1)  # background-color
    attr2 = "color: {}".format(color2)

    return pd.DataFrame(
        np.where(np.array(data.shape[1] * [is_significant]).T, attr1, attr2),
        index=data.index,
        columns=data.columns,
    )


def generate_tables(dist_mean: pd.DataFrame, dist_std: pd.DataFrame, sigma_th: float):
    """
    Takes the raw results from the CAT routine function and turns them into nice table per cluster.

    Parameters
    ----------
    redundant_distance_mean_df : pandas.core.frame.DataFrame
        The dataframe contains a matrix of the mean distance from all clusters to all other clusters
    redundant_distance_std_df : pandas.core.frame.DataFrame
        The dataframe contains a matrix of the std on all distance from all clusters to all other clusters

    Returns
    -------
    dict
        A dictionary containing tables for each cluster. The dictionary contains pandas dataframes and styled dataframes.
        The two representation can be access in the following way;

        tables[dataset_name_from][dataset_name_to]["example_cluster"]['dataframe']
        tables[dataset_name_from][dataset_name_to]["example_cluster"]['styled']

        In the example above, we get the distances from "example_cluster" in the dataset with name: "dataset_name_from", compared
        to all clusters in the dataset with the name: dataset_name_to.
    """
    tables = {name.split("_")[0]: {} for name in dist_mean.columns}
    pairwise = list(itertools.product(tables.keys(), repeat=2))

    for ds1_name, ds2_name in pairwise:
        tables[ds1_name][ds2_name] = {}

        filt_cols = dist_mean.columns[dist_mean.columns.str.startswith(ds1_name)]
        for cluster in filt_cols:
            filt_rows = dist_mean.columns[dist_mean.columns.str.startswith(ds2_name)]
            filt_rows = [i for i in filt_rows if cluster != i]

            dist_per_cluster = dist_mean.loc[filt_rows, filt_cols][cluster]
            std_per_cluster = dist_std.loc[filt_rows, filt_cols][cluster]
            dist_per_cluster_sorted = dist_per_cluster[np.argsort(dist_per_cluster)]
            std_per_cluster_sorted = std_per_cluster[np.argsort(dist_per_cluster)]
            df_per_cluster = pd.DataFrame(
                [dist_per_cluster_sorted, std_per_cluster_sorted]
            ).T
            df_per_cluster.columns = ["dist mean", "dist std"]
            df_per_cluster.columns.name = cluster

            df_per_cluster["diff to closest"] = [
                i - df_per_cluster["dist mean"][0] for i in df_per_cluster["dist mean"]
            ]
            df_per_cluster["diff uncertainty"] = [
                np.sqrt(i ** 2 + df_per_cluster["dist std"][0] ** 2)
                for i in df_per_cluster["dist std"]
            ]
            df_per_cluster["diff sigma away"] = (
                df_per_cluster["diff to closest"] / df_per_cluster["diff uncertainty"]
            )
            df_per_cluster["diff sigma away p"] = scipy.stats.norm.sf(
                df_per_cluster["diff sigma away"]
            )

            df_per_cluster["significant"] = df_per_cluster["diff sigma away"] < sigma_th

            styled = df_per_cluster.style.apply(
                highlight_max,
                significant=df_per_cluster["significant"],
                color1="green",
                color2="darkorange",
                axis=None,
            )
            tables[ds1_name][ds2_name][cluster] = {}
            tables[ds1_name][ds2_name][cluster]["dataframe"] = df_per_cluster
            tables[ds1_name][ds2_name][cluster]["styled"] = styled

    return tables


def save_tables(tables: Dict[str, Any], output: str, dashboard: pd.DataFrame):
    Path(output).mkdir(parents=True, exist_ok=True)

    for ds_from in tables:
        for ds_to in tables[ds_from]:
            writer = pd.ExcelWriter(
                f'{output}/{ds_from}_{ds_to}_{dashboard.loc["distance", 1]}.xlsx',
                engine="xlsxwriter",
            )
            dashboard.to_excel(writer, sheet_name="Dashboard")
            for cluster in tables[ds_from][ds_to]:
                sheet_name = cluster.replace(" ", "_").replace(":", ".")
                df = tables[ds_from][ds_to][cluster]["dataframe"]
                df.to_excel(writer, sheet_name=sheet_name)

                workbook = writer.book
                data_format = workbook.add_format({"bg_color": "#00C7CE"})

                worksheet = writer.sheets[sheet_name]
                for idx, row in enumerate(df["significant"]):
                    if row:
                        worksheet.set_row(idx + 1, cell_format=data_format)
            writer.save()
