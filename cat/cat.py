import argparse
import logging
import sys
from typing import Any, Dict, List, Optional, Tuple

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from tqdm import tqdm

from .dataset import Dataset
from .report import generate_tables, save_tables
from .utils import get_nz_median, read_features, rename_dataset


def normalize(mat: np.ndarray, method="median"):
    if method == "median":
        nzm = get_nz_median(mat)
        # If entire row is 0, it is okay to divide by 1 (and not is otherwise median of 0)
        nzm[~np.isfinite(nzm)] = 1
        return mat / nzm

    logging.error(f"Normalization method {method} not implemented!")
    sys.exit(1)


def sample_cells(adata: anndata.AnnData, replace: bool = True):
    clusters = adata.obs.cat_cluster

    return np.concatenate(
        [
            np.random.choice(x, size=len(x), replace=replace)
            for x in clusters.groupby(clusters).groups.values()
        ]
    )


def internal_preprocessing(
    ds1: Dataset, ds2: Dataset, features_file: Optional[str], re_normalize: bool = True
) -> Tuple[Dataset, Dataset]:
    """The internal preprocessing of the CAT algorithm.
    ensures that data is in the right shape and format.
    """

    logging.info(f"Before {ds1.name}: {ds1.adata.shape}")
    logging.info(f"Before {ds2.name}: {ds2.adata.shape}")

    genes = set(ds1.adata.var_names) & set(ds2.adata.var_names)
    if features_file is not None:
        features = read_features(features_file)
        genes = genes & set(features)
        logging.info(f"Reading list of genes from {features_file} => {len(genes)}")

    if len(genes) == 0:
        logging.error(f"No common genes found ...")
        sys.exit(1)
    elif len(genes) < 20:
        logging.warning("Only <20 genes are common ...")
    else:
        ds1.adata = ds1.adata[:, list(genes)].copy()
        ds2.adata = ds2.adata[:, list(genes)].copy()

    logging.info(f"After {ds1.name}: {ds1.adata.shape}")
    logging.info(f"After {ds2.name}: {ds2.adata.shape}")

    if re_normalize:
        sc.pp.normalize_total(ds1.adata, target_sum=1)
        sc.pp.normalize_total(ds2.adata, target_sum=1)

    if not np.all(ds1.adata.var_names == ds2.adata.var_names):
        logging.error("Gene intersection between two datasets don't match!")
        sys.exit(1)

    # Make sure the normalization is good
    if not (
        np.all(np.abs(ds1.adata.X.sum(axis=1) - 1) < 0.001)
        or np.all(np.abs(ds2.adata.X.sum(axis=1) - 1) < 0.001)
    ):
        logging.error(
            "Datasets are not normalized properly. Make sure you normalized data in advance."
        )
        sys.exit(1)

    return ds1, ds2


def bootstrap(adata: anndata.AnnData, clusters: np.array, distance: str):
    subset = adata[sample_cells(adata, replace=True), :]

    # FIFTH - CLUSTER AVERAGES
    cluster_means = (
        (subset.to_df().groupby(adata.obs.cat_cluster.values).mean())
        .loc[clusters, :]
        .values
    )

    # SIXTH - CALCULATE DISTANCES FOR THIS BOOTSTRAP ITERATION
    return scipy.spatial.distance.pdist(cluster_means, metric=distance)


def compare(
    dataset1: anndata.AnnData,
    dataset2: anndata.AnnData,
    n_iterations: int = 1_000,
    features: Optional[List[str]] = None,
    distance: str = "euclidean",
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    CAT routine calculated the inter cluster distances.
    The order of dataset1 and dataset2 does not matter.

    Parameters
    ----------
    dataset1 : AnnData
        First dataset, must be scanpy anndata, with labels for
        each cluster in an obs variable "XXX". Should be normalized.
    dataset2 : AnnData
        Second dataset, must be scanpy anndata, with labels for
        each cluster in an obs variable "XXX". Should be normalized.
    n_iterations : int, optional
        Number of iterations in the bootstrap process, by default 'value'
    features:


    Returns
    -------
    Pandas DataFrame
        A table containing the distances of each cluster relative to every other cluster
    """

    adata = dataset1.concatenate(dataset2)
    adata.obs.cat_cluster = adata.obs.cat_cluster.astype("category")

    # SECOND - ENABLE GENE SUBSET / FEATURE SELECTION
    if features:
        adata = adata[:, features].copy()

    # THIRD - NORMALIZE BY MEDIAN GENE
    adata.X = normalize(adata.X, method="median")

    # FOURTH - BOOTSTRAP
    clusters = np.sort(adata.obs.cat_cluster.cat.categories)
    distances = np.array(
        [
            bootstrap(adata, clusters=clusters, distance=distance)
            for _ in tqdm(range(n_iterations))
        ]
    )

    # SEVENTH - GATHER RESULT - MEAN and STD
    dist_mean, dist_std = distances.mean(axis=0), distances.std(axis=0)

    dist_mean_df = pd.DataFrame(
        scipy.spatial.distance.squareform(dist_mean),
        index=clusters,
        columns=clusters,
    )

    dist_std_df = pd.DataFrame(
        scipy.spatial.distance.squareform(dist_std),
        index=clusters,
        columns=clusters,
    )

    return dist_mean_df, dist_std_df


def run(args: argparse.Namespace):
    if not (args.ds1 or args.ds1_cluster or args.ds2 or args.ds2_cluster):
        logging.error("Two datasets with specified cluster column are required")
        sys.exit(1)

    ds1_name, ds2_name = rename_dataset([args.ds1_name, args.ds2_name])

    settings: Dict[str, Any] = {
        "dataset1": [args.ds1, ds1_name, args.ds1_cluster, args.ds1_genes],
        "dataset2": [args.ds2, ds2_name, args.ds2_cluster, args.ds2_genes],
        "features": args.features,
        "distance": args.distance,
        "sigma": args.sigma,
        "iterations": args.n_iter,
    }
    dashboard = pd.DataFrame(settings.items()).set_index(0)

    ds1 = Dataset(name=args.ds1_name, file=args.ds1)
    ds1.prepare(group_by=args.ds1_cluster, gene_symbol=args.ds1_genes)

    ds2 = Dataset(name=args.ds2_name, file=args.ds2)
    ds2.prepare(group_by=args.ds2_cluster, gene_symbol=args.ds2_genes)

    a, b = internal_preprocessing(ds1, ds2, features_file=args.features)

    dist_mean, dist_std = compare(
        a.adata,
        b.adata,
        n_iterations=args.n_iter,
        distance=args.distance,
    )

    logging.info(f"Saving results to {args.output}")
    tables = generate_tables(dist_mean, dist_std, args.sigma)
    save_tables(tables, output=args.output, dashboard=dashboard)
