import argparse
import logging
import sys
import warnings

import anndata
import numpy as np
import polars as pl
from rich.progress import track
from scipy.spatial.distance import pdist, squareform

from cat.constants import CLUSTER_FIELD, MIN_N_GENES, N_ITERATIONS
from cat.dataset import Dataset, DatasetDiff
from cat.report import generate_tables, save_tables
from cat.utils import normalize, read_features, rename_ds


def internal_preprocessing(
    ds1: Dataset, ds2: Dataset, features_file: str | None = None, normalize: bool = True
) -> tuple[Dataset, Dataset]:
    """The internal preprocessing of the CAT algorithm ensures that data is in the right shape and format.

    Parameters
    ----------
    ds1 : Dataset
        Source dataset
    ds2 : Dataset
        Target dataset
    features_file : str | None, optional
        List of genes, i.e. GO:TERM, by default None
    normalize : bool, optional
        Normalize data, by default True

    Returns
    -------
    tuple[Dataset, Dataset]
        CAT ready datasets
    """
    logging.info(f"Before {ds1.name}: {ds1.adata.shape}")
    logging.info(f"Before {ds2.name}: {ds2.adata.shape}")

    genes = set(ds1.adata.var_names) & set(ds2.adata.var_names)
    if features_file is not None:
        features = read_features(features_file)
        genes = genes & set(features)
        logging.info(f"Reading list of genes from {features_file} => {len(genes)}")
    genes = list(genes)

    if len(genes) == 0:
        logging.error("No common genes found ...")
        sys.exit(1)
    elif len(genes) < MIN_N_GENES:
        logging.warning(f"Only <{MIN_N_GENES} genes are common ...")
    else:
        ds1.adata = ds1.adata[:, genes].copy()
        ds2.adata = ds2.adata[:, genes].copy()

    logging.info(f"After {ds1.name}: {ds1.adata.shape}")
    logging.info(f"After {ds2.name}: {ds2.adata.shape}")

    if normalize:
        ds1.adata.X = ds1.adata.X / ds1.adata.X.sum(axis=1, keepdims=1)
        ds2.adata.X = ds2.adata.X / ds2.adata.X.sum(axis=1, keepdims=1)

    if not np.all(ds1.adata.var_names == ds2.adata.var_names):
        logging.error("Gene intersection between two datasets doesn't match!")
        sys.exit(1)

    # Make sure the normalization is good
    check_normalization = lambda x: np.all(np.abs(x.adata.X.sum(axis=1) - 1) < 0.001)
    if not (check_normalization(ds1) or check_normalization(ds2)):
        logging.error(
            "Datasets are not normalized properly. Make sure you normalized data in advance."
        )
        sys.exit(1)

    return ds1, ds2


def compare(
    ds1: Dataset,
    ds2: Dataset,
    n_iterations: int = N_ITERATIONS,
    features: list[str] | None = None,
    distance: str = "euclidean",
) -> DatasetDiff:
    """CAT routine calculated the inter cluster distances. The order of dataset1 and dataset2 does not matter.

    Parameters
    ----------
    dataset1 : Dataset
        First dataset, must contain anndata, with labels for
        each cluster in an obs variable "XXX". Should be normalized.
    dataset2 : Dataset
        Second dataset, must contain anndata, with labels for
        each cluster in an obs variable "XXX". Should be normalized.
    n_iterations : int, optional
        Number of iterations in the bootstrap process, by default 'value'
    features : list[str], optional
        Subset of gene names (i.e. GO:TERM)
    distance: str
        Euclidean distance, by default euclidean


    Returns
    -------
    Pandas DataFrame
        A table containing the distances of each cluster relative to every other cluster
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        adata = anndata.concat([ds1.adata, ds2.adata])
        adata.obs_names_make_unique()

    # SECOND - ENABLE GENE SUBSET / FEATURE SELECTION
    if features:
        adata = adata[:, features].copy()

    # THIRD - NORMALIZE BY MEDIAN GENE
    # Ignore Warning because nzm can return NaN, which is further fixed
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        df = pl.DataFrame(
            normalize(adata.X, method="median"), schema=adata.var_names.tolist()
        ).with_columns(pl.Series(name=CLUSTER_FIELD, values=adata.obs[CLUSTER_FIELD]))

    distances = []
    partitions: dict[str, pl.DataFrame] = df.partition_by(
        CLUSTER_FIELD, as_dict=True, include_key=False
    )
    cluster_names = list(sum(partitions.keys(), ()))

    for _ in track(range(n_iterations), description="Processing..."):
        # FIFTH - CLUSTER AVERAGES
        cluster_means = pl.concat(
            [
                x.sample(fraction=1.0, with_replacement=True).mean()
                for x in partitions.values()
            ]
        )

        # SIXTH - CALCULATE DISTANCES FOR THIS BOOTSTRAP ITERATION
        distances.append(pdist(cluster_means.to_numpy(), metric=distance))

    # SEVENTH - GATHER RESULT - MEAN and STD
    dist_mean = np.array(distances).mean(axis=0)
    dist_std = np.array(distances).std(axis=0)

    dist_mean_df = pl.DataFrame(
        squareform(dist_mean), schema=cluster_names
    ).with_columns(pl.Series(name=CLUSTER_FIELD, values=cluster_names))
    dist_std_df = pl.DataFrame(squareform(dist_std), schema=cluster_names).with_columns(
        pl.Series(name=CLUSTER_FIELD, values=cluster_names)
    )

    return DatasetDiff(ds1.name, ds2.name, dist_mean_df, dist_std_df)


def run(args: argparse.Namespace):
    """Run CAT inference

    Parameters
    ----------
    args : argparse.Namespace
        CLI arguments
    """
    if not (args.ds1 or args.ds1_cluster or args.ds2 or args.ds2_cluster):
        logging.error("Two datasets with specified cluster column are required")
        sys.exit(1)

    ds1_name, ds2_name = rename_ds([args.ds1_name, args.ds2_name])
    ds1 = Dataset(name=ds1_name, filename=args.ds1)
    ds1.prepare(group_by=args.ds1_cluster, gene_symbol_field=args.ds1_genes)

    ds2 = Dataset(name=ds2_name, filename=args.ds2)
    ds2.prepare(group_by=args.ds2_cluster, gene_symbol_field=args.ds2_genes)

    ds1, ds2 = internal_preprocessing(ds1, ds2, features_file=args.features)

    diff = compare(ds1, ds2, n_iterations=args.n_iter, distance=args.distance)

    logging.info(f"Saving results to {args.output}")
    tables = generate_tables(diff, args.sigma)

    save_tables(args, tables)
