import logging
import sys
from pathlib import Path
from typing import Literal

import numpy as np
import polars as pl
from numpy.typing import ArrayLike


def get_nz_mean(mat: ArrayLike) -> np.ndarray:
    """Calculate non-zero mean

    Parameters
    ----------
    mat : ArrayLike
        N-dimensional matrix

    Returns
    -------
    np.ndarray
        N-dimensional matrix
    """
    return np.apply_along_axis(lambda v: np.mean(v[np.nonzero(v)]), 0, mat)


def get_nz_median(mat: ArrayLike) -> np.ndarray:
    """Calculate non-zero median

    Parameters
    ----------
    mat : ArrayLike
        N-dimensional matrix

    Returns
    -------
    np.ndarray
        N-dimensional matrix
    """
    return np.apply_along_axis(lambda v: np.median(v[np.nonzero(v)]), 0, mat)


def normalize(
    mat: ArrayLike, method: Literal["median", "mean", "to_one"] = "median"
) -> np.ndarray:
    """Normalize count matrix

    Parameters
    ----------
    mat : ArrayLike
        Count matrix
    method : Literal["median", "mean"]
        Normalization methods, by default "median"

    Returns
    -------
    np.ndarray
        Normalized count matrix
    """
    if method == "median":
        norm = get_nz_median(mat)
    elif method == "mean":
        norm = get_nz_mean(mat)
    else:
        logging.error(f"Normalization method {method} not implemented!")
        sys.exit(1)

    # If entire row is 0, it is okay to divide by 1 (and not is otherwise median of 0)
    norm[~np.isfinite(norm)] = 1
    return mat / norm


def rename_ds(names: list[str]) -> list[str]:
    """Rename dataset names

    Parameters
    ----------
    names : list[str]
        Dataset names

    Returns
    -------
    list[str]
        Renamed dataset names
    """
    return [
        name.replace("(", "").replace(")", "").replace(".", "_").replace(" ", "")
        for name in names
    ]


def read_features(file: str) -> list[str]:
    """Read feature file containing list in TSV format (tab separated).

    Parameters
    ----------
    file : str
        Feature file

    Returns
    -------
    list[str]
        List of genes
    """
    if not Path(file).exists():
        logging.error(f"Provided file {file} not found!")
        sys.exit(1)

    return (
        pl.read_csv(file, has_header=False, separator="\t")
        .get_column("column_1")
        .cast(str)
        .str.to_lowercase()
        .to_list()
    )
