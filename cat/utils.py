import logging
import sys
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd


def init_logger(verbose: bool):
    logger = logging.getLogger(__name__)
    logger_format = (
        "[%(filename)s->%(funcName)s():%(lineno)s]%(levelname)s: %(message)s"
        if verbose
        else "[%(levelname)s]: %(message)s"
    )
    logging.basicConfig(format=logger_format, level=logging.INFO)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)


def get_nz_mean(mat: np.ndarray):
    return np.apply_along_axis(lambda v: np.mean(v[np.nonzero(v)]), 0, mat)


def get_nz_median(mat: np.ndarray):
    return np.apply_along_axis(lambda v: np.median(v[np.nonzero(v)]), 0, mat)


def read_features(file: str) -> List[str]:
    if not Path(file).exists():
        logging.error(f"Provided file {file} not found!")
        sys.exit(1)

    return pd.read_table(file, header=None)[0].str.lower().tolist()


def rename_dataset(names: List[str]):
    return [
        name.replace("(", "").replace(")", "").replace(".", "_").replace(" ", "")
        for name in names
    ]
