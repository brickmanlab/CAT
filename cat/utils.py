import logging
import sys
import time
from functools import wraps
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


def timeit(func):
    # stolen from: https://dev.to/kcdchennai/python-decorator-to-measure-execution-time-54hk
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f"Function {func.__name__}: Took {total_time:.4f} seconds")
        return result

    return timeit_wrapper
