import glob
import logging
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


def get_nz_mean(mat: np.ndarray):
    return np.apply_along_axis(lambda v: np.mean(v[np.nonzero(v)]), 0, mat)


def get_nz_median(mat: np.ndarray):
    return np.apply_along_axis(lambda v: np.median(v[np.nonzero(v)]), 0, mat)


def read_features(file: str) -> List[str]:
    if not Path(file).exists():
        logging.error(f"Provided file {file} not found!")
        sys.exit(1)

    return pd.read_table(file, header=None)[0].str.lower().tolist()
