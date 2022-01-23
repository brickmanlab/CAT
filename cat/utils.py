import glob
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


def get_nz_mean(mat: np.ndarray):
    return np.apply_along_axis(lambda v: np.mean(v[np.nonzero(v)]), 0, mat)


def get_nz_median(mat: np.ndarray):
    return np.apply_along_axis(lambda v: np.median(v[np.nonzero(v)]), 0, mat)


def load_pathways(path: str) -> Dict[str, np.ndarray]:
    if not Path(path).exists():
        raise f"Provided path `{path}` doesn't exists!"

    files: List[str] = glob.glob(f"{path}/*")
    genes = [pd.read_csv(file, header=None).values.flatten() for file in files]
    renamed_files: List[str] = [
        Path(file).with_suffix("").name.strip().lower().replace(" ", "_")
        for file in files
    ]

    pathways: Dict[str, np.ndarray] = {"all": np.array(["all"])}

    return {**pathways, **dict(zip(renamed_files, genes))}
