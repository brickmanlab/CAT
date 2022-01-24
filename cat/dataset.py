import logging
import sys
from cmath import log
from pathlib import Path
from typing import Any, Optional

import anndata
import numpy as np
import scipy


class Dataset:
    def __init__(self, name: str, file: str):

        if name == "" or file == "":
            logging.error(f"You forgot to provide all necessary parameters!")
            sys.exit(1)

        if not Path(file).exists():
            logging.error(f"Provided {file} doesn't exists!")
            sys.exit(1)

        self.adata = anndata.read_h5ad(file)
        self.name = name
        self.cluster = "cat_cluster"

    def _fix_metadata(self, group_by: str):
        if group_by == "":
            logging.error(
                "You forgot to specify `group_by` param for cluster comparisons."
            )
            sys.exit(1)

        if group_by not in self.adata.obs.columns:
            logging.error(f"Defined column {group_by} not found in the dataset.")
            sys.exit(1)

        self.adata.obs[self.cluster] = (
            self.name + "_" + self.adata.obs[group_by].astype(str)
        )

    def _fix_genes(self, gene_symbol: str):
        if gene_symbol is not None and gene_symbol in self.adata.var_keys():
            self.adata.var.index = self.adata.var[gene_symbol].to_numpy()

        self.adata.var_names = self.adata.var_names.str.lower()
        self.adata.var_names_make_unique()

        # check if they are really gene symbols and not Ensembl IDs
        if self.adata.var_names[0].startswith("ENS"):
            logging.error("`var_names` should contain gene symbols!")
            sys.exit(1)

    def _filter_genes(self, gene_type: str, pattern: Any):
        to_filter = self.adata.var_names.str.startswith(pattern)
        n_genes: int = np.sum(to_filter)
        if n_genes > 0:
            logging.info(f"Removing {n_genes} {gene_type} genes")
            self.adata = self.adata[:, ~to_filter].copy()

    def _save(self, save_path: str):
        if not Path(save_path).exists():
            Path(save_path).mkdir(parents=True, exist_ok=True)

        filename: str = f"{save_path}/{self.name}.h5ad"
        logging.info(f"Saving processed dataset into {filename}")
        self.adata.write(filename)

    def prepare(
        self, group_by: str, save_path: str = "./tmp", gene_symbol: Optional[str] = None
    ) -> anndata.AnnData:

        logging.info(f"Preprocessing {self.name}")

        if scipy.sparse.issparse(self.adata.X):
            self.adata.X = self.adata.X.todense()

        self._fix_metadata(group_by=group_by)
        self._fix_genes(gene_symbol=gene_symbol)

        self._filter_genes(gene_type="mitochondrial", pattern="mt-")
        self._filter_genes(gene_type="ribosomal", pattern=("rps", "rpl"))
        self._filter_genes(gene_type="spike", pattern="ercc-")

        # self._save(save_path=save_path)
