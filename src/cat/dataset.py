import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import anndata
import numpy as np
import polars as pl
import scipy

from cat.constants import CLUSTER_FIELD, DELIMITER


@dataclass
class Dataset:
    """Class describing scRNA-seq dataset"""

    filename: str
    name: str
    cluster_field: str = CLUSTER_FIELD

    def __post_init__(self):
        if self.name == "" or self.filename == "":
            logging.error("You forgot to provide all necessary parameters!")
            sys.exit(1)

        if not Path(self.filename).exists():
            logging.error(f"Provided {self.filename} doesn't exists!")
            sys.exit(1)

        self.adata = anndata.read_h5ad(self.filename)

    def _fix_metadata(self, group_by: str):
        """Helper function to fix metadata information and add `cluster_field` to the `obs`

        Parameters
        ----------
        group_by : str
            Column in `obs`
        """
        if group_by == "":
            logging.error(
                "You forgot to specify `group_by` param for cluster comparisons."
            )
            sys.exit(1)

        if group_by not in self.adata.obs.columns:
            logging.error(f"Defined column {group_by} not found in the dataset.")
            sys.exit(1)

        self.adata.obs[self.cluster_field] = (
            self.name + DELIMITER + self.adata.obs[group_by].astype(str)
        )

    def _fix_genes(self, gene_symbol_field: str):
        """Set correct column from `var` to gene symbol.

        Parameters
        ----------
        gene_symbol_field : str
            Column in `var`
        """
        if gene_symbol_field is not None and gene_symbol_field in self.adata.var_keys():
            self.adata.var.index = self.adata.var[gene_symbol_field].to_numpy()

        self.adata.var_names = self.adata.var_names.str.lower()
        self.adata.var_names_make_unique()

        # check if they are really gene symbols and not Ensembl IDs
        if self.adata.var_names[0].startswith("ENS"):
            logging.error("`var_names` should contain gene symbols!")
            sys.exit(1)

    def _filter_genes(self, gene_type: str, pattern: str):
        """Filter genes based on pattern

        Parameters
        ----------
        gene_type : str
            Gene type name (i.e. mitochondrial genes)
        pattern : str
            String/Regex pattern
        """
        to_filter = self.adata.var_names.str.match(pattern)
        n_genes: int = np.sum(to_filter)
        if n_genes > 0:
            logging.info(f"Removing {n_genes} {gene_type} genes")
            self.adata = self.adata[:, ~to_filter].copy()

    def _save(self, save_path: str):
        """Save dataset into anndata format

        Parameters
        ----------
        save_path : str
            File path
        """
        if not Path(save_path).exists():
            Path(save_path).mkdir(parents=True, exist_ok=True)

        filename: str = f"{save_path}/{self.name}.h5ad"
        logging.info(f"Saving processed dataset into {filename}")
        self.adata.write(filename)

    def prepare(
        self,
        group_by: str,
        gene_symbol_field: str | None = None,
        save_path: str | None = None,
    ) -> None:
        """Prepares dataset for CAT analysis

        Parameters
        ----------
        group_by : str
            Column in `obs`
        gene_symbol_field : str | None, optional
            Column in `var` where genes are gene symbols, by default None
        save_path : str, optional
            Path to safe processed dataset in anndata format, by default "./tmp"
        """
        logging.info(f"Preprocessing {self.name}")

        if scipy.sparse.issparse(self.adata.X):
            self.adata.X = self.adata.X.todense()

        self._fix_metadata(group_by=group_by)
        self._fix_genes(gene_symbol_field=gene_symbol_field)

        self._filter_genes(gene_type="mitochondrial", pattern="mt-")
        self._filter_genes(gene_type="ribosomal", pattern="rp[s,l]")
        self._filter_genes(gene_type="spike", pattern="ercc-")

        if save_path:
            self._save(save_path=save_path)


@dataclass
class DatasetDiff:
    """Pairwise dataset container"""

    ds1_name: str
    ds2_name: str
    mean: pl.DataFrame
    std: pl.DataFrame
