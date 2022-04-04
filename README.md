# Cluster Alignment Tool (CAT)

## Install with pip

```bash
pip install cat-python
```

## Install from source

```bash
git clone https://github.com/brickmanlab/CAT.git && cd CAT
conda create --name cat python=3.7
pip install -e .
```

## How to run

```bash
$ catcli \
    --ds1 ds1.h5ad \
    --ds1_name DS1 \
    --ds1_cluster seurat_clusters \
    --ds2 ds2.h5ad \
    --ds2_name DS2 \
    --ds2_cluster seurat_clusters \
    --output ./results/ds1-vs-ds2

# generate sankey plot
$ Rscript ./CAT/scripts/sankey.R \
  --excel ./results/ds1-vs-ds2/ds1_ds2_euclidean.xlsx \
  --output ./results/ds1-vs-ds2/
```

## Help

```bash
$ conda activate cat
$ catcli --help
usage: catcli [-h] [--ds1 DS1] [--ds1_name DS1_NAME]
              [--ds1_cluster DS1_CLUSTER] [--ds1_genes DS1_GENES] [--ds2 DS2]
              [--ds2_name DS2_NAME] [--ds2_cluster DS2_CLUSTER]
              [--ds2_genes DS2_GENES] [--features FEATURES] [--output OUTPUT]
              [--distance DISTANCE] [--sigma SIGMA] [--n_iter N_ITER]
              [--format {excel,html}] [--verbose] [--version]

Cluster Alignment Tool (CAT)

optional arguments:
  -h, --help            show this help message and exit
  --ds1 DS1             Processed dataset (h5/h5ad)
  --ds1_name DS1_NAME   Dataset name
  --ds1_cluster DS1_CLUSTER
                        Column name for comparison
  --ds1_genes DS1_GENES
                        Gene column, using `index` as default
  --ds2 DS2             Processed dataset (h5/h5ad)
  --ds2_name DS2_NAME   Dataset name
  --ds2_cluster DS2_CLUSTER
                        Column name for comparison
  --ds2_genes DS2_GENES
                        Gene column, using `index` as default
  --features FEATURES   File containing list of genes on new lines
  --output OUTPUT       Output location
  --distance DISTANCE   Distance measurement
  --sigma SIGMA         Sigma cutoff (1.6 => p-value: 0.05)
  --n_iter N_ITER       Number of bootstraps, default 1,000
  --format {excel,html}
                        Report output format
  --verbose             Verbose mode
  --version             show program's version number and exit
```
