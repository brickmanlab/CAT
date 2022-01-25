#!/bin/bash

set -e 

catcli \
    --ds1 ./datasets/06_log.h5ad \
    --ds1_name org \
    --ds1_cluster seurat_clusters \
    --ds2 ./datasets/13_log_counts.h5ad \
    --ds2_name vitro \
    --ds2_cluster seurat_clusters \
    --output ./results/org-vs-vitro \
    --n_iter 1000 \
    --distance euclidean

catcli \
    --ds1 ./datasets/06_log.h5ad \
    --ds1_name org \
    --ds1_cluster seurat_clusters \
    --ds2 ./datasets/02_rothova_log.h5ad \
    --ds2_name rot \
    --ds2_cluster cluster_names \
    --output ./results/org-vs-rot \
    --n_iter 1000 \
    --distance euclidean

catcli \
    --ds1 ./datasets/Kat_subset_log.h5ad \
    --ds1_name kat \
    --ds1_cluster CellType \
    --ds2 ./datasets/02_rothova_log.h5ad \
    --ds2_name rot \
    --ds2_cluster cluster_names \
    --output ./results/kat-vs-rot \
    --n_iter 1000 \
    --distance euclidean

catcli \
    --ds1 ./datasets/06_log.h5ad \
    --ds1_name org \
    --ds1_cluster seurat_clusters \
    --ds2 ./datasets/Kat_subset_log.h5ad \
    --ds2_name kat \
    --ds2_cluster CellType \
    --output ./results/org-vs-kat \
    --n_iter 1000 \
    --distance euclidean
