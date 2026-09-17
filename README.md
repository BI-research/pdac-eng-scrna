# PDAC ENG (CD105) single-cell RNA-seq analysis

This repository contains R code used for downstream analysis of publicly available single-cell RNA-seq data in the PDAC ENG (CD105) study.

The repository focuses on comparisons between **ENG-positive (CD105+)** and **ENG-negative (CD105−) cancer-associated fibroblasts (CAFs)**. Upstream single-cell RNA-seq preprocessing, integration, clustering, and cell-type annotation are not included.


## Repository structure

```text
scripts/
  01_ENG_classification.R
  02_DEG_naive_CAF.R
  03_DEG_treated_CAF.R
  04_GSEA_naive_CAF.R
  05_GSEA_treated_CAF.R

data/
  processed_naive_CAF.rds
  processed_treated_CAF.rds
```

## Main R packages

- Seurat
- dplyr
- ggplot2
- ggrepel
- presto
- fgsea
- msigdbr
- jsonlite
- tibble
