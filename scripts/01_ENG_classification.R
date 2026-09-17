# ENG (CD105) classification in CAFs
# Cells with detectable ENG expression (> 0) are classified as ENG_pos;
# cells with ENG expression == 0 are classified as ENG_neg.
#
# Input: processed Seurat objects containing CAFs and a `celltype` metadata column.

library(Seurat)
library(dplyr)

naive_CAF <- readRDS("data/processed_naive_CAF.rds")
treated_CAF <- readRDS("data/processed_treated_CAF.rds")

classify_ENG <- function(seurat_obj) {
  seurat_obj$ENG_status <- ifelse(
    FetchData(seurat_obj, vars = "ENG")[, "ENG"] > 0,
    "ENG_pos",
    "ENG_neg"
  )
  seurat_obj
}

naive_CAF <- classify_ENG(naive_CAF)
treated_CAF <- classify_ENG(treated_CAF)

cat("Naive CAFs\n")
print(table(naive_CAF$celltype, naive_CAF$ENG_status))
cat("\nTreated CAFs\n")
print(table(treated_CAF$celltype, treated_CAF$ENG_status))
