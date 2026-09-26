# GSE164983 LUAD scRNA-seq
# 01_load_data.R
# Load processed 10x H5 matrices, create Seurat objects, add metadata, and merge samples.

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(hdf5r)

# Run scripts from the repository root.
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

sample1_file <- "data/GSM5024081_591146_filtered_feature_bc_matrix.h5"
sample2_file <- "data/GSM5024082_614658_filtered_feature_bc_matrix.h5"

if (!all(file.exists(sample1_file, sample2_file))) {
  stop("Processed GEO H5 files are missing. Download the two GSE164983 filtered_feature_bc_matrix.h5 files and place them in data/.")
}

sample1_counts <- Read10X_h5(sample1_file, use.names = TRUE, unique.features = TRUE)
sample2_counts <- Read10X_h5(sample2_file, use.names = TRUE, unique.features = TRUE)

sample1 <- CreateSeuratObject(counts = sample1_counts, project = "LUAD_591146", min.cells = 3, min.features = 200)
sample2 <- CreateSeuratObject(counts = sample2_counts, project = "LUAD_614658", min.cells = 3, min.features = 200)

sample1$sample_id <- "591146"
sample2$sample_id <- "614658"
sample1$condition <- "LUAD"
sample2$condition <- "LUAD"
sample1$dataset <- "GSE164983"
sample2$dataset <- "GSE164983"

cat("Sample 591146:", ncol(sample1), "cells;", nrow(sample1), "genes\n")
cat("Sample 614658:", ncol(sample2), "cells;", nrow(sample2), "genes\n")

luad <- merge(x = sample1, y = sample2, add.cell.ids = c("591146", "614658"), project = "GSE164983_LUAD")

cat("Total cells before QC:", ncol(luad), "\n")
cat("Genes in merged object:", nrow(luad), "\n")
print(table(luad$sample_id))

saveRDS(luad, "results/GSE164983_merged_preQC.rds")
