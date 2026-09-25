# GSE164983 LUAD scRNA-seq
# 03_doublet_detection.R
# Detect doublets separately within each 10x sample using scDblFinder.

library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(dplyr)

project_dir <- "C:/Users/Admin/Downloads/scrna project"
setwd(project_dir)

dir.create("results", showWarnings = FALSE)

luad_qc <- readRDS("results/GSE164983_QC_filtered.rds")

# Each 10x capture is processed independently for doublet detection.
luad_split <- SplitObject(luad_qc, split.by = "sample_id")

sce_591146 <- as.SingleCellExperiment(luad_split[["591146"]])
sce_614658 <- as.SingleCellExperiment(luad_split[["614658"]])

set.seed(123)
sce_591146 <- scDblFinder(sce_591146)

set.seed(123)
sce_614658 <- scDblFinder(sce_614658)

# Transfer classifications and scores back to Seurat metadata.
luad_split[["591146"]]$doublet_class <- sce_591146$scDblFinder.class
luad_split[["591146"]]$doublet_score <- sce_591146$scDblFinder.score

luad_split[["614658"]]$doublet_class <- sce_614658$scDblFinder.class
luad_split[["614658"]]$doublet_score <- sce_614658$scDblFinder.score

doublet_summary <- data.frame(
  sample_id = c("591146", "614658"),
  total_cells = c(ncol(sce_591146), ncol(sce_614658)),
  singlets = c(
    sum(sce_591146$scDblFinder.class == "singlet"),
    sum(sce_614658$scDblFinder.class == "singlet")
  ),
  predicted_doublets = c(
    sum(sce_591146$scDblFinder.class == "doublet"),
    sum(sce_614658$scDblFinder.class == "doublet")
  )
)

doublet_summary$doublet_percent <- round(
  100 * doublet_summary$predicted_doublets / doublet_summary$total_cells,
  2
)

print(doublet_summary)
write.csv(doublet_summary, "results/doublet_summary.csv", row.names = FALSE)

# Retain predicted singlets.
sample1_clean <- subset(
  luad_split[["591146"]],
  subset = doublet_class == "singlet"
)

sample2_clean <- subset(
  luad_split[["614658"]],
  subset = doublet_class == "singlet"
)

luad_clean <- merge(
  sample1_clean,
  y = sample2_clean,
  project = "GSE164983_LUAD_clean"
)

cat("QC-passed cells before doublet removal:", ncol(luad_qc), "\n")
cat("Predicted doublets removed:", sum(doublet_summary$predicted_doublets), "\n")
cat("Final singlets retained:", ncol(luad_clean), "\n")
print(table(luad_clean$sample_id))

saveRDS(luad_clean, "results/GSE164983_singlets.rds")
