# GSE164983 LUAD scRNA-seq
# 02_quality_control.R
# Calculate QC metrics, inspect distributions, and apply dataset-informed filtering.

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

luad <- readRDS("results/GSE164983_merged_preQC.rds")
luad[["percent.mt"]] <- PercentageFeatureSet(luad, pattern = "^MT-")

qc_summary <- luad@meta.data %>% group_by(sample_id) %>%
  summarise(cells=n(), median_genes=median(nFeature_RNA), mean_genes=mean(nFeature_RNA),
            median_counts=median(nCount_RNA), mean_counts=mean(nCount_RNA),
            median_mt=median(percent.mt), mean_mt=mean(percent.mt),
            min_genes=min(nFeature_RNA), max_genes=max(nFeature_RNA),
            min_mt=min(percent.mt), max_mt=max(percent.mt), .groups="drop")
print(qc_summary)
write.csv(qc_summary, "results/QC_summary_before_filtering.csv", row.names=FALSE)

qc_plot <- VlnPlot(luad, features=c("nFeature_RNA","nCount_RNA","percent.mt"),
                   group.by="sample_id", ncol=3, pt.size=0.05)
ggsave("figures/01_QC_before_filtering.png", qc_plot, width=12, height=5, dpi=300, bg="white")

plot1 <- FeatureScatter(luad, feature1="nCount_RNA", feature2="percent.mt")
plot2 <- FeatureScatter(luad, feature1="nCount_RNA", feature2="nFeature_RNA")
ggsave("figures/02_QC_scatter_before_filtering.png", plot1 + plot2, width=11, height=5, dpi=300, bg="white")

# Dataset-informed thresholds selected after examining sample-specific distributions.
luad$QC_status <- ifelse(luad$nFeature_RNA >= 300 & luad$nFeature_RNA <= 6000 & luad$percent.mt <= 20, "Pass", "Fail")

qc_retention <- luad@meta.data %>% group_by(sample_id) %>%
  summarise(before=n(), passed=sum(QC_status=="Pass"), removed=sum(QC_status=="Fail"),
            retention_percent=round(100*passed/before,2), .groups="drop")
print(qc_retention)
write.csv(qc_retention, "results/QC_cell_retention.csv", row.names=FALSE)

luad_qc <- subset(luad, subset=nFeature_RNA >= 300 & nFeature_RNA <= 6000 & percent.mt <= 20)
cat("Cells before QC:", ncol(luad), "\n")
cat("Cells after QC:", ncol(luad_qc), "\n")
cat("Cells removed:", ncol(luad)-ncol(luad_qc), "\n")
cat("Overall retention:", round(100*ncol(luad_qc)/ncol(luad),2), "%\n")

qc_after <- VlnPlot(luad_qc, features=c("nFeature_RNA","nCount_RNA","percent.mt"),
                    group.by="sample_id", ncol=3, pt.size=0.05)
ggsave("figures/03_QC_after_filtering.png", qc_after, width=12, height=5, dpi=300, bg="white")

saveRDS(luad_qc, "results/GSE164983_QC_filtered.rds")
capture.output(sessionInfo(), file="results/sessionInfo.txt")
