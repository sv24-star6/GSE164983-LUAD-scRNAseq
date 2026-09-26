# GSE164983 LUAD scRNA-seq
# 04 - Normalization, dimensionality reduction and clustering

library(Seurat)
library(ggplot2)

luad_clean <- readRDS("results/GSE164983_singlets.rds")

luad_clean <- NormalizeData(luad_clean, normalization.method = "LogNormalize", scale.factor = 10000)
luad_clean <- FindVariableFeatures(luad_clean, selection.method = "vst", nfeatures = 2000)
luad_clean <- ScaleData(luad_clean, features = VariableFeatures(luad_clean))
luad_clean <- RunPCA(luad_clean, features = VariableFeatures(luad_clean), npcs = 50)

# Sample-associated structure was retained rather than forcing integration because
# patient biology and technical sample effects cannot be reliably separated with n = 2.
luad_clean <- FindNeighbors(luad_clean, dims = 1:30)
luad_clean <- FindClusters(luad_clean, resolution = 0.5, random.seed = 123)
luad_clean <- RunUMAP(luad_clean, dims = 1:30, seed.use = 123)

p_sample <- DimPlot(luad_clean, reduction = "umap", group.by = "sample_id", pt.size = 0.3) +
  ggtitle("LUAD single-cell landscape by sample") +
  theme_classic()

ggsave("figures/04_UMAP_by_sample.png", p_sample, width = 8, height = 6, dpi = 300, bg = "white")

# Seurat v5: join expression layers before marker testing.
luad_clean <- JoinLayers(luad_clean)
saveRDS(luad_clean, "results/GSE164983_clustered_unintegrated.rds")
