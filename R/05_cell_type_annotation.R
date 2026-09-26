# GSE164983 LUAD scRNA-seq
# 05 - Cell-type annotation

library(Seurat)
library(dplyr)
library(ggplot2)

luad_clean <- readRDS("results/GSE164983_clustered_unintegrated.rds")
Idents(luad_clean) <- "seurat_clusters"

markers <- FindAllMarkers(luad_clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "results/cluster_marker_genes.csv", row.names = FALSE)

cluster_annotations <- c(
  "0"="T cells","1"="Activated CD4 T","2"="Cytotoxic T","3"="Epithelial",
  "4"="T cells","5"="Epithelial","6"="CD8 tissue-resident T","7"="B-like immune",
  "8"="Treg","9"="NK","10"="Cytotoxic NK","11"="Mast","12"="Plasma",
  "13"="B cells","14"="Treg","15"="Gamma-delta T/NK-like","16"="Epithelial",
  "17"="Macrophage/Monocyte","18"="Gamma-delta T","19"="Interferon-responsive",
  "20"="Endothelial","21"="Neutrophil","22"="Fibroblast","23"="pDC"
)

cell_type_vector <- unname(cluster_annotations[as.character(luad_clean$seurat_clusters)])
stopifnot(length(cell_type_vector) == ncol(luad_clean), sum(is.na(cell_type_vector)) == 0)
luad_clean$cell_type <- cell_type_vector
luad_clean$cell_type_broad <- luad_clean$cell_type

luad_clean$cell_type_broad[luad_clean$cell_type %in% c("T cells","Activated CD4 T")] <- "CD4/Conventional T"
luad_clean$cell_type_broad[luad_clean$cell_type %in% c("Cytotoxic T","CD8 tissue-resident T","Gamma-delta T","Gamma-delta T/NK-like")] <- "Cytotoxic T"
luad_clean$cell_type_broad[luad_clean$cell_type %in% c("NK","Cytotoxic NK")] <- "NK"
luad_clean$cell_type_broad[luad_clean$cell_type %in% c("B cells","B-like immune")] <- "B cells"
luad_clean$cell_type_broad[luad_clean$cell_type == "Interferon-responsive"] <- "Interferon-responsive T"

print(table(luad_clean$cell_type_broad))
print(table(luad_clean$cell_type_broad, luad_clean$sample_id))

p_celltype <- DimPlot(luad_clean, reduction="umap", group.by="cell_type_broad", label=TRUE, repel=TRUE, pt.size=0.25) +
  ggtitle("Cellular landscape of primary lung adenocarcinoma") + theme_classic() +
  theme(plot.title=element_text(face="bold", hjust=0.5), legend.title=element_blank())
ggsave("figures/05_UMAP_cell_types.png", p_celltype, width=10, height=7, dpi=300, bg="white")

broad_composition <- as.data.frame(table(sample_id=luad_clean$sample_id, cell_type_broad=luad_clean$cell_type_broad))
colnames(broad_composition) <- c("sample_id","cell_type_broad","cell_count")
broad_composition <- broad_composition %>% group_by(sample_id) %>%
  mutate(proportion = cell_count / sum(cell_count)) %>% ungroup()
write.csv(broad_composition, "results/cell_type_broad_composition_by_sample.csv", row.names=FALSE)

p_composition <- ggplot(broad_composition, aes(x=sample_id, y=proportion, fill=cell_type_broad)) +
  geom_col(width=0.75) + scale_y_continuous(labels=scales::percent) +
  labs(x="Sample", y="Cell proportion", fill="Cell type", title="Cellular composition of LUAD samples") +
  theme_classic(base_size=12) + theme(plot.title=element_text(face="bold", hjust=0.5))
ggsave("figures/06_cell_type_composition.png", p_composition, width=9, height=6, dpi=300, bg="white")

marker_panel <- c("CD3D","IL7R","CD8A","CCL5","FOXP3","IL2RA","NKG7","GNLY",
                  "CD79A","MS4A1","JCHAIN","MZB1","LYZ","C1QA","TPSAB1","KIT",
                  "EPCAM","KRT19","VWF","EMCN","DCN","SFRP2","FCGR3B","S100A8",
                  "CLEC4C","LILRA4","ISG15","MX1")
p_markers <- DotPlot(luad_clean, features=marker_panel, group.by="cell_type_broad") +
  RotatedAxis() + labs(x="Marker gene", y="Annotated cell type",
  title="Canonical marker expression across annotated cell populations") +
  theme_classic(base_size=11) + theme(plot.title=element_text(face="bold", hjust=0.5))
ggsave("figures/07_cell_type_marker_dotplot.png", p_markers, width=15, height=7, dpi=300, bg="white")

saveRDS(luad_clean, "results/GSE164983_annotated.rds")
