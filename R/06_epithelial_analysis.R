# GSE164983 LUAD scRNA-seq
# 06 - Epithelial-associated transcriptional programme

library(Seurat)
library(dplyr)

luad_clean <- readRDS("results/GSE164983_annotated.rds")
Idents(luad_clean) <- "cell_type_broad"
print(table(Idents(luad_clean)))

# Epithelial cells are compared with all other annotated cells.
# These are epithelial-associated genes, not tumour-vs-normal DE genes.
epithelial_markers <- FindMarkers(
  luad_clean, ident.1="Epithelial", only.pos=TRUE,
  min.pct=0.25, logfc.threshold=0.5
)
epithelial_markers$gene <- rownames(epithelial_markers)
epithelial_markers <- epithelial_markers %>% arrange(desc(avg_log2FC))

print(head(epithelial_markers[,c("gene","avg_log2FC","pct.1","pct.2","p_val_adj")],30))
cat("\nNumber of epithelial-associated markers:", nrow(epithelial_markers), "\n")
write.csv(epithelial_markers, "results/epithelial_marker_genes.csv", row.names=FALSE)

epithelial_pathway_genes <- epithelial_markers %>%
  filter(avg_log2FC >= 1.5, pct.1 >= 0.25, p_val_adj < 0.05) %>% pull(gene)
cat("Genes retained for pathway enrichment:", length(epithelial_pathway_genes), "\n")

write.table(epithelial_pathway_genes, "results/epithelial_pathway_gene_list.txt",
            row.names=FALSE, col.names=FALSE, quote=FALSE)
