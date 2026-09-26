# GSE164983 LUAD scRNA-seq
# 07 - Functional pathway enrichment

library(dplyr)
library(ggplot2)
library(gprofiler2)

epithelial_markers <- read.csv("results/epithelial_marker_genes.csv")
epithelial_pathway_genes <- epithelial_markers %>%
  filter(avg_log2FC >= 1.5, pct.1 >= 0.25, p_val_adj < 0.05) %>% pull(gene)
cat("Genes submitted for pathway enrichment:", length(epithelial_pathway_genes), "\n")

gost_results <- gost(
  query=epithelial_pathway_genes, organism="hsapiens",
  ordered_query=FALSE, multi_query=FALSE, significant=TRUE,
  correction_method="fdr", sources=c("GO:BP","REAC","KEGG")
)
pathway_results <- gost_results$result %>% arrange(p_value)

# Convert g:Profiler list columns to delimited strings for CSV export.
pathway_export <- pathway_results
list_cols <- sapply(pathway_export, is.list)
pathway_export[list_cols] <- lapply(pathway_export[list_cols],
  function(x) sapply(x, function(y) paste(y, collapse=";")))
write.csv(pathway_export, "results/epithelial_pathway_enrichment.csv", row.names=FALSE)

reactome_top20 <- pathway_results %>% filter(source=="REAC") %>% arrange(p_value) %>%
  select(term_name,p_value,intersection_size,term_size) %>% head(20)
kegg_top20 <- pathway_results %>% filter(source=="KEGG") %>% arrange(p_value) %>%
  select(term_name,p_value,intersection_size,term_size) %>% head(20)
print(as.data.frame(reactome_top20))
print(as.data.frame(kegg_top20))

# Representative significant pathways selected to summarize major biological
# themes while reducing redundancy; this is not an algorithmic top-10 ranking.
selected_pathways <- pathway_results %>%
  filter(term_name %in% c(
    "Cell junction organization","Extracellular matrix organization",
    "Non-integrin membrane-ECM interactions","RHO GTPase cycle",
    "Signaling by Receptor Tyrosine Kinases","Tight junction","Adherens junction",
    "Integrin signaling","Focal adhesion","ECM-receptor interaction"
  )) %>%
  mutate(log10FDR=-log10(p_value), pathway_label=paste0(term_name," (",source,")")) %>%
  arrange(log10FDR)

selected_export <- selected_pathways %>%
  select(source,term_name,p_value,intersection_size,term_size,log10FDR)
write.csv(selected_export, "results/selected_epithelial_pathways.csv", row.names=FALSE)

pathway_plot <- ggplot(selected_pathways,
  aes(x=reorder(pathway_label,log10FDR), y=log10FDR)) +
  geom_col() + coord_flip() +
  labs(x=NULL, y=expression(-log[10]("adjusted P-value")),
       title="Pathway enrichment in LUAD epithelial cells") +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold",size=14,hjust=0.5),
        axis.text.y=element_text(size=10),axis.title.x=element_text(size=11))

ggsave("figures/08_epithelial_pathway_enrichment.png", pathway_plot,
       width=10,height=6,dpi=300,bg="white")
ggsave("figures/08_epithelial_pathway_enrichment.pdf", pathway_plot,width=10,height=6)

# Interpretation: enrichment is relative to other cell populations, not normal lung.
