# Loading of gene set enrichment results
GSEA_result <- read.table("C:/Users/pteve/OneDrive/Opslag_voor_ResearchProjec/GSEA Welsh test/GB1_enrichment_matrix_h.tsv")
list_of_GSEA_results <- list.files("C:/Users/pteve/OneDrive/Opslag_voor_ResearchProjec/GSEA Welsh test")

# For loop to make heatmap out of GSEA_result
for (gene_enrichments in list_of_GSEA_results) {
  table_gene_enrichments <- read.table(sprintf("C:/Users/pteve/OneDrive/Opslag_voor_ResearchProjec/GSEA Welsh test/%s", gene_enrichments), row.names = 1, sep="\t", header=TRUE, fill=TRUE)
  names_of_row <- rownames(table_gene_enrichments)
  altered_names_of_row <- sapply(names_of_row, function(x){y <- strsplit(x, " -- ")[[1]][1]})
  rownames(table_gene_enrichments) <- altered_names_of_row
  matrix_table_gene_enrichments <- as.matrix(table_gene_enrichments)
  gene_enrichments_shortened <- gsub("_enrichment_matrix_h.tsv", "", gene_enrichments)
  pdf(sprintf("C:/Users/pteve/OneDrive/Opslag_voor_ResearchProjec/Heatmap Welsh T/%s_heatmaps_gene_enrichment", gene_enrichments_shortened))
  heatmap(matrix_table_gene_enrichments,)
  dev.off()
}