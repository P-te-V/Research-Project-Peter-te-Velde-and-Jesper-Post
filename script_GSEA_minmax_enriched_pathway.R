library(readr)
library(writexl)
#this example is for the GB1 sample. this can be changed per sample
GB1_GSEA_scores <- read.table("/home4/s4111419/Documents/GB1/GB1_enrichment_matrix_h.tsv", sep = '\t', header = TRUE, fill = TRUE)

cell_type_combinations <- colnames(GB1_GSEA_scores)

# Initialize lists to store results
max_enriched_pathways <- list()
min_enriched_pathways <- list()
pathway_max <- list()
pathway_min <- list()

for (combination in 1:length(cell_type_combinations)) {
  
  column_data <- as.numeric(GB1_GSEA_scores[, combination])
  
  # get the cell in the table which contains the highest positive and negative enrichment score
  max_index <- which.max(column_data)
  min_index <- which.min(column_data)
  
  # get out the rowname 
  max_enriched_pathway <- rownames(GB1_GSEA_scores)[max_index]
  min_enriched_pathway <- rownames(GB1_GSEA_scores)[min_index]
  
  # get the actual name of the geneset + link to the MSIGDB 
  max_first_column_value <- GB1_GSEA_scores[max_index, 1]
  min_first_column_value <- GB1_GSEA_scores[min_index, 1]
  
  # add the name to the list
  max_enriched_pathways[[cell_type_combinations[combination]]] <- max_enriched_pathway
  min_enriched_pathways[[cell_type_combinations[combination]]] <- min_enriched_pathway
  
  # add the name to the list to put in the table 
  pathway_max[[cell_type_combinations[combination]]] <- max_first_column_value
  pathway_min[[cell_type_combinations[combination]]] <- min_first_column_value
}
#  make a dataframe and transpose the rows and columns. this is a personal preference making the table easier to read
# + add the names to the dataframe 
min_max_GSEA <- data_frame(pathway_max, pathway_min)
min_max_GSEA <- t(min_max_GSEA)

# the [-1] is to remove the first column name, which we will not use
colnames(min_max_GSEA)<-cell_type_combinations[-1]
min_max_GSEA<-as.data.frame(min_max_GSEA)
min_max_GSEA_GB1 <- t(min_max_GSEA)
min_max_GSEA_GB1<-as.data.frame(min_max_GSEA_GB1)

# export the table as an excel file
write.xlsx(min_max_GSEA_GB1,"/home4/s4111419/Documents/GB1/GB1_minmax_table.xlsx", rowNames=TRUE)



