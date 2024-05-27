# Loading of Data
list_of_files =list.files("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Peter_Jesper/Data/tables_differential_expression_analysis/")

# Change the gene name to a entrezid so that GSEA can be run
change_gene_symbol_to_entrezid = function(file_name){
  de_scores = data.frame(fread(paste0("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Peter_Jesper/Data/tables_differential_expression_analysis/",file_name)), row.names = 1)
  de_scores_row_sd = apply(de_scores,1, sd)
  de_scores_column_sd = apply(de_scores,2, sd)
  de_scores = de_scores[which(de_scores_row_sd >0),which(de_scores_column_sd>0) ]
  
  gene_annotation_file = data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Peter_Jesper/Data/hgnc_gene_anno_merged.tsv"), row.names = 1)
  gene_annotation_file_unique_entrezids = names(which(table(gene_annotation_file$NCBI_id)>1))
  sum(is.na(gene_annotation_file$NCBI_id))
  
  
  gene_annotation_file = gene_annotation_file[which(!is.na(gene_annotation_file$NCBI_id)),]
  rownames(gene_annotation_file) = gene_annotation_file$HGNC_symbol
  
  common_genes = intersect(gene_annotation_file$HGNC_symbol, rownames(de_scores))
  
  gene_annotation_file = gene_annotation_file[common_genes, ]
  de_scores = de_scores[common_genes, ]
  
  rownames(de_scores) = gene_annotation_file$NCBI_id
  file_name_updated = paste0("Entrezid_",file_name)
  write.table(de_scores, file = paste0("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Peter_Jesper/Data/tables_differential_expression_analysis/",file_name_updated), sep = "\t", quote = FALSE)
}

sapply(list_of_files,change_gene_symbol_to_entrezid )