# Loading spatial transcriptomics data
list_of_all_the_st_samples <- list.files("C:/Users/pteve/OneDrive/Datasets Research Projects/")
setwd("C:/Users/pteve/OneDrive/Datasets Research Projects/")

# Processing of the data
for (ST_Sample in list_of_all_the_st_samples[ST_Sample]) {
  pre_sample_seurat <- readRDS(sprintf("C:/Users/pteve/OneDrive/Datasets Research Projects/%s", ST_Sample))
  ST_Sample_shortened <- gsub(".rds", "", ST_Sample) 
  raw_counts = as.matrix(pre_sample_seurat@assays$Spatial@counts) 
  raw_counts_pca = raw_counts[rowSums(raw_counts)!=0,]
  raw_counts_pca = Seurat::LogNormalize(raw_counts_pca)
  raw_counts_pca_sd = apply(raw_counts_pca,2,sd)
  raw_counts_pca = raw_counts_pca[,which(raw_counts_pca_sd>0)]
  pca_obj_princomp = princomp(raw_counts_pca, cor=T, fix_sign = F)  
  raw_counts_pc1removed = pca_obj_princomp$loadings[,-1] %*% t(pca_obj_princomp$scores[,-1])
  write.table(raw_counts_pc1removed, sprintf("C:/Users/pteve/OneDrive/Opslag_voor_ResearchProjec/normalized data/%s.tsv", ST_Sample_shortened), sep = "\t")
}

