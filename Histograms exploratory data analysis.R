# loading of samples
list_of_all_the_st_samples = list.files("C:/Users/pteve/OneDrive/Datasets Research Projects/")

# loop function for exploratory data analysis to save in tab separated file
for (ST_Sample in list_of_all_the_st_samples[ST_Sample]) {
  
  sample_seurat <- readRDS(sprintf("C:/Users/pteve/OneDrive/Datasets Research Projects/%s", ST_Sample))

  count_ST_Sample <- sample_seurat@assays$Spatial@counts
  count_ST_Sample_matrix  <- as.matrix(count_ST_Sample)
  ST_Sample_shortened <- gsub(".rds", "", ST_Sample)
  coordinates_ST_Sample <- sample_seurat@images[[ST_Sample_shortened]]@coordinates
  Mean_spots      <- colMeans(count_ST_Sample_matrix)
  Median_spots    <- apply(count_ST_Sample_matrix, 2, median)
  Max_spots       <- apply(count_ST_Sample_matrix,2,max)
  Min_spots       <- apply(count_ST_Sample_matrix,2,min)
  NA_spots        <- colSums(is.na(count_ST_Sample_matrix))
  temp_df <- data.frame(Mean_spots, Median_spots, Max_spots, Min_spots, NA_spots)
  write.table(temp_df_spots, sprintf("C:/Users/pteve/OneDrive/Opslag_voor_ResearchProjec/%s.tsv", ST_Sample_shortened), sep = "\t")
}

# loop function for exploratory data analysis to save in tab separated file
for (ST_Sample in list_of_all_the_st_samples[ST_Sample]) {
  sample_seurat <- readRDS(sprintf("C:/Users/pteve/OneDrive/Datasets Research Projects/%s", ST_Sample))
  count_ST_Sample <- sample_seurat@assays$Spatial@counts
  count_ST_Sample_matrix <- as.matrix(count_i)
  ST_Sample_shortened <- gsub(".rds", "", ST_Sample)
  coordinates_ST_Sample <- sample_seurat@images[[ST_Sample_shortened]]@coordinates
  Mean_genes      <- rowMeans(count_ST_Sample_matrix)
  Median_genes    <- apply(count_ST_Sample_matrix, 1, median)
  Max_genes       <- apply(count_ST_Sample_matrix,1,max)
  Min_genes       <- apply(count_ST_Sample_matrix,1,min)
  NA_genes        <- rowSums(is.na(count_ST_Sample_matrix))
  temp_df <- data.frame(Mean_genes, Median_genes, Max_genes, Min_genes, NA_genes)
  write.table(temp_df_genes, sprintf("C:/Users/pteve/OneDrive/Opslag_voor_ResearchProjec/x columns/%s.tsv", ST_Sample_shortened), sep = "\t")
}

# for the histograms of genes ------------------------------------------------------------------------------------
colnames(temp_df_genes) <- cbind("Min_genes", "Max_genes", "Mean_genes", "Median_genes", "NA_genes")
par(mfrow = c(2,5))
hists_per_sample <- list()
for (Gene_Measurement in colnames(temp_df_genes)) {
  temp_hist <- ggplot(temp_df_genes, aes(x = temp_df_genes[[Gene_Measurement]]))+
    geom_histogram(bins = 30, fill = "darkgreen", color = "black", alpha = 0.5)+
    labs(title = Gene_Measurement, 
         x = Gene_Measurement, 
         y = "Frequency")+
    theme_minimal()
  temp_grob <- ggplotGrob(temp_hist)  # Convert ggplot object to grob
  hists_per_sample <- c(hists_per_sample, list(temp_grob))
}

# combined_plot <- do.call(grid.arrange, c(hists_per_sample, ncol = 2))
hists <- append(hists, list(combined_plot))
pdf(sprintf("C:/Users/pteve/OneDrive/Datasets Research Projects/R codes/exploratory plots/%s_combined_plot.pdf", i_shortened))
print(do.call(grid.arrange, c(hists_per_sample, ncol = 2)))
dev.off()

# for the histograms of spots --------------------------
par(mfrow = c(2,5))
hists_per_sample_2 <- list()
for (Spot_Measurements in colnames(temp_df_spots)) {
  temp_hist_2 <- ggplot(temp_df_spots, aes(x = temp_df_spots[[Spot_Measurements]]))+
    geom_histogram(bins = 30, fill = "purple", color = "purple", alpha = 0.5)+
    labs(title = Spot_Measurements, 
         x = Spot_Measurements, 
         y = "Frequency")+
    theme_minimal()
  temp_grob_2 <- ggplotGrob(temp_hist_2)  # Convert ggplot object to grob
  hists_per_sample_2 <- c(hists_per_sample_2, list(temp_grob_2))
}
# combined_plot_2 <- do.call(grid.arrange, c(hists_per_sample_2, ncol = 2))
hists <- append(hists, list(combined_plot_2))
pdf(sprintf("C:/Users/jespe/OneDrive/RUG/bachelors_research_project/plots/%s_combined_plot.pdf_2", i_shortened))
print( do.call(grid.arrange, c(hists_per_sample_2, ncol = 2)))
dev.off()
