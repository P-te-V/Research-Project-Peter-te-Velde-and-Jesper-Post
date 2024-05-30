# loading of the refrence data
ref.data <- celldex::HumanPrimaryCellAtlasData()

# loading of the processed data
list_pre_processing <- list.files("/home4/s4617487/pre processed data")
setwd("/home4/s4617487/pre processed data")

# Cell type prediction and visualization
for(processed_data in list_pre_processing[processed_data]){
  import.data <- read.table(sprintf("/home4/s4617487/pre processed data/%s", processed_data))
  import.data1 <- t(import.data)
  processed_data_shortened <- gsub("_pre_processed.tsv", "", processed_data)
  
  Test1 <- SingleR(test = import.data1,
                   assay.type.test =1,
                   ref = ref.data,
                   labels = ref.data$label.main)
  
  df_assignedcelltypes <- data.frame(Test1[,4, drop = F])
  
  df_assignedcelltypes$pruned.labels <- factor(df_assignedcelltypes$pruned.labels,
                                               levels = names(sort(table(df_assignedcelltypes$pruned.labels), decreasing = TRUE)))
  
  write.table(df_assignedcelltypes, sprintf("/home4/s4617487/Table with cell types/%s_cell_type_table",processed_data_shortened), sep = "\t")
  
  pdf(sprintf("/home4/s4617487/PDF celltypes/%s_CellTypeInformationGraphs.pdf", i_shortened))
  plotScoreHeatmap(Test1)
  plotDeltaDistribution(Test1, ncol = 4, dots.on.top = F)
  ggplot(df_assignedcelltypes, aes(x = pruned.labels)) +
    geom_histogram(fill = "red", alpha = 0.5, stat = "count", color = "black") +
    labs(title = processed_data_shortened, x = "celltypes", y = "frequency") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()
}
