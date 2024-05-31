library(dplyr)
library(magrittr)
library(parallelly)
library(data.table)
library(parallel)

#load the samples, removal of files that are not needed
List_CellTypes <- list.files("/home4/s4111419/Documents/annotations")
List_CellTypes <-List_CellTypes[-c(22,23)]
List_ProcessedData <- list.files("/home4/s4111419/Documents/preprocessed_data")

# number of cores that is going to be used for the parallel analysis
n.cores <- availableCores()
ST
# first loop, through all the samples
for (ST_sample in List_ProcessedData[1:21]){
  
  # load the data from the documents in the Hábrók environment
  Test1 <- data.frame(fread(sprintf("/home4/s4111419/Documents/preprocessed_data/%s", ST_sample)), row.names = 1)
  cell_type_file = paste0(strsplit(ST_sample, "_")[[1]][1], "_cell_type_table")
  Test2 <- data.frame(fread(sprintf("/home4/s4111419/Documents/annotations/%s", cell_type_file)), row.names = 1)
  
  # merging the cell type annotations with the preprocessed data
  common_spot_names = intersect(rownames(Test1),rownames(Test2))
  Test1 = Test1[common_spot_names,]
  Test2 = as.data.frame(Test2[common_spot_names,])
  rownames(Test2) = common_spot_names
  colnames(Test2) = "pruned.labels"
  cell_type_gen_expr <- merge(Test1, Test2, by.x = 0, by.y = 0)
  cell_type_gen_expr <- cell_type_gen_expr %>% relocate(pruned.labels)
  ST_sample_shortened <- paste0(gsub(".tsv", "", ST_sample),"_DE_score")
  
  # setting all the variables
  # all cell types in one vector
  all_cell_types <- unique(cell_type_gen_expr$pruned.labels)
  # remove the NAs
  all_cell_types <- all_cell_types[!is.na(all_cell_types)]
  # remove the names of the first two columns
  column_name <-names(cell_type_gen_expr)[c(-1,-2)]
  # number of cell type combinations
  number_of_genes <- length(column_name)
  number_of_celltype_combinations = length(all_cell_types)*(length(all_cell_types)-1)/2
  count = 1
  # create an empty matrix for the scores
  scores_combined = matrix(NA, number_of_genes, number_of_celltype_combinations)
  # generating the row_names of the data frame from the column (gene) names of the preprocessed data file
  rownames(scores_combined) = column_name
  colnames(scores_combined) =  character(number_of_celltype_combinations)
  cl <- makeCluster(n.cores)
  
  # loop for the cell_type_cell_type combinations, two loops are used to make the combinations

   for (cell_type_1 in 1:(length(all_cell_types)-1)){
      for (cell_type_2 in (cell_type_1+1):length(all_cell_types)){
        # assign both cell types
          cell_type1 <- all_cell_types[cell_type_1]
          cell_type2 <- all_cell_types[cell_type_2]
          # data for both cell types
          data1 <- cell_type_gen_expr[cell_type_gen_expr$pruned.labels == cell_type1, column_name]
          data2 <- cell_type_gen_expr[cell_type_gen_expr$pruned.labels == cell_type2, column_name]
          #loop for each individual gene, per cell_type_cell_type combination
          gene_loop_de_score <- function(genes) {
          for(genes in column_name){
            gene <- genes
            #  the actual test which compares the gene expression between the cells 
            res = wilcox.test(data1[,genes], data2[,genes])
            # transformation of the p-value the wilcox.test yields
            val <- (-log10(res$p.value))*(sign(median(data1[,genes], na.rm = TRUE)- median(data2[,genes], na.rm = TRUE)))
           return(val)
          }
            # start up the cluster to let it run in parallel
          clusterEvalQ(cl, {
            library(data.table)
            library(dplyr)
          })
          clusterExport(cl, c("gene_loop_de_score", "median", "wilcox.test", "data1", "data2"))
          results <- parLapply(cl, column_name, gene_loop_de_score)

          
          scores_combined[,count] = unlist(results)
          colnames(scores_combined)[count]<- paste(cell_type1, cell_type2, sep = "_")
        
        count = count + 1
        
        # to check if the script runs
        print(ST_sample)
        print(cell_type_1)
        print(cell_type_2)
        print(proc.time())
    }
   }
  # saving the generated data in a tsv.file
  write.table(scores_combined, sprintf("/home4/s4111419/Documents/tables_differential_expression_analysis/%s.tsv", ST_sample_shortened), sep = "\t")
  
  stop(cl)
  
}


