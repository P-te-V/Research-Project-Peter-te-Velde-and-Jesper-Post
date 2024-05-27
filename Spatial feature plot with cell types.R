# Loading of the cell types into table 
CellTypeTable <- read.table("C:/Users/pteve/OneDrive/Opslag_voor_ResearchProjec/Table cell type identification/GB1_cell_type_table")
CellTypeTable$pruned.labels[is.na(CellTypeTable$pruned.labels)] = "unknown"
CellTypeTable$numericvalue <- as.numeric(factor(unlist(CellTypeTable$pruned.labels)))
CellTypeTable$sampleID <- rownames(CellTypeTable) 

# Chaning the cell types table in to a wideformat
wideformat_celltypetable <- dcast(CellTypeTable, sampleID ~ pruned.labels)
rownames(wideformat_celltypetable) <- wideformat_celltypetable$sampleID
wideformat_celltypetable[!is.na(wideformat_celltypetable)] = 1
wideformat_celltypetable[is.na(wideformat_celltypetable)] = 0
wideformat_celltypetable$sampleID <- NULL
wideformat_celltypetable[] <- lapply(wideformat_celltypetable, as.numeric)

# Loading of the spatial data and adding the cell type data as metadata
SeuratData <- readRDS("C:/Users/pteve/OneDrive/Opslag_voor_ResearchProjec/Original Data/CRC3.rds")
SeuratData <- AddMetaData(SeuratData, metadata = wideformat_celltypetable)

# Making a spatial feature plot out of the new mata data
plotlist <- lapply(colnames(wideformat_celltypetable[,c(3,4,9,12,13,14,15,16,22)]), function(Cell_Type){
  p <- SpatialFeaturePlot(SeuratData, 
                     features = Cell_Type,
                     image.alpha = 0.3,
                     stroke = NA,
                     alpha = 1)&ggtitle(Cell_Type)
  p + theme(legend.position = "none")
  })
plot_grid(plotlist = plotlist, ncol = 3)

