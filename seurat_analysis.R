library(SingleCellExperiment)
library(Seurat)
library(SeuratDisk)

## Converting  h5ad file to seurat readble file
# This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory
Convert("submission_210120.h5ad", "seurat.h5seurat")

# reading converted file
skin_data <- LoadH5Seurat("seurat.h5seurat")


### adding individual ids to metadata
full_metadata <- read.delim(file = "E-MTAB-8142_sample_table.txt", sep = "\t", stringsAsFactors = F)

sample_to_patient <- setNames(object = full_metadata$Characteristics.individual., nm = full_metadata$Source.Name)

skin_data_seurat@meta.data$ind_id <- sample_to_patient[as.character(skin_data_seurat@meta.data$sample_id)]

### log normalizing data
skin_data_seurat <- NormalizeData(skin_data, normalization.method = "LogNormalize", scale.factor = 10000)
### detecting variable features to continue with reduced data
skin_data_seurat <- FindVariableFeatures(skin_data_seurat)
### scaling(centrilizing) dta
skin_data_seurat <- ScaleData(skin_data_seurat)
### running dimensional reduction with PCA
skin_data_seurat <- RunPCA(skin_data_seurat)
### choosing number of dimensions with the help of elbow plot
ElbowPlot(skin_data_seurat, ndims = 50)
skin_data_seurat <- FindNeighbors(skin_data_seurat, dims = 1:30)
### detecting clusters
skin_data_seurat <- FindClusters(skin_data_seurat, resolution = 0.8, verbose = FALSE)
### runing UMAP
skin_data_seurat <- RunUMAP(skin_data_seurat, dims = 1:30)
### plotting clustering results
DimPlot(skin_data_seurat, label = TRUE, 
        group.by = "seurat_clusters"
        ) + ggtitle("")



#### generating metacells for som analysis ####

### cluster labels celltype X Louvain X patient ###
labels <- paste0( gsub("_", ".", skin_data_seurat$full_clustering), "_c", skin_data_seurat@meta.data$seurat_clusters, "_", skin_data_seurat@meta.data$ind_id )
names(labels) <- colnames(skin_data_seurat)
table(labels)

### remove meta-cells with less than 5 cells ###
tab <- table(labels)
par(mar=c(11,4,1,1))
barplot( sort( tab[ grep("^F|Peric|Schwann",names(tab)) ] ), las=2, log="y" )
labels <- labels[ which( labels %in% names( which( table(labels)>=5 ) ) ) ]

### subclustering – determine number of subclusters ###
labels.clusterNo <- ceiling( sort( table(labels) ) / 20 )
par(mar=c(13,4,1,1))
barplot(labels.clusterNo,las=2)
labels.clusterNo.interest <- labels.clusterNo[ grep("^F|Peric|Schwann",names(labels.clusterNo)) ]
labels.clusterNo[which(labels.clusterNo>2)] <- 2
labels.clusterNo[ grep("^F|Peric|Schwann",names(labels.clusterNo)) ] <- labels.clusterNo.interest
par(mar=c(13,4,1,1))
barplot(labels.clusterNo,las=2)
o <- order( sapply(strsplit(names(labels.clusterNo), "_" ), "[", 1 ) )
labels.clusterNo <- labels.clusterNo[o]
o <- order( sapply( strsplit(names(labels.clusterNo), "_" ), function(x) as.numeric(substr(x[2],2,nchar(x[2]))) ) )
labels.clusterNo <- labels.clusterNo[o]

### subclustering - do it! - part 1 ###
metacell.labels <- rep(NA,ncol(skin_data_seurat)) 
names(metacell.labels) <- colnames(skin_data_seurat)
metacell.labels[names(labels)] <- labels
skin_data_seurat[["metacellLabels"]] <- metacell.labels
skin_data_seurat[["cellInMetacell"]] <- !is.na(metacell.labels)
metacell.data <- matrix(NA,nrow(skin_data_seurat),0,dimnames=list(rownames(skin_data_seurat),c()))


### subclustering - do it! – part 2 ###
pb <- txtProgressBar(min = 0, max = length(labels.clusterNo),style=3)

for( x in names(labels.clusterNo) ) {
  mc.cells <- names(which(metacell.labels==x))
  expr <- data.matrix( skin_data_seurat@assays$RNA@data[ , mc.cells ] )
  km <- kmeans(t(expr),centers=labels.clusterNo[x])
  lab <- paste0( x, "_k", seq(max(km$cluster)) )
  expr <- t(km$centers)
  colnames(expr) <- lab
  metacell.data <- cbind(metacell.data, expr)
  skin_data_seurat$metacellLabels[mc.cells] <- paste0( skin_data_seurat$metacellLabels[mc.cells], " x", km$cluster )
  setTxtProgressBar( pb, pb$getVal()+1 )
}
pb$kill()


saveRDS(list(seurat, metacell.data), file = "rds/metacell_seurat.rds")
