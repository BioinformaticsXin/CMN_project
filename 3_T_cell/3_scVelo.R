library(Seurat)
library(SeuratDisk)
library(RColorBrewer)
library(loomR)

################# scVelo analysis ###########################################################
out_dir <- "/3_T_cell/3_scVelo"
loom_combined <- readRDS("/2_Melanocyte/2_velocyto/loom_combined.rds")
T_cell <- readRDS("/3_T_cell/t_cell_annotation.rds")  
inter_cell <- intersect(colnames(T_cell), colnames(loom_combined))
loom_combined <- subset(loom_combined, cells = inter_cell)
T_cell <- subset(T_cell, cells = inter_cell)

#add barcode
T_cell.loom <- subset(loom_combined, cells = colnames(T_cell))
barcode <- data.frame(barcode=rownames(T_cell.loom@meta.data))
rownames(barcode) <- rownames(T_cell.loom@meta.data)
T_cell.loom <- AddMetaData(T_cell.loom, metadata = barcode)
DefaultAssay(T_cell.loom) <- "ambiguous"

setwd("/3_T_cell/3_scVelo")
SaveH5Seurat(T_cell.loom, filename = "T_cell.loom.h5Seurat", overwrite = TRUE)
Convert("T_cell.loom.h5Seurat", dest = "h5ad", overwrite = TRUE)

#######################################################################################
################################T_cell########################################
setwd("/3_T_cell/3_scVelo")
# save metadata table:
T_cell$barcode <- colnames(T_cell)
T_cell$UMAP_1 <- T_cell@reductions$umap@cell.embeddings[,1]
T_cell$UMAP_2 <- T_cell@reductions$umap@cell.embeddings[,2]
write.csv(T_cell@meta.data, file='T_cell_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
T_cell_counts_matrix <- GetAssayData(T_cell, assay='RNA', slot='counts')
writeMM(T_cell_counts_matrix, file='T_cell_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(T_cell@reductions$pca@cell.embeddings, file='T_cell_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(T_cell_counts_matrix)),file='T_cell_gene_names.csv',
  quote=F,row.names=F,col.names=F
)
