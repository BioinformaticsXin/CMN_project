library(Seurat)
library(SeuratDisk)
library(RColorBrewer)
library(loomR)

#################scVelo analysis###########################################################
out_dir <- "/5_Endothelial_cell/2_scVelo"
loom_combined <- readRDS("/2_Melanocyte/2_velocyto/loom_combined.rds")
Endothelial_cell <- readRDS("/5_Endothelial_cell/endothelial_cell_annotation.rds")
Endothelial_cell <- subset(Endothelial_cell, idents=c("venous ECs", "capillary ECs", "arterial ECs"))
Endothelial_cell@meta.data$UMAP_1 <- Endothelial_cell@reductions$umap@cell.embeddings[,1]
Endothelial_cell@meta.data$UMAP_2 <- Endothelial_cell@reductions$umap@cell.embeddings[,2]

Endothelial_cell <- subset(Endothelial_cell, UMAP_1 > -10)

inter_cell <- intersect(colnames(Endothelial_cell), colnames(loom_combined))
loom_combined <- subset(loom_combined, cells = inter_cell)
Endothelial_cell <- subset(Endothelial_cell, cells = inter_cell)

#add barcode
Endothelial_cell.loom <- subset(loom_combined, cells = colnames(Endothelial_cell))
barcode <- data.frame(barcode=rownames(Endothelial_cell.loom@meta.data))
rownames(barcode) <- rownames(Endothelial_cell.loom@meta.data)
Endothelial_cell.loom <- AddMetaData(Endothelial_cell.loom, metadata = barcode)
DefaultAssay(Endothelial_cell.loom) <- "ambiguous"

setwd(paste0("/5_Endothelial_cell/2_scVelo"))
SaveH5Seurat(Endothelial_cell.loom, filename = "Endothelial_cell.loom.h5Seurat", overwrite = TRUE)
Convert("Endothelial_cell.loom.h5Seurat", dest = "h5ad", overwrite = TRUE)

#######################################################################################
################################Endothelial_cell########################################

# save metadata table:
Endothelial_cell$barcode <- colnames(Endothelial_cell)
Endothelial_cell$UMAP_1 <- Endothelial_cell@reductions$umap@cell.embeddings[,1]
Endothelial_cell$UMAP_2 <- Endothelial_cell@reductions$umap@cell.embeddings[,2]
write.csv(Endothelial_cell@meta.data, file='Endothelial_cell_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
Endothelial_cell_counts_matrix <- GetAssayData(Endothelial_cell, assay='RNA', slot='counts')
writeMM(Endothelial_cell_counts_matrix, file='Endothelial_cell_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(Endothelial_cell@reductions$pca@cell.embeddings, file='Endothelial_cell_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(Endothelial_cell_counts_matrix)),file='Endothelial_cell_gene_names.csv',
  quote=F,row.names=F,col.names=F
)
