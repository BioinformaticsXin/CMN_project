library(Seurat)
library(SeuratDisk)
library(RColorBrewer)
library(loomR)

#################scVelo analysis###########################################################
out_dir <- "/4_Myeloid/3_scVelo"
loom_combined <- readRDS("/2_Melanocyte/2_velocyto/loom_combined.rds")
Myeloid_cell <- readRDS("/4_Myeloid/Myeloid_cell_annotation.rds")
inter_cell <- intersect(colnames(Myeloid_cell), colnames(loom_combined))
loom_combined <- subset(loom_combined, cells = inter_cell)
Myeloid_cell <- subset(Myeloid_cell, cells = inter_cell)

#add barcode
Myeloid_cell.loom <- subset(loom_combined, cells = colnames(Myeloid_cell))
barcode <- data.frame(barcode=rownames(Myeloid_cell.loom@meta.data))
rownames(barcode) <- rownames(Myeloid_cell.loom@meta.data)
Myeloid_cell.loom <- AddMetaData(Myeloid_cell.loom, metadata = barcode)
DefaultAssay(Myeloid_cell.loom) <- "ambiguous"

setwd("/4_Myeloid/3_scVelo")
SaveH5Seurat(Myeloid_cell.loom, filename = "Myeloid_cell.loom.h5Seurat", overwrite = TRUE)
Convert("Myeloid_cell.loom.h5Seurat", dest = "h5ad", overwrite = TRUE)

#######################################################################################
################################Myeloid_cell########################################
setwd("/4_Myeloid/3_scVelo")
# save metadata table:
Myeloid_cell$barcode <- colnames(Myeloid_cell)
Myeloid_cell$UMAP_1 <- Myeloid_cell@reductions$umap@cell.embeddings[,1]
Myeloid_cell$UMAP_2 <- Myeloid_cell@reductions$umap@cell.embeddings[,2]
write.csv(Myeloid_cell@meta.data, file='Myeloid_cell_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
Myeloid_cell_counts_matrix <- GetAssayData(Myeloid_cell, assay='RNA', slot='counts')
writeMM(Myeloid_cell_counts_matrix, file='Myeloid_cell_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(Myeloid_cell@reductions$pca@cell.embeddings, file='Myeloid_cell_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(Myeloid_cell_counts_matrix)),file='Myeloid_cell_gene_names.csv',
  quote=F,row.names=F,col.names=F
)
