#Melanocyte
library(Seurat)
library(loomR)
library(SeuratDisk)
library(RColorBrewer)

#################scVelo###########################################################
out_dir <- "/2_Melanocyte/2_scVelo"
loom_combined <- readRDS("/2_Melanocyte/2_velocyto/loom_combined.rds")
Melanocyte <- readRDS("/2_Melanocyte/Melanocyte_sub_merge.rds")
inter_cell <- intersect(colnames(Melanocyte), colnames(loom_combined))
loom_combined <- subset(loom_combined, cells = inter_cell)
Melanocyte <- subset(Melanocyte, cells = inter_cell)

#all Melanocyte
Melanocyte.loom <- subset(loom_combined, cells = colnames(Melanocyte))
barcode <- data.frame(barcode=rownames(Melanocyte.loom@meta.data))
rownames(barcode) <- rownames(Melanocyte.loom@meta.data)
Melanocyte.loom <- AddMetaData(Melanocyte.loom, metadata = barcode)
DefaultAssay(Melanocyte.loom) <- "ambiguous"

setwd("/2_Melanocyte/2_scVelo/")
SaveH5Seurat(Melanocyte.loom, filename = "Melanocyte_loom.h5Seurat", overwrite = TRUE)
Convert("Melanocyte_loom.h5Seurat", dest = "h5ad", overwrite = TRUE)

#######################################################################################
################################Melanocyte########################################
# save metadata table:
Melanocyte$barcode <- colnames(Melanocyte)
Melanocyte$UMAP_1 <- Melanocyte@reductions$umap@cell.embeddings[,1]
Melanocyte$UMAP_2 <- Melanocyte@reductions$umap@cell.embeddings[,2]
write.csv(Melanocyte@meta.data, file='Melanocyte_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
Melanocyte_counts_matrix <- GetAssayData(Melanocyte, assay='RNA', slot='counts')
writeMM(Melanocyte_counts_matrix, file='Melanocyte_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(Melanocyte@reductions$pca@cell.embeddings, file='Melanocyte_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(Melanocyte_counts_matrix)),file='Melanocyte_gene_names.csv',
  quote=F,row.names=F,col.names=F
)
