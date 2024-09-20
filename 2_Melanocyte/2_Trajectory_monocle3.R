#Melanocyte
library(RColorBrewer)
library(SeuratWrappers)
library(Seurat)
library(monocle3)
library(SingleCellExperiment)
library(VGAM)

Melanocyte <- readRDS("/2_Melanocyte/Melanocyte_sub_merge.rds")
Idents(object = Melanocyte) <- "RNA_snn_res.0.3"

downsample_ident <- "RNA_snn_res.0.3"
gene_annotation <- as.data.frame(rownames(Melanocyte@reductions[["pca"]]@feature.loadings), row.names = rownames(Melanocyte@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"
expression_matrix <- Melanocyte@assays$RNA@counts
expression_matrix <- expression_matrix[rownames(Melanocyte@reductions[["pca"]]@feature.loadings), ]
Melanocyte_cds <- new_cell_data_set(expression_data=expression_matrix, cell_metadata = Melanocyte@meta.data, gene_metadata = gene_annotation)
set.seed(42)
recreate.partition <- c(rep(1, length(Melanocyte_cds@colData@rownames)))
names(recreate.partition) <- Melanocyte_cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
Melanocyte_cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

list_cluster <- Melanocyte@meta.data[[downsample_ident]]
names(list_cluster) <- Melanocyte@assays[["RNA"]]@data@Dimnames[[2]]
Melanocyte_cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
### Could be a space-holder, but essentially fills out louvain parameters
Melanocyte_cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

reducedDims(Melanocyte_cds)[["UMAP"]] <- Melanocyte@reductions[["umap"]]@cell.embeddings
Melanocyte_cds@preprocess_aux$gene_loadings <- Melanocyte@reductions[["pca"]]@feature.loadings[,1:20]
# Learn graph, this step usually takes a significant period of time for larger samples
Melanocyte_cds <- learn_graph(Melanocyte_cds, use_partition = TRUE)
Melanocyte_cds <- cluster_cells(cds = Melanocyte_cds, reduction_method = "UMAP")


pdf("/2_Melanocyte/2_Trajectory/monocle3/monocle3.with_label.pdf")
print(plot_cells(Melanocyte_cds, reduction_method = "UMAP", label_groups_by_cluster=FALSE,  color_cells_by = "RNA_snn_res.0.3", label_principal_points = TRUE))
dev.off()
pdf(paste0(outfile_prefix, ".pdf"))
print(plot_cells(Melanocyte_cds, reduction_method = "UMAP", label_groups_by_cluster=FALSE,  color_cells_by = "RNA_snn_res.0.3", label_principal_points = FALSE))
dev.off()

saveRDS(Melanocyte_cds, file = "/2_Melanocyte/2_Trajectory/monocle3/Melanocyte_cds.rds")

Melanocyte_cds <- monocle3::order_cells(Melanocyte_cds)
pdf("/2_Melanocyte/2_Trajectory/monocle3/pseudotime.pdf", height=6, width=7)
plot_cells(Melanocyte_cds,
   color_cells_by = "pseudotime",
   label_cell_groups=FALSE,
   label_leaves=FALSE,
   label_branch_points=FALSE,
   cell_size=0.8,
   alpha=0.8,
   graph_label_size=1.5,
   trajectory_graph_color="#d1afe2")
dev.off()
saveRDS(Melanocyte_cds, file = "/2_Melanocyte/2_Trajectory/monocle3/Melanocyte_cds_ordercells.rds")
