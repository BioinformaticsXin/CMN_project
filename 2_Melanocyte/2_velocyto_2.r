#Melanocyte
library(RColorBrewer)
library(velocyto.R)
library(Seurat)
library(SeuratWrappers)

D1_D5 <- readRDS("/1_Cell_annotation/5_annotation/D1_D5.rds")

# Convert each loom object to a Seurat object and rename them accordingly
loom_list <- list()
for(type in c("D5_Nevus", "D5_Malignant_Nevus", "D5_Melanoma", "D1", "D2", "D3", "D4")){
	D1_D5_sub <- subset(D1_D5, subset = SampleID == type)
	loom <- read.loom.matrices(paste0("/2_Melanocyte/2_velocyto/", type, ".loom"))
	loom <- as.Seurat(x = loom)

	loom_name <- data.frame(loom_name = D1_D5_sub@meta.data$Cell_Index)
	loom_name$loom_name <- paste0(paste0("cellranger:", gsub("-1", "", loom_name$loom_name)), "x")
	rownames(loom_name) <- rownames(D1_D5_sub@meta.data)
	D1_D5_sub <- AddMetaData(D1_D5_sub, metadata = loom_name)

	inter_cell <- intersect(D1_D5_sub@meta.data$loom_name, colnames(loom))
	inter_cell <- D1_D5_sub@meta.data[D1_D5_sub@meta.data$loom_name %in% inter_cell, ]

	loom <- subset(loom, cells = inter_cell$loom_name)
	loom <- RenameCells(loom, new.names = rownames(inter_cell), for.merge = T)

	loom_list[[type]] <- loom
}

loom_combined <- purrr::reduce(loom_list, merge)
saveRDS(loom_combined, file="/2_Melanocyte/2_velocyto/loom_combined.rds")

#################velocity###########################################################
loom_combined <- readRDS("/2_Melanocyte/2_velocyto/loom_combined.rds")
Melanocyte <- readRDS("/2_Melanocyte/Melanocyte_final.rds")

library(velocyto.R)
emat <- loom_combined$spliced
nmat <- loom_combined$unspliced
com_cells <- intersect(colnames(emat), colnames(Melanocyte))
print(paste("Totally,", length(com_cells), "common cells were remained ..."))
print(paste("seurat cells:", length(colnames(Melanocyte))))
print(paste("loom cells:", length(colnames(emat))))
emat <- emat[, com_cells]
nmat <- nmat[, com_cells]
seurat_obj_clust_sub <- subset(Melanocyte, cells = com_cells)

cluster_label <- seurat_obj_clust_sub@meta.data
cluster.label <- cluster_label[, "RNA_snn_res.0.3"]
names(cluster.label) <- rownames(cluster_label)

print("filtering genes by cluster expression ...")
emat <- filter.genes.by.cluster.expression(emat, cluster.label, min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat, cluster.label, min.max.cluster.average = 0.05)

data_emb <- seurat_obj_clust_sub@reductions$umap@cell.embeddings

seurat_obj_clust_sub.dist <- as.dist(1-armaCor(t(seurat_obj_clust_sub@reductions$umap@cell.embeddings)))

cluster_colors=c(ggsci::pal_nejm("default",alpha = 1)(8), ggsci::pal_npg("nrc",alpha = 1)(10))
cell.colors <- cluster_colors[cluster.label]
names(cell.colors) <- names(cluster.label)

print("estimating gene relative velocity ...")
rvel.cd <- gene.relative.velocity.estimates(emat, nmat, deltaT=1, kCells=25, cell.dist=seurat_obj_clust_sub.dist, fit.quantile=0.02, n.cores=48)

print("printing gene relative velocity ...")
pdf("/2_Melanocyte/2_velocyto/umap.pdf")
data_visual <- show.velocity.on.embedding.cor(emb=data_emb, vel=rvel.cd, n = 100, scale = "log", cell.colors = ac(x = cell.colors, alpha = 0.5), 
cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
do.par = FALSE, cell.border.alpha = 0.1, n.cores=48, return.details=TRUE)
dev.off()

trajectory_velocity_list <- list(emat=emat, nmat=nmat, emb=data_emb, cell.colors=cell.colors, cell.dist=seurat_obj_clust_sub.dist, rvel.cd=rvel.cd, data_visual=data_visual)
saveRDS(trajectory_velocity_list, file="/2_Melanocyte/2_velocyto/trajectory_velocity_list.rds")

