library(Seurat)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(doMC)
library(openxlsx)
library(RColorBrewer)

D1_D5 <- readRDS("/1_Cell_annotation/4_harmony/D1_D5.rds")
D1_D5 <- FindNeighbors(D1_D5, reduction = "harmony")
D1_D5 <- FindClusters(D1_D5, resolution = 0.12)
markers <- FindAllMarkers(D1_D5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_sig <- markers[which(markers$avg_log2FC > 0.25 & markers$p_val_adj < 0.01),]
openxlsx::write.xlsx(markers_sig, file= "/1_Cell_annotation/5_annotation/MarkerGene_0.12_sig.xlsx", overwrite = TRUE)


cluster_to_cell <- list("0" = "Melanocyte", "1" = "T cell / NK", "2" = "Endothelial cell", "3" = "iCAF", "4" = "Keratinocyte", "5" = "myo-CAF", "6" = "Myeloid cell", "7" = "Mast cell", "8" = "Mitotic Melanocyte cell")
cell_color <- c(ggsci::pal_npg("nrc",alpha = 1)(10), ggsci::pal_nejm("default",alpha = 1)(8))

marker_gene <- c("PMEL","MLANA","TYRP1","DCT","CD3D","CD8A","CD3E","CD3G","NKG7","VWF","CDH5","DCN","COL1A1","COL6A2","PDGFRA","PDGFRB","ACTA2","RGS5","KRT1","KRT5","LYZ","CD14","CD163","CPA3","IL1RL1","TOP2A","MKI67")
marker_gene <- unique(unlist(marker_gene))
new.cluster.ids <- unlist(as.character(cluster_to_cell))
names(new.cluster.ids) <- names(cluster_to_cell)
D1_D5 <- RenameIdents(D1_D5, new.cluster.ids)
D1_D5@meta.data$cell.type <- Idents(D1_D5)

p1 <- DimPlot(D1_D5, reduction = "umap", label = TRUE, pt.size = .1, cols = cell_color) + coord_fixed()
ggsave(p1, file = "/1_Cell_annotation/5_annotation/umap.pdf", width=8, height=6)
p <- DotPlot(D1_D5, features = marker_gene) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + scale_color_viridis(option="magma") + 
     guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(p1, file = "/1_Cell_annotation/5_annotation/DotPlot.pdf", width=11, height=4)

p1 <- DimPlot(object = D1_D5, reduction = "umap", label = TRUE, pt.size = .1, group.by = "cell.type", raster = FALSE, cols = cell_color)
ggsave(filename="/1_Cell_annotation/5_annotation/umap_CellType.pdf",plot=p1,width=7.8,height=6)
p1 <- DimPlot(object = D1_D5, reduction = "umap", label = FALSE, pt.size = .1, group.by = "SampleType", raster = FALSE, cols = brewer.pal(5,"Paired"))
ggsave(filename="/1_Cell_annotation/5_annotation/umap_SampleType.pdf",plot=p1,width=7,height=6)
p1 <- DimPlot(object = D1_D5, reduction = "umap", label = FALSE, pt.size = .1, group.by = "SampleID", raster = FALSE, cols = c(brewer.pal(12,"Paired"), brewer.pal(3,"Set1")))
ggsave(filename="/1_Cell_annotation/5_annotation/umap_SampleID.pdf",plot=p1,width=7.8,height=6)
saveRDS(D1_D5, file = "/1_Cell_annotation/5_annotation/D1_D5.rds")


