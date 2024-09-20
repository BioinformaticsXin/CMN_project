library(Seurat)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

#Melanocyte
project_dir <- "/HSCR/hengya_work/lixin/Project/Congenital_melanoma_nevus_single_cell"
D1_D5 <- readRDS("/1_Cell_annotation/5_annotation/D1_D5.rds")
Melanocyte <- subset(D1_D5, idents = c("Melanocyte", "Mitotic Melanocyte cell"))
Melanocyte <- ScaleData(Melanocyte, features = rownames(Melanocyte))
Melanocyte <- RunPCA(Melanocyte, features = VariableFeatures(object = Melanocyte))
Melanocyte <- RunHarmony(Melanocyte, c("SampleID","Tech"))
Melanocyte <- RunUMAP(Melanocyte, reduction = "harmony", dims = 1:30)
Melanocyte <- FindNeighbors(Melanocyte, reduction = "harmony")
Melanocyte <- FindClusters(Melanocyte, resolution = 0.3)
markers <- FindAllMarkers(Melanocyte, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_sig <- markers[which(markers$avg_log2FC > 0.25 & markers$p_val_adj < 0.01),]
openxlsx::write.xlsx(markers_sig, file= "/2_Melanocyte/MarkerGene_0.3_sig.xlsx", overwrite = TRUE)
saveRDS(Melanocyte, file="/2_Melanocyte/Melanocyte_final.rds")


Idents(Melanocyte) <- "RNA_snn_res.0.3"
p <- FeaturePlot(Melanocyte, reduction="umap",features = c("DCT","PMEL","TYRP1"),min.cutoff = 0, max.cutoff = 4,ncol=3,cols = c("#DCDCDC","red"))
ggsave(p, file = "/2_Melanocyte/featureplot_check_marker.pdf", width=18,height=6)

p1 <- DimPlot(Melanocyte, reduction = "umap", group.by = "SampleType", pt.size = .5, cols = c("#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF"), label = FALSE)
ggsave(p1, file = "/2_Melanocyte/umap_SampleType.pdf", width=7.3, height=6)
p1 <- DimPlot(Melanocyte, reduction = "umap", group.by = "SampleID", pt.size = .5, cols = c(brewer.pal(12,"Paired"), brewer.pal(3,"Set1")), label = FALSE)
ggsave(p1, file = "/2_Melanocyte/umap_SampleID.pdf", width=7, height=5.5)

p1 <- DimPlot(Melanocyte, reduction = "umap", pt.size = .5, ncol=5, cols = c("#C92523", "#5061B0", "#50B33A", "#8A239B", "#E68310", "#F0FF18", "#8E5425", "#E767BA"), label = FALSE)
ggsave(p1, file = "/2_Melanocyte/umap.pdf", width=6.5, height=6)

#Subcluster merging
Idents(Melanocyte) <- "RNA_snn_res.0.3"
cluster_to_cell <- list("1" = "cluster1", "7" = "cluster1", "0" = "cluster2", "3" = "cluster2", "5" = "cluster2", "4" = "cluster3", "2" = "cluster4", "6" = "cluster4")
new.cluster.ids <- unlist(as.character(cluster_to_cell))
names(new.cluster.ids) <- names(cluster_to_cell)
Melanocyte <- RenameIdents(Melanocyte, new.cluster.ids)
Melanocyte@meta.data$RNA_snn_res.0.3 <- Idents(Melanocyte)

p1 <- DimPlot(Melanocyte, reduction = "umap", pt.size = .5, ncol=5, cols = ggsci::pal_nejm("default",alpha = 1)(8), label = FALSE)
ggsave(p1, file = "/2_Melanocyte/umap_submerge.pdf", width=7, height=6)

saveRDS(Melanocyte, file = "/2_Melanocyte/Melanocyte_sub_merge.rds")
