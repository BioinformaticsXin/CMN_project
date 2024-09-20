library(Seurat)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(dplyr)

# load data
T_cell <- readRDS("/3_T_cell/t_cell_annotation.rds")  
Idents(T_cell) <- factor(x = Idents(T_cell), levels = c("Treg", "Effector CD4+ T cell", "Naive CD4+ T cell", "Effector GZMB+ CD8+ T cell", "Effector GZMK+ CD8+ T cell", "NK"))

# findmarker
markers <- FindAllMarkers(T_cell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

top_genes <- markers %>% select(gene,cluster)
top_genes$cluster <- paste0('cluster',top_genes$cluster)
unstack(top_genes) %>% openxlsx::write.xlsx(file= "/3_T_cell/1_visualization/top_genes.xlsx", overwrite = TRUE)

# visualization
cell_color <- brewer.pal(7,"Accent")
p1 <- DimPlot(object = T_cell, reduction = "umap", label = TRUE, pt.size = .5, group.by = "cell.type", raster = FALSE, cols = cell_color)
ggsave(filename="/3_T_cell/1_visualization/umap.pdf",plot=p1,width=8,height=6)
p1 <- DimPlot(object = T_cell, reduction = "umap", label = FALSE, pt.size = .5, group.by = "SampleType", raster = FALSE, cols = brewer.pal(5,"Paired"))
ggsave(filename="/3_T_cell/1_visualization/umap_SampleType.pdf",plot=p1,width=7.3,height=6)
p1 <- DimPlot(object = T_cell, reduction = "umap", label = FALSE, pt.size = .5, group.by = "SampleID", raster = FALSE, cols = c(brewer.pal(12,"Paired"), brewer.pal(3,"Set1")))
ggsave(filename="/3_T_cell/1_visualization/umap_SampleID.pdf",plot=p1,width=7.5,height=6)
marker_gene <- c("CD3E", "CD3D", "CD4", "CD8A", "CD8B", "CCR7", "SELL", "FOXP3", "TIGIT", "CTLA4", "GZMK", "GZMB", "NKG7", "GNLY", "IL7R", "PRDM1")
p <- FeaturePlot(T_cell, pt.size = 1, features = marker_gene, cols = c("#f8f9df", "#f2814b", "#741f2b"), ncol = 4)
ggsave(p, file="/3_T_cell/1_visualization/feature.pdf", width=23, height=18)

# the proportion and number of cells in each SampleID
# the proportion
T_cell@meta.data$cell.type <- Idents(T_cell)
samples <- unique(T_cell@meta.data$SampleID)
Cell_sta <- c()
for(i in 1:length(samples)){
   this.meta <- T_cell@meta.data[which(T_cell@meta.data$SampleID == samples[i]),]
   this.sta <- data.frame(table(this.meta$cell.type), as.numeric(table(this.meta$cell.type)/sum(table(this.meta$cell.type))), SampleID = rep(samples[i], length(table(this.meta$cell.type))), stringsAsFactors=FALSE)
   Cell_sta <- rbind(Cell_sta, this.sta)
}
colnames(Cell_sta) <- c("Cell", "Number", "Proportion", "SampleID")

Cell_sta[,1] <- as.character(Cell_sta[,1])
Cell_sta$Cell <- factor(Cell_sta$Cell, levels=c("Treg", "Naive CD4+ T cell", "Effector CD4+ T cell", "Effector GZMB+ CD8+ T cell", "Effector GZMK+ CD8+ T cell", "NK"))

out_dir <- "/3_T_cell/1_visualization/"
Cell_sta$SampleID <- factor(Cell_sta$SampleID, levels=rev(samples))
p1 <- ggplot(Cell_sta, aes(x=SampleID, y=Proportion,fill=Cell))+
geom_bar(stat="identity")+scale_fill_manual(values=c("Treg"=cell_color[1], "Naive CD4+ T cell"=cell_color[2], 
"Effector CD4+ T cell"=cell_color[3], "Effector GZMB+ CD8+ T cell"=cell_color[4], 
"Effector GZMK+ CD8+ T cell"=cell_color[5], "NK"=cell_color[6]))+theme_classic()+coord_flip();
ggsave(p1, file = paste0(out_dir, "/bar_cell_pro.pdf"), width=8,height=7)

# the absolute number of cells in each sample
Cell_sta <- data.frame(SampleID = names(table(T_cell@meta.data$SampleID)), Cellnum = as.numeric(table(T_cell@meta.data$SampleID)), stringsAsFactors=FALSE)
Cell_sta$SampleID <- factor(Cell_sta$SampleID, levels=rev(c("D1_Melanoma","D1_Malignant_Nevus","D1_Nevus",
"D2_Melanoma","D2_Malignant_Nevus","D2_Nevus","D3_Melanoma","D3_Malignant_Nevus","D4_Melanoma","D4_Normal","D5_Melanoma",
"D5_Malignant_Nevus","D5_Nevus")))

p1 <- ggplot(Cell_sta, aes(x=SampleID, y=Cellnum))+geom_bar(stat="identity", color="grey")+theme_classic()+coord_flip();
ggsave(p1, file = paste0(out_dir, "/bar_cell_num.pdf"), width=6,height=7)

# the proportion of different cell types in each SampleType
T_cell@meta.data$cell.type <- Idents(T_cell)
cells <- as.character(unique(T_cell@meta.data$cell.type))
Cell_sta <- c()
for(i in 1:length(cells)){
   this.meta <- T_cell@meta.data[which(T_cell@meta.data$cell.type == cells[i]),]
   this.sta <- data.frame(table(this.meta$SampleType), as.numeric(table(this.meta$SampleType)/sum(table(this.meta$SampleType))), Cell = rep(cells[i], length(table(this.meta$SampleType))), stringsAsFactors=FALSE)
   Cell_sta <- rbind(Cell_sta, this.sta)
}
colnames(Cell_sta) <- c("SampleType", "Number", "Proportion", "Cell")
p1 <- ggplot(Cell_sta, aes(x=Cell, y=Proportion,fill=SampleType))+
geom_bar(stat="identity")+scale_fill_manual(values=brewer.pal(5,"Paired"))+theme_classic()+coord_flip();
ggsave(p1, file = paste0(out_dir, "/bar_sampletype_pro.pdf"), width=8,height=5)


# the absolute number of cells in each cell.type
Cell_sta <- data.frame(Cell = names(table(T_cell@meta.data$cell.type)), Cellnum = as.numeric(table(T_cell@meta.data$cell.type)), stringsAsFactors=FALSE)

p1 <- ggplot(Cell_sta, aes(x=Cell, y=Cellnum))+geom_bar(stat="identity", color="grey")+theme_classic()+coord_flip();
ggsave(p1, file = paste0(out_dir, "/bar_sampletype_num.pdf"), width=6,height=5)
