library(Seurat)
library(ggplot2)
library(ggthemes)

#functional annotation plot
C0 <- openxlsx::read.xlsx("/5_Endothelial_cell/4_capillary_EC_function/cluster0/metascape_result.xlsx", sheet = 2)
C1 <- openxlsx::read.xlsx("/5_Endothelial_cell/4_capillary_EC_function/cluster1/metascape_result.xlsx", sheet = 2)
C2 <- openxlsx::read.xlsx("/5_Endothelial_cell/4_capillary_EC_function/cluster2/metascape_result.xlsx", sheet = 2)
C3 <- openxlsx::read.xlsx("/5_Endothelial_cell/4_capillary_EC_function/cluster3/metascape_result.xlsx", sheet = 2)
C4 <- openxlsx::read.xlsx("/5_Endothelial_cell/4_capillary_EC_function/cluster4/metascape_result.xlsx", sheet = 2)

C0 <- C0[grepl("_Summary", C0[,1]),]
C1 <- C1[grepl("_Summary", C1[,1]),]
C2 <- C2[grepl("_Summary", C2[,1]),]
C3 <- C3[grepl("_Summary", C3[,1]),]
C4 <- C4[grepl("_Summary", C4[,1]),]

colnames(C0)[5] <- "log_p"
colnames(C1)[5] <- "log_p"
colnames(C2)[5] <- "log_p"
colnames(C3)[5] <- "log_p"
colnames(C4)[5] <- "log_p"

C0$Description <- factor(C0$Description, levels = rev(C0$Description))
p <- ggplot(C0[1:6,], aes(x = Description, y = -log_p, fill = -log_p)) + geom_bar(stat="identity", width = 0.8) + coord_flip() +
scale_fill_gradient(low = '#b1dde5', high = '#4cbbd2') + theme_classic()
ggsave(p, file = "/5_Endothelial_cell/4_capillary_EC_function/function_C0.pdf", width = 5, height = 2)


C1$Description <- factor(C1$Description, levels = rev(C1$Description))
p <- ggplot(C1[1:6,], aes(x = Description, y = -log_p, fill = -log_p)) + geom_bar(stat="identity", width = 0.8) + coord_flip() +
scale_fill_gradient(low = "#ecc8cc", high = "#e9949e") + theme_classic()
ggsave(p, file = "/5_Endothelial_cell/4_capillary_EC_function/function_C1.pdf", width = 4.5, height = 2)


C2$Description <- factor(C2$Description, levels = rev(C2$Description))
p <- ggplot(C2[1:6,], aes(x = Description, y = -log_p, fill = -log_p)) + geom_bar(stat="identity", width = 0.8) + coord_flip() +
scale_fill_gradient(low = "#cde0eb", high = "#9fc6db") + theme_classic()
ggsave(p, file = "/5_Endothelial_cell/4_capillary_EC_function/function_C2.pdf", width = 4.5, height = 2)


C3$Description <- factor(C3$Description, levels = rev(C3$Description))
p <- ggplot(C3[1:6,], aes(x = Description, y = -log_p, fill = -log_p)) + geom_bar(stat="identity", width = 0.8) + coord_flip() +
scale_fill_gradient(low = "#f5e6ec", high = "#efb4cc") + theme_classic()
ggsave(p, file = "/5_Endothelial_cell/4_capillary_EC_function/function_C3.pdf", width = 4.5, height = 2)


C4$Description <- factor(C4$Description, levels = rev(C4$Description))
p <- ggplot(C4[1:6,], aes(x = Description, y = -log_p, fill = -log_p)) + geom_bar(stat="identity", width = 0.8) + coord_flip() +
scale_fill_gradient(low = "#e2e0e0", high = "#c4aca6") + theme_classic()
ggsave(p, file = "/5_Endothelial_cell/4_capillary_EC_function/function_C4.pdf", width = 4.2, height = 2)


#plot heatmap
library(ComplexHeatmap)
library(circlize)

capillary_EC <- readRDS("/5_Endothelial_cell/4_capillary_EC_function/capillary_EC.rds")

C0 <- openxlsx::read.xlsx("/5_Endothelial_cell/4_capillary_EC_function/cluster0/metascape_result.xlsx", sheet = 2)
C1 <- openxlsx::read.xlsx("/5_Endothelial_cell/4_capillary_EC_function/cluster1/metascape_result.xlsx", sheet = 2)
C2 <- openxlsx::read.xlsx("/5_Endothelial_cell/4_capillary_EC_function/cluster2/metascape_result.xlsx", sheet = 2)
C3 <- openxlsx::read.xlsx("/5_Endothelial_cell/4_capillary_EC_function/cluster3/metascape_result.xlsx", sheet = 2)
C4 <- openxlsx::read.xlsx("/5_Endothelial_cell/4_capillary_EC_function/cluster4/metascape_result.xlsx", sheet = 2)

C0 <- C0[grepl("_Summary", C0[,1]),]
C1 <- C1[grepl("_Summary", C1[,1]),]
C2 <- C2[grepl("_Summary", C2[,1]),]
C3 <- C3[grepl("_Summary", C3[,1]),]
C4 <- C4[grepl("_Summary", C4[,1]),]

colnames(C0)[5] <- "log_p"
colnames(C1)[5] <- "log_p"
colnames(C2)[5] <- "log_p"
colnames(C3)[5] <- "log_p"
colnames(C4)[5] <- "log_p"

C0_gene <- unlist(strsplit(C0[C0$Description %in% "angiogenesis", "Symbols"], split = ","))
C1_gene <- unlist(strsplit(C1[C1$Description %in% "VEGFA-VEGFR2 signaling", "Symbols"], split = ","))
C2_gene <- unlist(strsplit(C2[C2$Description %in% "blood vessel morphogenesis", "Symbols"], split = ","))
C3_gene <- unlist(strsplit(C3[C3$Description %in% "Platelet activation, signaling and aggregation", "Symbols"], split = ","))
C4_gene <- unlist(strsplit(C4[C4$Description %in% "Cytokine Signaling in Immune system", "Symbols"], split = ","))

plot_heatmap <- capillary_EC[["RNA"]]@scale.data
col_anno <- data.frame(ID = rownames(capillary_EC@meta.data), Cluster = paste0("C", capillary_EC@meta.data[,c("RNA_snn_res.0.8")]))
rownames(col_anno) <- rownames(capillary_EC@meta.data)
col_anno$Cluster <- factor(col_anno$Cluster, levels = paste0("C",0:4))

col_anno <- col_anno[order(col_anno$Cluster),]
plot_heatmap <- plot_heatmap[,rownames(col_anno)]

color <- c("#4cbbd2", "#e9949e", "#9fc6db", "#efb4cc", "#c4aca6")
names(color) <- paste0("C",0:4)

row_anno <- data.frame(gene = c(C0_gene, C1_gene, C2_gene, C3_gene, C4_gene), cluster = c(rep("C0", length(C0_gene)), rep("C1", length(C1_gene)), rep("C2", length(C2_gene)), rep("C3", length(C3_gene)), rep("C4", length(C4_gene))))
row_anno$cluster <- factor(row_anno$cluster, levels = c(paste0("C",0:4)))

row_anno <- row_anno[!duplicated(row_anno$gene),]
rownames(row_anno) <- row_anno$gene

plot_heatmap[which(plot_heatmap <= -2)] <- -2
plot_heatmap[which(plot_heatmap >= 2)] <- 2

col_fun = colorRamp2(c(min(plot_heatmap[row_anno$gene,]), 0, max(plot_heatmap[row_anno$gene,])), c("#499414", "#f6f6f5","#b00577"))

pdf("/5_Endothelial_cell/4_capillary_EC_function/capillary_EC_heatmap.pdf")
Heatmap(plot_heatmap[row_anno$gene,], name = "mat",
   top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color), labels = levels(col_anno$Cluster), 
   labels_gp = gpar(col = "white", fontsize = 10))), 
   left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = color), labels = levels(row_anno$cluster), 
   labels_gp = gpar(col = "white", fontsize = 10))),  col = col_fun,
   cluster_columns = FALSE, cluster_rows = FALSE, column_split = col_anno$Cluster, row_split = row_anno$cluster,
   show_column_names = FALSE, show_row_names = TRUE, use_raster = TRUE)
dev.off()
