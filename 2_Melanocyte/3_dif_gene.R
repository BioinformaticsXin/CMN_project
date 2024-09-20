library(Seurat)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

#不同类型样本的差异表达
Melanocyte <- readRDS("/2_Melanocyte/Melanocyte_sub_merge.rds")
cluster1 <- Dif_cluster_enrich(seurat_object = Melanocyte, clustera = "cluster1", clusterb = NULL, group.by = NULL, out_dir = paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster1"), OrgDb = "org.Hs.eg.db", sig_gene_p = 0.05, sig_gene_fd = 1, organism = "hsa", cut_p = "p_val_adj")
cluster2 <- Dif_cluster_enrich(seurat_object = Melanocyte, clustera = "cluster2", clusterb = NULL, group.by = NULL, out_dir = paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster2"), OrgDb = "org.Hs.eg.db", sig_gene_p = 0.05, sig_gene_fd = 1, organism = "hsa", cut_p = "p_val_adj")
cluster3 <- Dif_cluster_enrich(seurat_object = Melanocyte, clustera = "cluster3", clusterb = NULL, group.by = NULL, out_dir = paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster3"), OrgDb = "org.Hs.eg.db", sig_gene_p = 0.05, sig_gene_fd = 1, organism = "hsa", cut_p = "p_val_adj")
cluster4 <- Dif_cluster_enrich(seurat_object = Melanocyte, clustera = "cluster4", clusterb = NULL, group.by = NULL, out_dir = paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster4"), OrgDb = "org.Hs.eg.db", sig_gene_p = 0.05, sig_gene_fd = 1, organism = "hsa", cut_p = "p_val_adj")

cluster2_1 <-  Dif_cluster_enrich(seurat_object = Melanocyte, clustera = "cluster2", clusterb = "cluster1", group.by = NULL, out_dir = paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster2_cluster1"), OrgDb = "org.Hs.eg.db", sig_gene_p = 0.05, sig_gene_fd = 1, organism = "hsa", cut_p = "p_val_adj")
cluster3_1 <-  Dif_cluster_enrich(seurat_object = Melanocyte, clustera = "cluster3", clusterb = "cluster1", group.by = NULL, out_dir = paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster3_cluster1"), OrgDb = "org.Hs.eg.db", sig_gene_p = 0.05, sig_gene_fd = 1, organism = "hsa", cut_p = "p_val_adj")
cluster4_1 <-  Dif_cluster_enrich(seurat_object = Melanocyte, clustera = "cluster4", clusterb = "cluster1", group.by = NULL, out_dir = paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster4_cluster1"), OrgDb = "org.Hs.eg.db", sig_gene_p = 0.05, sig_gene_fd = 1, organism = "hsa", cut_p = "p_val_adj")
cluster3_2 <-  Dif_cluster_enrich(seurat_object = Melanocyte, clustera = "cluster3", clusterb = "cluster2", group.by = NULL, out_dir = paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster3_cluster2"), OrgDb = "org.Hs.eg.db", sig_gene_p = 0.05, sig_gene_fd = 1, organism = "hsa", cut_p = "p_val_adj")
cluster4_3 <-  Dif_cluster_enrich(seurat_object = Melanocyte, clustera = "cluster4", clusterb = "cluster3", group.by = NULL, out_dir = paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster4_cluster3"), OrgDb = "org.Hs.eg.db", sig_gene_p = 0.05, sig_gene_fd = 1, organism = "hsa", cut_p = "p_val_adj")
cluster2_4 <-  Dif_cluster_enrich(seurat_object = Melanocyte, clustera = "cluster2", clusterb = "cluster4", group.by = NULL, out_dir = paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster2_cluster4"), OrgDb = "org.Hs.eg.db", sig_gene_p = 0.05, sig_gene_fd = 1, organism = "hsa", cut_p = "p_val_adj")

cluster1 <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster1/dif_gene.xlsx"))
cluster2 <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster2/dif_gene.xlsx"))
cluster3 <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster3/dif_gene.xlsx"))
cluster4 <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster4/dif_gene.xlsx"))
cluster2_1_up <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster2_cluster1/dif_gene.xlsx"), sheet = 1)
cluster2_1_down <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster2_cluster1/dif_gene.xlsx"), sheet = 2)
cluster3_1_up <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster3_cluster1/dif_gene.xlsx"), sheet = 1)
cluster3_1_down <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster3_cluster1/dif_gene.xlsx"), sheet = 2)
cluster4_1_up <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster4_cluster1/dif_gene.xlsx"), sheet = 1)
cluster4_1_down <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster4_cluster1/dif_gene.xlsx"), sheet = 2)
cluster3_2_up <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster3_cluster2/dif_gene.xlsx"), sheet = 1)
cluster3_2_down <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster3_cluster2/dif_gene.xlsx"), sheet = 2)
cluster4_3_up <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster4_cluster3/dif_gene.xlsx"), sheet = 1)
cluster4_3_down <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster4_cluster3/dif_gene.xlsx"), sheet = 2)
cluster2_4_up <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster2_cluster4/dif_gene.xlsx"), sheet = 1)
cluster2_4_down <- openxlsx::read.xlsx(paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/cluster2_cluster4/dif_gene.xlsx"), sheet = 2)


#cluster1
cluster1 <- cluster1[cluster1$SYMBOL %in% intersect(intersect(cluster2_1_down$SYMBOL, cluster3_1_down$SYMBOL), cluster4_1_down$SYMBOL),]
openxlsx::write.xlsx(cluster1, file = paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/function/cluster1.xlsx"))

#cluster2
cluster2 <- cluster2[cluster2$SYMBOL %in% intersect(intersect(cluster2_1_up$SYMBOL, cluster3_2_down$SYMBOL), cluster2_4_up$SYMBOL),]
openxlsx::write.xlsx(cluster2, file = paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/function/cluster2.xlsx"))

#cluster3
cluster3 <- cluster3[cluster3$SYMBOL %in% intersect(intersect(cluster3_1_up$SYMBOL, cluster3_2_up$SYMBOL), cluster4_3_down$SYMBOL),]
openxlsx::write.xlsx(cluster3, file = paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/function/cluster3.xlsx"))

#cluster4
cluster4 <- cluster4[cluster4$SYMBOL %in% intersect(intersect(cluster4_1_up$SYMBOL, cluster2_4_down$SYMBOL), cluster4_3_up$SYMBOL),]
openxlsx::write.xlsx(cluster4, file = paste0(project_dir, "/Res/2_Melanocyte/Dif/Cluster/function/cluster4.xlsx"))

p <- FeaturePlot(Melanocyte, reduction="umap", features = c("SPP1", "PRAME"), min.cutoff = 0, max.cutoff = 4, ncol=2, cols = c("#2888bd","#95c7dd","#ffffff","#e09786","#bb2438"))
ggsave(p, file = paste0(project_dir, "/Res/2_Melanocyte/Dif/SPP1_PRAME.pdf"), width=12, height=5)

exp_matrix <- Melanocyte[["RNA"]]@scale.data
SOD3_frame <- data.frame(cell_id = colnames(exp_matrix), SOD3 = unlist(exp_matrix["SOD3",]), AQP1 = unlist(exp_matrix["AQP1",]),
IFNG = exp_matrix["IFNG", ], PRKCD = exp_matrix["PRKCD", ], IRF1 = exp_matrix["IRF1", ],
IRF9 = exp_matrix["IRF9", ],
SOD3_group = rep("SOD3_high", ncol(exp_matrix)), IFNG_group = rep("IFNG_high", ncol(exp_matrix)),
AQP1_group = rep("AQP1_high", ncol(exp_matrix)), PRKCD_group = rep("PRKCD_high", ncol(exp_matrix)),
IRF1_group = rep("IRF1_high", ncol(exp_matrix)))

SOD3_frame[which(SOD3_frame$SOD3 < mean(SOD3_frame$SOD3)),]$SOD3_group <- "SOD3_low"
SOD3_frame[which(SOD3_frame$IFNG < mean(SOD3_frame$IFNG)),]$IFNG_group <- "IFNG_low"
SOD3_frame[which(SOD3_frame$AQP1 < mean(SOD3_frame$AQP1)),]$AQP1_group <- "AQP1_low"
SOD3_frame[which(SOD3_frame$PRKCD < mean(SOD3_frame$PRKCD)),]$PRKCD_group <- "PRKCD_low"
SOD3_frame[which(SOD3_frame$IRF1 < mean(SOD3_frame$IRF1)),]$IRF1_group <- "IRF1_low"
save(SOD3_frame, file="/Res/2_Melanocyte/Dif/SOD3_frame.Rdata")




################################### dif gene enrichment
library(Seurat)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)


Melanocyte <- readRDS(paste0(project_dir, "/Res/2_Melanocyte/Melanocyte_sub_merge.rds"))
cluster1 <- openxlsx::read.xlsx("/2_Melanocyte/Dif/Cluster/function/cluster1/all.tqghqosju/metascape_result.xlsx")[,c(1,7)] #Metascape
cluster3 <- openxlsx::read.xlsx("/2_Melanocyte/Dif/Cluster/function/cluster3/all.t9n1b0ywk/metascape_result.xlsx")[,c(1,7,8,11,12,13)]
cluster4 <- openxlsx::read.xlsx("/2_Melanocyte/Dif/Cluster/function/cluster4/all.toxyum_t0/metascape_result.xlsx")[,c(1,7,9,12,22)]

cluster1 <- cluster1[which(cluster1[,2] == 1),]
cluster3 <- cluster3[which(apply(matrix(as.numeric(unlist(cluster3[,-1])), ncol = ncol(cluster3)-1), 1, sum) > 0),]
cluster4 <- cluster4[which(apply(matrix(as.numeric(unlist(cluster4[,-1])), ncol = ncol(cluster4)-1), 1, sum) > 0),]

#不同类型样本的差异表达
library(ggplot2)
library(viridis)
p <- DotPlot(Melanocyte, features = rev(c(cluster1$MyList, cluster3$MyList))) + geom_point(aes(size=pct.exp), color = NA, shape = 21, colour="black", stroke=0.5) +
scale_colour_gradient2(low="#5bacdb", mid="#ffffff", high="#ed97b3") + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5)) + coord_flip() + theme_few()
ggsave(p, file = "/2_Melanocyte/Dif/Cluster/function/DotPlot1.pdf", width=4, height=6)

p <- DotPlot(Melanocyte, features = rev(cluster4$MyList)) + geom_point(aes(size=pct.exp), color = NA, shape = 21, colour="black", stroke=0.5) +
scale_colour_gradient2(low="#5bacdb", mid="#ffffff", high="#ed97b3") + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5)) + coord_flip() + theme_few()
ggsave(p, file = "/2_Melanocyte/Dif/Cluster/function/DotPlot2.pdf", width=4, height=7)


################################cluster1 cluster3#################################
cluster1_frame <- reshape2::melt(cluster1, id.vars=c("MyList"),variable.name="Term",value.name="edge")
cluster3_frame <- reshape2::melt(cluster3, id.vars=c("MyList"),variable.name="Term",value.name="edge")
cluster1_3 <- reshape2::dcast(rbind(cluster1_frame, cluster3_frame), MyList~Term)
rownames(cluster1_3) <- cluster1_3[,1]
cluster1_3 <- cluster1_3[,-1]
cluster1_3[is.na(cluster1_3)] <- 0
cluster1_3_matrix <- matrix(as.numeric(unlist(cluster1_3)), ncol = ncol(cluster1_3))
colnames(cluster1_3_matrix) <- colnames(cluster1_3)
rownames(cluster1_3_matrix) <- rownames(cluster1_3)
cluster1_3_matrix <- cluster1_3_matrix[c(cluster1$MyList, cluster3$MyList),]

library(pheatmap)
pdf("/2_Melanocyte/Dif/Cluster/function/Heatmap1.pdf", height = 12, width = 3)
  pheatmap(cluster1_3_matrix, color = c("#f7f9ff", "#6b8bc9"), breaks = c(-0.1, 0, 1), border_color = "#ffffff", cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()

################################cluster4#################################
rownames(cluster4) <- cluster4[,1]
cluster4 <- cluster4[,-1]
cluster4_matrix <- matrix(as.numeric(unlist(cluster4)), ncol = ncol(cluster4))
colnames(cluster4_matrix) <- colnames(cluster4)
rownames(cluster4_matrix) <- rownames(cluster4)
pdf("/2_Melanocyte/Dif/Cluster/function/Heatmap2.pdf", height = 15, width = 2.8)
  pheatmap(cluster4_matrix, color = c("#f7f9ff", "#6b8bc9"), breaks = c(-0.1, 0, 1), border_color = "#ffffff", cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()


###################################SPP1 PRAME##################################
Melanocyte_scale <- as.matrix(Melanocyte[["RNA"]]@scale.data[c("SPP1","PRAME"),])
Melanocyte_scale <- data.table::melt(Melanocyte_scale)
colnames(Melanocyte_scale) <- c("Gene","Cell","value")
cell_cluster <- data.frame(Cell = colnames(Melanocyte), Cell_cluster = Melanocyte@meta.data$RNA_snn_res.0.3)
Melanocyte_scale <- merge(Melanocyte_scale, cell_cluster, by = "Cell")
Melanocyte_scale$Cell_cluster <- factor(Melanocyte_scale$Cell_cluster, levels=paste0("cluster",1:4))

library(ggplot2)
library(ggthemes)
library(ggpubr)
p <- ggplot(Melanocyte_scale,aes(x = Cell_cluster, y = value, fill = Cell_cluster)) +
  geom_violin(position = position_dodge(0.7),alpha = 0.5,
              width = 1.4,trim = T,
              color = NA) +
  geom_boxplot(width = .2,show.legend = F,
               position = position_dodge(0.7),
               color = 'grey20',alpha = 0.5,
               outlier.color = 'grey50') +
  theme_few(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
        legend.position = 'top') +
  scale_fill_manual(values = c('cluster1'='#993a29','cluster2'='#3d4da0','cluster3'='#c08427','cluster4'='#2e7c40'), name = '')+
  facet_wrap(. ~ Gene, ncol=2, scales = "free")+
  stat_compare_means(label = "p.signif", comparisons=list(c("cluster1","cluster4"),c("cluster2","cluster4"),c("cluster3","cluster4")))
ggsave(p, file="/2_Melanocyte/Dif/Cluster/function/SPP1_PRAME.pdf", width=8, height=7)
