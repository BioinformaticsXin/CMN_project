library(Seurat)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

# Calculate the relationship between genes and pseudotime using linear regression
Melanocyte <- readRDS("/2_Melanocyte/Melanocyte_sub_merge.rds")
Idents(object = Melanocyte) <- "RNA_snn_res.0.3"
pseudotime <- read.csv("/2_Melanocyte/2_Trajectory/monocle3/pseudotime.csv")
cluster4_1 <- openxlsx::read.xlsx("/2_Melanocyte/Dif/Cluster/cluster4_cluster1/dif_gene.xlsx")
dif_gene <- cluster4_1$SYMBOL

Melanocyte_explore <- Melanocyte[["RNA"]]@scale.data[dif_gene,]
linear_regression <- apply(Melanocyte_explore, 1, function(x){
   x <- as.numeric(unlist(x))
   lr <- lm(pseudotime$pseudotime ~ x)
   return(data.frame(p = as.numeric(summary(lr)$coefficients[,4][2]), r_squared = summary(lr)$r.squared))
})

linear_regression <- unlist(linear_regression)
linear_regression <- data.frame(gene = rownames(Melanocyte_explore), p = linear_regression[(1:nrow(Melanocyte_explore))*2-1], r_squared = linear_regression[(1:nrow(Melanocyte_explore))*2])
linear_regression <- data.frame(fdr = p.adjust(linear_regression$p), linear_regression)
saveRDS(linear_regression, file = "/2_Melanocyte/3_lineage_gene/linear_regression.rds")

###################################permutation############################################################
library(foreach)
library(doParallel)
registerDoParallel(cores = 10)

permutation_p <- foreach(i=1:1000, .combine="cbind") %dopar% {
   apply(Melanocyte_explore, 1, function(x){
   x <- as.numeric(unlist(x))
   this.pseudotime <- sample(pseudotime$pseudotime, size = length(pseudotime$pseudotime), replace = FALSE)
   lr <- lm(this.pseudotime ~ x)
   return(as.numeric(summary(lr)$coefficients[,4][2]))
  })
}
saveRDS(permutation_p, file = "/2_Melanocyte/3_lineage_gene/permutation_p.rds")


################### Select lineage genes based on permutation ##############################
obs_p <- readRDS("/2_Melanocyte/3_lineage_gene/linear_regression.rds")

expec_p <- apply(cbind(obs_p$p, permutation_p), 1, function(x){
   this.p <- (length(which(x[-1] < x[1]))+1)/1001
   return(this.p)
})

obs_p <- data.frame(obs_p, expec_p = expec_p)
rownames(obs_p) <- obs_p$gene
openxlsx::write.xlsx(obs_p, file = "/2_Melanocyte/3_lineage_gene/obs_p.xlsx", overwrite = TRUE)


###############################################################################
# Based on these final selected genes, plot a heatmap to visualize the differential expression fold change
###############################################################################
library(ComplexHeatmap)
library(circlize)
enrich_result <- openxlsx::read.xlsx("/2_Melanocyte/3_lineage_gene/all.t3rarsbj_/metascape_result.xlsx")
enrich_result <- enrich_result[,c(1, 7, 9, 11, 12, 13, 22, 23, 26)]
enrich_result <- reshape2::melt(enrich_result, id.vars=c("MyList"),variable.name="Term",value.name="edge")
enrich_result <- enrich_result[which(enrich_result$edge == "1.0"),]
enrich_result <- unique(enrich_result[,1:2],)

# Sort expression values by pseudotime
Melanocyte_scale <- Melanocyte[["RNA"]]@scale.data
pseudotime <- read.csv("/2_Melanocyte/2_Trajectory/monocle3/pseudotime.csv")
pseudotime <- pseudotime[order(pseudotime$pseudotime),]
Melanocyte_scale <- Melanocyte_scale[,pseudotime$X]
Melanocyte_scale <- Melanocyte_scale[rownames(Melanocyte_scale) %in% enrich_result$MyList,]

heatmap_exp <- Melanocyte_scale
heatmap_exp[which(heatmap_exp < -4)] <- -4
heatmap_exp[which(heatmap_exp > 4)] <- 4

obs_p <- openxlsx::read.xlsx("/2_Melanocyte/3_lineage_gene/obs_p.xlsx")
cluster4_cluster1 <- openxlsx::read.xlsx("/2_Melanocyte/Dif/Cluster/cluster4_cluster1/dif_gene.xlsx")
rownames(cluster4_cluster1) <- cluster4_cluster1$SYMBOL
colnames(cluster4_cluster1)[1] <- "gene"
right_anno <- merge(obs_p, cluster4_cluster1, by = "gene")
rownames(right_anno) <- right_anno$gene
num <- table(enrich_result$Term)

set.seed(123)
col_1 = colorRamp2(seq(min(right_anno$r_squared), max(right_anno$r_squared), length.out = 50), colorRampPalette(c("#ffffff", "#1f389b"))(50))
col_2 = colorRamp2(seq(min(right_anno$avg_log2FC), max(right_anno$avg_log2FC), length.out = 50), colorRampPalette(c("#ffffff", "#0b9387"))(50))
col_pseu = colorRamp2(seq(min(pseudotime$pseudotime), max(pseudotime$pseudotime), length.out = 50), colorRampPalette(c("#ffffff", "#eb39b3"))(50))

this.right_anno <- right_anno[enrich_result[which(enrich_result$Term == "R-HSA-163210.Formation.of.ATP.by.chemiosmotic.coupling"),]$MyList,]
ht1 <- Heatmap(heatmap_exp[enrich_result[which(enrich_result$Term == "R-HSA-163210.Formation.of.ATP.by.chemiosmotic.coupling"),]$MyList,], cluster_columns = FALSE,
      use_raster = TRUE, show_column_names = FALSE, show_row_names = FALSE, col = colorRampPalette(c("#3161a9", "#34b7f8", "#d6d3b3", "#e94527", "#a31b1b"))(100),
      top_annotation = HeatmapAnnotation(pseudotime = pseudotime$pseudotime, cell_cluster = pseudotime$RNA_snn_res.0.3,
      col = list(pseudotime = col_pseu, cell_cluster = c("cluster1"="#db6657","cluster2"="#977dac","cluster3"="#b1ccea","cluster4"="#4876b4"))),
      left_annotation = rowAnnotation(Term = rep("Formation of ATP by chemiosmotic coupling", num[1]), col = list(Term = c("Formation of ATP by chemiosmotic coupling"="#d1342b"))),
      right_annotation = rowAnnotation(r_squared = this.right_anno[enrich_result[which(enrich_result$Term == "R-HSA-163210.Formation.of.ATP.by.chemiosmotic.coupling"),]$MyList,]$r_squared, 
      avg_log2FC = this.right_anno[enrich_result[which(enrich_result$Term == "R-HSA-163210.Formation.of.ATP.by.chemiosmotic.coupling"),]$MyList,]$avg_log2FC,
      link = anno_mark(at = which(this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)), 
      labels = rownames(this.right_anno)[this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)], 
      labels_gp = gpar(fontsize = 6), padding = unit(1, "mm")), 
      col = list(r_squared = col_1, avg_log2FC = col_2), show_annotation_name = FALSE), show_heatmap_legend = TRUE,
      heatmap_legend_param = list(title = "value"))

this.right_anno <- right_anno[enrich_result[which(enrich_result$Term == "R-HSA-9006934.Signaling.by.Receptor.Tyrosine"),]$MyList,]
ht2 <- Heatmap(heatmap_exp[enrich_result[which(enrich_result$Term == "R-HSA-9006934.Signaling.by.Receptor.Tyrosine"),]$MyList,], cluster_columns = FALSE,
      use_raster = TRUE, show_column_names = FALSE, show_row_names = FALSE, col = colorRampPalette(c("#3161a9", "#34b7f8", "#d6d3b3", "#e94527", "#a31b1b"))(100),
      left_annotation = rowAnnotation(Term = rep("Signaling by Receptor Tyrosine", num[2]), col = list(Term = c("Signaling by Receptor Tyrosine"="#67ac56"))),
      right_annotation = rowAnnotation(r_squared = this.right_anno[enrich_result[which(enrich_result$Term == "R-HSA-9006934.Signaling.by.Receptor.Tyrosine"),]$MyList,]$r_squared, 
      avg_log2FC = this.right_anno[enrich_result[which(enrich_result$Term == "R-HSA-9006934.Signaling.by.Receptor.Tyrosine"),]$MyList,]$avg_log2FC,
      link = anno_mark(at = which(this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)), 
      labels = rownames(this.right_anno)[this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)], 
      labels_gp = gpar(fontsize = 6), padding = unit(1, "mm")), 
      col = list(r_squared = col_1, avg_log2FC = col_2), show_annotation_name = FALSE), show_heatmap_legend = FALSE)

this.right_anno <- right_anno[enrich_result[which(enrich_result$Term == "R-HSA-194315.Signaling.by.Rho.GTPases"),]$MyList,]
ht3 <- Heatmap(heatmap_exp[enrich_result[which(enrich_result$Term == "R-HSA-194315.Signaling.by.Rho.GTPases"),]$MyList,], cluster_columns = FALSE,
      use_raster = TRUE, show_column_names = FALSE, show_row_names = FALSE, col = colorRampPalette(c("#3161a9", "#34b7f8", "#d6d3b3", "#e94527", "#a31b1b"))(100),
      left_annotation = rowAnnotation(Term = rep("Signaling by Rho GTPases", num[3]), col = list(Term = c("Signaling by Rho GTPases"="#ef8433"))),
      right_annotation = rowAnnotation(r_squared = this.right_anno[enrich_result[which(enrich_result$Term == "R-HSA-194315.Signaling.by.Rho.GTPases"),]$MyList,]$r_squared, 
      avg_log2FC = this.right_anno[enrich_result[which(enrich_result$Term == "R-HSA-194315.Signaling.by.Rho.GTPases"),]$MyList,]$avg_log2FC,
      link = anno_mark(at = which(this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)), 
      labels = rownames(this.right_anno)[this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)], 
      labels_gp = gpar(fontsize = 6), padding = unit(1, "mm")), 
      col = list(r_squared = col_1, avg_log2FC = col_2), show_annotation_name = FALSE), show_heatmap_legend = FALSE)

this.right_anno <- right_anno[enrich_result[which(enrich_result$Term == "GO:0050678.regulation.of.epithelial.cell.proliferation"),]$MyList,]
ht4 <- Heatmap(heatmap_exp[enrich_result[which(enrich_result$Term == "GO:0050678.regulation.of.epithelial.cell.proliferation"),]$MyList,], cluster_columns = FALSE,
      use_raster = TRUE, show_column_names = FALSE, show_row_names = FALSE, col = colorRampPalette(c("#3161a9", "#34b7f8", "#d6d3b3", "#e94527", "#a31b1b"))(100),
      top_annotation = HeatmapAnnotation(pseudotime = pseudotime$pseudotime, cell_cluster = pseudotime$RNA_snn_res.0.3,
      col = list(pseudotime = col_pseu, cell_cluster = c("cluster1"="#db6657","cluster2"="#977dac","cluster3"="#b1ccea","cluster4"="#4876b4"))),
      left_annotation = rowAnnotation(Term = rep("regulation of epithelial cell proliferation", num[4]), col = list(Term = c("regulation of epithelial cell proliferation"="#fcea5a"))),
      right_annotation = rowAnnotation(r_squared = this.right_anno[enrich_result[which(enrich_result$Term == "GO:0050678.regulation.of.epithelial.cell.proliferation"),]$MyList,]$r_squared, 
      avg_log2FC = this.right_anno[enrich_result[which(enrich_result$Term == "GO:0050678.regulation.of.epithelial.cell.proliferation"),]$MyList,]$avg_log2FC,
      link = anno_mark(at = which(this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)), 
      labels = rownames(this.right_anno)[this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)], 
      labels_gp = gpar(fontsize = 6), padding = unit(1, "mm")), 
      col = list(r_squared = col_1, avg_log2FC = col_2), show_annotation_name = FALSE), show_heatmap_legend = TRUE,
      heatmap_legend_param = list(title = "value"))

this.right_anno <- right_anno[enrich_result[which(enrich_result$Term == "WP3888.VEGFA-VEGFR2.signaling.pathway"),]$MyList,]
ht5 <- Heatmap(heatmap_exp[enrich_result[which(enrich_result$Term == "WP3888.VEGFA-VEGFR2.signaling.pathway"),]$MyList,], cluster_columns = FALSE,
      use_raster = TRUE, show_column_names = FALSE, show_row_names = FALSE, col = colorRampPalette(c("#3161a9", "#34b7f8", "#d6d3b3", "#e94527", "#a31b1b"))(100),
      left_annotation = rowAnnotation(Term = rep("VEGFA-VEGFR2 signaling pathway", num[5]), col = list(Term = c("VEGFA-VEGFR2 signaling pathway"="#9c5a32"))),
      right_annotation = rowAnnotation(r_squared = this.right_anno[enrich_result[which(enrich_result$Term == "WP3888.VEGFA-VEGFR2.signaling.pathway"),]$MyList,]$r_squared, 
      avg_log2FC = this.right_anno[enrich_result[which(enrich_result$Term == "WP3888.VEGFA-VEGFR2.signaling.pathway"),]$MyList,]$avg_log2FC,
      link = anno_mark(at = which(this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)), 
      labels = rownames(this.right_anno)[this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)], 
      labels_gp = gpar(fontsize = 6), padding = unit(1, "mm")), 
      col = list(r_squared = col_1, avg_log2FC = col_2), show_annotation_name = FALSE), show_heatmap_legend = FALSE)

this.right_anno <- right_anno[enrich_result[which(enrich_result$Term == "GO:0002520.immune.system.development"),]$MyList,]
ht6 <- Heatmap(heatmap_exp[enrich_result[which(enrich_result$Term == "GO:0002520.immune.system.development"),]$MyList,], cluster_columns = FALSE,
      use_raster = TRUE, show_column_names = FALSE, show_row_names = FALSE, col = colorRampPalette(c("#3161a9", "#34b7f8", "#d6d3b3", "#e94527", "#a31b1b"))(100),
      left_annotation = rowAnnotation(Term = rep("immune system development", num[6]), col = list(Term = c("immune system development"="#bcdb78"))),
      right_annotation = rowAnnotation(r_squared = this.right_anno[enrich_result[which(enrich_result$Term == "GO:0002520.immune.system.development"),]$MyList,]$r_squared, 
      avg_log2FC = this.right_anno[enrich_result[which(enrich_result$Term == "GO:0002520.immune.system.development"),]$MyList,]$avg_log2FC,
      link = anno_mark(at = which(this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)), 
      labels = rownames(this.right_anno)[this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)], 
      labels_gp = gpar(fontsize = 6), padding = unit(1, "mm")), 
      col = list(r_squared = col_1, avg_log2FC = col_2), show_annotation_name = FALSE), show_heatmap_legend = FALSE)

this.right_anno <- right_anno[enrich_result[which(enrich_result$Term == "GO:0070848.response.to.growth.factor"),]$MyList,]
ht7 <- Heatmap(heatmap_exp[enrich_result[which(enrich_result$Term == "GO:0070848.response.to.growth.factor"),]$MyList,], cluster_columns = FALSE,
      use_raster = TRUE, show_column_names = FALSE, show_row_names = FALSE, col = colorRampPalette(c("#3161a9", "#34b7f8", "#d6d3b3", "#e94527", "#a31b1b"))(100),
      left_annotation = rowAnnotation(Term = rep("response to growth factor", num[7]), col = list(Term = c("response to growth factor"="#f5d0e4"))),
      right_annotation = rowAnnotation(r_squared = this.right_anno[enrich_result[which(enrich_result$Term == "GO:0070848.response.to.growth.factor"),]$MyList,]$r_squared, 
      avg_log2FC = this.right_anno[enrich_result[which(enrich_result$Term == "GO:0070848.response.to.growth.factor"),]$MyList,]$avg_log2FC,
      link = anno_mark(at = which(this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)), 
      labels = rownames(this.right_anno)[this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)], 
      labels_gp = gpar(fontsize = 6), padding = unit(1, "mm")), 
      col = list(r_squared = col_1, avg_log2FC = col_2), show_annotation_name = FALSE), show_heatmap_legend = FALSE)

this.right_anno <- right_anno[enrich_result[which(enrich_result$Term == "R-HSA-6798695.Neutrophil.degranulation"),]$MyList,]
ht8 <- Heatmap(heatmap_exp[enrich_result[which(enrich_result$Term == "R-HSA-6798695.Neutrophil.degranulation"),]$MyList,], cluster_columns = FALSE,
      use_raster = TRUE, show_column_names = FALSE, show_row_names = FALSE, col = colorRampPalette(c("#3161a9", "#34b7f8", "#d6d3b3", "#e94527", "#a31b1b"))(100),
      left_annotation = rowAnnotation(Term = rep("Neutrophil degranulation", num[8]), col = list(Term = c("Neutrophil degranulation"="#d2e9c9"))),
      right_annotation = rowAnnotation(r_squared = this.right_anno[enrich_result[which(enrich_result$Term == "R-HSA-6798695.Neutrophil.degranulation"),]$MyList,]$r_squared, 
      avg_log2FC = this.right_anno[enrich_result[which(enrich_result$Term == "R-HSA-6798695.Neutrophil.degranulation"),]$MyList,]$avg_log2FC,
      link = anno_mark(at = which(this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)), 
      labels = rownames(this.right_anno)[this.right_anno$r_squared > quantile(this.right_anno$r_squared, 0.7)], 
      labels_gp = gpar(fontsize = 6), padding = unit(1, "mm")), 
      col = list(r_squared = col_1, avg_log2FC = col_2), show_annotation_name = FALSE), show_heatmap_legend = FALSE)

ht_list = ht1 %v% ht2 %v% ht3 %v% ht4 %v% ht5 %v% ht6 %v% ht7 %v% ht8
pdf("/2_Melanocyte/3_lineage_gene/heatmap.pdf", width = 12, height = 16)
  draw(ht_list,column_km = 1)
dev.off()


ht_list = ht1 %v% ht2 %v% ht3
pdf("/2_Melanocyte/3_lineage_gene/heatmap1.pdf", width = 12, height = 8)
  draw(ht_list,column_km = 1)
dev.off()

ht_list = ht4 %v% ht5 %v% ht6 %v% ht7 %v% ht8
pdf("/2_Melanocyte/3_lineage_gene/heatmap2.pdf", width = 10, height = 8)
  draw(ht_list,column_km = 1)
dev.off()


################# Determine at which stage lineage genes begin to upregulate ######################
enrich_result <- openxlsx::read.xlsx("/2_Melanocyte/3_lineage_gene/all.t3rarsbj_/metascape_result.xlsx")
enrich_result <- enrich_result[,c(1, 7, 9, 11, 12, 13, 22, 23, 26)]
enrich_result <- reshape2::melt(enrich_result, id.vars=c("MyList"),variable.name="Term",value.name="edge")
enrich_result <- enrich_result[which(enrich_result$edge == "1.0"),]
enrich_result <- unique(enrich_result[,1:2],)

obs_p <- openxlsx::read.xlsx("/2_Melanocyte/3_lineage_gene/obs_p.xlsx")

library(Seurat)
Melanocyte <- readRDS("/2_Melanocyte/Melanocyte_sub_merge.rds")
Idents(object = Melanocyte) <- "RNA_snn_res.0.3"

cluster2_vs_1 <- FindMarkers(Melanocyte, ident.1 = "cluster2", ident.2 = "cluster1", features = unique(obs_p$gene),
min.cells.group = 1, min.cells.feature = 1, min.pct = 0, logfc.threshold = 0, only.pos = FALSE)[unique(obs_p$gene),]

cluster3_vs_1 <- FindMarkers(Melanocyte, ident.1 = "cluster3", ident.2 = "cluster1", features = unique(obs_p$gene),
min.cells.group = 1, min.cells.feature = 1, min.pct = 0, logfc.threshold = 0, only.pos = FALSE)[unique(obs_p$gene),]

cluster4_vs_1 <- FindMarkers(Melanocyte, ident.1 = "cluster4", ident.2 = "cluster1", features = unique(obs_p$gene),
min.cells.group = 1, min.cells.feature = 1, min.pct = 0, logfc.threshold = 0, only.pos = FALSE)[unique(obs_p$gene),]

avg_log2FC <- data.frame(cluster2 = cluster2_vs_1$avg_log2FC, cluster3 = cluster3_vs_1$avg_log2FC, cluster4 = cluster4_vs_1$avg_log2FC)
rownames(avg_log2FC) <- unique(obs_p$gene)
sig_p <- data.frame(cluster2 = cluster2_vs_1$p_val_adj, cluster3 = cluster3_vs_1$p_val_adj, cluster4 = cluster4_vs_1$p_val_adj)
avg <- AverageExpression(Melanocyte, features = unique(obs_p$gene), slot = "scale.data")$RNA

up_in_cluster2 <- rownames(avg_log2FC[which(avg_log2FC$cluster2 > 1 & avg_log2FC$cluster2 > 1 & sig_p$cluster2 < 0.05 & sig_p$cluster3 < 0.05),])
up_in_cluster3 <- setdiff(rownames(avg_log2FC[which(avg_log2FC$cluster3 > 1 & sig_p$cluster3 < 0.05),]), up_in_cluster2)
up_in_cluster4 <- setdiff(obs_p$gene, c(up_in_cluster2, up_in_cluster3))
up_cluster <- list(up_in_cluster2=up_in_cluster2, up_in_cluster3=up_in_cluster3, up_in_cluster4=up_in_cluster4)

saveRDS(up_cluster, file = "/2_Melanocyte/3_lineage_gene/up_cluster.rds")
openxlsx::write.xlsx(up_cluster, file = "/2_Melanocyte/3_lineage_gene/up_cluster.xlsx")

up_cluster1 <- list(up_in_cluster2 = cluster2_vs_1[up_in_cluster2,], 
up_in_cluster3 = cluster3_vs_1[up_in_cluster3,],
up_in_cluster4 = cluster4_vs_1[up_in_cluster4,])
openxlsx::write.xlsx(up_cluster1, row.names = TRUE, file = "/2_Melanocyte/3_lineage_gene/up_cluster1.xlsx")

p <- StackedVlnPlot(obj = Melanocyte, features = up_in_cluster2, plot.margin = unit(c(-0.75, 0, -0.75, 0)),cols=c(ggsci::pal_npg("nrc",alpha = 1)(10), ggsci::pal_nejm("default",alpha = 1)(8)),ncol=1, test_sign=NULL)
ggsave(p, file="/2_Melanocyte/3_lineage_gene/up_in_cluster2_VlnPlot.pdf", width=3, height=15)


up_cluster <- readRDS("/2_Melanocyte/3_lineage_gene/up_cluster.rds")
up_cluster <- sapply(up_cluster, function(x){
      x <- x[x %in% enrich_result[,1]]
      return(x)
})
up_in_cluster2 <- up_cluster$up_in_cluster2
up_in_cluster3 <- up_cluster$up_in_cluster3
up_in_cluster4 <- up_cluster$up_in_cluster4

avg_log2FC <- avg_log2FC[rownames(avg_log2FC) %in% enrich_result[,1],]
avg <- avg[rownames(avg) %in% enrich_result[,1],]

library(ComplexHeatmap)
library(circlize)
col_1 = colorRamp2(seq(min(avg_log2FC), max(avg_log2FC), length.out = 50), colorRampPalette(c("#ffffff", "#c42d53"))(50))
col_2 = colorRamp2(seq(min(avg), max(avg), length.out = 50), colorRampPalette(c("#ffffff", "#2278c3"))(50))

ht1 <- Heatmap(avg_log2FC[rownames(avg_log2FC) %in% up_in_cluster2, ], cluster_columns = FALSE,
      use_raster = FALSE, show_column_names = FALSE, show_row_names = FALSE, col = col_1,
      left_annotation = rowAnnotation(Term = rep("up_in_cluster2", length(up_in_cluster2)), col = list(Term = c("up_in_cluster2"="#2d6ca4"))), show_heatmap_legend = TRUE, heatmap_legend_param = list(title = "avg_log2FC vs cluster1"))
ht2 <- Heatmap(avg_log2FC[rownames(avg_log2FC) %in% up_in_cluster3, ], cluster_columns = FALSE,
      use_raster = FALSE, show_column_names = FALSE, show_row_names = FALSE, col = col_1,
      left_annotation = rowAnnotation(Term = rep("up_in_cluster3", length(up_in_cluster3)), col = list(Term = c("up_in_cluster3"="#cd8641"))), show_heatmap_legend = FALSE)
ht3 <- Heatmap(avg_log2FC[rownames(avg_log2FC) %in% up_in_cluster4, ], cluster_columns = FALSE,
      use_raster = FALSE, show_column_names = FALSE, show_row_names = FALSE, col = col_1,
      left_annotation = rowAnnotation(Term = rep("up_in_cluster4", length(up_in_cluster4)), col = list(Term = c("up_in_cluster4"="#3d7d51"))), show_heatmap_legend = FALSE)

pdf("/2_Melanocyte/3_lineage_gene/temp.pdf", width = 6, height = 12)
ht1_draw <- draw(ht1)
dev.off()

pdf("/2_Melanocyte/3_lineage_gene/temp.pdf", width = 6, height = 12)
ht2_draw <- draw(ht2)
dev.off()

pdf("/2_Melanocyte/3_lineage_gene/temp.pdf", width = 6, height = 12)
ht3_draw <- draw(ht3)
dev.off()

ht1_order <- row_order(ht1_draw)
ht2_order <- row_order(ht2_draw)
ht3_order <- row_order(ht3_draw)

avg <- as.data.frame(avg)
ht4 <- Heatmap(avg[rownames(avg) %in% up_in_cluster2, ][ht1_order,], cluster_columns = FALSE, cluster_rows = FALSE,
      use_raster = FALSE, show_column_names = FALSE, show_row_names = FALSE, col = col_2,
      show_heatmap_legend = TRUE, heatmap_legend_param = list(title = "avg_value"),
      right_annotation = rowAnnotation(link = anno_mark(at = which(avg[rownames(avg) %in% up_in_cluster2, ][ht1_order,]$cluster2 > quantile(avg[rownames(avg) %in% up_in_cluster2, ][ht1_order,]$cluster2, 0.7)), 
      labels = rownames(avg[rownames(avg) %in% up_in_cluster2, ])[avg[rownames(avg) %in% up_in_cluster2, ][ht1_order,]$cluster2 > quantile(avg[rownames(avg) %in% up_in_cluster2, ][ht1_order,]$cluster2, 0.7)], 
      labels_gp = gpar(fontsize = 6), padding = unit(1, "mm"))))

ht5 <- Heatmap(avg[rownames(avg) %in% up_in_cluster3, ][ht2_order,], cluster_columns = FALSE, cluster_rows = FALSE,
      use_raster = FALSE, show_column_names = FALSE, show_row_names = FALSE, col = col_2,
      show_heatmap_legend = FALSE, heatmap_legend_param = list(title = "avg_value"),
      right_annotation = rowAnnotation(link = anno_mark(at = which(avg[rownames(avg) %in% up_in_cluster3, ][ht2_order,]$cluster3 > quantile(avg[rownames(avg) %in% up_in_cluster3, ][ht2_order,]$cluster3, 0.7)), 
      labels = rownames(avg[rownames(avg) %in% up_in_cluster3, ])[avg[rownames(avg) %in% up_in_cluster3, ][ht2_order,]$cluster3 > quantile(avg[rownames(avg) %in% up_in_cluster3, ][ht2_order,]$cluster3, 0.7)], 
      labels_gp = gpar(fontsize = 6), padding = unit(1, "mm"))))

ht6 <- Heatmap(avg[rownames(avg) %in% up_in_cluster4, ][ht3_order,], cluster_columns = FALSE, cluster_rows = FALSE,
      use_raster = FALSE, show_column_names = FALSE, show_row_names = FALSE, col = col_2,
      show_heatmap_legend = FALSE, heatmap_legend_param = list(title = "avg_value"),
      right_annotation = rowAnnotation(link = anno_mark(at = which(avg[rownames(avg) %in% up_in_cluster4, ][ht3_order,]$cluster4 > quantile(avg[rownames(avg) %in% up_in_cluster4, ][ht3_order,]$cluster4, 0.7)), 
      labels = rownames(avg[rownames(avg) %in% up_in_cluster4, ])[avg[rownames(avg) %in% up_in_cluster4, ][ht3_order,]$cluster4 > quantile(avg[rownames(avg) %in% up_in_cluster4, ][ht3_order,]$cluster4, 0.7)], 
      labels_gp = gpar(fontsize = 6), padding = unit(1, "mm"))))

ht_list = ht1 %v% ht2 %v% ht3
pdf("/2_Melanocyte/3_lineage_gene/avg_log2FC.pdf", width = 5, height = 12)
  draw(ht_list,column_km = 1)
dev.off()

ht_list = ht4 %v% ht5 %v% ht6
pdf("/2_Melanocyte/3_lineage_gene/avg.pdf", width = 3, height = 12)
  draw(ht_list,column_km = 1)
dev.off()


############# Pathways enriched for upregulated genes #############
enrich_result[,2] <- as.character(enrich_result[,2])
enrich_result[,2] <- gsub("\\.", " ", enrich_result[,2])
enrich_result[,2] <- gsub("R-HSA-163210 ","",enrich_result[,2])
enrich_result[,2] <- gsub("R-HSA-9006934 ","",enrich_result[,2])
enrich_result[,2] <- gsub("R-HSA-194315 ","",enrich_result[,2])
enrich_result[,2] <- gsub("GO:0050678 ","",enrich_result[,2])
enrich_result[,2] <- gsub("GO:0002520 ","",enrich_result[,2])
enrich_result[,2] <- gsub("GO:0070848 ","",enrich_result[,2])
enrich_result[,2] <- gsub("R-HSA-6798695 ","",enrich_result[,2])
enrich_result[,2] <- gsub("WP3888 ","",enrich_result[,2])

terms <- unique(enrich_result[,2])
fisher_re <- c()
for(i in 1:length(terms)){
      #up in cluster2
      enrich_matrix <- as.matrix(data.frame(up_in_cluster = c(length(unique(intersect(up_in_cluster2, enrich_result[which(enrich_result[,2] == terms[i]),1]))), length(unique(intersect(up_in_cluster2, enrich_result[which(enrich_result[,2] != terms[i]),1])))),
      up_in_other = c(length(unique(intersect(c(up_in_cluster3, up_in_cluster4), enrich_result[which(enrich_result[,2] == terms[i]),1]))), length(unique(intersect(c(up_in_cluster3, up_in_cluster4), enrich_result[which(enrich_result[,2] != terms[i]),1]))))))
      this.test <- fisher.test(enrich_matrix)
      this.fisher_re <- c(this.test$estimate, this.test$conf.int[1], this.test$conf.int[2], this.test$p.value, "up_in_cluster2", terms[i])
      fisher_re <- rbind(fisher_re, this.fisher_re)
      
      #up in cluster3
      enrich_matrix <- as.matrix(data.frame(up_in_cluster = c(length(unique(intersect(up_in_cluster3, enrich_result[which(enrich_result[,2] == terms[i]),1]))), length(unique(intersect(up_in_cluster3, enrich_result[which(enrich_result[,2] != terms[i]),1])))),
      up_in_other = c(length(unique(intersect(c(up_in_cluster2, up_in_cluster4), enrich_result[which(enrich_result[,2] == terms[i]),1]))), length(unique(intersect(c(up_in_cluster2, up_in_cluster4), enrich_result[which(enrich_result[,2] != terms[i]),1]))))))
      this.test <- fisher.test(enrich_matrix)
      this.fisher_re <- c(this.test$estimate, this.test$conf.int[1], this.test$conf.int[2], this.test$p.value, "up_in_cluster3", terms[i])
      fisher_re <- rbind(fisher_re, this.fisher_re)

      #up in cluster4
      enrich_matrix <- as.matrix(data.frame(up_in_cluster = c(length(unique(intersect(up_in_cluster4, enrich_result[which(enrich_result[,2] == terms[i]),1]))), length(unique(intersect(up_in_cluster4, enrich_result[which(enrich_result[,2] != terms[i]),1])))),
      up_in_other = c(length(unique(intersect(c(up_in_cluster2, up_in_cluster3), enrich_result[which(enrich_result[,2] == terms[i]),1]))), length(unique(intersect(c(up_in_cluster2, up_in_cluster3), enrich_result[which(enrich_result[,2] != terms[i]),1]))))))
      this.test <- fisher.test(enrich_matrix)
      this.fisher_re <- c(this.test$estimate, this.test$conf.int[1], this.test$conf.int[2], this.test$p.value, "up_in_cluster4", terms[i])
      fisher_re <- rbind(fisher_re, this.fisher_re)
}
fisher_re <- as.data.frame(fisher_re)
colnames(fisher_re) <- c("odds_ratio", "conf.int1", "conf.int2", "p", "class", "term")
fisher_re[,1] <- as.numeric(fisher_re[,1])
fisher_re[,2] <- as.numeric(fisher_re[,2])
fisher_re[,3] <- as.numeric(fisher_re[,3])
fisher_re[,4] <- as.numeric(fisher_re[,4])

library(ggplot2)
library(ggthemes)
fisher_re$class <- factor(fisher_re$class, levels=c("up_in_cluster2","up_in_cluster3","up_in_cluster4"))
p <- ggplot(fisher_re, aes(x = class, y = term)) + geom_tile(aes(fill = odds_ratio)) + 
            scale_fill_gradient2(low="#5bacdb", mid="#ffffff", high="#ed97b3", midpoint = 1, # <-- look at this
            limit = c(0, 2.3)) + geom_point(aes(colour="white", size=-log10(p)), shape = 21,colour = "black") + 
            scale_size_continuous(range=c(2,10)) + theme_few() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(p, filename = "/2_Melanocyte/3_lineage_gene/fisher_re.pdf", width = 6, height = 5.5)


########################### Calculate the specificity of lineage genes in melanocytes ###################
obs_p <- openxlsx::read.xlsx("/2_Melanocyte/3_lineage_gene/obs_p.xlsx")
obs_p <- obs_p[rev(order(obs_p$r_squared)),]

D1_D5 <- readRDS("/1_Cell_annotation/5_annotation/D1_D5_final.rds")
D1_D5@meta.data[which(D1_D5@meta.data$cell.type == "Mitotic Melanocyte cell"), ]$cell.type <- "Melanocyte"
D1_D5 <- SetIdent(D1_D5, value = "cell.type")
cluster_mel_vs_other <- FindMarkers(D1_D5, ident.1 = "Melanocyte", ident.2 = NULL, features = obs_p$gene,
min.cells.group = 1, min.cells.feature = 1, min.pct = 0, logfc.threshold = 0, only.pos = FALSE)[obs_p$gene,]
cluster_mel_vs_other_sig <- cluster_mel_vs_other[which(cluster_mel_vs_other$p_val_adj < 0.05 & cluster_mel_vs_other$avg_log2FC > 1),]

saveRDS(cluster_mel_vs_other_sig, file = "/2_Melanocyte/3_lineage_gene/cluster_mel_vs_other_sig.rds")
