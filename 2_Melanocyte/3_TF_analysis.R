library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(foreach)
library(AUCell)
library(ComplexHeatmap)
library(circlize)

#TF_analysis
Melanocyte <- readRDS("/2_Melanocyte/Melanocyte_sub_merge.rds")
cellInfo <- Melanocyte@meta.data[,c("SampleType", "RNA_snn_res.0.3")]
colnames(cellInfo)[2] <- "CellType"
saveRDS(cellInfo, file = "/2_Melanocyte/3_lineage_gene/network/SCENIC/int/cellInfo.Rds")

setwd("/2_Melanocyte/3_lineage_gene/network/SCENIC")
org <- "hgnc"
dbDir <- "/SCENIC/cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "Melanocyte TF analysis" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=300)
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# Co-expression network
exprMat <- as.matrix(Melanocyte@assays$RNA@counts)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,minCountsPerGene=3*.01*ncol(exprMat),minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)

##GENIE3
exprMat_filtered <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered, scenicOptions, nParts = 10)

exprMat_log <- log2(exprMat+1)
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_log)

#cell-type specific regulators
scenicOptions <- readRDS("int/scenicOptions.Rds")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
cellInfo <- readRDS("int/cellInfo.Rds")
scenicOptions <- readRDS("int/scenicOptions.Rds")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
saveRDS(rss, file = "output/rss.Rds")

binaryRegulonActivity <- readRDS("int/4.1_binaryRegulonActivity.Rds")
binaryRegulonActivity <- binaryRegulonActivity[!grepl("_extended", rownames(binaryRegulonActivity)),]


rss <- rss[!grepl("_extended", rownames(rss)),]
rss <- as.data.frame(rss)
cluster1_RSS <- data.frame(Cluster = rep("cluster1", nrow(rss)), TF = rownames(rss[rev(order(rss$cluster1)),]), RSS = rss[rev(order(rss$cluster1)),]$cluster1, Regulons = 1:nrow(rss), stringsAsFactors = FALSE)
cluster2_RSS <- data.frame(Cluster = rep("cluster2", nrow(rss)), TF = rownames(rss[rev(order(rss$cluster2)),]), RSS = rss[rev(order(rss$cluster2)),]$cluster2, Regulons = 1:nrow(rss), stringsAsFactors = FALSE)
cluster3_RSS <- data.frame(Cluster = rep("cluster3", nrow(rss)), TF = rownames(rss[rev(order(rss$cluster3)),]), RSS = rss[rev(order(rss$cluster3)),]$cluster3, Regulons = 1:nrow(rss), stringsAsFactors = FALSE)
cluster4_RSS <- data.frame(Cluster = rep("cluster4", nrow(rss)), TF = rownames(rss[rev(order(rss$cluster4)),]), RSS = rss[rev(order(rss$cluster4)),]$cluster4, Regulons = 1:nrow(rss), stringsAsFactors = FALSE)
RSS <- rbind(cluster1_RSS, cluster2_RSS, cluster3_RSS, cluster4_RSS)
openxlsx::write.xlsx(RSS, file = "RSS.xlsx")

RSS$label <- RSS$TF
RSS[which(RSS$Regulons >15),]$label <- NA

colnames(RSS)[3] <- "rss"
p <- ggplot(RSS,aes(x=Regulons, y=rss, fill = rss, color = rss))+
  geom_point(size=3.2,stroke = 0.5,shape = 21)+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#143962","#2888bd","#95c7dd","#ffffff","#e09786","#bb2438","#660022"))+
  scale_fill_gradientn(values = seq(0,1,0.2),
                        colors = alpha(c("#143962","#2888bd","#95c7dd","#ffffff","#e09786","#bb2438","#660022"), 0.3))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_text_repel(aes(label=label), color = "#313131", size = 3, segment.color="#313131", show.legend=FALSE)+
  facet_wrap(~Cluster, scale = "free", ncol = 4)
ggsave(p, file="/2_Melanocyte/3_lineage_gene/network/SCENIC/RSS_plot.pdf", width=9, height=4.5)


###################binaryRegulonActivity###############
TF_plot <- unique(RSS[which(RSS$Regulons <= 10),]$TF)
TF_omit <- c("DDIT3 (19g)", "STAT3 (48g)", "LHX2 (142g)", "ETV1 (4576g)", "LUZP1 (58g)", "JUNB (224g)", "USF2 (452g)", 
             "GTF2F1 (2934g)", "REL (108g)", "JUN (246g)", "FOSB (78g)", "FOS (452g)", "FOSL1 (475g)", "CEBPD (79g)", "TGIF1 (319g)")
TF_plot <- TF_plot[!(TF_plot %in% TF_omit)]

Melanocyte <- readRDS("/2_Melanocyte/Melanocyte_sub_merge.rds")
pseudotime <- read.csv("/2_Melanocyte/2_Trajectory/monocle3/pseudotime.csv")
pseudotime <- pseudotime[order(pseudotime$pseudotime),]

cellInfo <- pseudotime[,c("SampleType", "RNA_snn_res.0.3", "pseudotime")]
colnames(cellInfo) <- c("SampleType", "Cluster", "Pseudotime")
rownames(cellInfo) <- pseudotime[,1]
binaryRegulonActivity <- binaryRegulonActivity[,rownames(cellInfo)]

col_pseu = colorRamp2(seq(min(pseudotime$pseudotime), max(pseudotime$pseudotime), length.out = 50), colorRampPalette(c("#ffffff", "#eb39b3"))(50))

pdf("/2_Melanocyte/3_lineage_gene/network/SCENIC/binaryRegulonActivity_heatmap.pdf", width = 7, height = 6)
Heatmap(binaryRegulonActivity[TF_plot,], cluster_columns = FALSE, use_raster = TRUE, 
        show_column_names = FALSE, show_row_names = TRUE, col = colorRampPalette(c("#ffffff", "#000000"))(10),
        top_annotation = HeatmapAnnotation(SampleType = cellInfo$SampleType, Cluster = cellInfo$Cluster,
        col = list(SampleType = c("Malignant Nevus"="#7876B1FF", "Melanoma"="#6F99ADFF", "Nevus"="#FFDC91FF", "Normal"="#EE4C97FF"), 
        Cluster = c("cluster1"="#db6657","cluster2"="#977dac","cluster3"="#b1ccea","cluster4"="#4876b4"))),
        show_heatmap_legend = TRUE, row_names_gp = gpar(fontsize = 7),
        heatmap_legend_param = list(title = "value"))
dev.off()

