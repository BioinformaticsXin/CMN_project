library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(Seurat)
library(SeuratDisk)

Endothelial_cell <- readRDS("/5_Endothelial cell/endothelial_cell_annotation.rds")

# subclustering of capillary_EC
capillary_EC <- subset(Endothelial_cell, idents = "capillary ECs")
all.genes <- rownames(capillary_EC)
capillary_EC <- ScaleData(capillary_EC, features = all.genes)
capillary_EC <- RunPCA(capillary_EC, features = VariableFeatures(object = capillary_EC))
capillary_EC <- RunHarmony(seurat_object = capillary_EC, group.by.vars = c("SampleID", "Tech"))
capillary_EC <- FindNeighbors(capillary_EC, FindNeighbor_reduction = "harmony")
capillary_EC <- FindClusters(seurat_object, resolution = 0.8, algorithm = 1, method = "matrix")

# score of Stalk/Tip
Tip_Stalk <- openxlsx::read.xlsx("/5_Endothelial_cell/3_sub_capillary_EC/tipMarker.xlsx")
Tip <- na.omit(as.character(Tip_Stalk[,1]))
Stalk <- na.omit(as.character(Tip_Stalk[,2]))
Tip_Stalk <- list(Tip=as.character(Tip), Stalk=as.character(Stalk))

# GSVA
library(GSVA)
gsva_res <- gsva(expr=as.matrix(capillary_EC[["RNA"]]@scale.data), gset.idx.list=Tip_Stalk, 
method="gsva", kcdf="Gaussian", ssgsea.norm=TRUE, parallel.sz=50)
saveRDS(gsva_res, file = "/5_Endothelial_cell/3_sub_capillary_EC/Tip_Stalk_score.rds")

#AUCell
library(AUCell)
library(clusterProfiler)
library(ggplot2)
library(Seurat)
library(Matrix)
#Build the rankings
Seurat_matrix <- as.matrix(capillary_EC@assays$RNA@data)
cells_rankings <- AUCell_buildRankings(Seurat_matrix, plotStats = FALSE) 

#Calculate the Area Under the Curve (AUC)
cells_AUC <- AUCell_calcAUC(Tip_Stalk, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
saveRDS(cells_AUC, file= "/5_Endothelial_cell/3_sub_capillary_EC/Tip_Stalk_AUC.rds")

####################
capillary_EC <- AddMetaData(capillary_EC, metadata = as.data.frame(t(cells_AUC@assays@data$AUC)))
saveRDS(capillary_EC, file= "/5_Endothelial_cell/3_sub_capillary_EC/capillary_EC.rds")


### umap of capillary_EC
capillary_EC@meta.data$UMAP_1 <- capillary_EC@reductions$umap@cell.embeddings[,1]
capillary_EC@meta.data$UMAP_2 <- capillary_EC@reductions$umap@cell.embeddings[,2]

p <- ggplot(capillary_EC@meta.data, aes(color = RNA_snn_res.0.8, x = UMAP_1, y = UMAP_2)) + 
geom_point(shape = 16, size = 1.5, alpha = .5) +
theme_minimal() +
scale_color_manual(values=c("0"="#4cbbd2", "1"="#e9949e", "2"="#9fc6db", "3"="#efb4cc", "4"="#c4aca6"))
ggsave(filename= "/5_Endothelial_cell/3_sub_capillary_EC/umap.pdf",plot=p,width=4.4,height=2.5)


# violin plot of tip_stalk score
Tip_Stalk <- capillary_EC$Stalk / capillary_EC$Tip
capillary_EC$Tip_Stalk <- Tip_Stalk

plot_frame <- capillary_EC@meta.data
plot_frame$RNA_snn_res.0.8 <- paste0("Cluster", plot_frame$RNA_snn_res.0.8)

plot_frame$RNA_snn_res.0.8 <- factor(plot_frame$RNA_snn_res.0.8, levels = c(paste0("Cluster", 0:4)))
plot_frame[which(plot_frame$Tip_Stalk > 2),]$Tip_Stalk <- 2

p1 <- ggplot(plot_frame, aes(RNA_snn_res.0.8, Tip_Stalk, fill = RNA_snn_res.0.8))+
  geom_violin() + 
  scale_fill_manual(values = c("Cluster0"="#4cbbd2", "Cluster1"="#e9949e", "Cluster2"="#9fc6db", "Cluster3"="#efb4cc", "Cluster4"="#c4aca6")) +
  stat_summary(fun= mean, geom = "point",
               shape = 19, size = 2, color = "black")+
  theme_classic()+
  geom_hline(aes(yintercept=1), colour="#565354", linetype="dashed")+
  xlab('')+
  ylab('Stalk / Tip')+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y = element_text(size = 12,
                                   face="bold"),
        legend.title=element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))

ggsave(p1, file = '/5_Endothelial_cell/3_sub_capillary_EC/violin_tip_stalk.pdf', width = 3, height = 2)


# subcluster_enrichment
#######################################################################################
#################################calculate Pearson Residuals###########################
#######################################################################################
# create the table of cell numbers
library(dplyr)
table_count <- data.frame(capillary_EC@meta.data %>% group_by(SampleType, RNA_snn_res.0.8) %>% summarise(n = n()))
table_count <- reshape2::dcast(table_count, SampleType~RNA_snn_res.0.8)
rownames(table_count) <- table_count[,1]
table_count <- table_count[,-1]
table_count[is.na(table_count)] <- 0
colnames(table_count) <- paste0("cluster", colnames(table_count))
table_count <- table_count[-4,]

# Calculate Pearson Residuals according to the table_count
pearson_residuals <- matrix(0, nrow = 3, ncol = 5)
colnames(pearson_residuals) <- colnames(table_count)
rownames(pearson_residuals) <- rownames(table_count)
for(i in 1:nrow(table_count)){
    for(j in 1:ncol(table_count)){
        #print(paste(i,j,sep=" "))
        this.expected <- sum(table_count[i,]) * sum(table_count[,j])/sum(table_count)
        this.observation <- table_count[i,j]
        this.r <- (this.observation-this.expected)/sqrt(this.expected)
        pearson_residuals[i,j] <- this.r
    }
}

pearson_residuals <- reshape2::melt(pearson_residuals, id.vars=c("SampleType"),variable.name="Cluster",value.name="pearson_residuals")
colnames(pearson_residuals)[1:2] <- c("SampleType", "Cluster")

##########################################################################################
###################################Cell Fraction##########################################
##########################################################################################
cell_fraction <- table_count
for(i in 1:ncol(cell_fraction)){
    cell_fraction[,i] <- cell_fraction[,i]/sum(cell_fraction[,i])
}
cell_fraction <- reshape2::melt(data.frame(SampleType = rownames(cell_fraction), cell_fraction), value.name="cell_fraction")
cell_fraction[,2] <- as.character(cell_fraction[,2])
cell_fraction[,2] <- gsub("\\.", " ", cell_fraction[,2])
colnames(cell_fraction)[1:2] <- c("SampleType", "Cluster")

enrichment_result <- merge(pearson_residuals, cell_fraction)

###########################################################################################
#######################################visualization#######################################
###########################################################################################
library(ggplot2)
library(ggthemes)
enrichment_result$Cluster <- as.character(enrichment_result$Cluster)
enrichment_result$SampleType <- factor(enrichment_result$SampleType, levels=c("Nevus","Malignant Nevus","Melanoma"))
p <- ggplot(enrichment_result, aes(x = Cluster, y = SampleType)) + geom_tile(aes(fill = pearson_residuals)) + 
scale_fill_gradient2(low = "#64bfe7", mid = "white", high = "#de7479") + geom_point(aes(colour="white", size=cell_fraction), shape = 21,colour = "black") + 
scale_size_continuous(range=c(2,10)) + theme_few() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(p, filename =  "/5_Endothelial_cell/3_sub_capillary_EC/cell_enrich.pdf", width = 6.2, height = 2.8)

