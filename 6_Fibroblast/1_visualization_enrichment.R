library(Seurat)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(dplyr)
library(patchwork)

CAF <- readRDS("/6_CAF/CAF.rds")

# umap
CAF@meta.data$UMAP_1 <- CAF@reductions$umap@cell.embeddings[,1]
CAF@meta.data$UMAP_2 <- CAF@reductions$umap@cell.embeddings[,2]
p <- ggplot(CAF@meta.data, aes(color = cell.type, x = UMAP_1, y = UMAP_2)) + 
geom_point(shape = 16, size = 1.5, alpha = .5) +
theme_minimal() +
scale_color_manual(values=c("metaCAF - FGFBP2"="#ce5d7b", "angiCAF - CCN5"="#ae7b82", "TactCAF - SLPI"="#4165a7", 
"angiCAF - TPD52"="#cf8392", "infCAF - PLAU"="#d7dfab", "angiCAF - COCH"="#659bb0", "TdiffCAF - ZFP36"="#765b9a",
"apoCAF - NAMPT"="#6c8940", "ecmCAF - NPTX2"="#2f4862", "ecmCAF - COL11A1"="#d8bb98", "antCAF - COL6A5"="#8f3165",
"angiCAF - APCDD1"="#ba9314", "antCAF - C3"="#33c8e6", "antCAF - CCL19"="#646567", "CAF - AHNAK"="#f78b07"))
ggsave(filename="/6_CAF/1_visuliztion/umap.pdf",plot=p,width=5,height=3)

p <- ggplot(CAF@meta.data, aes(color = SampleType, x = UMAP_1, y = UMAP_2)) + 
geom_point(shape = 16, size = 1.5, alpha = .5) +
theme_minimal() +
scale_color_manual(values=c("#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF"))
ggsave(filename="/6_CAF/1_visuliztion/umap_SampleType.pdf",plot=p,width=4.7,height=3)

# subcluster_enrichment
#######################################################################################
#################################calculate Pearson Residuals###########################
#######################################################################################
# create the table of cell numbers
library(dplyr)
table_count <- data.frame(CAF@meta.data %>% group_by(SampleType, cell.type) %>% summarise(n = n()))
table_count <- reshape2::dcast(table_count, SampleType~cell.type)
rownames(table_count) <- table_count[,1]
table_count <- table_count[,-1]
table_count[is.na(table_count)] <- 0

# Calculate Pearson Residuals according to the table_count
pearson_residuals <- matrix(0, nrow = 4, ncol = 15)
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
cell_fraction[,2] <- gsub("\\.\\.\\.", " - ", cell_fraction[,2])
colnames(cell_fraction)[1:2] <- c("SampleType", "Cluster")

enrichment_result <- merge(pearson_residuals, cell_fraction)

###########################################################################################
#######################################visualization#######################################
###########################################################################################
library(ggplot2)
library(ggthemes)
enrichment_result$Cluster <- as.character(enrichment_result$Cluster)
enrichment_result$SampleType <- factor(enrichment_result$SampleType, levels=c("Nevus","Malignant Nevus","Melanoma"))
p <- ggplot(enrichment_result[which(enrichment_result$SampleType != "Normal"),], aes(x = Cluster, y = SampleType)) + geom_tile(aes(fill = pearson_residuals)) + 
scale_fill_gradient2(low = "#64bfe7", mid = "white", high = "#de7479") + geom_point(aes(colour="white", size=cell_fraction), shape = 21,colour = "black") + 
scale_size_continuous(range=c(2,10)) + theme_few() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(p, filename =  "/6_CAF/1_visuliztion/cell_enrich.pdf", width = 10.5, height = 3)



############featureplot##############
CAF <- readRDS("/6_CAF/1_visuliztion/CAF.rds")
p <- FeaturePlot(CAF, pt.size = 1, features = c("DCN", "COL1A1", "COL6A2", "PDGFRA", "PDGFRB", "ACTA2", "RGS5"), cols = c("#f8f9df", "#f2814b", "#741f2b"), ncol = 7)
ggsave(p, file= "/6_CAF/1_visuliztion/CAF_feature.pdf", width=37, height=6)



