library(Seurat)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(dplyr)

Myeloid_cell <- readRDS("/4_Myeloid/Myeloid_cell_annotation.rds")
Myeloid_color = c("CD209+ SPP1- TAMs"="#FCB8B9", "MDSCs"="#A7BBD6", "CD1C DCs"="#BA94C7", "CD209- SPP1+ TAMs"="#ECBA7D", "CD209- SPP1low TAMs"="#EEDFCA", "Langerhans Cells"="#F6F195")
p1 <- DimPlot(Myeloid_cell, reduction = "umap", group.by = "SampleType", pt.size = 1.2, cols = c("#FCB8B9", "#A7BBD6", "#BA94C7", "#ECBA7D", "#EEDFCA", "#F6F195"), label = FALSE)
ggsave(p1, file ="/4_Myeloid/1_visualization/umap.pdf", width=7.2, height=5)

features <- c("CD163", "CD209", "MRC1", "CCL18", "CD14", "SPP1", "S100A8", "S100A9", "MMP19", "CLEC10A", "CD1C", "CD207")
plots <- lapply(features, function(gene) {
  VlnPlot(Myeloid_cell, features = gene, pt.size = 0, cols = Myeloid_color) + 
  theme(legend.position = "none", plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"))
})
p <- wrap_plots(plots, ncol = 1)  # Stack the plots vertically with 1 column
ggsave(filename ="/4_Myeloid/1_visualization/VlnPlot.pdf", plot = p, width = 3.5, height = 5)

# subcluster_enrichment
#### Enrichment of each cluster of myeloid cells in different SampleTypes ("Nevus","Malignant Nevus","Melanoma") 

#######################################################################################
#################################calculate Pearson Residuals###########################
#######################################################################################
# create the table of cell numbers
library(dplyr)
table_count <- data.frame(Myeloid_cell@meta.data %>% group_by(SampleType, cell.type) %>% summarise(n = n()))
table_count <- reshape2::dcast(table_count, SampleType~cell.type)
rownames(table_count) <- table_count[,1]
table_count <- table_count[,-1]
table_count[is.na(table_count)] <- 0

# Calculate Pearson Residuals according to the table_count
pearson_residuals <- matrix(0, nrow = 4, ncol = 6)
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
cell_fraction[which(cell_fraction[,2] == "CD209..SPP1..TAMs"),2] <- "CD209+ SPP1- TAMs"
cell_fraction[which(cell_fraction[,2] == "CD209..SPP1..TAMs.1"),2] <- "CD209- SPP1+ TAMs"
cell_fraction[which(cell_fraction[,2] == "CD209..SPP1low.TAMs"),2] <- "CD209- SPP1low TAMs"
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
p <- ggplot(enrichment_result[which(enrichment_result$SampleType != "Normal"),], aes(x = Cluster, y = SampleType)) + geom_tile(aes(fill = pearson_residuals)) + 
scale_fill_gradient2(low = "#64bfe7", mid = "white", high = "#de7479") + geom_point(aes(colour="white", size=cell_fraction), shape = 21,colour = "black") + 
scale_size_continuous(range=c(2,10)) + theme_few() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + coord_flip()
ggsave(p, filename ="/4_Myeloid/1_visualization/cell_enrich.pdf", width = 5.2, height = 4.2)

