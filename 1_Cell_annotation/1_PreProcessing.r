library(dplyr)
library(Seurat)
#for BD
Tag_Sample_list <- list(D1 = data.frame(Tag = c("SampleTag07_hs","SampleTag12_hs","SampleTag05_hs"), SampleType = c("Nevus","Malignant Nevus","Melanoma")),
                        D2 = data.frame(Tag = c("SampleTag01_hs","SampleTag02_hs","SampleTag03_hs"), SampleType = c("Nevus","Malignant Nevus","Melanoma")),
                        D3 = data.frame(Tag = c("SampleTag02_hs","SampleTag01_hs"), SampleType = c("Malignant Nevus","Melanoma")),
                        D4 = data.frame(Tag = c("SampleTag04_hs","SampleTag08_hs"), SampleType = c("Normal","Melanoma")))

object_list <- list()
for(Donor_id in names(Tag_Sample_list)){
    RSEC <- read.table(paste0("/run_raw_data/", Donor_id, "_RSEC_MolsPerCell.csv"), sep=",", check.names=F, header=T, row.names=1, stringsAsFactors=F, quote = "")
    Tag <- read.table(paste0("/run_raw_data/", Donor_id, "_Sample_Tag_Calls.csv"), sep=",", check.names=F, header=T, stringsAsFactors=F)
    Tag <- Tag[Tag[,3] %in% Tag_Sample_list[[Donor_id]]$Tag,]
    RSEC <- RSEC[rownames(RSEC) %in% as.character(Tag[,1]),]
    Tag <- Tag[as.character(Tag[,1]) %in% rownames(RSEC),]

    for(i in 1:dim(Tag_Sample_list[[Donor_id]])[1]){
        Tag[which(Tag[,2] == as.character(Tag_Sample_list[[Donor_id]][i,1])),3] <- as.character(Tag_Sample_list[[Donor_id]][i,2])
    }
    colnames(Tag)[3] <- "SampleType"

    RSEC <- RSEC[as.character(Tag[,1]),]
    Tag[,1] <- paste(rep(Donor_id,dim(Tag)[1]),Tag[,1],Tag[,2],sep="_")
    rownames(RSEC) <- paste(rep(Donor_id,dim(Tag)[1]),rownames(RSEC),Tag[,2],sep="_")

    Tag <- data.frame(Donor = rep(Donor_id, dim(Tag)[1]), Cell_Index = Tag$Cell_Index, SampleType = Tag$SampleType, Tech = rep("BD", dim(Tag)[1]),stringsAsFactors=FALSE)
    Tag <- cbind(Tag, SampleID = paste(Tag$Donor, Tag$SampleType, sep="_"), stringsAsFactors = FALSE)
    Tag$SampleID <- gsub(" ", "_", Tag$SampleID)
    rownames(Tag) <- Tag$Cell_Index
    object_list[[Donor_id]] <- CreateSeuratObject(counts = t(RSEC), project = Donor_id, min.cells = 3, min.features = 200)
    object_list[[Donor_id]] <- AddMetaData(object_list[[Donor_id]], metadata = Tag)
}


#for 10X
for(SampleType in c("Nevus", "Malignant Nevus", "Melanoma")){
    Object <- Read10X(data.dir = paste0(SampleType, "/outs/filtered_gene_bc_matrices/"))
    Object <- CreateSeuratObject(counts = Object, project = "D5", min.cells = 3, min.features = 200)

    Tag <- data.frame(Donor = rep("D5", dim(Object)[2]), Cell_Index = colnames(Object), SampleType = rep(SampleType, dim(Object)[2]), Tech = rep("10X", dim(Object)[2]),stringsAsFactors=FALSE)
    rownames(Tag) <- Tag$Cell_Index
    Tag <- cbind(Tag, SampleID = paste(Tag$Donor, Tag$SampleType, sep="_"), stringsAsFactors = FALSE)
    Tag$SampleID <- gsub(" ", "_", Tag$SampleID)
    object_list[[paste0("D5_", SampleType)]] <- AddMetaData(Object, metadata = Tag)
}


D1_D4 <- purrr::reduce(object_list[-grep("D5",names(object_list))], merge, do.normalize = FALSE)
D5 <- purrr::reduce(object_list[grep("D5",names(object_list))], merge, do.normalize = FALSE)
D1_D5 <- purrr::reduce(object_list, merge, do.normalize = FALSE)
saveRDS(D1_D4, file="/1_Cell_annotation/1_Create_Object/D1_D4.rds")
saveRDS(D5, file="/1_Cell_annotation/1_Create_Object/D5.rds")
saveRDS(D1_D5, file="/1_Cell_annotation/1_Create_Object/D1_D5.rds")

#QC
D1_D5[["percent.mt"]] <- PercentageFeatureSet(D1_D5, pattern = "^MT-")

house_keeping_gene <- openxlsx::read.xlsx("house_keeping_gene.xlsx", startRow = 2)
house_keeping_gene <- unique(house_keeping_gene[,1])
house_keeping_object <- D1_D5[house_keeping_gene, ]
D1_D5[["nCount_house_keeping"]] <- Matrix::colSums(as.matrix(D1_D5[house_keeping_gene, ]@assays[["RNA"]]@counts))
D1_D5[["nhouse_keeping"]] <- apply(as.matrix(D1_D5[house_keeping_gene, ]@assays[["RNA"]]@counts), 2, function(x){return(length(which(x > 0)))})
p <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "nCount_house_keeping", "nhouse_keeping"), ncol = 3, pt.size = 0, group.by = "SampleID")
ggsave(p, file="/1_Cell_annotation/2_filter_cell/VlnPlot_no_point.pdf", width=15, height=10)
plot1 <- FeatureScatter(D1_D5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(D1_D5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p <- CombinePlots(plots = list(plot1, plot2))
ggsave(p, file="/1_Cell_annotation/2_filter_cell/ScatterPlot.pdf", width=10, height=4)


meta_frame <- D1_D5@meta.data
meta_frame <- meta_frame[,c("nCount_RNA", "nFeature_RNA", "SampleID", "percent.mt", "nhouse_keeping")]
meta_frame <- data.frame(id = rownames(meta_frame), meta_frame, type = rep("single_cell", nrow(meta_frame)), stringsAsFactors=FALSE)

meta_frame <- reshape2::melt(meta_frame, id.vars=c("id", "SampleID", "type"), variable.name="feature", value.name="value")
library(tidyverse)
outlier1 <- as.data.frame(meta_frame[meta_frame$feature %in% c("nhouse_keeping"),] %>% group_by(SampleID, feature) %>% summarise(line = median(value)-3*mad(value)))
outlier2 <- as.data.frame(meta_frame[meta_frame$feature %in% c("nCount_RNA", "nFeature_RNA"),] %>% group_by(SampleID, feature) %>% summarise(line = median(value)+3*mad(value)))
outlier3 <- as.data.frame(meta_frame[meta_frame$feature %in% c("percent.mt"),] %>% group_by(SampleID, feature) %>% summarise(line = median(value)+3*mad(value)))
outlier3[outlier3$SampleID %in% c(paste("D", 1:4, rep(c("_Malignant_Nevus", "_Melanoma", "_Nevus", "_Normal"), each=4), sep="")), 3] <- 25
outlier3[outlier3$SampleID %in% c(paste("D", 5, c("_Malignant_Nevus", "_Melanoma", "_Nevus"), sep="")), 3] <- 10
outlier <- rbind(outlier1, outlier2, outlier3)

library(ggplot2)
library(ggthemes)
p <- ggplot(meta_frame, aes(x=SampleID, y=value)) + geom_violin(aes(fill=SampleID), trim = FALSE, position = position_dodge(0.9)) + 
  geom_boxplot(aes(fill=SampleID),width = 0.15,position = position_dodge(0.9))+
  scale_fill_manual(values = c(ggsci::pal_npg("nrc",alpha = 1)(8), ggsci::pal_nejm("default",alpha = 1)(8)))+ 
  facet_wrap(feature ~ SampleID, scales = "free", ncol = 13) + coord_cartesian(ylim = c(0, NA)) +
  geom_hline(data=outlier, aes(yintercept=line), linetype="dashed", color = "#353333") +
  theme_few() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, filename = "/1_Cell_annotation/2_filter_cell/qc/D1_D5/qc.pdf", width=18, height=18)

#Filter cells using nhouse_keeping, nCount_RNA, nFeature_RNA, and percent.mt.
D1_D5_list <- list()
for(i in 1:length(unique(outlier$SampleID))){
   this.sample <- subset(D1_D5, subset = SampleID == unique(outlier$SampleID)[i])
   this.sample <- subset(this.sample, subset = nhouse_keeping > outlier[which(outlier$SampleID == unique(outlier$SampleID)[i] & outlier$feature == "nhouse_keeping"), 3])
   this.sample <- subset(this.sample, subset = nCount_RNA < outlier[which(outlier$SampleID == unique(outlier$SampleID)[i] & outlier$feature == "nCount_RNA"), 3])
   this.sample <- subset(this.sample, subset = nFeature_RNA < outlier[which(outlier$SampleID == unique(outlier$SampleID)[i] & outlier$feature == "nFeature_RNA"), 3])
   this.sample <- subset(this.sample, subset = percent.mt < outlier[which(outlier$SampleID == unique(outlier$SampleID)[i] & outlier$feature == "percent.mt"), 3])
   D1_D5_list <- c(D1_D5_list, list(this.sample))
}
#merge sample after filter cell
D1_D5_filter <- purrr::reduce(D1_D5_list, merge)



#Remove doublets
library(DropletUtils)
tmp_dir <- "/1_Cell_annotation/2_filter_cell/qc/temp_file"
dir.create(tmp_dir)
# Note: When working with data from multiple samples, run Scrublet on each sample separately. 
object_10X <- D1_D5_filter@assays$RNA@counts
sample_list <- unique(D1_D5_filter@meta.data$SampleID)
for(i in 1:length(sample_list)){
	this_10X <- object_10X[, colnames(object_10X) %in% rownames(D1_D5_filter@meta.data[which(D1_D5_filter@meta.data$SampleID == sample_list[i]),])]
	barcodes <- colnames(this_10X)
	gene.id <- rownames(this_10X)
	gene.symbol <- rownames(this_10X)
	write_dir <- paste(tmp_dir, "/", sample_list[i], "/", sep = "")
	#print(write_dir)
	write10xCounts(path=write_dir, this_10X, version = "2", barcodes = barcodes, gene.id = gene.id, gene.symbol == gene.symbol)
}

for(i in 1:length(sample_list)){
	input_dir <- paste(tmp_dir, sample_list[i], sep = "/")
	out <- Scrublet_py(sample_id = sample_list[i], input_dir = input_dir, out_dir = "/1_Cell_annotation/2_filter_cell/qc/")
}
unlink(tmp_dir, recursive = TRUE)



doublet_threshold <- list("M1_870_BCG" = 0.25, "M2_870" = 0.25, "M3_CON" = 0.25)
tmp_dir <- "/1_Cell_annotation/2_filter_cell/temp_file"
dir.create(tmp_dir)

object_10X <- D1_D5_filter@assays$RNA@counts
sample_list <- unique(D1_D5_filter@meta.data$SampleID)
for(i in 1:length(sample_list)){
	this_10X <- object_10X[, colnames(object_10X) %in% rownames(D1_D5_filter@meta.data[which(D1_D5_filter@meta.data$SampleID == sample_list[i]),])]
	barcodes <- colnames(this_10X)
	gene.id <- rownames(this_10X)
	gene.symbol <- rownames(this_10X)
	write_dir <- paste(tmp_dir, "/", sample_list[i], "/", sep = "")
	#print(write_dir)
	write10xCounts(path=write_dir, this_10X,
	version = "2", barcodes = barcodes, gene.id = gene.id, gene.symbol == gene.symbol)
}

meta_data <- D1_D5_filter@meta.data
for(i in 1:length(sample_list)){
	input_dir <- paste(tmp_dir, sample_list[i], sep = "/")
	print(paste0("Calculating ", sample_list[i]))
	print(paste0("doublet_threshold: ",doublet_threshold))
	out <- Scrublet_py(sample_id = sample_list[i], input_dir = input_dir, out_dir = "/1_Cell_annotation/2_filter_cell/D1_D5/")
}

# Filter
D1_D5_filter <- subset(D1_D5_filter, subset = predicted_doublets == "False")
mt.genes <- rownames(D1_D5_filter)[grep("^MT-",rownames(D1_D5_filter))]
rp.genes <- rownames(D1_D5_filter)[grep("^RP[SL]",rownames(D1_D5_filter))]
D1_D5_filter[["percent.rp"]] <- PercentageFeatureSet(D1_D5_filter, pattern = "^RP[SL]")
D1_D5_filter <- subset(D1_D5_filter, subset = nFeature_RNA > 200)
D1_D5_filter <- subset(D1_D5_filter, subset = nCount_RNA < 20000)
D1_D5_filter <- subset(D1_D5_filter, subset = percent.mt < 10)
D1_D5_filter <- D1_D5_filter[!rownames(D1_D5_filter) %in% mt.genes,]
D1_D5_filter <- D1_D5_filter[!rownames(D1_D5_filter) %in% rp.genes,]

saveRDS(D1_D5_filter,file="/1_Cell_annotation/2_filter_cell/D1_D5.rds")


# Scale
library(dplyr)
library(Seurat)
library(ggplot2)
library(doMC)
registerDoMC(300)
D1_D5_filter <- NormalizeData(D1_D5_filter, normalization.method = "LogNormalize", scale.factor = 10000)
D1_D5_filter <- CellCycleScoring(object = D1_D5_filter, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
D1_D5_filter <- FindVariableFeatures(D1_D5_filter, selection.method = "vst", nfeatures = 2000)
D1_D5_filter <- ScaleData(D1_D5_filter, features = rownames(D1_D5_filter))
D1_D5_filter <- RunPCA(D1_D5_filter, features = VariableFeatures(object = D1_D5_filter))
D1_D5_filter <- FindNeighbors(D1_D5_filter, dims = 1:30)
D1_D5_filter <- RunUMAP(D1_D5_filter, dims = 1:30, label = T)
saveRDS(D1_D5_filter,file="/1_Cell_annotation/3_scale/D1_D5.rds")


# Batch Remove
library(harmony)
library(doMC)
registerDoMC(30)
D1_D5 <- RunHarmony(D1_D5_filter, c("SampleID","Tech"))
D1_D5 <- RunUMAP(D1_D5, reduction = "harmony", dims = 1:30)
saveRDS(D1_D5, file = "/1_Cell_annotation/4_harmony/D1_D5.rds")


