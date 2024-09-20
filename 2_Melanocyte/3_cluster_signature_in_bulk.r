library(dplyr)
library(ggplot2)
library(ggthemes)


#Melanocyte
Melanocyte <- readRDS("/2_Melanocyte/Melanocyte_sub_merge.rds")


########### TCGA
load("/TCGA/TCGA-SKCM/data_exprs_FPKM.RData")
SKCM_exp <- data_exprs_FPKM
load("/TCGA/TCGA-SKCM/data_clinical.RData")
TCGA_SKCM_clinical <- data_clinical

gene_probe <- read.table("/TCGA/gencode.v22.annotation.gene.probeMap", header = T, sep = "\t", stringsAsFactors = F)
gene_probe[,1] <- gsub("\\..*", "", gene_probe[,1])

inter_genes <- intersect(rownames(SKCM_exp), gene_probe[,1])
gene_probe <- gene_probe[gene_probe[,1] %in% inter_genes,]
SKCM_exp <- SKCM_exp[rownames(SKCM_exp) %in% inter_genes,]
rownames(gene_probe) <- gene_probe[,1]
rownames(SKCM_exp) <- gene_probe[rownames(SKCM_exp),]$gene
colnames(SKCM_exp) <- substring(colnames(SKCM_exp), 1, 12)

sample_frame <- TCGA_SKCM_clinical[,c("bcr_patient_barcode", "submitted_tumor_location")]
sample_frame <- sample_frame[which(sample_frame[,2] != ""),]
sample_frame <- sample_frame[which(sample_frame[,2] != "Regional Cutaneous or Subcutaneous Tissue (includes satellite and in-transit metastasis)"),]
sample_frame[,2] <- as.character(sample_frame[,2])
colnames(sample_frame) <- c("sample_id", "sample_type")


##############Directly examine the scores of specifically expressed genes in each sample##############
#cluster1
cluster1 <- openxlsx::read.xlsx("/2_Melanocyte/Dif/Cluster/function/cluster1.xlsx")
cluster1 <- cluster1$SYMBOL
#cluster2
cluster2 <- openxlsx::read.xlsx("/2_Melanocyte/Dif/Cluster/function/cluster2.xlsx")
cluster2 <- cluster2$SYMBOL
#cluster3
cluster3 <- openxlsx::read.xlsx("/2_Melanocyte/Dif/Cluster/function/cluster3.xlsx")
cluster3 <- cluster3$SYMBOL
#cluster4
cluster4 <- openxlsx::read.xlsx("/2_Melanocyte/Dif/Cluster/function/cluster4.xlsx")
cluster4 <- cluster4$SYMBOL

#计算GSVA得分
library(GSVA)
gsva_res <- gsva(expr=as.matrix(SKCM_exp), gset.idx.list=list(cluster1=cluster1, cluster2=cluster2, cluster3=cluster3, cluster4=cluster4), 
                 method="gsva", kcdf="Gaussian", ssgsea.norm=TRUE, parallel.sz=50)

gsva_res <- reshape2::melt(gsva_res, id.vars=c("cluster"),variable.name="sample_id",value.name="score")
colnames(gsva_res)[1:2] <- c("cluster", "sample_id")
gsva_res <- merge(gsva_res, sample_frame, by = "sample_id")
data.frame(gsva_res %>% group_by(sample_type, cluster) %>% summarise(median = median(score)))

gsva_res$cluster <- factor(gsva_res$cluster, levels = c("cluster1", "cluster2", "cluster3", "cluster4"))
gsva_res$sample_type <- factor(gsva_res$sample_type, levels = c("Primary Tumor", "Distant Metastasis", "Regional Lymph Node"))
save(gsva_res, file="/2_Melanocyte/3_cluster_signature_in_bulk/TCGA/gsva_res.Rdata")



########### GSE98394
SKCM_exp <- read.table("/GSE98394/GSE98394_expression.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
gene_probe <- read.table("/TCGA/gencode.v22.annotation.gene.probeMap", header = T, sep = "\t", stringsAsFactors = F)
gene_probe[,1] <- gsub("\\..*", "", gene_probe[,1])

inter_genes <- intersect(rownames(SKCM_exp), gene_probe[,1])
gene_probe <- gene_probe[gene_probe[,1] %in% inter_genes,]
SKCM_exp <- SKCM_exp[rownames(SKCM_exp) %in% inter_genes,]
rownames(gene_probe) <- gene_probe[,1]
gene_probe <- gene_probe[!duplicated(gene_probe$gene),]
SKCM_exp <- SKCM_exp[rownames(gene_probe),]
rownames(SKCM_exp) <- gene_probe[rownames(SKCM_exp),]$gene
write.table(SKCM_exp, file = "/GSE98394/GSE98394_exp_matrix.txt", col.names = T, row.names = T, sep = "\t", quote = F)

sample_frame <- read.table("/GSE98394/GSE98394_series_matrix.txt", skip = 51, nrow = 102, quote = "\"", fill = TRUE)
sample_frame <- data.frame(sampleid = unlist(sample_frame[2,-1]), sample_type = unlist(sample_frame[1,-1]), sample_id = unlist(sample_frame[1,-1]), stringsAsFactors = F)
sample_frame[grepl("common acquired nevus",sample_frame[,2]), 2] <- "nevus"
sample_frame[grepl("primary melanoma",sample_frame[,2]), 2] <- "primary"

id <- function(x){
  new_id <- unlist(strsplit(x, split = "\\("))[2]
  new_id <- gsub("\\)", "", new_id)
  return(new_id)
}

sample_frame[,3] <- sapply(sample_frame[,3], id)
#########################################################
gsva_res <- gsva(expr=as.matrix(SKCM_exp), gset.idx.list=list(cluster1=cluster1, cluster2=cluster2, cluster3=cluster3, cluster4=cluster4), 
                 method="gsva", kcdf="Gaussian", ssgsea.norm=TRUE, parallel.sz=50)

gsva_res <- reshape2::melt(gsva_res, id.vars=c("cluster"),variable.name="sample_id",value.name="score")
colnames(gsva_res)[1:2] <- c("cluster", "sample_id")
gsva_res <- merge(gsva_res, sample_frame, by = "sample_id")
data.frame(gsva_res %>% group_by(sample_type, cluster) %>% summarise(median = median(score)))

gsva_res$cluster <- factor(gsva_res$cluster, levels = c("cluster1", "cluster2", "cluster3", "cluster4"))
gsva_res$sample_type <- factor(gsva_res$sample_type, levels = c("nevus", "primary"))
save(gsva_res, file="/2_Melanocyte/3_cluster_signature_in_bulk/GSE98394/gsva_res.Rdata")


########### GSE98394
GPL96_57554 <- read.table("/GSE46517/GPL96-57554.txt", header = T, sep = "\t", stringsAsFactors = F, quote = "\"")
GPL96_57554 <- GPL96_57554[,c("ID", "Gene.Symbol")]
GPL96_57554 <- GPL96_57554[!grepl("///", GPL96_57554$Gene.Symbol),]
GPL96_57554 <- GPL96_57554[which(GPL96_57554[,2] != ""),]
id_gene <- split(GPL96_57554[,1], GPL96_57554[,2])
GSE46517_series <- read.table("/GSE46517/GSE46517_series_matrix.txt", header = T, sep = "\t", stringsAsFactors = F, quote = "\"", comment.char = "!")
rownames(GSE46517_series) <- GSE46517_series[,1]
GSE46517_series <- GSE46517_series[,-1]

exp_matrix <- sapply(id_gene, function(x){
   if(length(x) == 1){
      return(as.numeric(GSE46517_series[x,]))
   }else{
      this.exp <- GSE46517_series[x,]
      this.exp <- apply(this.exp, 1, mean)
      return(this.exp)
   }
})
exp_matrix <- Reduce(rbind, exp_matrix)
colnames(exp_matrix) <- colnames(GSE46517_series)
rownames(exp_matrix) <- names(id_gene)

write.table(exp_matrix, file = "/GSE46517/GSE46517_exp_matrix.txt", col.names = T, row.names = T, sep = "\t", quote = F)

GSE46517_clinical <- read.table("/GSE46517/GSE46517_series_matrix.txt", skip = 72, nrow = 92, quote = "\"", fill = TRUE)
GSE46517_clinical <- GSE46517_clinical[which(GSE46517_clinical[,1] == "!Sample_description"),-1]

exp_matrix <- read.table("/GSE46517/GSE46517_exp_matrix.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
sample_frame <- data.frame(sample_id = colnames(exp_matrix), sample_type = unlist(GSE46517_clinical), stringsAsFactors = FALSE)
sample_frame[grepl("Met",sample_frame[,2]), 2] <- "metastatic"
sample_frame[grepl("Pri",sample_frame[,2]), 2] <- "primary"
sample_frame[grepl("Nev",sample_frame[,2]), 2] <- "nevi"
sample_frame[grepl("Cont",sample_frame[,2]), 2] <- "normal"

gsva_res <- gsva(expr=as.matrix(exp_matrix), gset.idx.list=list(cluster1=cluster1, cluster2=cluster2, cluster3=cluster3, cluster4=cluster4), 
                 method="gsva", kcdf="Gaussian", ssgsea.norm=TRUE, parallel.sz=50)

gsva_res <- reshape2::melt(gsva_res, id.vars=c("cluster"),variable.name="sample_id",value.name="score")
colnames(gsva_res)[1:2] <- c("cluster", "sample_id")
gsva_res <- merge(gsva_res, sample_frame, by = "sample_id")
data.frame(gsva_res %>% group_by(sample_type, cluster) %>% summarise(median = median(score)))

gsva_res$cluster <- factor(gsva_res$cluster, levels = c("cluster1", "cluster2", "cluster3", "cluster4"))
gsva_res$sample_type <- factor(gsva_res$sample_type, levels = c("normal", "nevi", "primary", "metastatic"))
save(gsva_res, file="/2_Melanocyte/3_cluster_signature_in_bulk/GSE46517/gsva_res.Rdata")
