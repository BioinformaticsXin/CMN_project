library(Seurat)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(msigdbr)
library(fgsea)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(GO.db)
library(KEGG.db)
library(openxlsx)
library(dplyr)
library(pryr)

Melanocyte <- readRDS("/2_Melanocyte/Melanocyte_sub_merge.rds")
H_set <- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets <- H_set %>% split(x = .$gene_symbol, f = .$gs_name)

hallmarkAnno <- read.table("/Cancer_hallmark/hallmarkAnno.txt",header=F,sep="\t",stringsAsFactors=F)
hallmarks <- hallmarkAnno[,3]
names(hallmarks) <- hallmarkAnno[,1]
hallmarks <- sapply(hallmarks,function(x){x <- unlist(strsplit(x,split=",")); return(x)})
names(hallmarks) <- gsub(" ","_",names(hallmarks))

hallmarks <- sapply(hallmarks,function(x){
   id_symbol <- bitr(x, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb = "org.Hs.eg.db")
   return(id_symbol$SYMBOL)})

Hset_fgseaRes <- list()
hallmarks_fgseaRes <- list()
for(group in c("cluster1", "cluster2", "cluster3", "cluster4")){
	cmarkers <- FindMarkers(Melanocyte, ident.1 = group, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0)
	cmarkers$genes = rownames(cmarkers)
	cluster.genes <- cmarkers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC)
	ranks_cluster <- cluster.genes$avg_log2FC
	names(ranks_cluster) <- rownames(cluster.genes)

	Hset_fgseaRes[[group]] <- fgsea(pathways = fgsea_sets, stats = ranks_cluster, minSize=15, maxSize=2000, scoreType = "pos")
	hallmarks_fgseaRes[[group]] <- fgsea(pathways = hallmarks, stats = ranks_cluster, minSize=15, maxSize=2000, scoreType = "pos")
}

#merge
########### H set
Hset_result <- rbind(data.frame(Object=rep("cluster1", nrow(Hset_fgseaRes[["cluster1"]])), Hset_fgseaRes[["cluster1"]], stringsAsFactors = FALSE),
                     data.frame(Object=rep("cluster2", nrow(Hset_fgseaRes[["cluster2"]])), Hset_fgseaRes[["cluster2"]], stringsAsFactors = FALSE),
                     data.frame(Object=rep("cluster3", nrow(Hset_fgseaRes[["cluster3"]])), Hset_fgseaRes[["cluster3"]], stringsAsFactors = FALSE),
                     data.frame(Object=rep("cluster4", nrow(Hset_fgseaRes[["cluster4"]])), Hset_fgseaRes[["cluster4"]], stringsAsFactors = FALSE))
Hset_result <- cbind(Hset_result, sig=rep("sig", nrow(Hset_result)), stringsAsFactors=F)
Hset_result[which(Hset_result$padj >= 0.05), ]$sig <- "non-sig"
Hset_result[,2] <- gsub("HALLMARK_", "", Hset_result[,2])
openxlsx::write.xlsx(Hset_result, file = "/2_Melanocyte/Dif/hallmark/Hset_result.xlsx")


########### cancer hallmark
hallmark_result <- rbind(data.frame(Object=rep("cluster1", nrow(hallmarks_fgseaRes[["cluster1"]])), hallmarks_fgseaRes[["cluster1"]], stringsAsFactors = FALSE),
                         data.frame(Object=rep("cluster2", nrow(hallmarks_fgseaRes[["cluster2"]])), hallmarks_fgseaRes[["cluster2"]], stringsAsFactors = FALSE),
                         data.frame(Object=rep("cluster3", nrow(hallmarks_fgseaRes[["cluster3"]])), hallmarks_fgseaRes[["cluster3"]], stringsAsFactors = FALSE),
                         data.frame(Object=rep("cluster4", nrow(hallmarks_fgseaRes[["cluster4"]])), hallmarks_fgseaRes[["cluster4"]], stringsAsFactors = FALSE))
hallmark_result <- cbind(hallmark_result, sig=rep("sig", nrow(hallmark_result)), stringsAsFactors=F)
hallmark_result[which(hallmark_result$padj >= 0.05), ]$sig <- "non-sig"
openxlsx::write.xlsx(hallmark_result, file = "/2_Melanocyte/Dif/hallmark/cancer_hallmark_result.xlsx")
