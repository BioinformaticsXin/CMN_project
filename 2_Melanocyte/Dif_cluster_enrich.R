#' @TODO 计算两个细胞群之间的差异基因，并使用差异基因进行富集分析
#' @param seurat_object 输入的seurat对象
#' @param clustera cluster1的名字
#' @param clusterb cluster2的名字
#' @param group.by 进行差异基因的分组信息
#' @param OrgDb 物种注释信息，可以为"org.Hs.eg.db"或"org.Mm.eg.db"
#' @param organism 功能富集时用的物种，现在支持"mmu"和"hsa"
#' @param cut_p 使用哪一列p值进行差异基因的筛选，可选择"p_val"和"p_val_adj"
#' @param out_dir 输出结果的路径
#' @returnType
#' @return 
#' 
#' author LX HLJ



Dif_cluster_enrich <- function(seurat_object = NULL, clustera = NULL,clusterb = NULL,group.by = NULL, out_dir = NULL, OrgDb = NULL, organism = NULL, cut_p = NULL, sig_gene_p = NULL, sig_gene_fd = NULL){
  library(Seurat)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(GO.db)
  library(KEGG.db)
  library(openxlsx)
  library(dplyr)
  library(pryr)

  #比较cluster0和cluster1的差异表达基因
  print("Calculating dif genes...")
  dge.cluster <- FindMarkers(seurat_object,ident.1 = clustera,ident.2 = clusterb, group.by = group.by, min.cells.group = 1, min.cells.feature = 1, min.pct = 0, logfc.threshold = 0, only.pos = FALSE)
  dge.cluster <- dge.cluster[rev(order(dge.cluster$avg_log2FC)),]
  dge.cluster <- data.frame(SYMBOL=rownames(dge.cluster), dge.cluster, stringsAsFactors=F)

  #symbol to entrezid 转换
  gene.tx <- bitr(rownames(dge.cluster), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb = OrgDb)
  dge.cluster <- left_join(dge.cluster, gene.tx, by="SYMBOL")

  #分为上下调基因，以及所有差异基因集合
  if(cut_p == "p_val"){
    sig_dge.cluster_up <- subset(dge.cluster, p_val < sig_gene_p & avg_log2FC > sig_gene_fd)
    sig_dge.cluster_down <- subset(dge.cluster, p_val < sig_gene_p & avg_log2FC < -sig_gene_fd)
    sig_dge.cluster_all <- subset(dge.cluster, p_val < sig_gene_p & abs(avg_log2FC) > sig_gene_fd)
  }
  if(cut_p == "p_val_adj"){
    sig_dge.cluster_up <- subset(dge.cluster, p_val_adj < sig_gene_p & avg_log2FC > sig_gene_fd)
    sig_dge.cluster_down <- subset(dge.cluster, p_val_adj < sig_gene_p & avg_log2FC < -sig_gene_fd)
    sig_dge.cluster_all <- subset(dge.cluster, p_val_adj < sig_gene_p & abs(avg_log2FC) > sig_gene_fd)
  }
  sig_dge.cluster <- list("up" = sig_dge.cluster_up, "down" = sig_dge.cluster_down)
  #输出差异基因结果
  openxlsx::write.xlsx(sig_dge.cluster, file = paste0(out_dir, "/dif_gene.xlsx"), overwrite = TRUE)

  #生成GSEA注释的基因文件
  dge_List <- dge.cluster$avg_log2FC;
  names(dge_List) <- dge.cluster$ENTREZID
  dge_List <- sort(dge_List ,decreasing = T)

  #当基因不足做功能富集时，使用enrich
  setClass("enrich",slots=list(result="data.frame"))

  #差异基因GO富集分析
  if(length(sig_dge.cluster_up$SYMBOL) > 1){
  print("GO enrichment analysis for up-regulated genes...")
  ego_up <- enrichGO(gene          = sig_dge.cluster_up$SYMBOL,
                      universe     = dge.cluster$SYMBOL,
                      OrgDb         = OrgDb,
                      keyType       = 'SYMBOL',
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05)
  }else{
    ego_up <- new("enrich",result=data.frame(ONTOLOGY=c(),ID=c(),Description=c(),GeneRatio=c(),BgRatio=c(),pvalue=c(),p.adjust=c(),qvalue=c(),geneID=c(),Count=c()))
  }

  if(length(sig_dge.cluster_down$SYMBOL) > 1){
  print("GO enrichment analysis for down-regulated genes...")
  ego_down <- enrichGO(gene         = sig_dge.cluster_down$SYMBOL,
                      universe     = dge.cluster$SYMBOL,
                      OrgDb         = OrgDb,
                      keyType       = 'SYMBOL',
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05)
  }else{
    ego_down <- new("enrich",result=data.frame(ONTOLOGY=c(),ID=c(),Description=c(),GeneRatio=c(),BgRatio=c(),pvalue=c(),p.adjust=c(),qvalue=c(),geneID=c(),Count=c()))
  }

  if(length(sig_dge.cluster_all$SYMBOL) > 1){
  print("GO enrichment analysis for all dif genes...")
  ego_all <- enrichGO(gene         = sig_dge.cluster_all$SYMBOL,
                      universe     = dge.cluster$SYMBOL,
                      OrgDb         = OrgDb,
                      keyType       = 'SYMBOL',
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05)
  }else{
    ego_all <- new("enrich",result=data.frame(ONTOLOGY=c(),ID=c(),Description=c(),GeneRatio=c(),BgRatio=c(),pvalue=c(),p.adjust=c(),qvalue=c(),geneID=c(),Count=c()))
  }

  
  #差异基因KEGG分析
  if(length(sig_dge.cluster_up$ENTREZID[!is.na(sig_dge.cluster_up$ENTREZID)]) > 1){
  print("KEGG enrichment analysis for up-regulated genes...")
  ekegg_up <- enrichKEGG(gene = sig_dge.cluster_up$ENTREZID, universe = dge.cluster$ENTREZID, organism = organism, pvalueCutoff = 0.05, use_internal_data = T)
  ekegg_up <- setReadable(ekegg_up, OrgDb = OrgDb, keyType="ENTREZID")
  }else{
    ekegg_up <- new("enrich",result=data.frame(ID=c(),Description=c(),GeneRatio=c(),BgRatio=c(),pvalue=c(),p.adjust=c(),qvalue=c(),geneID=c(),Count=c()))
  }

  if(length(sig_dge.cluster_down$ENTREZID[!is.na(sig_dge.cluster_down$ENTREZID)]) > 1){
  print("KEGG enrichment analysis for down-regulated genes...")
  ekegg_down <- enrichKEGG(gene = sig_dge.cluster_down$ENTREZID, universe = dge.cluster$ENTREZID, organism = organism, pvalueCutoff = 0.05, use_internal_data = T)
  ekegg_down <- setReadable(ekegg_down, OrgDb, keyType="ENTREZID")
  }else{
    ekegg_down <- new("enrich",result=data.frame(ID=c(),Description=c(),GeneRatio=c(),BgRatio=c(),pvalue=c(),p.adjust=c(),qvalue=c(),geneID=c(),Count=c()))
  }

  if(length(sig_dge.cluster_all$ENTREZID[!is.na(sig_dge.cluster_all$ENTREZID)]) > 1){
  print("KEGG enrichment analysis for all dif genes...")
  ekegg_all <- enrichKEGG(gene = sig_dge.cluster_all$ENTREZID, universe = dge.cluster$ENTREZID, organism = organism, pvalueCutoff = 0.05, use_internal_data = T)
  ekegg_all <- setReadable(ekegg_all, OrgDb, keyType="ENTREZID")

  #GSEA KEGG
  gsea_kegg <- gseKEGG(dge_List, organism=organism, keyType="ncbi-geneid",
                    minGSSize=10, maxGSSize=500, use_internal_data=T, eps=0,
                    pvalueCutoff=0.05, pAdjustMethod = "BH")
  }else{
    ekegg_all <- new("enrich",result=data.frame(ID=c(),Description=c(),GeneRatio=c(),BgRatio=c(),pvalue=c(),p.adjust=c(),qvalue=c(),geneID=c(),Count=c()))
    gsea_kegg <- new("enrich",result=data.frame(ID=c(),Description=c(),setSize=c(),enrichmentScore=c(),NES=c(),pvalue=c(),p.adjust=c(),qvalues=c(),rank=c(),leading_edge=c(),core_enrichment=c()))
  }
  
  sheet.list <- list("GO_up" = ego_up@result[which(ego_up@result$pvalue < 0.05),], "GO_down" = ego_down@result[which(ego_down@result$pvalue < 0.05),], 
  "GO_all" = ego_all@result[which(ego_all@result$pvalue < 0.05),], "KEGG_up" = ekegg_up@result[which(ekegg_up@result$pvalue < 0.05),], 
  "KEGG_down" = ekegg_down@result[which(ekegg_down@result$pvalue < 0.05),], "KEGG_all" = ekegg_all@result[which(ekegg_all@result$pvalue < 0.05),],
  "GSEA_KEGG" = gsea_kegg)
  openxlsx::write.xlsx(sheet.list, file = paste0(out_dir, "/enrich_result.xlsx"), overwrite = TRUE)
  return(sheet.list)
}