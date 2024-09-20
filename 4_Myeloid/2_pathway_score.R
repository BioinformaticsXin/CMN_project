library(Seurat)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(dplyr)
library(KEGG.db)
library(clusterProfiler)
library(org.Hs.eg.db)

# load data
Myeloid_cell <- readRDS("/4_Myeloid/Myeloid_cell_annotation.rds")

# hsa04668	TNF signaling pathway
# hsa04657	IL-17 signaling pathway
# hsa04064	NF-kappa B signaling pathway
# hsa04060	Cytokine-cytokine receptor interaction
# hsa04621	NOD-like receptor signaling pathway

pathway <- data.frame(id = c("hsa04668","hsa04657","hsa04064","hsa04060","hsa04621","hsa04145","hsa04370","hsa04666"), 
name = c("TNF signaling pathway","IL-17 signaling pathway","NF-kappa B signaling pathway","Cytokine-cytokine receptor interaction","NOD-like receptor signaling pathway","Phagosome","VEGF signaling pathway","Fc gamma R-mediated phagocytosis"), stringsAsFactors = FALSE)

# Obtain pathway genes and convert ENTREZID to SYMBOL
KEGG.pathway <- as.list(KEGG.db::KEGGPATHID2EXTID)
KEGG.pathway <- KEGG.pathway[pathway[,1]]
KEGG.pathway.symbol <- lapply(KEGG.pathway, function(x){
    ID_SYMBOL <- bitr(x, fromType="ENTREZID", toType="SYMBOL", OrgDb = "org.Hs.eg.db")
    return(ID_SYMBOL$SYMBOL)
})
names(KEGG.pathway.symbol) <- pathway[,2]

##############################################################################################################################################
#add Angiogenesis_signature, Phagocytosis_signature, from "A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells"
##############################################################################################################################################
Angiogenesis_signature <- c("CCND2","CCNE1","CD44","CXCR4","E2F3","EDN1","EZH2","FGF18","FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2","MMP9","NOTCH1","PDGFA","PTK2","SPP1","STC1","TNFAIP6","TYMP","VAV2","VCAN","VEGFA") 
Phagocytosis_signature <- c("MRC1","CD163","MERTK","C1QB")

KEGG.pathway.symbol[["Angiogenesis_signature"]] <- Angiogenesis_signature
KEGG.pathway.symbol[["Phagocytosis_signature"]] <- Phagocytosis_signature

# GSVA score
library(GSVA)
gsva_res <- gsva(expr=as.matrix(Myeloid_cell[["RNA"]]@scale.data), gset.idx.list=KEGG.pathway.symbol, 
method="gsva", kcdf="Gaussian", ssgsea.norm=TRUE, parallel.sz=50)

library(data.table)
gsva_res_frame <- melt(gsva_res)
Cell_frame <- data.frame(Cell_id = rownames(Myeloid_cell@meta.data), Cell_type = Myeloid_cell@meta.data$cell.type, SampleType = Myeloid_cell@meta.data$SampleType, Sample_ID = Myeloid_cell@meta.data$orig.ident, stringsAsFactors=FALSE)
colnames(gsva_res_frame) <- c("Pathway", "Cell_id", "score")
gsva_res_frame <- merge(gsva_res_frame, Cell_frame, by = "Cell_id")
gsva_res_frame <- gsva_res_frame[which(gsva_res_frame$SampleType != "Normal"),]
gsva_res_frame$SampleType <- factor(gsva_res_frame$SampleType, levels=c("Nevus", "Malignant Nevus", "Melanoma"))
saveRDS(gsva_res_frame, file="/4_Myeloid/2_pathway_score/gsva_res_frame.rds")


# Phagocytosis_signature score
library("ggplot2")
library("ggthemes")
library("ggpubr")
library("ggsci")
library(forcats)
###########################################
# Angiogenesis_signature score
fac <- with(gsva_res_frame[(gsva_res_frame$Pathway %in% "Angiogenesis_signature"),], reorder(Cell_type, score, median, order = TRUE))
gsva_res_frame$Cell_type <- factor(gsva_res_frame$Cell_type, levels = rev(levels(fac)))
p <- ggplot(gsva_res_frame[(gsva_res_frame$Pathway %in% "Angiogenesis_signature"),],aes(x = Cell_type, y = score, fill = Cell_type, color = Cell_type)) +
  # boxplot
    geom_boxplot(width = .8,show.legend = F,
                position = position_dodge(0.9),
                alpha = 0.5,
                outlier.color = 'grey50') +
    geom_point(position=position_jitterdodge(jitter.width = 0.5)) +
    # theme
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
          legend.position = 'top') +
    # color
    scale_fill_manual(values = c("CD209+ SPP1- TAMs"="#FCB8B9","CD209- SPP1+ TAMs"="#ECBA7D","CD1C DCs"="#BA94C7","MDSCs"="#A7BBD6","CD209- SPP1low TAMs"="#EEDFCA",
    "Langerhans Cells"="#F6F195")) +
    scale_color_manual(values = c("CD209+ SPP1- TAMs"="#FCB8B9","CD209- SPP1+ TAMs"="#ECBA7D","CD1C DCs"="#BA94C7","MDSCs"="#A7BBD6","CD209- SPP1low TAMs"="#EEDFCA",
    "Langerhans Cells"="#F6F195"))
ggsave(p, file= "/4_Myeloid/2_pathway_score/Angiogenesis_signature_score.pdf", width=5, height=6)

# Angiogenesis_signature score of myeloid clusters in different SampleTypes 
p <- ggplot(gsva_res_frame[(gsva_res_frame$Pathway %in% "Angiogenesis_signature") & (gsva_res_frame$Cell_type %in% c("CD209+ SPP1- TAMs", "CD209- SPP1+ TAMs", "CD209- SPP1low TAMs")),],aes(x = SampleType, y = score, fill = SampleType)) +
  # boxplot
    geom_boxplot(width = .5,show.legend = F,
                position = position_dodge(0.9),
                color = 'grey20',alpha = 0.5,
                outlier.color = 'grey50') +
    # theme
    theme_few(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
          legend.position = 'top') +
    # color
    scale_fill_manual(values = c('Melanoma'='#6F99ADFF','Malignant Nevus'='#7876B1FF','Nevus'='#FFDC91FF'), name = '')+
    facet_wrap(. ~ Cell_type, ncol=3, scales = "free")+
    stat_compare_means(label = "p.signif", comparisons=list(c("Nevus","Malignant Nevus"),c("Nevus","Melanoma"),c("Malignant Nevus","Melanoma")))
ggsave(p, file= "/4_Myeloid/2_pathway_score/Angiogenesis_signature.pdf", width=7, height=3.8)

###########################################
#Phagocytosis_signature score 
fac <- with(gsva_res_frame[(gsva_res_frame$Pathway %in% "Phagocytosis_signature"),], reorder(Cell_type, score, median, order = TRUE))
gsva_res_frame$Cell_type <- factor(gsva_res_frame$Cell_type, levels = rev(levels(fac)))
p <- ggplot(gsva_res_frame[(gsva_res_frame$Pathway %in% "Phagocytosis_signature"),],aes(x = Cell_type, y = score, fill = Cell_type, color = Cell_type)) +
  # boxplot
    geom_boxplot(width = .8,show.legend = F,
                position = position_dodge(0.9),
                alpha = 0.5,
                outlier.color = 'grey50') +
    geom_point(position=position_jitterdodge(jitter.width = 0.5)) +
    # theme
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
          legend.position = 'top') +
    # color
    scale_fill_manual(values = c("CD209+ SPP1- TAMs"="#FCB8B9","CD209- SPP1+ TAMs"="#ECBA7D","CD1C DCs"="#BA94C7","MDSCs"="#A7BBD6","CD209- SPP1low TAMs"="#EEDFCA",
    "Langerhans Cells"="#F6F195")) +
    scale_color_manual(values = c("CD209+ SPP1- TAMs"="#FCB8B9","CD209- SPP1+ TAMs"="#ECBA7D","CD1C DCs"="#BA94C7","MDSCs"="#A7BBD6","CD209- SPP1low TAMs"="#EEDFCA",
    "Langerhans Cells"="#F6F195"))
ggsave(p, file= "/4_Myeloid/2_pathway_score/Phagocytosis_signature_score.pdf", width=5, height=6)

#Phagocytosis_signature score of myeloid clusters in different SampleTypes 
p <- ggplot(gsva_res_frame[(gsva_res_frame$Pathway %in% "Phagocytosis_signature") & (gsva_res_frame$Cell_type %in% c("CD209+ SPP1- TAMs", "CD209- SPP1+ TAMs", "CD209- SPP1low TAMs")),],aes(x = SampleType, y = score, fill = SampleType)) +
  # boxplot
    geom_boxplot(width = .5,show.legend = F,
                position = position_dodge(0.9),
                color = 'grey20',alpha = 0.5,
                outlier.color = 'grey50') +
    # theme
    theme_few(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
          legend.position = 'top') +
    # color
    scale_fill_manual(values = c('Melanoma'='#6F99ADFF','Malignant Nevus'='#7876B1FF','Nevus'='#FFDC91FF'), name = '')+
    facet_wrap(. ~ Cell_type, ncol=3, scales = "free")+
    stat_compare_means(label = "p.signif", comparisons=list(c("Nevus","Malignant Nevus"),c("Nevus","Melanoma"),c("Malignant Nevus","Melanoma")))
ggsave(p, file= "/4_Myeloid/2_pathway_score/Phagocytosis_signature.pdf", width=7, height=3.8)



