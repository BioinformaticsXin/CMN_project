library(ggpubr)
library(tidyverse)
library(viridis)
library(ggplot)
library(forcats)

CAF <- readRDS("/6_CAF/CAF.rds")
all.genes <- rownames(CAF)
CAF <- ScaleData(CAF, features = all.genes)

#caculate Angiogenesis score
library(qusage)
Angiogenesis <- read.gmt("/HALLMARK_ANGIOGENESIS.v7.5.1.gmt")

library(GSVA)
gsva_res <- gsva(expr=as.matrix(CAF[["RNA"]]@scale.data), gset.idx.list=Angiogenesis, 
method="gsva", kcdf="Gaussian", ssgsea.norm=TRUE, parallel.sz=50)
gsva_res <- data.frame(cell_id = colnames(gsva_res), Angiogenesis_score = gsva_res[1,], stringsAsFactors = FALSE)

#
CAF_frame <- data.frame(cell_id = rownames(CAF@meta.data), CAF@meta.data, stringsAsFactors = FALSE)

CAF_frame <- merge(CAF_frame, gsva_res, by = "cell_id")
saveRDS(CAF_frame, file = "/6_CAF/Angiogenesis_score.rds")

CAF_frame <- readRDS("/6_CAF/Angiogenesis_score.rds")


fac <- with(CAF_frame, reorder(cell.type, Angiogenesis_score, median, order = TRUE))
CAF_frame$cell.type <- factor(CAF_frame$cell.type, levels = rev(levels(fac)))

#######################compare Angiogenesis score between ACTA2+ MYH11- fibroblasts and ACTA2+ MYH11+ fibroblasts####################################
p <- ggplot(CAF_frame, aes(x = cell.type, y = Angiogenesis_score, fill = cell.type)) +
  # violin
  geom_violin(position = position_dodge(0.9),alpha = 0.5,
              width = 1, trim = T,
              color = NA) +
  # boxplot
  geom_boxplot(width = .2,show.legend = F,
               position = position_dodge(0.9),
               color = 'grey20',alpha = 0.5,
               outlier.color = 'grey50') +
  # theme
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
        legend.position = 'top') +
  # color
  scale_fill_manual(values = c("metaCAF - FGFBP2"="#ce5d7b", "angiCAF - CCN5"="#ae7b82", "TactCAF - SLPI"="#4165a7", 
"angiCAF - TPD52"="#cf8392", "infCAF - PLAU"="#d7dfab", "angiCAF - COCH"="#659bb0", "TdiffCAF - ZFP36"="#765b9a",
"apoCAF - NAMPT"="#6c8940", "ecmCAF - NPTX2"="#2f4862", "ecmCAF - COL11A1"="#d8bb98", "antCAF - COL6A5"="#8f3165",
"angiCAF - APCDD1"="#ba9314", "antCAF - C3"="#33c8e6", "antCAF - CCL19"="#646567", "CAF - AHNAK"="#f78b07"), name = '')
ggsave(p, file="/6_CAF/Angiogenesis_celltype.pdf", width=10, height=6)



#######################compare among sampletype####################################
CAF_frame <- CAF_frame[!(CAF_frame$SampleType %in% c("Normal")),]
CAF_frame$SampleType <- factor(CAF_frame$SampleType, levels=c("Nevus", "Malignant Nevus", "Melanoma"))
p <- ggplot(CAF_frame,aes(x = SampleType, y = Angiogenesis_score, fill = SampleType)) +
  # violin
  geom_violin(position = position_dodge(0.9),alpha = 0.5,
              width = 1, trim = T,
              color = NA) +
  # boxplot
  geom_boxplot(width = .2,show.legend = F,
               position = position_dodge(0.9),
               color = 'grey20',alpha = 0.5,
               outlier.color = 'grey50') +
  # theme
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
        legend.position = 'top') +
  # color
  scale_fill_manual(values = c("Malignant Nevus"="#7876B1", "Melanoma"="#6F99AD", "Nevus"="#FFDC91"), name = '')+
  stat_compare_means(label = "p.signif", comparisons=list(c("Nevus","Melanoma"), c("Malignant Nevus","Melanoma"))) +
  facet_wrap(~cell.type, scale = "free", ncol = 15)
ggsave(p, file="/6_CAF/Angiogenesis_sampletype.pdf", width=33, height=6)
