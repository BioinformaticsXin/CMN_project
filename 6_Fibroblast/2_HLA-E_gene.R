library(ggpubr)
library(ggplot2)
library(ggthemes)
library(Seurat)

CAF <- readRDS("/6_CAF/CAF.rds")

Idents(CAF) <- "SampleType"
CAF <- subset(CAF, subset = SampleType != "Normal")

all.genes <- rownames(CAF)
CAF <- ScaleData(CAF, features = all.genes)

this.exp <-reshape2::melt(CAF[["RNA"]]@scale.data[c("DCN", "HLA-E"),], id.vars=c("gene"),variable.name="cell",value.name="value")
colnames(this.exp) <- c("gene", "cell", "value")
this.meta <- data.frame(cell=rownames(CAF@meta.data), SampleType=CAF@meta.data$SampleType, CellType=CAF@meta.data$cell.type, stringsAsFactors = FALSE)
this.exp <- merge(this.exp, this.meta, by = "cell")
this.exp$SampleType <- factor(this.exp$SampleType, levels=c("Nevus", "Malignant Nevus", "Melanoma")) 

p <- ggplot(this.exp[this.exp$gene %in% "HLA-E", ],aes(x = SampleType, y = value, fill = SampleType)) +
# boxplot
geom_boxplot(width = .6,show.legend = F,
            position = position_dodge(0.9),
            color = 'grey20',alpha = 0.5,
            outlier.color = 'grey50') +
# theme
theme_classic(base_size = 16) +
theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
      legend.position = 'top') +
# color
scale_fill_manual(values = c('Malignant Nevus'="#7876B1FF",'Melanoma'='#6F99ADFF','Nevus'='#FFDC91FF'), name = '')+
stat_compare_means(label = "p.signif", comparisons=list(c("Nevus","Malignant Nevus"),c("Malignant Nevus","Melanoma"),c("Nevus","Melanoma")))+
facet_wrap(~CellType, ncol=16)
ggsave(p, file= "/6_CAF/2_HLA-E/HLA-E_CAF.pdf", width=26, height=4.3)
