library(ggpubr)
library(ggplot2)
library(ggthemes)
library(Seurat)

# load data
T_cell <- readRDS("/3_T_cell/t_cell_annotation.rds")

############################################ calculations ###########################################################
### The exhaustion score was defined as the average expression of CXCL13, HAVCR2, PDCD1, TIGIT, LAG3, CTLA4, LAYN, RBPJ, VCAM1, GZMB, TOX and MYO7A, while effector memory score as the mean expression of effector- or memory-related genes of PRF1, IFNG, CCL4, HLA-DQA1, GZMK, GZMA, GZMH, CD44, DUSP2, KLRB1, KLRD1 and CTSW.
exhaustion_gene <- c("CXCL13", "HAVCR2", "PDCD1", "TIGIT", "LAG3", "CTLA4", "LAYN", "RBPJ", "VCAM1", "TOX", "MYO7A")
effector_memory_gene <- c("PRF1", "IFNG", "CCL4", "HLA-DQA1", "GZMK", "GZMA", "GZMH", "CD44", "DUSP2", "KLRB1", "KLRD1", "CTSW")
naive_gene <- c("CCR7", "TCF7", "LEF1", "SELL")
cytotoxicity_gene <- c("PRF1", "IFNG", "GNLY", "NKG7", "GZMA", "GZMH", "KLRK1", "KLRB1", "KLRD1", "CTSW", "CST7")

scaled_data <- T_cell[["RNA"]]@scale.data
x <- scaled_data[, 2]  # Extracts the second column
x <- unlist(x)  # Unlist turns the column into a simple vector
names(x) <- rownames(scaled_data)  # Assign gene names to the vector
naive_gene_expression <- x[naive_gene]  # Extract expression values for the naive genes
T_cell_naive <- mean(naive_gene_expression)  # Compute the mean expression for the naive genes

# The other scores of exhaustion_gene, effector_memory_gene and cytotoxicity_gene were calculated in the same way.
T_cell_score <- data.frame(SampleType = T_cell@meta.data$SampleType, CellType = T_cell@meta.data$cell.type, 
                           exhaustion_score = T_cell_exhaustion, 
                           naive_score = T_cell_naive, 
                           cytotoxicity_score = T_cell_cytotoxicity,
                           effector_memory_score = T_cell_effector_memory, 
                           stringsAsFactors = FALSE)
rownames(T_cell_score) <- rownames(T_cell@meta.data)
T_cell_score$SampleType <- factor(T_cell_score$SampleType, levels=c("Normal", "Nevus", "Malignant Nevus", "Melanoma")) 

############################################ plot ###########################################################
### T_cell_naive_score
p <- ggplot(T_cell_score[T_cell_score$CellType %in% c("effector GZMK+ CD8+ T", "effector GZMB+ CD8+ T"),],aes(x = SampleType, y = naive_score, fill = SampleType)) +
  geom_violin(position = position_dodge(0.9),alpha = 0.5,
              width = 1.2,trim = T,
              color = NA) +
  geom_boxplot(width = .2,show.legend = F,
               position = position_dodge(0.9),
               color = 'grey20',alpha = 0.5,
               outlier.color = 'grey50') +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
        legend.position = 'top') +
  scale_fill_manual(values = c('Malignant Nevus'="#7876B1FF",'Melanoma'='#6F99ADFF','Nevus'='#FFDC91FF'), name = '')+
  facet_wrap(. ~ CellType, ncol=5, scales = "free")+
  stat_compare_means(label = "p.signif", comparisons=list(c("Nevus","Malignant Nevus"),c("Nevus","Melanoma"),c("Malignant Nevus","Melanoma")))
ggsave(p, file= "/3_T_cell/2_score/CD8_naive_score.pdf", width=7, height=7)

### The plots of exhaustion and cytotoxicity scores among different sample types in effector GZMK+ CD8+ T cells are like the codes above.
### The expression of HLA-E in T cells among different sample types. 
this.exp <-reshape2::melt(T_cell[["RNA"]]@scale.data[c("HLA-E", "CD274"),], id.vars=c("gene"),variable.name="cell",value.name="value")
colnames(this.exp) <- c("gene", "cell", "value")
this.meta <- data.frame(cell=rownames(T_cell@meta.data), SampleType=T_cell@meta.data$SampleType, CellType=T_cell@meta.data$cell.type, stringsAsFactors = FALSE)
this.exp <- merge(this.exp, this.meta, by = "cell")
this.exp$SampleType <- factor(this.exp$SampleType, levels=c("Nevus", "Malignant Nevus", "Melanoma")) 


p <- ggplot(this.exp[this.exp$gene %in% "HLA-E" & !(this.exp$CellType %in% "NK"), ],aes(x = SampleType, y = value, fill = SampleType)) +
geom_boxplot(width = .6,show.legend = F,
            position = position_dodge(0.9),
            color = 'grey20',alpha = 0.5,
            outlier.color = 'grey50') +
theme_classic(base_size = 16) +
theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
      legend.position = 'top') +
scale_fill_manual(values = c('Malignant Nevus'="#7876B1FF",'Melanoma'='#6F99ADFF','Nevus'='#FFDC91FF'), name = '')+
stat_compare_means(label = "p.signif", comparisons=list(c("Nevus","Malignant Nevus"),c("Malignant Nevus","Melanoma"),c("Nevus","Melanoma")))+
facet_wrap(~CellType, ncol=5)
ggsave(p, file= "/3_T_cell/2_score/HLA-E_T.pdf", width=10, height=5)
