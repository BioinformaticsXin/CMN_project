library(CytoTRACE)
library(RColorBrewer)

CAF <- readRDS("/6_CAF/CAF.rds")
Idents(object = CAF) <- "cell.type"

IDs <- unique(CAF@meta.data$SampleID)

datasets <- list()

for(i in 1:length(IDs)){
    this.sample <- subset(CAF, subset = SampleID == IDs[i])
    this.sample <- as.matrix(this.sample[["RNA"]]@counts)
    datasets <- c(datasets, list(this.sample))
}

results <- iCytoTRACE(datasets, ncores = 200)
# results <- iCytoTRACE(as.matrix(CAF[["RNA"]]@counts), ncores = 200)
saveRDS(results, file="/6_CAF/4_cytotrace/cytotrace.rds")



CAF_pheno <- as.character(CAF@meta.data$cell.type)
names(CAF_pheno) <- rownames(CAF@meta.data)

results <- readRDS("/6_CAF/4_cytotrace/cytotrace.rds")
CytoTRACE_result <- data.frame(cluster = CAF@meta.data$cell.type, stem_score = results$CytoTRACE[rownames(CAF@meta.data)], stringsAsFactors = FALSE)
fac <- with(CytoTRACE_result, reorder(cluster, stem_score, median, order = TRUE))
CytoTRACE_result$cluster <- factor(CytoTRACE_result$cluster, levels = rev(levels(fac)))
p <- ggplot(CytoTRACE_result,aes(x = cluster, y = stem_score, fill = cluster, color = cluster)) +
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
    scale_fill_manual(values = c("metaCAF - FGFBP2"="#ce5d7b", "angiCAF - CCN5"="#ae7b82", "TactCAF - SLPI"="#4165a7", 
"angiCAF - TPD52"="#cf8392", "infCAF - PLAU"="#d7dfab", "angiCAF - COCH"="#659bb0", "TdiffCAF - ZFP36"="#765b9a",
"apoCAF - NAMPT"="#6c8940", "ecmCAF - NPTX2"="#2f4862", "ecmCAF - COL11A1"="#d8bb98", "antCAF - COL6A5"="#8f3165",
"angiCAF - APCDD1"="#ba9314", "antCAF - C3"="#33c8e6", "antCAF - CCL19"="#646567", "CAF - AHNAK"="#f78b07")) +
    scale_color_manual(values = c("metaCAF - FGFBP2"="#ce5d7b", "angiCAF - CCN5"="#ae7b82", "TactCAF - SLPI"="#4165a7", 
"angiCAF - TPD52"="#cf8392", "infCAF - PLAU"="#d7dfab", "angiCAF - COCH"="#659bb0", "TdiffCAF - ZFP36"="#765b9a",
"apoCAF - NAMPT"="#6c8940", "ecmCAF - NPTX2"="#2f4862", "ecmCAF - COL11A1"="#d8bb98", "antCAF - COL6A5"="#8f3165",
"angiCAF - APCDD1"="#ba9314", "antCAF - C3"="#33c8e6", "antCAF - CCL19"="#646567", "CAF - AHNAK"="#f78b07")) + ylab("Predicted ordering by CytoTRACE")
ggsave(p, file="/6_CAF/4_cytotrace/cytotrace_median.pdf", width=8, height=5)

