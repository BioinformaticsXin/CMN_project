#Melanocyte
library(RColorBrewer)
library(CytoTRACE)


Melanocyte <- readRDS("/2_Melanocyte/Melanocyte_sub_merge.rds")
Idents(object = Melanocyte) <- "RNA_snn_res.0.3"

#############################################################
IDs <- unique(Melanocyte@meta.data$SampleID)
datasets <- list()

for(i in 1:length(IDs)){
    this.sample <- subset(Melanocyte, subset = SampleID == IDs[i])
    this.sample <- as.matrix(this.sample[["RNA"]]@counts)
    datasets <- c(datasets, list(this.sample))
}

results <- iCytoTRACE(datasets, ncores = 200)
saveRDS(results, file="/2_Melanocyte/2_cytotrace/cytotrace.rds")

Melanocyte_pheno <- as.character(Melanocyte@meta.data$RNA_snn_res.0.3)
names(Melanocyte_pheno) <- rownames(Melanocyte@meta.data)

cell_color <- c("#20854EFF", "#E18727FF", "#BC3C29FF", "#0072B5FF")
plotCytoTRACE(results, phenotype = Melanocyte_pheno, colors = cell_color, gene = "PMEL", outputDir = "/2_Melanocyte/2_cytotrace/")

results <- readRDS("/2_Melanocyte/2_cytotrace/cytotrace.rds")
CytoTRACE_result <- data.frame(cluster = Melanocyte@meta.data$RNA_snn_res.0.3, stem_score = results$CytoTRACE[rownames(Melanocyte@meta.data)], stringsAsFactors = FALSE)
fac <- with(CytoTRACE_result, reorder(cluster, stem_score, mean, order = TRUE))
CytoTRACE_result$cluster <- factor(CytoTRACE_result$cluster, levels = rev(levels(fac)))
p <- ggplot(CytoTRACE_result,aes(x = cluster, y = stem_score, fill = cluster, color = cluster)) +
    geom_boxplot(width = .8,show.legend = F,
                position = position_dodge(0.9),
                alpha = 0.5,
                outlier.color = 'grey50') +
    geom_point(position=position_jitterdodge(jitter.width = 0.5)) +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
        legend.position = 'top') +
    scale_fill_manual(values = c("cluster1"="#BC3C29FF","cluster2"="#0072B5FF","cluster3"="#E18727FF","cluster4"="#20854EFF")) +
    scale_color_manual(values = c("cluster1"="#BC3C29FF","cluster2"="#0072B5FF","cluster3"="#E18727FF","cluster4"="#20854EFF"))
ggsave(p, file="/2_Melanocyte/2_cytotrace/cytotrace_mean.pdf", width=7, height=6)

