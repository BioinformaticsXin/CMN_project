library(Seurat)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(dplyr)
library(Matrix)
library(monocle)

Melanocyte <- readRDS("/2_Melanocyte/Melanocyte_final.rds")
Malignant_gene <- openxlsx::read.xlsx("/2_Melanocyte/Dif/Malignant_level/Malignant/dif_gene.xlsx")
Malignant_potential_gene <- openxlsx::read.xlsx("/2_Melanocyte/Dif/Malignant_level/Malignant-potential/dif_gene.xlsx")
Normal_gene <- openxlsx::read.xlsx("/2_Melanocyte/Dif/Malignant_level/Normal/dif_gene.xlsx")

######################################Malignant_gene#############################################
Melanocyte_Malignant <- data.table::melt(Melanocyte[["RNA"]]@scale.data[Malignant_gene$SYMBOL,])
colnames(Melanocyte_Malignant) <- c("Gene", "Cell_index", "Value")
Melanocyte_Malignant$Gene <- factor(Melanocyte_Malignant$Gene, levels = Malignant_gene$SYMBOL)

Melanocyte_meta.data <- data.frame(Cell_index = rownames(Melanocyte@meta.data), Cell_cluster = Melanocyte@meta.data$RNA_snn_res.0.3, stringsAsFactors = FALSE)
Melanocyte_Malignant <- merge(Melanocyte_Malignant, Melanocyte_meta.data)

# Sort by Monocle2 pseudotime
for(i in 1:length(unique(Melanocyte_Malignant$Gene))){
    Melanocyte_Malignant$Cell_cluster <- factor(Melanocyte_Malignant$Cell_cluster, levels = c(7,1,0,5,3,4,2,6))

    p <- ggplot(Melanocyte_Malignant[(Melanocyte_Malignant$Gene %in% unique(Melanocyte_Malignant$Gene)[i]),],aes(x = Cell_cluster, y = Value, fill = Cell_cluster, color = Cell_cluster)) +
        geom_boxplot(width = .8,show.legend = F,
                    position = position_dodge(0.9),
                    alpha = 0.5,
                    outlier.color = 'grey50') +
        geom_point(position=position_jitterdodge(jitter.width = 0.5)) +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
            legend.position = 'top') +
        scale_fill_manual(values = c("0"=ggsci::pal_nejm("default",alpha = 1)(8)[1],"1"=ggsci::pal_nejm("default",alpha = 1)(8)[2],"2"=ggsci::pal_nejm("default",alpha = 1)(8)[3],
        "3"=ggsci::pal_nejm("default",alpha = 1)(8)[4], "4"=ggsci::pal_nejm("default",alpha = 1)(8)[5],"5"=ggsci::pal_nejm("default",alpha = 1)(8)[6],
        "6"=ggsci::pal_nejm("default",alpha = 1)(8)[7])) +
        scale_color_manual(values = c("0"=ggsci::pal_nejm("default",alpha = 1)(8)[1],"1"=ggsci::pal_nejm("default",alpha = 1)(8)[2],"2"=ggsci::pal_nejm("default",alpha = 1)(8)[3],
        "3"=ggsci::pal_nejm("default",alpha = 1)(8)[4], "4"=ggsci::pal_nejm("default",alpha = 1)(8)[5],"5"=ggsci::pal_nejm("default",alpha = 1)(8)[6],
        "6"=ggsci::pal_nejm("default",alpha = 1)(8)[7]))
    ggsave(p, file=paste0("/2_Melanocyte/dif_gene_along_Trajectory/Malignant/", unique(Melanocyte_Malignant$Gene)[i], ".pdf"), width=10, height=5)
}

######################################Malignant_potential_gene#############################################
Melanocyte_Malignant_potential <- data.table::melt(Melanocyte[["RNA"]]@scale.data[Malignant_potential_gene$SYMBOL,])
colnames(Melanocyte_Malignant_potential) <- c("Gene", "Cell_index", "Value")
Melanocyte_Malignant_potential$Gene <- factor(Melanocyte_Malignant_potential$Gene, levels = Malignant_potential_gene$SYMBOL)

Melanocyte_meta.data <- data.frame(Cell_index = rownames(Melanocyte@meta.data), Cell_cluster = Melanocyte@meta.data$RNA_snn_res.0.3, stringsAsFactors = FALSE)
Melanocyte_Malignant_potential <- merge(Melanocyte_Malignant_potential, Melanocyte_meta.data)

for(i in 1:length(unique(Melanocyte_Malignant_potential$Gene))){
    Melanocyte_Malignant_potential$Cell_cluster <- factor(Melanocyte_Malignant_potential$Cell_cluster, levels = c(7,1,0,5,3,4,2,6))

    p <- ggplot(Melanocyte_Malignant_potential[(Melanocyte_Malignant_potential$Gene %in% unique(Melanocyte_Malignant_potential$Gene)[i]),],aes(x = Cell_cluster, y = Value, fill = Cell_cluster, color = Cell_cluster)) +
        geom_boxplot(width = .8,show.legend = F,
                    position = position_dodge(0.9),
                    alpha = 0.5,
                    outlier.color = 'grey50') +
        geom_point(position=position_jitterdodge(jitter.width = 0.5)) +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
            legend.position = 'top') +
        scale_fill_manual(values = c("0"=ggsci::pal_nejm("default",alpha = 1)(8)[1],"1"=ggsci::pal_nejm("default",alpha = 1)(8)[2],"2"=ggsci::pal_nejm("default",alpha = 1)(8)[3],
        "3"=ggsci::pal_nejm("default",alpha = 1)(8)[4], "4"=ggsci::pal_nejm("default",alpha = 1)(8)[5],"5"=ggsci::pal_nejm("default",alpha = 1)(8)[6],
        "6"=ggsci::pal_nejm("default",alpha = 1)(8)[7])) +
        scale_color_manual(values = c("0"=ggsci::pal_nejm("default",alpha = 1)(8)[1],"1"=ggsci::pal_nejm("default",alpha = 1)(8)[2],"2"=ggsci::pal_nejm("default",alpha = 1)(8)[3],
        "3"=ggsci::pal_nejm("default",alpha = 1)(8)[4], "4"=ggsci::pal_nejm("default",alpha = 1)(8)[5],"5"=ggsci::pal_nejm("default",alpha = 1)(8)[6],
        "6"=ggsci::pal_nejm("default",alpha = 1)(8)[7]))
    ggsave(p, file=paste0("/2_Melanocyte/dif_gene_along_Trajectory/Malignant_potential/", unique(Melanocyte_Malignant_potential$Gene)[i], ".pdf"), width=10, height=5)
}

######################################Normal_gene#############################################
Melanocyte_Normal <- data.table::melt(Melanocyte[["RNA"]]@scale.data[Normal_gene$SYMBOL,])
colnames(Melanocyte_Normal) <- c("Gene", "Cell_index", "Value")
Melanocyte_Normal$Gene <- factor(Melanocyte_Normal$Gene, levels = Normal_gene$SYMBOL)

Melanocyte_meta.data <- data.frame(Cell_index = rownames(Melanocyte@meta.data), Cell_cluster = Melanocyte@meta.data$RNA_snn_res.0.3, stringsAsFactors = FALSE)
Melanocyte_Normal <- merge(Melanocyte_Normal, Melanocyte_meta.data)

for(i in 1:length(unique(Melanocyte_Normal$Gene))){
    Melanocyte_Normal$Cell_cluster <- factor(Melanocyte_Normal$Cell_cluster, levels = c(7,1,0,5,3,4,2,6))

    p <- ggplot(Melanocyte_Normal[(Melanocyte_Normal$Gene %in% unique(Melanocyte_Normal$Gene)[i]),],aes(x = Cell_cluster, y = Value, fill = Cell_cluster, color = Cell_cluster)) +
        geom_boxplot(width = .8,show.legend = F,
                    position = position_dodge(0.9),
                    alpha = 0.5,
                    outlier.color = 'grey50') +
        geom_point(position=position_jitterdodge(jitter.width = 0.5)) +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
            legend.position = 'top') +
        scale_fill_manual(values = c("0"=ggsci::pal_nejm("default",alpha = 1)(8)[1],"1"=ggsci::pal_nejm("default",alpha = 1)(8)[2],"2"=ggsci::pal_nejm("default",alpha = 1)(8)[3],
        "3"=ggsci::pal_nejm("default",alpha = 1)(8)[4], "4"=ggsci::pal_nejm("default",alpha = 1)(8)[5],"5"=ggsci::pal_nejm("default",alpha = 1)(8)[6],
        "6"=ggsci::pal_nejm("default",alpha = 1)(8)[7])) +
        scale_color_manual(values = c("0"=ggsci::pal_nejm("default",alpha = 1)(8)[1],"1"=ggsci::pal_nejm("default",alpha = 1)(8)[2],"2"=ggsci::pal_nejm("default",alpha = 1)(8)[3],
        "3"=ggsci::pal_nejm("default",alpha = 1)(8)[4], "4"=ggsci::pal_nejm("default",alpha = 1)(8)[5],"5"=ggsci::pal_nejm("default",alpha = 1)(8)[6],
        "6"=ggsci::pal_nejm("default",alpha = 1)(8)[7]))
    ggsave(p, file=paste0("/2_Melanocyte/dif_gene_along_Trajectory/Normal/", unique(Melanocyte_Normal$Gene)[i], ".pdf"), width=10, height=5)
}
