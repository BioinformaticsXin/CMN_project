library(RColorBrewer)
library(Seurat)
library(dplyr)
library(Matrix)
library(monocle)

#Melanocyte
Melanocyte <- readRDS("/2_Melanocyte/Melanocyte_sub_merge.rds")
Idents(object = Melanocyte) <- "RNA_snn_res.0.3"

data <- as(as.matrix(Melanocyte@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = Melanocyte@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
cds <- newCellDataSet(data, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size());
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds)))
var.genes <- VariableFeatures(Melanocyte)[1:1500]
cds <- setOrderingFilter(cds, var.genes)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
write.csv(pData(cds), paste0(out_dir,"/pseudotime.csv"))
Melanocyte_cds.variable <- cds
saveRDS(Melanocyte_cds.variable, file="/2_Melanocyte/2_Trajectory/variable/Melanocyte_cds.variable.rds")

#State trajectory distribution map
plot1 <- plot_cell_trajectory(Melanocyte_cds.variable, color_by = "State", cell_size = 1, show_branch_points=FALSE)+theme(legend.position = "​right") 
##Cluster trajectory distribution map
plot2 <- plot_cell_trajectory(Melanocyte_cds.variable, color_by = "SampleType", cell_size = 1, show_branch_points=FALSE)+scale_colour_manual(values = brewer.pal(5,"Paired"))+theme(legend.position = "​right")
##Pseudotime trajectory map
##合并作图
plotc <- plot1|plot2
ggsave("/2_Melanocyte/2_Trajectory/variable/Melanocyte_Combination.pdf", plot = plotc, width = 10, height = 5)

#Calculate differential genes among the three branches
BEAM_res <- BEAM(Melanocyte_cds.variable, branch_point = 2, cores = 20)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
saveRDS(BEAM_res, file="/2_Melanocyte/2_Trajectory/variable/BEAM_res.rds")
