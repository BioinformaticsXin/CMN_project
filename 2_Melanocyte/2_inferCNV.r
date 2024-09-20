# inferCNV data preparation
library(Seurat)
project_dir <- "/HSCR/hengya_work/lixin/Project/Congenital_melanoma_nevus_single_cell"
Melanocyte <- readRDS("/2_Melanocyte/Melanocyte_sub_merge.rds")
Idents(Melanocyte) <- "Donor"

# gene
gene <- read.table("Homo_sapiens.GRCh38.99_gene_pos.txt",header=F,sep="\t",stringsAsFactors=F);
gene <- unique(gene,);
gene <- gene[gene[,1] %in% rownames(Melanocyte),]
write.table(gene, '/2_Melanocyte/1_Infer_CNV/inferCNV_input/gene.txt', sep = '\t', row.names = F, col.names = F, quote = F)

for(donorid in unique(Melanocyte$Donor)){
	subMelanocyte <- subset(Melanocyte, idents = donorid)
	# count
	counts <- as.matrix(GetAssayData(object = subMelanocyte, slot = "counts"))[gene[,1],]
	# annotation
	annotation <- data.frame(id=rownames(subMelanocyte@meta.data),SampleType=subMelanocyte@meta.data$SampleType,stringsAsFactors=F)
	rownames(annotation) <- annotation[,1]
	annotation <- annotation[colnames(counts),]
	# write
	write.table(annotation, paste0('/2_Melanocyte/1_Infer_CNV/inferCNV_input/', donorid, '_annotation.txt'), sep = '\t', row.names = F, col.names = F, quote = F)
	write.table(counts, paste0('/2_Melanocyte/1_Infer_CNV/inferCNV_input/', donorid, '_counts.txt'), sep = '\t', row.names = T, col.names = T, quote = F)
}


# run infercnv
# D1
library(infercnv);
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="/2_Melanocyte/1_Infer_CNV/inferCNV_input/D1_counts.txt", 
                                    annotations_file="/2_Melanocyte/1_Infer_CNV/inferCNV_input/D1_annotation.txt", delim="\t", 
                                    gene_order_file="/2_Melanocyte/1_Infer_CNV/inferCNV_input/gene.txt", 
                                    ref_group_names=c("Nevus","Malignant Nevus"))
infercnv_re = infercnv::run(infercnv_obj, cutoff=0.1, out_dir="/2_Melanocyte/1_Infer_CNV/D1", cluster_by_groups=F, k_obs_groups=2,denoise=F, HMM=F,output_format="pdf", num_threads=300)

# D2
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="/2_Melanocyte/1_Infer_CNV/inferCNV_input/D2_counts.txt", 
                                    annotations_file="/2_Melanocyte/1_Infer_CNV/inferCNV_input/D2_annotation.txt", delim="\t", 
                                    gene_order_file="/2_Melanocyte/1_Infer_CNV/inferCNV_input/gene.txt", 
                                    ref_group_names=c("Nevus","Malignant Nevus"))
infercnv_re = infercnv::run(infercnv_obj, cutoff=0.1, out_dir="/2_Melanocyte/1_Infer_CNV/D2", cluster_by_groups=F, k_obs_groups=2,denoise=F, HMM=F,output_format="pdf", num_threads=300)

# D5
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="/2_Melanocyte/1_Infer_CNV/inferCNV_input/D5_counts.txt", 
                                    annotations_file="/2_Melanocyte/1_Infer_CNV/inferCNV_input/D5_annotation.txt", delim="\t", 
                                    gene_order_file="/2_Melanocyte/1_Infer_CNV/inferCNV_input/gene.txt", 
                                    ref_group_names=c("Nevus","Malignant Nevus"))
infercnv_re = infercnv::run(infercnv_obj, cutoff=0.1, out_dir="/2_Melanocyte/1_Infer_CNV/D5_4", cluster_by_groups=F, k_obs_groups=4,denoise=F, HMM=F,output_format="pdf", num_threads=300)

