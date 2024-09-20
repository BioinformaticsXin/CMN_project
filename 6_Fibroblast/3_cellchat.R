library(dplyr)
library(Seurat)
library(purrr)

CAF <- readRDS("/6_CAF/CAF.rds")
T_cell <- readRDS("/3_T_cell/t_cell_annotation.rds")  
T_cell <- subset(T_cell, idents = "NK")

### merge
CAF <- Seurat_merge(seurat_object_list <- list(CAF=CAF, T_cell=T_cell))
Idents(CAF) <- "cell.type"
### CellChat analysis

# load the required libraries
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(openxlsx)
options(stringsAsFactors = FALSE)

library(doMC)
registerDoMC(cores)


seurat_object = CAF
seurat_annotations = "cell.type"
CellChatDB = "CellChatDB.human"
cores = 20
PPIdata = "PPI.human"

data.input  <- seurat_object@assays$RNA@data
# create a dataframe consisting of the cell labels
command <- paste0("identity <- data.frame(group =seurat_object$",seurat_annotations,", row.names = names(seurat_object$",seurat_annotations,"))")
eval(parse(text = command))

# create an object of CellChat
cellchat <- createCellChat(data.input)

# add the information of  metadata to the CellChat object
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity

# import receptor database
command <- paste0("CellChatDB <- ", CellChatDB)
eval(parse(text = command))

if(!is.null(search_subsetDB)){
    # use a subset of CellChatDB for cell-cell communication analysis
    CellChatDB.use <- subsetDB(CellChatDB, search = search_subsetDB) # use Secreted Signaling
    # set the used database in the object
    cellchat@DB <- CellChatDB.use
}else{
    # use all CellChatDB for cell-cell communication analysis
    CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    # set the used database in the object
    cellchat@DB <- CellChatDB.use
}

# preprocessing
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


# Interaction inference

if(raw.use){
    cellchat <- computeCommunProb(cellchat, raw.use = raw.use)

}else{
    print("Using the precompiled Protein-protein-Interactions as a priori network information...")
    command <- paste0("cellchat <- projectData(cellchat, ", PPIdata, ")")
    eval(parse(text = command))
    cellchat <- computeCommunProb(cellchat, raw.use = raw.use)
}

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)


cellchat <- computeCommunProbPathway(cellchat)

# Computational integrated cellular communication network
cellchat <- aggregateNet(cellchat)

# Visualization of cellular communication networks
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle_plot <- function(x){
    if(x == "count"){
        netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    }
    if(x == "weight"){
        netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    }
}

out_dir = "/6_CAF/3_cellchat/"
pdf(paste0(out_dir, "/netVisual.pdf"), onefile=TRUE)
par(mfrow = c(1,2), xpd=TRUE)
#split.screen(c(1,2))
sapply(c("count","weight"), function(x) netVisual_circle_plot(x))
dev.off()

pdf(paste0(out_dir, "/netVisual_split.pdf"), onefile=TRUE)
mat <- cellchat@net$weight
par(mfrow = c(3,ceiling(length(unique(identity$group))/3)), xpd=TRUE)
sapply(1:nrow(mat), function(i){
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
})
dev.off()

# save data
saveRDS(cellchat, file="/6_CAF/3_cellchat/cellchat.rds")


# plot
pathways.show = "HLA-E"  # MHC-I
pdf("/6_CAF/3_cellchat/HLA-E.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = c(1,2,3,4,5),  color.use = c("#4cbbd2", "#e9949e", "#9fc6db", "#efb4cc", "#c4aca6", "#db6657", "#977dac", "#b1ccea", "#4876b4"))
dev.off()

