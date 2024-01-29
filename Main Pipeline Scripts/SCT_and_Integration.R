##############################
#                            #
#     SCT and Integration    #
# Peters - Single Cell Liver #
#     Erica DePasquale       #
#        30 Aug 2023         #
#                            #
##############################

# Load libraries
library(Seurat)
library(patchwork)
library(gdata)
library(sctransform, lib.loc="~")

############################
#                          #
# Create Integrated Object #
#                          #
############################

setwd("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR")

# Load objects
expt_list=c("ALP-0003-BX1",
            "ALP-0003-BX2",
            "ALP-0003-BX3",
            "ALP-0003-BX4",
            "ALP-0003-BX5",
            "ALP-00011-BX2",
            "ALP-00017",
            "ALP-00023-BX1",
            "ALP-00023-BX2",
            "ALP-00023-BX3",
            "ALP-00030-BX1",
            "ALP-00030-BX4",
            "ALP-00033-BX2",
            "ALP-00033-BX3",
            "ALP-00036",
            "ALP-00036-BX2",
            "ALP-00036-BX3",
            "ALP-00041-TXP",
            "ALP-00048",
            "ALP-00066",
            "ALP-00069",
            "ALP-00078",
            "ALP-00017-BX2",
            "ALP-00033-BX4",
            "ALP-00036-BX4",
            "ALP-00036-BX5",
            "ALP-00042-BX2",
            "ALP-00042-TXP",
            "ALP-00045-BX1",
            "ALP-00045-TXP")

ifnb.list=NULL
for(i in expt_list){
  print(i)
  x=readRDS(paste0("QC/", i, "/", i, "_filtered_matrices/SCTransform_v2/Seurat.rds"))
  ifnb.list=c(ifnb.list, x)
}
remove(x)

resolution=0.5
dir.create(paste0("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_", resolution, "_SCT_08.30.23"))
setwd(paste0("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_", resolution, "_SCT_08.30.23"))

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, normalization.method = "SCT")

# this command creates an 'integrated' data assay
immune.combined2 <- IntegrateData(anchorset = immune.anchors, k.weight=50, normalization.method = "SCT")

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined2) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined2 <- RunPCA(immune.combined2, npcs = 30, verbose = FALSE)
immune.combined2 <- RunUMAP(immune.combined2, reduction = "pca", dims = 1:30)
immune.combined2 <- RunTSNE(immune.combined2, reduction = "pca", dims = 1:30, check_duplicates = FALSE)
immune.combined2 <- FindNeighbors(immune.combined2, reduction = "pca", dims = 1:30)
immune.combined2 <- FindClusters(immune.combined2, resolution = as.numeric(resolution))


############################
#                          #
#        Make Plots        #
#                          #
############################

# Visualization
p1 <- DimPlot(immune.combined2, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(immune.combined2, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(immune.combined2, reduction = "tsne", group.by = "orig.ident")
p4 <- DimPlot(immune.combined2, reduction = "tsne", label = TRUE, repel = TRUE)

pdf(file = "UMAP_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
  p1 + p2
dev.off()

pdf(file = "UMAP_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
  DimPlot(immune.combined2, reduction = "umap", split.by = "orig.ident")
dev.off()

pdf(file = "TSNE_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
  p3 + p4
dev.off()

pdf(file = "TSNE_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
  DimPlot(immune.combined2, reduction = "tsne", split.by = "orig.ident")
dev.off()


##############################
#                            #
# Score Cells for Signatures #
#                            #
##############################

source("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/200116_FunctionsGeneral.R")

signatures <- read.xls("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/220420_signatures_top150.xlsx", header = T, sheet = 1)
cexsize <- round(max(c(0.3, 0.5-nrow(immune.combined2@assays[["SCT"]]@data)/40000)),2)

# Color by signature score
message("\nCalculating signature scores")
# Signatures
signatures <- lapply(as.list(signatures), function(y) y[!is.na(y)])
signatures <- lapply(signatures, intersect, rownames(immune.combined2@assays[["SCT"]]@data))
# Average gene expression for scoreSignature
CM.mean <- rowMeans(as.matrix(immune.combined2@assays[["SCT"]]@data))
signScore <- lapply(names(signatures), function(g) {
  message(g)
  scoreSignature(CM = as.matrix(immune.combined2@assays[["SCT"]]@data), signatures = signatures[[g]], CM.mean = CM.mean, verbose = T)
})
names(signScore) <- names(signatures)

# Plot all signatures - UMAP
pdf(file = "UMAP_Liver_30_Seurat.pdf", width = 10, height = 10)
par(mar=c(4, 4, 4, 4))

for (n in names(signScore)) {
  mycol <- colItay(signScore[[n]])
  plot(immune.combined2@reductions[["umap"]]@cell.embeddings, pch = 16, cex = cexsize, col = mycol, main = n)
}
dev.off()

# Plot all signatures - TSNE
pdf(file = "TSNE_Liver_30_Seurat.pdf", width = 10, height = 10)
par(mar=c(4, 4, 4, 4))

for (n in names(signScore)) {
  mycol <- colItay(signScore[[n]])
  plot(immune.combined2@reductions[["tsne"]]@cell.embeddings, pch = 16, cex = cexsize, col = mycol, main = n)
}
dev.off()


##################################
#                                #
# Save New Plots and Save Object #
#                                #
##################################

# Plot Seurat UMAP
pdf(file = "UMAP_Liver_30_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
  DimPlot(object = immune.combined2, reduction = "umap",label= TRUE)
dev.off()

# Plot Seurat TSNE
pdf(file = "TSNE_Liver_30_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
  DimPlot(object = immune.combined2, reduction = "tsne",label= TRUE)
dev.off()

saveRDS(immune.combined2, file = "Seurat_Liver_30.rds")

temp=as.data.frame(cbind(immune.combined2@meta.data[["orig.ident"]],immune.combined2@meta.data[[paste0("integrated_snn_res.", resolution)]]))
write.table(table(temp), "Seurat_Liver_30_table.txt", sep="\t", quote=F)