####################
#                  #
# Explore Genes in #
#   Supp. Table 5  #
#     Feng 2019    #
#    12 Feb 2024   # 
# Erica DePasquale #
#                  #
####################

# Load libraries
#library(sctransform, lib.loc="~")
library(sctransform)
library(Seurat)
library(patchwork)
library(readxl)
library(ggplot2)
library(cowplot)

# Read in Seurat object
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/")
M <- readRDS("Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")
DefaultAssay(M)<-"RNA"

dir.create("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Feng2019_SuppTable5")
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Feng2019_SuppTable5")

###############
# Score cells #
###############
source("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Scripts_Used/Peters-single-cell/T Cell Figure Scripts/AccessoryFunctions.R")

signatures <- read_excel("/Volumes/GI-Informatics/DePasquale/Signature_Files/12Feb2024_FengSupp5.xlsx", sheet = 1)


####
# SCT data slot
####
cexsize <- round(max(c(0.3, 0.5-nrow(M@assays[["SCT"]]@data)/40000)),2)

# Color by signature score
message("\nCalculating signature scores")
# Signatures
signatures <- lapply(as.list(signatures), function(y) y[!is.na(y)])
signatures <- lapply(signatures, intersect, rownames(M@assays[["SCT"]]@data))
# Average gene expression for scoreSignature
CM.mean <- rowMeans(as.matrix(M@assays[["SCT"]]@data))
signScore <- lapply(names(signatures), function(g) {
  message(g)
  scoreSignature(CM = as.matrix(M@assays[["SCT"]]@data), signatures = signatures[[g]], CM.mean = CM.mean, verbose = T)
})
names(signScore) <- names(signatures)

# Plot all signatures - TSNE
pdf(file = paste0("TSNE_Liver_Feng2019_SCT.pdf"), width = 10, height = 10)
par(mar=c(4, 4, 4, 4))
  for (n in names(signScore)) {
    mycol <- colItay(signScore[[n]])
    plot(M@reductions[["tsne"]]@cell.embeddings, pch = 16, cex = cexsize, col = mycol, main = n)
  }
dev.off()

####
# RNA data slot
####
M <- NormalizeData(M, assay="RNA")
all.genes <- rownames(M)
M <- ScaleData(M, features = all.genes, assay="RNA")

cexsize <- round(max(c(0.3, 0.5-nrow(M@assays[["RNA"]]@data)/40000)),2)

# Color by signature score
message("\nCalculating signature scores")
# Signatures
signatures <- lapply(as.list(signatures), function(y) y[!is.na(y)])
signatures <- lapply(signatures, intersect, rownames(M@assays[["RNA"]]@data))
# Average gene expression for scoreSignature
CM.mean <- rowMeans(as.matrix(M@assays[["RNA"]]@data))
signScore <- lapply(names(signatures), function(g) {
  message(g)
  scoreSignature(CM = as.matrix(M@assays[["RNA"]]@data), signatures = signatures[[g]], CM.mean = CM.mean, verbose = T)
})
names(signScore) <- names(signatures)

# Plot all signatures - TSNE
pdf(file = paste0("TSNE_Liver_Feng2019_RNA.pdf"), width = 10, height = 10)
par(mar=c(4, 4, 4, 4))
for (n in names(signScore)) {
  mycol <- colItay(signScore[[n]])
  plot(M@reductions[["tsne"]]@cell.embeddings, pch = 16, cex = cexsize, col = mycol, main = n)
}
dev.off()


######################
# Average expression #
######################
results<-AverageExpression(M, return.seurat = T, features=signatures[[1]])
results2=as.data.frame(results@assays[["RNA"]]@layers[["data"]])
results3<-AverageExpression(M, return.seurat = F, features=signatures[[1]])
row.names(results2)=results3[["RNA"]]@Dimnames[[1]]
colnames(results2)=results3[["RNA"]]@Dimnames[[2]]

write.table(results2, "AverageExpression_Feng2019genes.txt", sep="\t", quote=F)