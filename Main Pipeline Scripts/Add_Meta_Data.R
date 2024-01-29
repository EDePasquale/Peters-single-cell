##############################
#                            #
#        Add Meta Data       #
# Peters - Single Cell Liver #
#     Erica DePasquale       #
#        31 Aug 2023         #
#                            #
##############################

# Load libraries
library(Seurat)
library(dplyr)
library(plyr)

# Read in new Seurat object
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23")
M<-readRDS("Seurat_Liver_30_subcluster_names.rds")

# Read in metadata file
metadata=read.table("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/sample_metadata.txt", sep="\t", header=T)
metadata$Diagnosis[which(metadata$Diagnosis=="")]<-NA # to fix blank space issues
metadata=metadata[metadata$Sample %in% unique(M@meta.data[["orig.ident"]]),]

# Add metadata to the Seurat object
M@meta.data[["Diagnosis"]] <- mapvalues(M@meta.data[["orig.ident"]], from=metadata$Sample, to=metadata$Diagnosis)
M@meta.data[["Condition"]] <- mapvalues(M@meta.data[["orig.ident"]], from=metadata$Sample, to=metadata$Condition)
M@meta.data[["ACR"]] <- mapvalues(M@meta.data[["orig.ident"]], from=metadata$Sample, to=metadata$ACR)
M@meta.data[["SteroidResistant"]] <- mapvalues(M@meta.data[["orig.ident"]], from=metadata$Sample, to=metadata$SteroidResistant)
M@meta.data[["SampleType"]] <- mapvalues(M@meta.data[["orig.ident"]], from=metadata$Sample, to=metadata$SampleType)
M@meta.data[["ACRType"]] <- mapvalues(M@meta.data[["orig.ident"]], from=metadata$Sample, to=metadata$Type)

# Plot UMAPs
pdf(file = "UMAP_metadata_Sample.pdf", width = 11, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "umap", pt.size=0.15, group.by="orig.ident", shuffle=T))
dev.off()

pdf(file = "UMAP_metadata_Diagnosis.pdf", width = 11, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "umap", pt.size=0.15, group.by="Diagnosis", shuffle=T))
dev.off()

pdf(file = "UMAP_metadata_Condition.pdf", width = 9, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "umap", pt.size=0.15, group.by="Condition", shuffle=T))
dev.off()

pdf(file = "UMAP_metadata_ACR.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "umap", pt.size=0.15, group.by="ACR", shuffle=T))
dev.off()

pdf(file = "UMAP_metadata_SteroidResistant.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "umap", pt.size=0.15, group.by="SteroidResistant", shuffle=T))
dev.off()

pdf(file = "UMAP_metadata_SampleType.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "umap", pt.size=0.15, group.by="SampleType", shuffle=T))
dev.off()

pdf(file = "UMAP_metadata_ACRType.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "umap", pt.size=0.15, group.by="ACRType", shuffle=T))
dev.off()

# Plot TSNEs
pdf(file = "TSNE_metadata_Sample.pdf", width = 11, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "tsne", pt.size=0.15, group.by="orig.ident", shuffle=T))
dev.off()

pdf(file = "TSNE_metadata_Diagnosis.pdf", width = 11, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "tsne", pt.size=0.15, group.by="Diagnosis", shuffle=T))
dev.off()

pdf(file = "TSNE_metadata_Condition.pdf", width = 9, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "tsne", pt.size=0.15, group.by="Condition", shuffle=T))
dev.off()

pdf(file = "TSNE_metadata_ACR.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "tsne", pt.size=0.15, group.by="ACR", shuffle=T))
dev.off()

pdf(file = "TSNE_metadata_SteroidResistant.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "tsne", pt.size=0.15, group.by="SteroidResistant", shuffle=T))
dev.off()

pdf(file = "TSNE_metadata_SampleType.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "tsne", pt.size=0.15, group.by="SampleType", shuffle=T))
dev.off()

pdf(file = "TSNE_metadata_ACRType.pdf", width = 9.5, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "tsne", pt.size=0.15, group.by="ACRType", shuffle=T))
dev.off()

# Save new object
saveRDS(M, file = "Seurat_Liver_30_subcluster_names_meta.rds")