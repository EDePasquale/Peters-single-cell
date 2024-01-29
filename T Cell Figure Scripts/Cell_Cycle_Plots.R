# Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(Seurat)

# Read in object
M <- readRDS("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")
M <- PrepSCTFindMarkers(M)
ET_markers <- FindMarkers(M, assay="SCT", ident.1 = "CD8 Effector T", ident.2="CD8 Effector Memory T")
EMT_markers <- FindMarkers(M, assay="SCT", ident.1 = "CD8 Effector Memory T", ident.2="CD8 Effector T")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
M <- CellCycleScoring(M, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Subset to T-cells only
Idents(object = M) <- "cluster_names_new"
Clusters_Keep=c("CD8 Effector Memory T", "CD4 naive T", "CD8 TRM Activated", "Gamma Delta T", "MAIT", "CD8 Effector T", "Cycling T")
M_T=subset(x = M, idents = Clusters_Keep)

# Make TSNE plot
M_T <- RunPCA(M_T, npcs = 30, verbose = FALSE)
M_T <- RunUMAP(M_T, reduction = "pca", dims = 1:30)
M_T <- RunTSNE(M_T, reduction = "pca", dims = 1:30)

P1<-FeaturePlot(M_T, reduction="tsne", features = "G2M.Score")
P2<-FeaturePlot(M_T, reduction="tsne", features = "S.Score")

pdf(file = "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/TSNE_Cell_Cycle.pdf", width = 5.5, height = 10)
par(mar=c(2, 2, 2, 2))
P1 + P2
dev.off()