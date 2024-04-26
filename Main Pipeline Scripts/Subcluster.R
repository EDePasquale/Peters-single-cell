##############################
#                            #
#         Subcluster         #
# Peters - Single Cell Liver #
#     Erica DePasquale       #
#        31 Aug 2023         #
#                            #
##############################

# Load libraries
library(Seurat)
library(gdata)
library(plyr)
library(ggplot2)

#################
#               #
# Subcluster 12 # # endothelial cells
#               #
#################

M = readRDS("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Seurat_Liver_30.rds")

# Split cluster 12
M_12=subset(x = M, idents = c(12))

# Subcluster
M_12 <- RunPCA(M_12, npcs = 30, verbose = FALSE)
M_12 <- RunUMAP(M_12, reduction = "pca", dims = 1:30)
M_12 <- RunTSNE(M_12, reduction = "pca", dims = 1:30)
M_12 <- FindNeighbors(M_12, reduction = "pca", dims = 1:30)
M_12 <- FindClusters(M_12, resolution = 0.5)

dir.create("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_12")
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_12")

# Visualization
p1 <- DimPlot(M_12, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(M_12, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(M_12, reduction = "tsne", group.by = "orig.ident")
p4 <- DimPlot(M_12, reduction = "tsne", label = TRUE, repel = TRUE)

pdf(file = "UMAP_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p1 + p2
dev.off()

pdf(file = "UMAP_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_12, reduction = "umap", split.by = "orig.ident")
dev.off()

pdf(file = "TSNE_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p3 + p4
dev.off()

pdf(file = "TSNE_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_12, reduction = "tsne", split.by = "orig.ident")
dev.off()

# Plot Seurat UMAP
pdf(file = "UMAP_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_12, reduction = "umap",label= TRUE)
dev.off()

# Plot Seurat TSNE
pdf(file = "TSNE_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_12, reduction = "tsne",label= TRUE)
dev.off()

# Add new split clusters to full Seurat object
temp=cbind(as.numeric(as.character(M@active.ident)),as.numeric(as.character(M@active.ident)))
row.names(temp)=M@assays[["SCT"]]@data@Dimnames[[2]]
new_clusters=mapvalues(as.numeric(as.character(M_12@active.ident)), from=c(0,1,2), to=c(12,16,17))
temp[names(M_12@active.ident),2]<-new_clusters
M@meta.data[["subcluster_12"]] <- temp[,2]


################
#              #
# Subcluster 5 # # unknown T, kupffer 5, CD8 TRM activated
#              #
################

# Split cluster 5
M_5=subset(x = M, idents = c(5))

# Subcluster
M_5 <- RunPCA(M_5, npcs = 30, verbose = FALSE)
M_5 <- RunUMAP(M_5, reduction = "pca", dims = 1:30)
M_5 <- RunTSNE(M_5, reduction = "pca", dims = 1:30)
M_5 <- FindNeighbors(M_5, reduction = "pca", dims = 1:30)
M_5 <- FindClusters(M_5, resolution = 0.2)

dir.create("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_5")
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_5")

# Visualization
p1 <- DimPlot(M_5, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(M_5, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(M_5, reduction = "tsne", group.by = "orig.ident")
p4 <- DimPlot(M_5, reduction = "tsne", label = TRUE, repel = TRUE)

pdf(file = "UMAP_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p1 + p2
dev.off()

pdf(file = "UMAP_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_5, reduction = "umap", split.by = "orig.ident")
dev.off()

pdf(file = "TSNE_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p3 + p4
dev.off()

pdf(file = "TSNE_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_5, reduction = "tsne", split.by = "orig.ident")
dev.off()

# Plot Seurat UMAP
pdf(file = "UMAP_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_5, reduction = "umap",label= TRUE)
dev.off()

# Plot Seurat TSNE
pdf(file = "TSNE_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_5, reduction = "tsne",label= TRUE)
dev.off()

# Add new split clusters to full Seurat object
temp=cbind(as.numeric(as.character(M@meta.data[["subcluster_12"]])),as.numeric(as.character(M@meta.data[["subcluster_12"]])))
row.names(temp)=M@assays[["SCT"]]@data@Dimnames[[2]]
new_clusters=mapvalues(as.numeric(as.character(M_5@active.ident)), from=c(0,1,2), to=c(100,5,18))
temp[names(M_5@active.ident),2]<-new_clusters
M@meta.data[["subcluster_5"]] <- temp[,2]


################
#              #
# Subcluster 3 # # B cells
#              #
################

# Split cluster 3
M_3=subset(x = M, idents = c(3))

# Subcluster
M_3 <- RunPCA(M_3, npcs = 30, verbose = FALSE)
M_3 <- RunUMAP(M_3, reduction = "pca", dims = 1:30)
M_3 <- RunTSNE(M_3, reduction = "pca", dims = 1:30)
M_3 <- FindNeighbors(M_3, reduction = "pca", dims = 1:30)
M_3 <- FindClusters(M_3, resolution = 0.2)

dir.create("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_3")
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_3")

# Visualization
p1 <- DimPlot(M_3, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(M_3, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(M_3, reduction = "tsne", group.by = "orig.ident")
p4 <- DimPlot(M_3, reduction = "tsne", label = TRUE, repel = TRUE)

pdf(file = "UMAP_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p1 + p2
dev.off()

pdf(file = "UMAP_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_3, reduction = "umap", split.by = "orig.ident")
dev.off()

pdf(file = "TSNE_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p3 + p4
dev.off()

pdf(file = "TSNE_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_3, reduction = "tsne", split.by = "orig.ident")
dev.off()

# Plot Seurat UMAP
pdf(file = "UMAP_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_3, reduction = "umap",label= TRUE)
dev.off()

# Plot Seurat TSNE
pdf(file = "TSNE_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_3, reduction = "tsne",label= TRUE)
dev.off()

# Add new split clusters to full Seurat object
temp=cbind(as.numeric(as.character(M@meta.data[["subcluster_5"]])),as.numeric(as.character(M@meta.data[["subcluster_5"]])))
row.names(temp)=M@assays[["SCT"]]@data@Dimnames[[2]]
new_clusters=mapvalues(as.numeric(as.character(M_3@active.ident)), from=c(0,1,2), to=c(3,3,20))
temp[names(M_3@active.ident),2]<-new_clusters
M@meta.data[["subcluster_3"]] <- temp[,2]


################
#              #
# Subcluster 6 # # MAIT and CD4 naive T
#              #
################

# Split cluster 6
M_6=subset(x = M, idents = c(6))

# Subcluster
M_6 <- RunPCA(M_6, npcs = 30, verbose = FALSE)
M_6 <- RunUMAP(M_6, reduction = "pca", dims = 1:30)
M_6 <- RunTSNE(M_6, reduction = "pca", dims = 1:30)
M_6 <- FindNeighbors(M_6, reduction = "pca", dims = 1:30)
M_6 <- FindClusters(M_6, resolution = 0.2)

dir.create("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_6")
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_6")

# Visualization
p1 <- DimPlot(M_6, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(M_6, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(M_6, reduction = "tsne", group.by = "orig.ident")
p4 <- DimPlot(M_6, reduction = "tsne", label = TRUE, repel = TRUE)

pdf(file = "UMAP_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p1 + p2
dev.off()

pdf(file = "UMAP_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_6, reduction = "umap", split.by = "orig.ident")
dev.off()

pdf(file = "TSNE_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p3 + p4
dev.off()

pdf(file = "TSNE_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_6, reduction = "tsne", split.by = "orig.ident")
dev.off()

# Plot Seurat UMAP
pdf(file = "UMAP_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_6, reduction = "umap",label= TRUE)
dev.off()

# Plot Seurat TSNE
pdf(file = "TSNE_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_6, reduction = "tsne",label= TRUE)
dev.off()

# Add new split clusters to full Seurat object
temp=cbind(as.numeric(as.character(M@meta.data[["subcluster_3"]])),as.numeric(as.character(M@meta.data[["subcluster_3"]])))
row.names(temp)=M@assays[["SCT"]]@data@Dimnames[[2]]
new_clusters=mapvalues(as.numeric(as.character(M_6@active.ident)), from=c(0,1,2), to=c(6,21,21))
temp[names(M_6@active.ident),2]<-new_clusters
M@meta.data[["subcluster_6"]] <- temp[,2]


##################
#                #
# Subcluster 7/8 # # NK1 and CD8Effector T
#                #
##################

# Split cluster 7_8
M_7_8=subset(x = M, idents = c(7,8))

# Subcluster
M_7_8 <- RunPCA(M_7_8, npcs = 30, verbose = FALSE)
M_7_8 <- RunUMAP(M_7_8, reduction = "pca", dims = 1:30)
M_7_8 <- RunTSNE(M_7_8, reduction = "pca", dims = 1:30)
M_7_8 <- FindNeighbors(M_7_8, reduction = "pca", dims = 1:30)
M_7_8 <- FindClusters(M_7_8, resolution = 0.3)

dir.create("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_7_8")
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_7_8")

# Visualization
p1 <- DimPlot(M_7_8, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(M_7_8, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(M_7_8, reduction = "tsne", group.by = "orig.ident")
p4 <- DimPlot(M_7_8, reduction = "tsne", label = TRUE, repel = TRUE)

pdf(file = "UMAP_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p1 + p2
dev.off()

pdf(file = "UMAP_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_7_8, reduction = "umap", split.by = "orig.ident")
dev.off()

pdf(file = "TSNE_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p3 + p4
dev.off()

pdf(file = "TSNE_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_7_8, reduction = "tsne", split.by = "orig.ident")
dev.off()

# Plot Seurat UMAP
pdf(file = "UMAP_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_7_8, reduction = "umap",label= TRUE)
dev.off()

# Plot Seurat TSNE
pdf(file = "TSNE_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_7_8, reduction = "tsne",label= TRUE)
dev.off()

# Add new split clusters to full Seurat object
temp=cbind(as.numeric(as.character(M@meta.data[["subcluster_6"]])),as.numeric(as.character(M@meta.data[["subcluster_6"]])))
row.names(temp)=M@assays[["SCT"]]@data@Dimnames[[2]]
new_clusters=mapvalues(as.numeric(as.character(M_7_8@active.ident)), from=c(0,1,2), to=c(7,7,8))
temp[names(M_7_8@active.ident),2]<-new_clusters
M@meta.data[["subcluster_7_8"]] <- temp[,2]


######################
#                    #
# Subcluster Kupffer # # kupffer clusters
#                    #
######################

# Split clusters 2,4,9,15,18
M<-SetIdent(M, value="subcluster_7_8")
M_K=subset(x = M, idents = c(2,4,9,15,18))

# Subcluster
M_K <- RunPCA(M_K, npcs = 30, verbose = FALSE)
M_K <- RunUMAP(M_K, reduction = "pca", dims = 1:30)
M_K <- RunTSNE(M_K, reduction = "pca", dims = 1:30)
M_K <- FindNeighbors(M_K, reduction = "pca", dims = 1:30)
M_K <- FindClusters(M_K, resolution = 0.2)

dir.create("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_K")
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_K")

# Visualization
p1 <- DimPlot(M_K, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(M_K, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(M_K, reduction = "tsne", group.by = "orig.ident")
p4 <- DimPlot(M_K, reduction = "tsne", label = TRUE, repel = TRUE)

pdf(file = "UMAP_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p1 + p2
dev.off()

pdf(file = "UMAP_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_K, reduction = "umap", split.by = "orig.ident")
dev.off()

pdf(file = "TSNE_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p3 + p4
dev.off()

pdf(file = "TSNE_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_K, reduction = "tsne", split.by = "orig.ident")
dev.off()

# Plot Seurat UMAP
pdf(file = "UMAP_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_K, reduction = "umap",label= TRUE)
dev.off()

# Plot Seurat TSNE
pdf(file = "TSNE_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_K, reduction = "tsne",label= TRUE)
dev.off()

# Add new split clusters to full Seurat object
temp=cbind(as.numeric(as.character(M@meta.data[["subcluster_7_8"]])),as.numeric(as.character(M@meta.data[["subcluster_7_8"]])))
row.names(temp)=M@assays[["SCT"]]@data@Dimnames[[2]]
new_clusters=mapvalues(as.numeric(as.character(M_K@active.ident)), from=c(0,1,2,3,4,5,6,7), to=c(2,4,9,9,15,18,19,9))
temp[names(M_K@active.ident),2]<-new_clusters
M@meta.data[["subcluster_K"]] <- temp[,2]


##################
#                #
# Subcluster 9/2 # # more kupffer splits (refine boundaries)
#                #
##################

# Split cluster 9_2
M<-SetIdent(M, value="subcluster_K")
M_9_2=subset(x = M, idents = c(9,2))

# Subcluster
M_9_2 <- RunPCA(M_9_2, npcs = 30, verbose = FALSE)
M_9_2 <- RunUMAP(M_9_2, reduction = "pca", dims = 1:30)
M_9_2 <- RunTSNE(M_9_2, reduction = "pca", dims = 1:30)
M_9_2 <- FindNeighbors(M_9_2, reduction = "pca", dims = 1:30)
M_9_2 <- FindClusters(M_9_2, resolution = 0.5)

dir.create("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_9_2")
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_9_2")

# Visualization
p1 <- DimPlot(M_9_2, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(M_9_2, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(M_9_2, reduction = "tsne", group.by = "orig.ident")
p4 <- DimPlot(M_9_2, reduction = "tsne", label = TRUE, repel = TRUE)

pdf(file = "UMAP_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p1 + p2
dev.off()

pdf(file = "UMAP_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_9_2, reduction = "umap", split.by = "orig.ident")
dev.off()

pdf(file = "TSNE_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p3 + p4
dev.off()

pdf(file = "TSNE_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_9_2, reduction = "tsne", split.by = "orig.ident")
dev.off()

# Plot Seurat UMAP
pdf(file = "UMAP_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_9_2, reduction = "umap",label= TRUE)
dev.off()

# Plot Seurat TSNE
pdf(file = "TSNE_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_9_2, reduction = "tsne",label= TRUE)
dev.off()

# Add new split clusters to full Seurat object
temp=cbind(as.numeric(as.character(M@meta.data[["subcluster_K"]])),as.numeric(as.character(M@meta.data[["subcluster_K"]])))
row.names(temp)=M@assays[["SCT"]]@data@Dimnames[[2]]
new_clusters=mapvalues(as.numeric(as.character(M_9_2@active.ident)), from=c(0,1,2,3,4,5,6), to=c(2,9,9,2,9,2,2))
temp[names(M_9_2@active.ident),2]<-new_clusters
M@meta.data[["subcluster_9_2"]] <- temp[,2]
DimPlot(M, reduction="tsne", group.by = "subcluster_9_2", label=T)


##################
#                #
# Subcluster CD4 # # split CD4 into naive and activated (and remove 6 contaminating B cells)
#                #
##################

# Split cluster 21
M<-SetIdent(M, value="subcluster_9_2")
M_21=subset(x = M, idents = c(21))

# Subcluster
M_21=RunPCA(M_21, npcs = 30, verbose = FALSE)
M_21=RunUMAP(M_21, reduction = "pca", dims = 1:30)
M_21=RunTSNE(M_21, reduction = "pca", dims = 1:30)
M_21=FindNeighbors(M_21, reduction = "pca", dims = 1:30)
M_21=FindClusters(M_21, resolution = 0.5)

DefaultAssay(M_21)<-"RNA"

X=M_21[which(row.names(M_21)=="MS4A1"),]
to_remove=names(X@assays$RNA@data[1,][which(X@assays$RNA@data[1,] > 2)]) # remove 6 cells with "high" MS4A1
to_keep=setdiff(Cells(M_21), to_remove)
M_21_sub=subset(M_21, cells=to_keep)

DefaultAssay(M_21)<-"integrated"
DefaultAssay(M_21_sub)<-"integrated"

M_21_sub=FindNeighbors(M_21_sub, reduction = "pca", dims = 1:30)
M_21_sub=FindClusters(M_21_sub, resolution = 0.5)

dir.create("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_21")
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Subcluster_21")

# Visualization
p1 <- DimPlot(M_21_sub, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(M_21_sub, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(M_21_sub, reduction = "tsne", group.by = "orig.ident")
p4 <- DimPlot(M_21_sub, reduction = "tsne", label = TRUE, repel = TRUE)

pdf(file = "UMAP_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p1 + p2
dev.off()

pdf(file = "UMAP_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_21_sub, reduction = "umap", split.by = "orig.ident")
dev.off()

pdf(file = "TSNE_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p3 + p4
dev.off()

pdf(file = "TSNE_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_21_sub, reduction = "tsne", split.by = "orig.ident")
dev.off()

# Plot Seurat UMAP
pdf(file = "UMAP_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_21_sub, reduction = "umap",label= TRUE)
dev.off()

# Plot Seurat TSNE
pdf(file = "TSNE_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_21_sub, reduction = "tsne",label= TRUE)
dev.off()

# Add new split clusters to full Seurat object
temp=cbind(as.numeric(as.character(M@meta.data[["subcluster_9_2"]])),as.numeric(as.character(M@meta.data[["subcluster_9_2"]])))
row.names(temp)=M@assays[["SCT"]]@data@Dimnames[[2]]
new_clusters=mapvalues(as.numeric(as.character(M_21_sub@active.ident)), from=c(0,1,2), to=c(21,22,23))
temp[names(M_21_sub@active.ident),2]<-new_clusters
M@meta.data[["subcluster_21"]] <- temp[,2]
DimPlot(M, reduction="tsne", group.by = "subcluster_21", label=T)
Idents(M)<-"subcluster_21"


############################
#                          #
# Remove Unknown T Cluster # # This cluster has a high ribosomal mRNA count, so we decided to move it
#                          #
############################

M=subset(M, idents = setdiff(unique(M@active.ident), "100"))

dir.create("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/No_Unknown")
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/No_Unknown")

# Visualization
p1 <- DimPlot(M, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(M, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(M, reduction = "tsne", group.by = "orig.ident")
p4 <- DimPlot(M, reduction = "tsne", label = TRUE, repel = TRUE)

pdf(file = "UMAP_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p1 + p2
dev.off()

pdf(file = "UMAP_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M, reduction = "umap", split.by = "orig.ident")
dev.off()

pdf(file = "TSNE_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p3 + p4
dev.off()

pdf(file = "TSNE_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M, reduction = "tsne", split.by = "orig.ident")
dev.off()

# Plot Seurat UMAP
pdf(file = "UMAP_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M, reduction = "umap",label= TRUE)
dev.off()

# Plot Seurat TSNE
pdf(file = "TSNE_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M, reduction = "tsne",label= TRUE)
dev.off()


###################
#                 #
# Save New Object #
#                 #
###################

# Save new object
saveRDS(M, file = "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Seurat_Liver_30_subcluster_CD4.rds")
