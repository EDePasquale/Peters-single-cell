##############################
#                            #
#      Add Cluster Names     #
# Peters - Single Cell Liver #
#     Erica DePasquale       #
#        2 Oct 2023         #
#                            #
##############################

# Load libraries
library(Seurat)

# Pull in colors
myColors=read.table("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Color_Guide.txt", sep = "\t", header=T, comment.char="*")
myColors[4,2]<-"Naïve B"
myColors_full=myColors[,2:3]
myColors_full=myColors_full[order(myColors_full$Cluster_Name),]
myColors_redu=myColors[,4:5]
myColors_redu=myColors_redu[order(myColors_redu$Cluster_Name_Redu),]
myColors_redu=unique(myColors_redu)

# Read in subclustered objects
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/")
M<-readRDS("Seurat_Liver_30_subcluster.rds")

# Add cluster names
cluster_names=c("CD8 Effector Memory T", #0
                "CD56bright NK", #1
                "Monocyte-Derived Macrophage", #2
                "Naïve B", #3
                "LST1+ Kupffer", #4
                "CD8 TRM Activated", #5
                "MAIT", #6
                "CD56dim NK", #7
                "CD8 Effector T", #8
                "IFI27+ Kupffer", #9
                "Cycling T", #10
                "Gamma Delta T", #11
                "LSEC/VEC", #12
                "Hepatocyte", #13
                "pDC", #14
                "C1Q+ Kupffer", #15
                "Cholangiocyte", #16
                "Plasma", #17
                "CD1C+ Kupffer", #18
                "PTPRC+ Kupffer", #19
                "Transitional B", #20
                "CD4 naive T") #21
M@meta.data[["cluster_names_new"]] <- mapvalues(M@meta.data[["subcluster_9_2"]], from=0:21, to=cluster_names)
Idents(M) = M@meta.data$cluster_names_new
DimPlot(M, reduction="umap", label=T)
DimPlot(M, reduction="tsne", label=T)

# Plot Seurat UMAP
pdf(file = "UMAP_Liver_25_Seurat_Groups_New.pdf", width = 12, height = 8)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M, reduction = "umap", group.by= "cluster_names_new", cols=myColors_full$Cluster_Color, label= TRUE, repel=T)
dev.off()

# Plot Seurat TSNE
pdf(file = "TSNE_Liver_25_Seurat_Groups_New.pdf", width = 12, height = 8)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M, reduction = "tsne",group.by= "cluster_names_new", cols=myColors_full$Cluster_Color, label= TRUE, repel=T)
dev.off()

# Reduced cluster names
cluster_names_redu=c("CD8T", #0
                     "NK", #1
                     "Monocyte-Derived Macrophage", #2
                     "B", #3
                     "Kupffer", #4
                     "CD8T", #5
                     "CD8T", #6
                     "NK", #7
                     "CD8T", #8
                     "Kupffer", #9
                     "CD8T", #10
                     "Gamma Delta T", #11
                     "LSEC/VEC", #12
                     "Hepatocyte", #13
                     "pDC", #14
                     "Kupffer", #15
                     "Cholangiocyte", #16
                     "Plasma", #17
                     "Kupffer", #18
                     "Kupffer", #19
                     "B", #20
                     "CD4T") #21
M@meta.data[["cluster_names_new_redu"]] <- mapvalues(M@meta.data[["subcluster_9_2"]], from=0:21, to=cluster_names_redu)

# Plot Seurat UMAP
pdf(file = "UMAP_Liver_25_Seurat_Groups_New_Redu.pdf", width = 10, height = 8)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M, reduction = "umap",group.by= "cluster_names_new_redu", cols=myColors_redu$Cluster_Color_Redu, label= TRUE, repel=T)
dev.off()

# Plot Seurat TSNE
pdf(file = "TSNE_Liver_25_Seurat_Groups_New_Redu.pdf", width = 10, height = 8)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M, reduction = "tsne",group.by= "cluster_names_new_redu", cols=myColors_redu$Cluster_Color_Redu, label= TRUE, repel=T)
dev.off()

saveRDS(M, file = "Seurat_Liver_30_subcluster_names.rds")
