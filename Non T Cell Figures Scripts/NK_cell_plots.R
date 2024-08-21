##############################
#                            #
#       NK Cell Plots        #
# Peters - Single Cell Liver #
#     Erica DePasquale       #
#        09 Oct 2023         #
#                            #
##############################

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(scales)
library(cowplot)

# Read in object
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/")
M <- readRDS("Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")

# Subset to NK cells only
Idents(object = M) <- "cluster_names_new"
Clusters_Keep=c("CD56dim NK", "CD56bright NK")
M_N=subset(x = M, idents = Clusters_Keep)

# Make TSNE plot
M_N <- RunPCA(M_N, npcs = 30, verbose = FALSE)
M_N <- RunUMAP(M_N, reduction = "pca", dims = 1:30)
M_N <- RunTSNE(M_N, reduction = "pca", dims = 1:30)

DimPlot(M_N, reduction="tsne", label=T, repel=T, raster=F)

# Pull in new colors
myColors=read.table("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Color_Guide.txt", sep = "\t", header=T, comment.char="*")
myColors[4,2]<-"Naïve B"
myColors=myColors[which(myColors$Cluster_Name %in% Clusters_Keep),]
row.names(myColors)=myColors$Cluster_Name
myColors=myColors[levels(M_N@active.ident),]

pdf(file = "TSNE_NK_Groups_New_Colors.pdf", width = 9.5, height = 8)
par(mar=c(2, 2, 2, 2))
DimPlot(M_N, reduction="tsne", label=T, repel=T, raster=F, cols=myColors$Cluster_Color, pt.size=1)
dev.off()

# Plot by ACR type
pdf(file = "TSNE_NK_ACRType.pdf", width = 9, height = 8)
par(mar=c(2, 2, 2, 2))
DimPlot(M_N, reduction="tsne", group.by="ACRType", label=F, repel=T, raster=F, pt.size=1)
dev.off()

# Prop plot by ACR type
data_long=as.data.frame(cbind(ACRType=M_N@meta.data[["ACRType"]], Cluster=M_N@meta.data[["cluster_names_new"]]))
data=as.data.frame.matrix(table(data_long))

# Read and wrangle data
data=t(data)
data <- sweep(data, 2, colSums(data), "/")*100

myColors=myColors[order(myColors$Cluster_Name, decreasing=F),]
my_colors=myColors$Cluster_Color
data=as.data.frame(cbind(data, colors=my_colors))
my_colors=data$colors
data=data[,-5]

# Make plot
pdf("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellTypeFrequencies_NK_byACRType_new_colors.pdf", width = 9, height = 9)
par(mar = c(8,4,4,16), xpd = T)
barplot(as.matrix(data[nrow(data):1,]), col = rev(my_colors), xaxt = "n", ylab = "Population frequency (%)", border = NA)
axis(side = 1, at = seq(1,ncol(data))*1.2-0.5, labels = colnames(data), las = 2)
legend(x = ncol(data)*1.2+0.5, y = 100, legend = row.names(data), fill = my_colors, bty = "n", border = NA)
dev.off()

write.table(data, "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellTypeFrequencies_NK_cell_byACRType_new_colors.txt", sep="\t", quote=F)

# Normalize data using log for visualization purposes
M_N <- NormalizeData(M_N, assay="RNA")
all.genes <- rownames(M_N)
M_N <- ScaleData(M_N, features = all.genes, assay="RNA")

myGenes=unique(c("NCAM1", "IL2RB", "EOMES", "CX3CR1", "GZMB", "CD160"))

pdf(file = "NK_Manual_RNA_Stacked_Violin_Names.pdf", width = 4, height = 6)
par(mar=c(2, 2, 2, 2))
VlnPlot(M_N, assay="RNA", group.by="cluster_names_new", stack=T, features=myGenes, flip=T, sort=T) + NoLegend()
dev.off()


pdf(file = "NK_Manual_RNA_DotPlot_Names_log.pdf", width = 8, height = 4)
par(mar=c(4, 4, 4, 4))
DotPlot(object=M_N, assay="RNA", features = myGenes, dot.scale=10, scale=F) + labs(y = NULL, x=NULL) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

p1<-DotPlot(object=M_N, assay="RNA", features = myGenes, dot.scale=10, scale=F) + labs(y = NULL, x=NULL) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
write.table(p1[["data"]], "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/NK_Manual_RNA_DotPlot_Names_log.txt", sep="\t", quote=F)

