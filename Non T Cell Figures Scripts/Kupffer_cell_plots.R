##############################
#                            #
#    Kupffer Cell Plots      #
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

# Subset to Kupffer cells only
Idents(object = M) <- "cluster_names_new"
Clusters_Keep=c("CD1C+ Kupffer", "C1Q+ Kupffer", "LST1+ Kupffer", "PTPRC+ Kupffer", "IFI27+ Kupffer", "Monocyte-Derived Macrophage")
M_K=subset(x = M, idents = Clusters_Keep)

# Make TSNE plot
M_K <- RunPCA(M_K, npcs = 30, verbose = FALSE)
M_K <- RunUMAP(M_K, reduction = "pca", dims = 1:30)
M_K <- RunTSNE(M_K, reduction = "pca", dims = 1:30)

# Pull in new colors
myColors=read.table("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Color_Guide.txt", sep = "\t", header=T, comment.char="*")
myColors[4,2]<-"Naïve B"
myColors=myColors[which(myColors$Cluster_Name %in% Clusters_Keep),]
row.names(myColors)=myColors$Cluster_Name
myColors=myColors[levels(M_K@active.ident),]

pdf(file = "TSNE_Kupffer_Groups_New_Colors.pdf", width = 10.5, height = 8)
par(mar=c(2, 2, 2, 2))
DimPlot(M_K, reduction="tsne", label=T, repel=T, raster=F, cols=myColors$Cluster_Color, pt.size=1)
dev.off()

# Plot by ACR type
pdf(file = "TSNE_Kupffer_ACRType.pdf", width = 9, height = 8)
par(mar=c(2, 2, 2, 2))
DimPlot(M_K, reduction="tsne", group.by="ACRType", label=F, repel=T, raster=F, pt.size=1)
dev.off()

# Prop plot by ACR type
data_long=as.data.frame(cbind(ACRType=M_K@meta.data[["ACRType"]], Cluster=M_K@meta.data[["cluster_names_new"]]))
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
pdf("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellTypeFrequencies_Kupffer_byACRType_new_colors.pdf", width = 9, height = 9)
par(mar = c(8,4,4,16), xpd = T)
barplot(as.matrix(data[nrow(data):1,]), col = rev(my_colors), xaxt = "n", ylab = "Population frequency (%)", border = NA)
axis(side = 1, at = seq(1,ncol(data))*1.2-0.5, labels = colnames(data), las = 2)
legend(x = ncol(data)*1.2+0.5, y = 100, legend = row.names(data), fill = my_colors, bty = "n", border = NA)
dev.off()

# Normalize data using log for visualization purposes
M_K <- NormalizeData(M_K, assay="RNA")
all.genes <- rownames(M_K)
M_K <- ScaleData(M_K, features = all.genes, assay="RNA")

#myGenes=unique(c("LILRB5", "CD5L", "MARCO", "HMOX1", "CD1C", "FCER1A", "CCR2", "FCGR2A", "CD68", "CD14"))
myGenes=unique(c("CD68", "CX3CR1", "ITGAM", "ADGRE1", "LY6G5C", "CCR2", "CSF1R", "CLEC4F", "MARCO", "C1QC", "IL18"))

pdf(file = "Kupffer_Manual_RNA_Stacked_Violin_Names.pdf", width = 11, height = 8)
par(mar=c(2, 2, 2, 2))
VlnPlot(M_K, assay="RNA", group.by="cluster_names_new", stack=T, features=myGenes, flip=T, sort=T) + NoLegend()
dev.off()


myGenes2=unique(c("CD68", "CD14", "LILRB5", "CD5L", "MARCO", "HMOX1", "CD1C", "FCER1A", "CCR2", "FCGR2A"))

pdf(file = "Kupffer_Manual_RNA_DotPlot_Names_log.pdf", width = 7.75, height = 4)
par(mar=c(4, 4, 4, 4))
DotPlot(object=M_K, assay="RNA", features = myGenes2, dot.scale=10, scale=F) + labs(y =NULL, x=NULL) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

DefaultAssay(M_K)<-"SCT"

pdf(file = "Kupffer_Manual_RNA_Markes_Anna.pdf", width = 11, height = 15)
par(mar=c(4, 4, 4, 4))
FeaturePlot(M_K, reduction="tsne", features=c("C1QA", "CD1C", "FTL", "LST1", "MALAT1", "S100A8"), slot="data", keep.scale="all")
dev.off()

pdf(file = "Kupffer_Manual_RNA_Markers_Erica.pdf", width = 11, height = 15)
par(mar=c(4, 4, 4, 4))
FeaturePlot(M_K, reduction="tsne", features=c("C1QA", "CD1C", "IFI27", "LST1", "PTPRC", "S100A8"), slot="data", keep.scale="all")
dev.off()
