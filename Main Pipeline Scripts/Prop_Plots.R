##############################
#                            #
#    Make Proportion Plot    #
# Peters - Single Cell Liver #
#     Erica DePasquale       #
#        03 Oct 2022         #
#                            #
##############################


########################################
#                                      #
# <><><><> FULL CLUSTER NAMES <><><><> #
#                                      #
########################################


##################
# Prop by Sample #
##################

# Load libraries
library(ggplot2)
library(Seurat)

# Read in object
M <- readRDS("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Seurat_Liver_30_subcluster_names_meta.rds")
data_long=as.data.frame(cbind(Sample=M@meta.data[["orig.ident"]], Cluster=M@meta.data[["cluster_names_new"]]))
data=as.data.frame.matrix(table(data_long))

# Read and wrangle data
data=t(data)
data <- sweep(data, 2, colSums(data), "/")*100

myColors=read.table("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Color_Guide.txt", sep = "\t", header=T, comment.char="*")
myColors[4,2]<-"Naïve B"
myColors[22,2]<-"CD4 Naïve T"
myColors=myColors[order(myColors$Cluster_Name, decreasing=F),]
my_colors=myColors$Cluster_Color
data=as.data.frame(cbind(data, colors=my_colors))
my_colors=data$colors
data=data[,-31]

# Make plot
pdf("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellTypeFrequencies_byOrigIdent.pdf", width = 13, height = 9)
par(mar = c(8,4,4,16), xpd = T)
  barplot(as.matrix(data[nrow(data):1,]), col = rev(my_colors), xaxt = "n", ylab = "Population frequency (%)", border = NA)
  axis(side = 1, at = seq(1,ncol(data))*1.2-0.5, labels = colnames(data), las = 2)
  legend(x = ncol(data)*1.2+0.5, y = 100, legend = row.names(data), fill = my_colors, bty = "n", border = NA)
dev.off()


###################
# Prop by ACRType #
###################

# Load libraries
library(ggplot2)
library(Seurat)

data_long=as.data.frame(cbind(ACRType=M@meta.data[["ACRType"]], Cluster=M@meta.data[["cluster_names_new"]]))
data=as.data.frame.matrix(table(data_long))

# Read and wrangle data
data=t(data)
data <- sweep(data, 2, colSums(data), "/")*100

myColors=read.table("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Color_Guide.txt", sep = "\t", header=T, comment.char="*")
myColors[4,2]<-"Naïve B"
myColors[22,2]<-"CD4 Naïve T"
myColors=myColors[order(myColors$Cluster_Name, decreasing=F),]
my_colors=myColors$Cluster_Color
data=as.data.frame(cbind(data, colors=my_colors))
my_colors=data$colors
data=data[,-5]

# Make plot
pdf("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellTypeFrequencies_byACRType.pdf", width = 9, height = 9)
par(mar = c(8,4,4,16), xpd = T)
  barplot(as.matrix(data[nrow(data):1,]), col = rev(my_colors), xaxt = "n", ylab = "Population frequency (%)", border = NA)
  axis(side = 1, at = seq(1,ncol(data))*1.2-0.5, labels = colnames(data), las = 2)
  legend(x = ncol(data)*1.2+0.5, y = 100, legend = row.names(data), fill = my_colors, bty = "n", border = NA)
dev.off()


########################################
#                                      #
# <><><><> REDU CLUSTER NAMES <><><><> #
#                                      #
########################################


##################
# Prop by Sample #
##################

# Load libraries
library(ggplot2)
library(Seurat)

data_long=as.data.frame(cbind(Sample=M@meta.data[["orig.ident"]], Cluster=M@meta.data[["cluster_names_new_redu"]]))
data=as.data.frame.matrix(table(data_long))

# Read and wrangle data
data=t(data)
data <- sweep(data, 2, colSums(data), "/")*100

my_colors=c("#E31A1C", "#FFD92F", "#33A02C", "#FDC086", "#A6CEE3", "#999999", "#386CB0", "#B15928", "#F781BF", "#6A3D9A", "#BC80BD", "#A6D854")
data=as.data.frame(cbind(data, colors=my_colors))
my_colors=data$colors
data=data[,-31]

# Make plot
pdf("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellTypeFrequencies_byOrigIdent_redu.pdf", width = 13, height = 9)
par(mar = c(8,4,4,16), xpd = T)
  barplot(as.matrix(data[nrow(data):1,]), col = rev(my_colors), xaxt = "n", ylab = "Population frequency (%)", border = NA)
  axis(side = 1, at = seq(1,ncol(data))*1.2-0.5, labels = colnames(data), las = 2)
  legend(x = ncol(data)*1.2+0.5, y = 100, legend = row.names(data), fill = my_colors, bty = "n", border = NA)
dev.off()


###################
# Prop by ACRType #
###################

# Load libraries
library(ggplot2)
library(Seurat)

data_long=as.data.frame(cbind(ACRType=M@meta.data[["ACRType"]], Cluster=M@meta.data[["cluster_names_new_redu"]]))
data=as.data.frame.matrix(table(data_long))

# Read and wrangle data
data=t(data)
data <- sweep(data, 2, colSums(data), "/")*100

my_colors=c("#E31A1C", "#FFD92F", "#33A02C", "#FDC086", "#A6CEE3", "#999999", "#386CB0", "#B15928", "#F781BF", "#6A3D9A", "#BC80BD", "#A6D854")
data=as.data.frame(cbind(data, colors=my_colors))
my_colors=data$colors
data=data[,-5]

# Make plot
pdf("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellTypeFrequencies_byACRType_redu.pdf", width = 9, height = 9)
par(mar = c(8,4,4,16), xpd = T)
barplot(as.matrix(data[nrow(data):1,]), col = rev(my_colors), xaxt = "n", ylab = "Population frequency (%)", border = NA)
axis(side = 1, at = seq(1,ncol(data))*1.2-0.5, labels = colnames(data), las = 2)
legend(x = ncol(data)*1.2+0.5, y = 100, legend = row.names(data), fill = my_colors, bty = "n", border = NA)
dev.off()
