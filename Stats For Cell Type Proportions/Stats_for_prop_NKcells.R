# Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(Seurat)

# Read in object
M <- readRDS("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")
data_long=as.data.frame(cbind(Sample=M@meta.data[["orig.ident"]], Cluster=M@meta.data[["cluster_names_new"]]))
data=as.data.frame.matrix(table(data_long))
data=cbind(sample=row.names(data), data)
data=data[order(data$sample),]

# Read and wrangle data
data=t(data)

x=as.data.frame(cbind(M@meta.data$orig.ident, M@meta.data$ACRType))
y=distinct(x)
y=y[order(y$V1),]
sum(y$V1!=colnames(data)) # 0, meaning they match

data=data[-1,]
rownamesdata=row.names(data)
data=apply(data, 2, as.numeric)
row.names(data)<-rownamesdata
data<-data[c(4, 5),]
data <- sweep(data, 2, colSums(data), "/")*100

#################
# CD56bright NK #
#################
z=cbind(y, data[1,]) #CD56bright NK
colnames(z)<-c("sample", "ACRType", "bright")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(bright ~ ACRType,  data = z)
# .y.   group1   group2                   p   p.adj p.format p.signif method  
# <chr>  <chr>    <chr>             <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 bright Donor    Not ACR           1         1 1.00     ns       Wilcoxon
# 2 bright Donor    Late ACR          0.850     1 0.85     ns       Wilcoxon
# 3 bright Donor    Resolved Late ACR 0.589     1 0.59     ns       Wilcoxon
# 4 bright Not ACR  Late ACR          0.654     1 0.65     ns       Wilcoxon
# 5 bright Not ACR  Resolved Late ACR 0.381     1 0.38     ns       Wilcoxon
# 6 bright Late ACR Resolved Late ACR 0.381     1 0.38     ns       Wilcoxon


##############
# CD56dim NK #
##############
z=cbind(y, data[2,]) #CD56dim NK
colnames(z)<-c("sample", "ACRType", "dim")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(dim ~ ACRType,  data = z)
# .y.   group1   group2                   p   p.adj p.format p.signif method  
# <chr> <chr>    <chr>             <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 dim   Donor    Not ACR           1         1 1.00     ns       Wilcoxon
# 2 dim   Donor    Late ACR          0.850     1 0.85     ns       Wilcoxon
# 3 dim   Donor    Resolved Late ACR 0.589     1 0.59     ns       Wilcoxon
# 4 dim   Not ACR  Late ACR          0.654     1 0.65     ns       Wilcoxon
# 5 dim   Not ACR  Resolved Late ACR 0.381     1 0.38     ns       Wilcoxon
# 6 dim   Late ACR Resolved Late ACR 0.381     1 0.38     ns       Wilcoxon
