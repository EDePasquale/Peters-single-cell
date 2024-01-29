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
data<-data[c(18, 22),]
data <- sweep(data, 2, colSums(data), "/")*100

###########
# Naïve B #
###########
z=cbind(y, data[1,]) #Naïve B
colnames(z)<-c("sample", "ACRType", "naiveB")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(naiveB ~ ACRType,  data = z)
# .y.   group1   group2                   p   p.adj p.format p.signif method  
# <chr>  <chr>    <chr>             <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 naiveB Donor    Not ACR           0.511     1 0.51     ns       Wilcoxon
# 2 naiveB Donor    Late ACR          0.427     1 0.43     ns       Wilcoxon
# 3 naiveB Donor    Resolved Late ACR 0.807     1 0.81     ns       Wilcoxon
# 4 naiveB Not ACR  Late ACR          1         1 1.00     ns       Wilcoxon
# 5 naiveB Not ACR  Resolved Late ACR 0.511     1 0.51     ns       Wilcoxon
# 6 naiveB Late ACR Resolved Late ACR 0.531     1 0.53     ns       Wilcoxon



##################
# Transitional B #
##################
z=cbind(y, data[2,]) #Transitional B
colnames(z)<-c("sample", "ACRType", "TransB")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(TransB ~ ACRType,  data = z)
# .y.   group1   group2                   p   p.adj p.format p.signif method  
# <chr>  <chr>    <chr>             <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 TransB Donor    Not ACR           0.511     1 0.51     ns       Wilcoxon
# 2 TransB Donor    Late ACR          0.427     1 0.43     ns       Wilcoxon
# 3 TransB Donor    Resolved Late ACR 0.807     1 0.81     ns       Wilcoxon
# 4 TransB Not ACR  Late ACR          1         1 1.00     ns       Wilcoxon
# 5 TransB Not ACR  Resolved Late ACR 0.511     1 0.51     ns       Wilcoxon
# 6 TransB Late ACR Resolved Late ACR 0.531     1 0.53     ns       Wilcoxon
