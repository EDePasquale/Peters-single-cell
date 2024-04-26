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
data<-data[c(3,4,5,8,9,10,12,13,18),]
data <- sweep(data, 2, colSums(data), "/")*100

###################
# CD4 activated T #
###################
z=cbind(y, data[1,]) #CD4 activated T
colnames(z)<-c("sample", "ACRType", "CD4activatedT")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(CD4activatedT ~ ACRType,  data = z)
# .y.           group1   group2                p p.adj p.format p.signif method  
# <chr>         <chr>    <chr>             <dbl> <dbl> <chr>    <chr>    <chr>   
# 1 CD4activatedT Donor    Not ACR           1         1 1.00     ns       Wilcoxon
# 2 CD4activatedT Donor    Late ACR          0.785     1 0.79     ns       Wilcoxon
# 3 CD4activatedT Donor    Resolved Late ACR 0.809     1 0.81     ns       Wilcoxon
# 4 CD4activatedT Not ACR  Late ACR          0.906     1 0.91     ns       Wilcoxon
# 5 CD4activatedT Not ACR  Resolved Late ACR 0.896     1 0.90     ns       Wilcoxon
# 6 CD4activatedT Late ACR Resolved Late ACR 0.785     1 0.79     ns       Wilcoxon

###############
# CD4 naive T #
###############
z=cbind(y, data[2,]) #CD4 naive T
colnames(z)<-c("sample", "ACRType", "CD4naiveT")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(CD4naiveT ~ ACRType,  data = z)
# .y.       group1   group2                 p p.adj p.format p.signif method  
# <chr>     <chr>    <chr>              <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 CD4naiveT Donor    Not ACR           0.381   1    0.381    ns       Wilcoxon
# 2 CD4naiveT Donor    Late ACR          0.235   1    0.235    ns       Wilcoxon
# 3 CD4naiveT Donor    Resolved Late ACR 0.0649  0.39 0.065    ns       Wilcoxon
# 4 CD4naiveT Not ACR  Late ACR          0.301   1    0.301    ns       Wilcoxon
# 5 CD4naiveT Not ACR  Resolved Late ACR 0.262   1    0.262    ns       Wilcoxon
# 6 CD4naiveT Late ACR Resolved Late ACR 0.340   1    0.340    ns       Wilcoxon

#############
# CD4 TCM #
#############
z=cbind(y, data[3,]) #CD4 TCM
colnames(z)<-c("sample", "ACRType", "CD4TCM")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(CD4TCM ~ ACRType,  data = z)
# .y.      group1   group2                 p p.adj p.format p.signif method  
# <chr>    <chr>    <chr>              <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 CD4TCM Donor    Not ACR           0.860   1    0.860    ns       Wilcoxon
# 2 CD4TCM Donor    Late ACR          0.0261  0.16 0.026    *        Wilcoxon
# 3 CD4TCM Donor    Resolved Late ACR 0.462   1    0.462    ns       Wilcoxon
# 4 CD4TCM Not ACR  Late ACR          0.143   0.72 0.143    ns       Wilcoxon
# 5 CD4TCM Not ACR  Resolved Late ACR 0.878   1    0.878    ns       Wilcoxon
# 6 CD4TCM Late ACR Resolved Late ACR 0.240   0.96 0.240    ns       Wilcoxon


#########################
# CD8 Effector Memory T #
#########################
z=cbind(y, data[4,]) #CD8 Effector Memory T
colnames(z)<-c("sample", "ACRType", "CD8EM")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(CD8EM ~ ACRType,  data = z)
# .y.   group1   group2                p p.adj p.format p.signif method  
# <chr> <chr>    <chr>             <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 CD8EM Donor    Not ACR           0.604  1    0.60     ns       Wilcoxon
# 2 CD8EM Donor    Late ACR          0.161  0.97 0.16     ns       Wilcoxon
# 3 CD8EM Donor    Resolved Late ACR 0.630  1    0.63     ns       Wilcoxon
# 4 CD8EM Not ACR  Late ACR          0.738  1    0.74     ns       Wilcoxon
# 5 CD8EM Not ACR  Resolved Late ACR 0.714  1    0.71     ns       Wilcoxon
# 6 CD8EM Late ACR Resolved Late ACR 0.559  1    0.56     ns       Wilcoxon


##################
# CD8 Effector T #
##################
z=cbind(y, data[5,]) #CD8 Effector T
colnames(z)<-c("sample", "ACRType", "CD8E")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(CD8E ~ ACRType,  data = z)
# .y.   group1   group2                p p.adj p.format p.signif method  
# <chr> <chr>    <chr>             <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 CD8E  Donor    Not ACR           0.381  1    0.38     ns       Wilcoxon
# 2 CD8E  Donor    Late ACR          0.205  0.92 0.21     ns       Wilcoxon
# 3 CD8E  Donor    Resolved Late ACR 0.818  1    0.82     ns       Wilcoxon
# 4 CD8E  Not ACR  Late ACR          0.738  1    0.74     ns       Wilcoxon
# 5 CD8E  Not ACR  Resolved Late ACR 0.167  0.92 0.17     ns       Wilcoxon
# 6 CD8E  Late ACR Resolved Late ACR 0.154  0.92 0.15     ns       Wilcoxon

#################
# CD8 NK-like T #
#################
z=cbind(y, data[6,]) #CD8 NK-like T
colnames(z)<-c("sample", "ACRType", "CD8NKlike")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(CD8NKlike ~ ACRType,  data = z)
# .y.       group1   group2                p p.adj p.format p.signif method  
# <chr>     <chr>    <chr>             <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 CD8NKlike Donor    Not ACR           0.604     1 0.60     ns       Wilcoxon
# 2 CD8NKlike Donor    Late ACR          0.470     1 0.47     ns       Wilcoxon
# 3 CD8NKlike Donor    Resolved Late ACR 0.310     1 0.31     ns       Wilcoxon
# 4 CD8NKlike Not ACR  Late ACR          0.574     1 0.57     ns       Wilcoxon
# 5 CD8NKlike Not ACR  Resolved Late ACR 0.905     1 0.90     ns       Wilcoxon
# 6 CD8NKlike Late ACR Resolved Late ACR 0.533     1 0.53     ns       Wilcoxon

#############
# Cycling T #
#############
z=cbind(y, data[7,]) #Cycling T
colnames(z)<-c("sample", "ACRType", "Cycling")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(Cycling ~ ACRType,  data = z)
# .y.   group1   group2                   p   p.adj p.format p.signif method  
# <chr>   <chr>    <chr>             <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 Cycling Donor    Not ACR           0.905     1 0.90     ns       Wilcoxon
# 2 Cycling Donor    Late ACR          0.235     1 0.23     ns       Wilcoxon
# 3 Cycling Donor    Resolved Late ACR 0.630     1 0.63     ns       Wilcoxon
# 4 Cycling Not ACR  Late ACR          0.301     1 0.30     ns       Wilcoxon
# 5 Cycling Not ACR  Resolved Late ACR 0.548     1 0.55     ns       Wilcoxon
# 6 Cycling Late ACR Resolved Late ACR 0.791     1 0.79     ns       Wilcoxon

#################
# Gamma Delta T #
#################
z=cbind(y, data[8,]) #Gamma Delta T
colnames(z)<-c("sample", "ACRType", "GDT")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(GDT ~ ACRType,  data = z)
# .y.   group1   group2                   p   p.adj p.format p.signif method  
# <chr> <chr>    <chr>               <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 GDT   Donor    Not ACR           0.0906  0.36  0.0906   ns       Wilcoxon
# 2 GDT   Donor    Late ACR          0.00446 0.027 0.0045   **       Wilcoxon
# 3 GDT   Donor    Resolved Late ACR 0.462   1     0.4624   ns       Wilcoxon
# 4 GDT   Not ACR  Late ACR          1       1     1.0000   ns       Wilcoxon
# 5 GDT   Not ACR  Resolved Late ACR 0.364   1     0.3642   ns       Wilcoxon
# 6 GDT   Late ACR Resolved Late ACR 0.0322  0.16  0.0322   *        Wilcoxon

cell_type="GDT"
z=cbind(y, data[cell_type,])
colnames(z)<-c("sample", "ACRType", cell_type)
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))
formula <- paste(cell_type, "~ ACRType")
print(compare_means(as.formula(formula),  data = z))
### Look to pick comparisons
my_comparisons <- list( c("Donor", "Not ACR"))
pdf(paste0("Stats_for_prop_",cell_type,".pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "ACRType", y = cell_type, color = "ACRType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="ACR Type", y = paste0("Proportion ",cell_type," (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle(cell_type)
dev.off()

########
# MAIT #
########
z=cbind(y, data[9,]) #MAIT
colnames(z)<-c("sample", "ACRType", "MAIT")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(MAIT ~ ACRType,  data = z)
# .y.   group1   group2                   p   p.adj p.format p.signif method  
# <chr> <chr>    <chr>              <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 MAIT  Donor    Not ACR           0.604   1    0.604    ns       Wilcoxon
# 2 MAIT  Donor    Late ACR          0.970   1    0.970    ns       Wilcoxon
# 3 MAIT  Donor    Resolved Late ACR 0.468   1    0.468    ns       Wilcoxon
# 4 MAIT  Not ACR  Late ACR          0.301   1    0.301    ns       Wilcoxon
# 5 MAIT  Not ACR  Resolved Late ACR 0.896   1    0.896    ns       Wilcoxon
# 6 MAIT  Late ACR Resolved Late ACR 0.0866  0.52 0.087    ns       Wilcoxon
