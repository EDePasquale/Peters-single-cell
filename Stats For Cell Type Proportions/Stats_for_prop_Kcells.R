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
data<-data[c(1,2,13,15,17,21),]
data <- sweep(data, 2, colSums(data), "/")*100

################
# C1Q+ Kupffer #
################
z=cbind(y, data[1,]) #C1Q+ Kupffer
colnames(z)<-c("sample", "ACRType", "C1Q")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(C1Q ~ ACRType,  data = z)
# .y.   group1   group2                   p   p.adj p.format p.signif method  
# <chr> <chr>    <chr>              <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 C1Q   Donor    Not ACR           0.0238  0.14 0.024    *        Wilcoxon
# 2 C1Q   Donor    Late ACR          0.0564  0.28 0.056    ns       Wilcoxon
# 3 C1Q   Donor    Resolved Late ACR 0.818   1    0.818    ns       Wilcoxon
# 4 C1Q   Not ACR  Late ACR          0.722   1    0.722    ns       Wilcoxon
# 5 C1Q   Not ACR  Resolved Late ACR 0.548   1    0.548    ns       Wilcoxon
# 6 C1Q   Late ACR Resolved Late ACR 0.508   1    0.508    ns       Wilcoxon


#################
# CD1C+ Kupffer #
#################
z=cbind(y, data[2,]) #CD1C+ Kupffer
colnames(z)<-c("sample", "ACRType", "CD1C")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(CD1C ~ ACRType,  data = z)
# .y.   group1   group2                   p   p.adj p.format p.signif method  
# <chr> <chr>    <chr>              <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 CD1C  Donor    Not ACR           0.0275 0.14  0.028    *        Wilcoxon
# 2 CD1C  Donor    Late ACR          0.0113 0.068 0.011    *        Wilcoxon
# 3 CD1C  Donor    Resolved Late ACR 0.0438 0.18  0.044    *        Wilcoxon
# 4 CD1C  Not ACR  Late ACR          0.635  1     0.635    ns       Wilcoxon
# 5 CD1C  Not ACR  Resolved Late ACR 0.604  1     0.604    ns       Wilcoxon
# 6 CD1C  Late ACR Resolved Late ACR 0.969  1     0.969    ns       Wilcoxon

cell_type="CD1C"
z=cbind(y, data[cell_type,])
colnames(z)<-c("sample", "ACRType", cell_type)
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))
formula <- paste(cell_type, "~ ACRType")
print(compare_means(as.formula(formula),  data = z))
### Look to pick comparisons
my_comparisons <- list( c("Donor", "Not ACR"), c("Donor", "Late ACR"))
pdf(paste0("Stats_for_prop_",cell_type,".pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "ACRType", y = cell_type, color = "ACRType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="ACR Type", y = paste0("Proportion ",cell_type," (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle(cell_type)
dev.off()

##################
# IFI27+ Kupffer #
##################
z=cbind(y, data[3,]) #IFI27+ Kupffer
colnames(z)<-c("sample", "ACRType", "IFI27")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(IFI27 ~ ACRType,  data = z)
# .y.   group1   group2                   p   p.adj p.format p.signif method  
# <chr> <chr>    <chr>              <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 IFI27 Donor    Not ACR           0.167   0.83 0.167    ns       Wilcoxon
# 2 IFI27 Donor    Late ACR          0.0322  0.19 0.032    *        Wilcoxon
# 3 IFI27 Donor    Resolved Late ACR 0.180   0.83 0.180    ns       Wilcoxon
# 4 IFI27 Not ACR  Late ACR          0.594   0.83 0.594    ns       Wilcoxon
# 5 IFI27 Not ACR  Resolved Late ACR 0.167   0.83 0.167    ns       Wilcoxon
# 6 IFI27 Late ACR Resolved Late ACR 0.330   0.83 0.330    ns       Wilcoxon
  

#################
# LST1+ Kupffer #
#################
z=cbind(y, data[4,]) #LST1+ Kupffer
colnames(z)<-c("sample", "ACRType", "LST1")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(LST1 ~ ACRType,  data = z)
# .y.   group1   group2                   p   p.adj p.format p.signif method  
# <chr> <chr>    <chr>               <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 LST1  Donor    Not ACR           0.0476  0.19  0.0476   *        Wilcoxon
# 2 LST1  Donor    Late ACR          0.00640 0.038 0.0064   **       Wilcoxon
# 3 LST1  Donor    Resolved Late ACR 0.0152  0.076 0.0152   *        Wilcoxon
# 4 LST1  Not ACR  Late ACR          0.722   1     0.7221   ns       Wilcoxon
# 5 LST1  Not ACR  Resolved Late ACR 0.795   1     0.7954   ns       Wilcoxon
# 6 LST1  Late ACR Resolved Late ACR 0.907   1     0.9070   ns       Wilcoxon


###############################
# Monocyte-Derived Macrophage #
###############################
z=cbind(y, data[5,]) #Monocyte-Derived Macrophage
colnames(z)<-c("sample", "ACRType", "MDM")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(MDM ~ ACRType,  data = z)
# .y.   group1   group2                   p   p.adj p.format p.signif method  
# <chr> <chr>    <chr>               <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 MDM   Donor    Not ACR           0.167    0.67 0.1667   ns       Wilcoxon
# 2 MDM   Donor    Late ACR          0.00837  0.05 0.0084   **       Wilcoxon
# 3 MDM   Donor    Resolved Late ACR 0.132    0.66 0.1320   ns       Wilcoxon
# 4 MDM   Not ACR  Late ACR          0.203    0.67 0.2034   ns       Wilcoxon
# 5 MDM   Not ACR  Resolved Late ACR 0.381    0.76 0.3810   ns       Wilcoxon
# 6 MDM   Late ACR Resolved Late ACR 0.791    0.79 0.7910   ns       Wilcoxon


##################
# PTPRC+ Kupffer #
##################
z=cbind(y, data[6,]) #PTPRC+ Kupffer
colnames(z)<-c("sample", "ACRType", "PTPRC")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

# Perform pairwise comparisons
compare_means(PTPRC ~ ACRType,  data = z)
# .y.   group1   group2                   p   p.adj p.format p.signif method  
# <chr> <chr>    <chr>              <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 PTPRC Donor    Not ACR           0.905   1    0.905    ns       Wilcoxon
# 2 PTPRC Donor    Late ACR          0.0428  0.26 0.043    *        Wilcoxon
# 3 PTPRC Donor    Resolved Late ACR 0.572   1    0.572    ns       Wilcoxon
# 4 PTPRC Not ACR  Late ACR          0.0660  0.33 0.066    ns       Wilcoxon
# 5 PTPRC Not ACR  Resolved Late ACR 1       1    1.000    ns       Wilcoxon
# 6 PTPRC Late ACR Resolved Late ACR 0.159   0.64 0.159    ns       Wilcoxon

