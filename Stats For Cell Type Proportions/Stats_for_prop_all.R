# Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(tidyr)

# Read in object
M <- readRDS("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")
data_long=as.data.frame(cbind(Sample=M@meta.data[["orig.ident"]], Cluster=M@meta.data[["cluster_names_new_redu"]]))
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
data=apply(data, 2, as.numeric)
data <- sweep(data, 2, colSums(data), "/")*100

z=cbind(y, data[3,]) #CD8T
colnames(z)<-c("sample", "ACRType", "CD8T")
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))

A=z[,2:3]
A %>%
  group_by(ACRType) %>%
  summarise(mean = mean(CD8T), sd = sd(CD8T), .groups = 'drop')


# Perform pairwise comparisons
compare_means(CD8T ~ ACRType,  data = z)
# .y.   group1   group2                   p   p.adj p.format p.signif method  
# <chr> <chr>    <chr>                <dbl>   <dbl> <chr>    <chr>    <chr>   
#   1 CD8T  Donor    Not ACR           0.905    0.9     0.90476  ns       Wilcoxon
# 2 CD8T  Donor    Late ACR          0.000147 0.00088 0.00015  ***      Wilcoxon
# 3 CD8T  Donor    Resolved Late ACR 0.00216  0.011   0.00216  **       Wilcoxon
# 4 CD8T  Not ACR  Late ACR          0.0392   0.16    0.03922  *        Wilcoxon
# 5 CD8T  Not ACR  Resolved Late ACR 0.0476   0.16    0.04762  *        Wilcoxon
# 6 CD8T  Late ACR Resolved Late ACR 0.205    0.41    0.20511  ns       Wilcoxon

# Visualize: Specify the comparisons you want
my_comparisons <- list( c("Donor", "Late ACR"), c("Donor", "Resolved Late ACR"))

pdf("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Stats_for_prop_CD8T.pdf", width = 6, height = 6)
par(mar = c(8,4,4,16), xpd = T)
  ggboxplot(z, x = "ACRType", y = "CD8T", color = "ACRType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="ACR Type", y = "Proportion CD8T (%)")
dev.off()

#########################

################
# ALL CLUSTERS #
################

# Read in object
M <- readRDS("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")
data_long=as.data.frame(cbind(Sample=M@meta.data[["orig.ident"]], Cluster=M@meta.data[["cluster_names_new_redu"]]))
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
rnames=row.names(data)
rnames2=gsub(" ", "", rnames)
rnames3=gsub("-", "", rnames2)
rnames4=gsub("/", "", rnames3)
data=apply(data, 2, as.numeric)
data <- sweep(data, 2, colSums(data), "/")*100
row.names(data)<-rnames4

for(i in row.names(data)){
  z=cbind(y, data[i,])
  colnames(z)<-c("sample", "ACRType", i)
  z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))
  formula <- paste(i, "~ ACRType")
  print(compare_means(as.formula(formula),  data = z))
}

dir.create("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Stat_Prop")
setwd("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Stat_Prop")

#B
cell_type="B"
z=cbind(y, data[cell_type,])
colnames(z)<-c("sample", "ACRType", cell_type)
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))
A=z[,2:3]
A %>%
  group_by(ACRType) %>%
  summarise(mean = mean(B), sd = sd(B), .groups = 'drop')
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

#CD4T
cell_type="CD4T"
z=cbind(y, data[cell_type,])
colnames(z)<-c("sample", "ACRType", cell_type)
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))
A=z[,2:3]
A %>%
  group_by(ACRType) %>%
  summarise(mean = mean(CD4T), sd = sd(CD4T), .groups = 'drop')
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

#CD8T
cell_type="CD8T"
z=cbind(y, data[cell_type,])
colnames(z)<-c("sample", "ACRType", cell_type)
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))
A=z[,2:3]
A %>%
  group_by(ACRType) %>%
  summarise(mean = mean(CD8T), sd = sd(CD8T), .groups = 'drop')
formula <- paste(cell_type, "~ ACRType")
print(compare_means(as.formula(formula),  data = z))
### Look to pick comparisons
my_comparisons <- list( c("Donor", "Late ACR"), c("Donor", "Resolved Late ACR"), c("Not ACR", "Late ACR"), c("Not ACR", "Resolved Late ACR"))
pdf(paste0("Stats_for_prop_",cell_type,".pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "ACRType", y = cell_type, color = "ACRType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="ACR Type", y = paste0("Proportion ",cell_type," (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle(cell_type)
dev.off()

#Cholangiocyte
cell_type="Cholangiocyte"
z=cbind(y, data[cell_type,])
colnames(z)<-c("sample", "ACRType", cell_type)
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))
formula <- paste(cell_type, "~ ACRType")
print(compare_means(as.formula(formula),  data = z))
### Look to pick comparisons
my_comparisons <- list()
pdf(paste0("Stats_for_prop_",cell_type,".pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "ACRType", y = cell_type, color = "ACRType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="ACR Type", y = paste0("Proportion ",cell_type," (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle(cell_type)
dev.off()

#GammaDeltaT
cell_type="GammaDeltaT"
z=cbind(y, data[cell_type,])
colnames(z)<-c("sample", "ACRType", cell_type)
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))
A=z[,2:3]
A %>%
  group_by(ACRType) %>%
  summarise(mean = mean(GammaDeltaT), sd = sd(GammaDeltaT), .groups = 'drop')
formula <- paste(cell_type, "~ ACRType")
print(compare_means(as.formula(formula),  data = z))
### Look to pick comparisons
my_comparisons <- list( c("Donor", "Late ACR"))
pdf(paste0("Stats_for_prop_",cell_type,".pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "ACRType", y = cell_type, color = "ACRType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="ACR Type", y = paste0("Proportion ",cell_type," (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle(cell_type)
dev.off()

#Hepatocyte
cell_type="Hepatocyte"
z=cbind(y, data[cell_type,])
colnames(z)<-c("sample", "ACRType", cell_type)
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))
formula <- paste(cell_type, "~ ACRType")
print(compare_means(as.formula(formula),  data = z))
### Look to pick comparisons
my_comparisons <- list()
pdf(paste0("Stats_for_prop_",cell_type,".pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "ACRType", y = cell_type, color = "ACRType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="ACR Type", y = paste0("Proportion ",cell_type," (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle(cell_type)
dev.off()

#Kupffer
cell_type="Kupffer"
z=cbind(y, data[cell_type,])
colnames(z)<-c("sample", "ACRType", cell_type)
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))
A=z[,2:3]
A %>%
  group_by(ACRType) %>%
  summarise(mean = mean(Kupffer), sd = sd(Kupffer), .groups = 'drop')
formula <- paste(cell_type, "~ ACRType")
print(compare_means(as.formula(formula),  data = z))
### Look to pick comparisons
my_comparisons <- list( c("Donor", "Not ACR"), c("Donor", "Late ACR"), c("Donor", "Resolved Late ACR"))
pdf(paste0("Stats_for_prop_",cell_type,".pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "ACRType", y = cell_type, color = "ACRType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="ACR Type", y = paste0("Proportion ",cell_type," (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle(cell_type)
dev.off()

#LSECVEC
cell_type="LSECVEC"
z=cbind(y, data[cell_type,])
colnames(z)<-c("sample", "ACRType", cell_type)
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))
formula <- paste(cell_type, "~ ACRType")
print(compare_means(as.formula(formula),  data = z))
### Look to pick comparisons
my_comparisons <- list()
pdf(paste0("Stats_for_prop_",cell_type,".pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "ACRType", y = cell_type, color = "ACRType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="ACR Type", y = paste0("Proportion ",cell_type," (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle(cell_type)
dev.off()

#MonocyteDeriveMacrophage
cell_type="MonocyteDerivedMacrophage"
z=cbind(y, data[cell_type,])
colnames(z)<-c("sample", "ACRType", cell_type)
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))
A=z[,2:3]
A %>%
  group_by(ACRType) %>%
  summarise(mean = mean(MonocyteDerivedMacrophage), sd = sd(MonocyteDerivedMacrophage), .groups = 'drop')
formula <- paste(cell_type, "~ ACRType")
print(compare_means(as.formula(formula),  data = z))
### Look to pick comparisons
my_comparisons <- list( c("Donor", "Not ACR"), c("Donor", "Late ACR"), c("Donor", "Resolved Late ACR"), c("Not ACR", "Late ACR"))
pdf(paste0("Stats_for_prop_",cell_type,".pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "ACRType", y = cell_type, color = "ACRType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="ACR Type", y = paste0("Proportion ",cell_type," (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle(cell_type)
dev.off()

#NK
cell_type="NK"
z=cbind(y, data[cell_type,])
colnames(z)<-c("sample", "ACRType", cell_type)
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))
A=z[,2:3]
A %>%
  group_by(ACRType) %>%
  summarise(mean = mean(NK), sd = sd(NK), .groups = 'drop')
formula <- paste(cell_type, "~ ACRType")
print(compare_means(as.formula(formula),  data = z))
### Look to pick comparisons
my_comparisons <- list()
pdf(paste0("Stats_for_prop_",cell_type,".pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "ACRType", y = cell_type, color = "ACRType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="ACR Type", y = paste0("Proportion ",cell_type," (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle(cell_type)
dev.off()

#pDC
cell_type="pDC"
z=cbind(y, data[cell_type,])
colnames(z)<-c("sample", "ACRType", cell_type)
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))
formula <- paste(cell_type, "~ ACRType")
print(compare_means(as.formula(formula),  data = z))
### Look to pick comparisons
my_comparisons <- list(c("Donor", "Late ACR"))
pdf(paste0("Stats_for_prop_",cell_type,".pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "ACRType", y = cell_type, color = "ACRType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="ACR Type", y = paste0("Proportion ",cell_type," (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle(cell_type)
dev.off()

#Plasma
cell_type="Plasma"
z=cbind(y, data[cell_type,])
colnames(z)<-c("sample", "ACRType", cell_type)
z$ACRType=factor(z$ACRType, levels = c("Donor", "Not ACR", "Late ACR", "Resolved Late ACR"))
formula <- paste(cell_type, "~ ACRType")
print(compare_means(as.formula(formula),  data = z))
### Look to pick comparisons
my_comparisons <- list( c("Donor", "Resolved Late ACR"), c("Late ACR", "Resolved Late ACR"))
pdf(paste0("Stats_for_prop_",cell_type,".pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "ACRType", y = cell_type, color = "ACRType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="ACR Type", y = paste0("Proportion ",cell_type," (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle(cell_type)
dev.off()

