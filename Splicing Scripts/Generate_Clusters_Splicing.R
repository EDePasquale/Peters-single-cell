# Prepare files for splicing analysis

# Load libraries
library(Seurat)

# Read in full object
M<-readRDS("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")
M@meta.data[["cluster_names_new_clean"]]<-gsub("-", "", gsub("/", "", gsub(" ", "", gsub("ï", "i", gsub("\\+", "", M@meta.data$cluster_names_new)))))

# Set working directory to where cluster files will be stored
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/splicing/cluster_files")

# Create table for all samples together
X=M@meta.data
Y=row.names(X)
Z=gsub("_[0-9]*_[0-9]*", "", Y)
A=cbind(Z, X$cluster_names_new_clean, X$cluster_names_new_clean)
write.table(A, "Clusters_all.txt", quote=F, row.names=F, col.names=F, sep="\t")

# Create tables for each sample
for(i in unique(M@meta.data$orig.ident)){
  X=M@meta.data[which(M@meta.data$orig.ident == i), ]
  Y=row.names(X)
  Z=gsub("_[0-9]*_[0-9]*", "", Y)
  A=cbind(Z, X$cluster_names_new_CD4sub_clean)
  write.table(A, paste0("Clusters_", i, ".txt"), quote=F, row.names=F, col.names=F, sep="\t")
}