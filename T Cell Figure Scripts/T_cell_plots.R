##############################
#                            #
#       T Cell Plots         #
# Peters - Single Cell Liver #
#     Erica DePasquale       #
#        03 Oct 2023         #
#                            #
##############################

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(scales)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(stringr)

# Read in object
setwd("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/")
M <- readRDS("Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")

# Subset to T-cells only
Idents(object = M) <- "cluster_names_new"
Clusters_Keep=c("CD8 Effector Memory T", "CD4 Naïve T", "CD4 Activated T", "CD4 TCM", "CD8 NK-like T", "Gamma Delta T", "MAIT", "CD8 Effector T", "Cycling T")
M_T=subset(x = M, idents = Clusters_Keep)

# Make TSNE plot
M_T <- RunPCA(M_T, npcs = 30, verbose = FALSE)
M_T <- RunUMAP(M_T, reduction = "pca", dims = 1:30)
M_T <- RunTSNE(M_T, reduction = "pca", dims = 1:30)

# Pull in new colors
myColors=read.table("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Color_Guide.txt", sep = "\t", header=T, comment.char="*")
myColors[4,2]<-"Naïve B"
myColors[22,2]<-"CD4 Naïve T"
myColors=myColors[which(myColors$Cluster_Name %in% Clusters_Keep),]
row.names(myColors)=myColors$Cluster_Name
myColors=myColors[levels(M_T@active.ident),]

pdf(file = "TSNE_Tcell_Groups_New_Colors.pdf", width = 9.5, height = 8)
par(mar=c(2, 2, 2, 2))
  DimPlot(M_T, reduction="tsne", label=T, repel=T, raster=F, cols=myColors$Cluster_Color, pt.size=1)
dev.off()

# Plot by ACR type
pdf(file = "TSNE_Tcell_ACRType.pdf", width = 9, height = 8)
par(mar=c(2, 2, 2, 2))
  DimPlot(M_T, reduction="tsne", group.by="ACRType", label=F, repel=T, raster=F, pt.size=1)
dev.off()

# Plot by Clone size
pdf(file = "TSNE_Tcell_Clonesize.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
FeaturePlot(M_T, reduction="tsne", features="clonesize", raster=F, pt.size=1, cols = c("lightgrey", "red"))
dev.off()

# Full clusters clone size
pdf(file = "TSNE_Clonesize.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
FeaturePlot(M, reduction="tsne", features ="clonesize", raster=F, cols=c("lightgrey", "red"))
dev.off()

# Prop plot by ACR type
data_long=as.data.frame(cbind(ACRType=M_T@meta.data[["ACRType"]], Cluster=M_T@meta.data[["cluster_names_new"]]))
data=as.data.frame.matrix(table(data_long))

# Read and wrangle data
data=t(data)
data <- sweep(data, 2, colSums(data), "/")*100

myColors=myColors[order(myColors$Cluster_Name, decreasing=F),]
my_colors=myColors$Cluster_Color
data=as.data.frame(cbind(data, colors=my_colors))
my_colors=data$colors
data=data[,-5]

write.table(data, "/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellTypeFrequencies_T_cell_byACRType_new_colors.txt", sep="\t", quote=F)

# Make plot
pdf("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellTypeFrequencies_T_cell_byACRType_new_colors.pdf", width = 9, height = 9)
par(mar = c(8,4,4,16), xpd = T)
  barplot(as.matrix(data[nrow(data):1,]), col = rev(my_colors), xaxt = "n", ylab = "Population frequency (%)", border = NA)
  axis(side = 1, at = seq(1,ncol(data))*1.2-0.5, labels = colnames(data), las = 2)
  legend(x = ncol(data)*1.2+0.5, y = 100, legend = row.names(data), fill = my_colors, bty = "n", border = NA)
dev.off()

# Normalize data using log for visualization purposes
M_T <- NormalizeData(M_T, assay="RNA")
all.genes <- rownames(M_T)
M_T <- ScaleData(M_T, features = all.genes, assay="RNA")

myGenes2=unique(c("CD3G", "CD4", "CCR7", "SELL", "CD8A", "TRAV1-2", "CD27", "S1PR1", "PRF1", "GZMB", "CD69", "NCAM1", "KLRB1", "KLRC1", "KLRC2", "KLRD1", "CXCR6", "TRDV1", "MKI67"))

# Reorder levels for dotplot
M_T@meta.data[["cluster_names_new"]]<-factor(M_T@meta.data[["cluster_names_new"]], levels=c("Cycling T", "Gamma Delta T", "CD8 NK-like T", "CD8 Effector T", "CD8 Effector Memory T", "MAIT", "CD4 TCM", "CD4 Activated T","CD4 Naïve T"))

pdf(file = "Tcell_Manual_RNA_DotPlot_Names_log.pdf", width = 8, height = 4)
par(mar=c(4, 4, 4, 4))
  DotPlot(object=M_T, assay="RNA", group.by="cluster_names_new", features = myGenes2, dot.scale=10, scale=F) + labs(y = NULL, x=NULL) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

p1<-DotPlot(object=M_T, assay="RNA", group.by="cluster_names_new", features = myGenes2, dot.scale=10, scale=F) + labs(y = NULL, x=NULL) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
write.table(p1[["data"]], "/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Tcell_Manual_RNA_DotPlot_Names_log.txt", sep="\t", quote=F)

# Shared clones heatmap
paste2 <- function(...,sep=", ") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  gsub(paste0("(^",sep,"|",sep,"$)"),"",
       gsub(paste0(sep,sep),sep,
            do.call(paste,c(L,list(sep=sep)))))
}

clones_list=NULL
expt_list=unique(M@meta.data[["orig.ident"]])
for(expt in expt_list){
  
  # Pull in TCR data files
  clonotypes <- read.csv(paste0("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/", expt, "/clonotypes.csv"), sep=',', header=T)
  filtered_contig <- read.csv(paste0("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/", expt, "/filtered_contig_annotations.csv"), sep=',', header=T)
  merged_TCR=left_join(filtered_contig, clonotypes, by=c("raw_clonotype_id" = "clonotype_id"))
  clones=clonotypes[order(clonotypes$frequency, decreasing=T),"cdr3s_aa"]
  X=str_count(clones, "TRA") # added to keep only those clones with 1 TRA and 1 TRB
  Y=str_count(clones, "TRB")
  Z=cbind(X,Y, keep=rep(0, length(X)))
  for(rrow in 1:nrow(Z)){
    if(Z[rrow,1]==1 && Z[rrow,2]==1){Z[rrow,3]<-1}
  }
  myInd=which(Z[,3]==1)
  clones<-clones[myInd]
  clones_list[[expt]]<-clones
}

# Make plots for shared clones
unique_clones=unique(as.character(unlist(clones_list)))
clones_matrix_wide=as.data.frame(matrix(nrow=30, ncol=length(unique_clones)))
row.names(clones_matrix_wide)=expt_list
colnames(clones_matrix_wide)=unique_clones
for(i in 1:length(expt_list)){
  clones_matrix_wide[i,unique(unlist(clones_list[i]))]<-1
}
clones_matrix_wide[is.na(clones_matrix_wide)] <- 0

# Sort
clones_matrix_wide_sums=rbind(clones_matrix_wide, colSums(clones_matrix_wide))
Y=as.matrix(clones_matrix_wide_sums[nrow(clones_matrix_wide_sums),])
Z=clones_matrix_wide_sums[,order(-Y)]
clones_matrix_wide_sums=clones_matrix_wide_sums[,order(-Y)]

# Heatmap
clones_matrix_wide_redu=clones_matrix_wide_sums[,which(clones_matrix_wide_sums[nrow(clones_matrix_wide_sums),] > 1)]
clones_matrix_wide_redu=clones_matrix_wide_redu[-nrow(clones_matrix_wide_redu),]
clones_matrix_wide_redu=clones_matrix_wide_redu[order(row.names(clones_matrix_wide_redu)), ]
clones_matrix_wide_redu=clones_matrix_wide_redu[c(7,8,9,10,11,1,2,3,4,5,6,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30),]
pheatmap(clones_matrix_wide_redu,
         cluster_rows = F,
         color=c("gray92", "red"),
         gaps_row=c(5,6,8,11,13,16,21,22,24,26,27,28,29),
         width=20,
         height = 11.5,
         filename = "/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/heatmap_all_shared_clones_updated.pdf")

pheatmap(t(clones_matrix_wide_redu),
         cluster_cols = F,
         color=c("gray92", "red"),
         gaps_col=c(5,6,8,11,13,16,21,22,24,26,27,28,29),
         width=11.5,
         height = 20,
         angle_col=45,
         filename = "/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/heatmap_all_shared_clones_long_updated.pdf")


# TRA TRB frequncy barplot
data_list=NULL
top_ten_list=NULL
for(expt in 1:length(expt_list)){
  
  # Set up path and read in Seurat object
  name <- expt_list[expt]
  path <- paste0("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/", name)
  M <- readRDS(paste0("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/QC/", name, "/", name, "_filtered_matrices/Seurat.rds"))
  
  # Pull in TCR data files
  clonotypes <- read.csv(paste0("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/", name, "/clonotypes.csv"), sep=',', header=T)
  filtered_contig <- read.csv(paste0("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/", name, "/filtered_contig_annotations.csv"), sep=',', header=T)
  merged_TCR=left_join(filtered_contig, clonotypes, by=c("raw_clonotype_id" = "clonotype_id"))
  
  # Add full TCR and BCR data to Seurat object
  M@misc[["TCR_full_table"]]<-merged_TCR
  
  # Add TCR data to each cell
  AB_table=as.data.frame(matrix(nrow=length(unique(merged_TCR$barcode)),ncol=2))
  colnames(AB_table)=c("TRA", "TRB")
  row.names(AB_table)=unique(merged_TCR$barcode)
  for(i in unique(merged_TCR$barcode)){
    temp=merged_TCR[which(merged_TCR$barcode==i),] #subset table to just unique cell barcode i
    temp=temp[,c(1,3,6,23,28)] #reduce to columns we care about, for ease of visualization when coding
    TRA=temp[which(temp$chain=="TRA"),] #pull out column(s) that are TRA
    TRB=temp[which(temp$chain=="TRB"),] #pull out column(s) that are TRB
    AB_table[i,1]=paste(TRA$cdr3, collapse=",")
    AB_table[i,2]=paste(TRB$cdr3, collapse=",")
  }
  AB_table=as.data.frame(cbind(barcodes=row.names(AB_table), AB_table))
  data_to_add=as.data.frame(cbind(barcodes_full=M@assays[["RNA"]]@data@Dimnames[[2]], barcodes=M@assays[["RNA"]]@data@Dimnames[[2]]))
  data_to_add[,2]=gsub("_.*", "", data_to_add[,2])
  data_to_add=left_join(data_to_add, AB_table)
  data_to_add[which(data_to_add[,3] == ""),3]<-NA #deal with blank spaces from the above paste
  data_to_add[which(data_to_add[,4] == ""),4]<-NA
  M@meta.data[["TRA"]]<-data_to_add$TRA
  M@meta.data[["TRB"]]<-data_to_add$TRB
  
  # Save data to make histogram of TRA and TRB presence
  hist_data=data_to_add[,3:4]
  row.names(hist_data)=data_to_add$barcodes_full
  hist_data[which(!is.na(hist_data$TRA)),1]<-1 #binarize
  hist_data[which(!is.na(hist_data$TRB)),2]<-1
  hist_data[which(is.na(hist_data$TRA)),1]<-0
  hist_data[which(is.na(hist_data$TRB)),2]<-0
  hist_data$TRA=as.numeric(hist_data$TRA)
  hist_data$TRB=as.numeric(hist_data$TRB)
  
  TRA_TRB_Freq=as.data.frame(table(hist_data))$Freq
  names(TRA_TRB_Freq)=c("None", "TRA Only", "TRB Only", "TRA and TRB")
  M@misc[["TRA_TRB_Freq"]]<-TRA_TRB_Freq
  data_list[[expt]]<-TRA_TRB_Freq
  
  top_ten=clonotypes[order(clonotypes$frequency, decreasing=T),"cdr3s_aa"][1:10]
  top_ten_list[[expt]]<-top_ten
  
  # Color UMAP by clone size
  clonesize_table=unique(merged_TCR[,c(1,32)])
  clonesize_table=left_join(data_to_add, clonesize_table, by=c("barcodes" = "barcode"))
  M@meta.data[["clonesize"]]<-clonesize_table$frequency
  clones_table=unique(merged_TCR[,c(1,34)])
  clones_table=left_join(data_to_add, clones_table, by=c("barcodes" = "barcode"))
  M@meta.data[["clones"]]<-clones_table$cdr3s_aa
  clones_top_ten<-M@meta.data[["clones"]]
  clones_top_ten[which(!clones_top_ten %in% top_ten)]<-NA
  M@meta.data[["clones_top_ten"]]<-clones_top_ten
}

names(data_list)=expt_list
data_matrix=do.call("rbind", data_list)
data_matrix=data_matrix[order(row.names(data_matrix)), ]

# Without the "None" group
pdf(file = "/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/barplot_TCRonly.pdf", width = 10, height = 8.5)
par(mar=c(10, 6, 4,4))
barplot(height = t(data_matrix)[2:4,], 
        beside = TRUE, 
        col = hue_pal()(3), 
        las=2,
        ylab="Barcodes",
        ylim=c(0,300))
legend("topright", legend=colnames(data_matrix)[2:4], fill = hue_pal()(3))
dev.off()

write.table(t(data_matrix)[2:4,], "/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/barplot_TCRonly.txt", sep="\t", quote=F)


#################
# Violin Plots

# Add new classifications to metadata
normal=c("ALP-00036", "ALP-00041-TXP", "ALP-00042-TXP", "ALP-00045-TXP", "ALP-00066", "ALP-00078")
notACR=c("ALP-00042-BX2", "ALP-00045-BX1")
indexRejection=c("ALP-0003-BX1", "ALP-00017", "ALP-00023-BX1", "ALP-00023-BX2", "ALP-00036-BX2")
recurrent=c("ALP-0003-BX2",  "ALP-00023-BX3", "ALP-00036-BX3", "ALP-00036-BX4", "ALP-00036-BX5")
resolved=c("ALP-0003-BX3", "ALP-0003-BX4", "ALP-0003-BX5", "ALP-00017-BX2")

comb=c(normal, notACR, indexRejection, recurrent, resolved)

Idents(object = M_T) <- "orig.ident"
M_T=subset(x = M_T, idents = comb)

temp=M_T@meta.data$orig.ident
temp[temp %in% normal]<-"Pre-TXP"
temp[temp %in% notACR]<-"No ACR"
temp[temp %in% indexRejection]<-"Index Rejection"
temp[temp %in% recurrent]<-"Recurrent"
temp[temp %in% resolved]<-"Resolved"
M_T@meta.data[["Classification"]]<-temp

Idents(object = M_T) <- "Classification"
Idents(object = M_T) <- factor(x = Idents(M_T), levels = c("Pre-TXP", "No ACR", "Index Rejection", "Recurrent", "Resolved"))

# Activated
myGenes_activated=unique(c("PRF1", "GZMB", "GZMK", "IFNG", "HLA-DRA", "CX3CR1", "TBX21", "MKI67", "PCNA"))
# Resident/Circulation
myGenes_resident=unique(c("ZNF683", "PRDM1", "CD69", "ITGAE", "CXCR6", "S1PR1", "S1PR5", "SELL", "TOX2"))
# Exhaustion
myGenes_exhausted=unique(c("TOX", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "NR4A1", "NKG7"))

expanded=rep("Not Expanded", nrow(M_T@meta.data))
expanded[which(M_T$clonesize>2)]<-"Expanded"
M_T[["expanded"]] <-expanded
table(M_T@meta.data$expanded)

# Reorder levels for split plot by expanded status
M_T@meta.data[["cluster_names_new"]]<-factor(M_T@meta.data[["cluster_names_new"]], levels=c("CD4 Naïve T", "CD4 Activated T","CD4 TCM","MAIT","CD8 Effector Memory T","CD8 Effector T","CD8 NK-like T", "Gamma Delta T", "Cycling T"))
colors_order=myColors[c("CD4 Naïve T", "CD4 Activated T","CD4 TCM","MAIT","CD8 Effector Memory T","CD8 Effector T","CD8 NK-like T", "Gamma Delta T", "Cycling T"),]
pdf(file = "Split_Tcell_by_expanded_clusterColors.pdf", width = 17, height = 8)
par(mar=c(2, 2, 2, 2))
  DimPlot(M_T, reduction="tsne", group.by = "cluster_names_new", split.by = "expanded", cols=colors_order$Cluster_Color, pt.size=1) + ggtitle(label=NULL)
dev.off()

# Make stacked barplot for split expanded status
data_long=as.data.frame(cbind(Expanded=M_T@meta.data[["expanded"]], Cluster=as.character(M_T@meta.data[["cluster_names_new"]])))
data=as.data.frame.matrix(table(data_long))

# Read and wrangle data
data=t(data)
data <- sweep(data, 2, colSums(data), "/")*100

myColors=myColors[order(myColors$Cluster_Name, decreasing=F),]
my_colors=myColors$Cluster_Color
data=as.data.frame(cbind(data, colors=my_colors))
my_colors=data$colors
data=data[,-3]

write.table(data, "/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellTypeFrequencies_T_cell_byexpanded_new_colors.txt", sep="\t", quote=F)


# Make plot
pdf("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellTypeFrequencies_T_cell_byexpanded_new_colors.pdf", width = 9, height = 9)
par(mar = c(8,4,4,16), xpd = T)
barplot(as.matrix(data[nrow(data):1,]), col = rev(my_colors), xaxt = "n", ylab = "Population frequency (%)", border = NA)
axis(side = 1, at = seq(1,ncol(data))*1.2-0.5, labels = colnames(data), las = 2)
legend(x = ncol(data)*1.2+0.5, y = 100, legend = row.names(data), fill = my_colors, bty = "n", border = NA)
dev.off()

pdf(file = "Tcell_sidebyside_activated_RNA_Violin_newClassifications_justSerial.pdf", width = 8, height = 9)
par(mar=c(2, 2, 2, 2))
  VlnPlot(M_T, assay="RNA", group.by="Classification", split.by="expanded", stack=T, features=myGenes_activated, flip=T, sort=F)
dev.off()

pdf(file = "Tcell_sidebyside_resident_RNA_Violin_newClassifications_justSerial.pdf", width = 8, height = 9)
par(mar=c(2, 2, 2, 2))
  VlnPlot(M_T, assay="RNA", group.by="Classification", split.by="expanded", stack=T, features=myGenes_resident, flip=T, sort=F)
dev.off()

pdf(file = "Tcell_sidebyside_exhausted_RNA_Violin_newClassifications_justSerial.pdf", width = 8, height = 9)
par(mar=c(2, 2, 2, 2))
  VlnPlot(M_T, assay="RNA", group.by="Classification", split.by="expanded", stack=T, features=myGenes_exhausted, flip=T, sort=F)
dev.off()

myGenes_new=c("GZMB", "HLA-DRA", "LAG3",  "TIGIT", "NKG7", "CD69", "CXCR6")
pdf(file = "Tcell_sidebyside_new_RNA_Violin_newClassifications_justSerial.pdf", width = 8, height = 9)
par(mar=c(2, 2, 2, 2))
  VlnPlot(M_T, assay="RNA", group.by="Classification", split.by="expanded", stack=T, features=myGenes_new, flip=T, sort=F)
dev.off()

Idents(M_T)<-"cluster_names_new"
M_T_CD8<-subset(M_T, idents=c("MAIT","CD8 Effector Memory T","CD8 Effector T","CD8 NK-like T", "Gamma Delta T", "Cycling T"))
Idents(M_T_CD8)<-"Classification"

pdf(file = "Tcell_sidebyside_new_RNA_Violin_newClassifications_justSerial_CD8.pdf", width = 8, height = 9)
par(mar=c(2, 2, 2, 2))
  VlnPlot(M_T_CD8, assay="RNA", group.by="Classification", split.by="expanded", stack=T, features=myGenes_new, flip=T, sort=F)
dev.off()

# for anna
Idents(object = M_T) <- "expanded"
M_T_exp=subset(x = M_T, idents = "Expanded")
M_T_Notexp=subset(x = M_T, idents = "Not Expanded")

Notexp_inds_2TRA=grep(",", M_T_Notexp$TRA)
Notexp_inds_1TRB=setdiff(setdiff(1:length(M_T_Notexp$TRB), grep(",", M_T_Notexp$TRB)), which(is.na(M_T_Notexp$TRB)))
Notexp_meta=M_T_Notexp@meta.data[intersect(Notexp_inds_2TRA, Notexp_inds_1TRB),]
nrow(Notexp_meta) #116
nrow(M_T_Notexp@meta.data) #2834
116/2834 #4.1% (was 4.4)
sum(!is.na(M_T_Notexp@meta.data$clonesize)) #1542
116/1542 #7.5% (was 7.6)

Exp_inds_2TRA=grep(",", M_T_exp$TRA)
Exp_inds_1TRB=setdiff(setdiff(1:length(M_T_exp$TRB), grep(",", M_T_exp$TRB)), which(is.na(M_T_exp$TRB)))
Exp_meta=M_T_exp@meta.data[intersect(Exp_inds_2TRA, Exp_inds_1TRB),]
nrow(Exp_meta) #49
nrow(M_T_exp@meta.data) #514
49/514 #9.5% (was 10.0)
sum(!is.na(M_T_exp@meta.data$clonesize)) #514
49/514 #9.5% (was 10.0)

