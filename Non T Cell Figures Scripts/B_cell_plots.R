##############################
#                            #
#       B Cell Plots         #
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
setwd("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/")
M <- readRDS("Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")

# Subset to B-cells only
Idents(object = M) <- "cluster_names_new"
Clusters_Keep=c("Naïve B", "Transitional B")
M_B=subset(x = M, idents = Clusters_Keep)

# Make TSNE plot
M_B <- RunPCA(M_B, npcs = 30, verbose = FALSE)
M_B <- RunUMAP(M_B, reduction = "pca", dims = 1:30)
M_B <- RunTSNE(M_B, reduction = "pca", dims = 1:30)

# Pull in new colors
myColors=read.table("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Color_Guide.txt", sep = "\t", header=T, comment.char="*")
myColors[4,2]<-"Naïve B"
myColors=myColors[which(myColors$Cluster_Name %in% Clusters_Keep),]
row.names(myColors)=myColors$Cluster_Name
myColors=myColors[levels(M_B@active.ident),]

pdf(file = "TSNE_Bcell_Groups_New_Colors.pdf", width = 9.5, height = 8)
par(mar=c(2, 2, 2, 2))
DimPlot(M_B, reduction="tsne", label=T, repel=T, raster=F, cols=myColors$Cluster_Color, pt.size=1)
dev.off()

# Plot by ACR type
pdf(file = "TSNE_Bcell_ACRType.pdf", width = 9, height = 8)
par(mar=c(2, 2, 2, 2))
DimPlot(M_B, reduction="tsne", group.by="ACRType", label=F, repel=T, raster=F, pt.size=1)
dev.off()

# Plot by Clone size
pdf(file = "TSNE_Bcell_Clonesize.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
FeaturePlot(M_B, reduction="tsne", features="clonesize_B_Krish", raster=F, pt.size=1, cols = c("lightgrey", "red"), order=T)
dev.off()

# B cell mutation
pdf(file = "TSNE_Bcell_Mutation.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
FeaturePlot(M_B, reduction="tsne", features ="mutation_level", raster=F, pt.size=1, cols=c("lightgrey", "blue"), order=T)
dev.off()

# Plot by Clone size
pdf(file = "TSNE_Bcell_Clonesize_Bigger.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
FeaturePlot(M_B, reduction="tsne", features="clonesize_B_Krish", raster=F, pt.size=2, cols = c("lightgrey", "red"), order=T)
dev.off()

# B cell mutation
pdf(file = "TSNE_Bcell_Mutation_Bigger.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
FeaturePlot(M_B, reduction="tsne", features ="mutation_level", raster=F, pt.size=2, cols=c("lightgrey", "blue"), order=T)
dev.off()

# All Clusters clone size
pdf(file = "TSNE_Clonesize_B.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
FeaturePlot(M, reduction="tsne", features ="clonesize_B_Krish", raster=F, cols=c("lightgrey", "red"))
dev.off()

# All Clusters mutation
pdf(file = "TSNE_Mutation.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
FeaturePlot(M, reduction="tsne", features ="mutation_level", raster=F, cols=c("lightgrey", "blue"))
dev.off()

# Prop plot by ACR type
data_long=as.data.frame(cbind(ACRType=M_B@meta.data[["ACRType"]], Cluster=M_B@meta.data[["cluster_names_new"]]))
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
pdf("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellTypeFrequencies_B_cell_byACRType_new_colors.pdf", width = 9, height = 9)
par(mar = c(8,4,4,16), xpd = T)
barplot(as.matrix(data[nrow(data):1,]), col = rev(my_colors), xaxt = "n", ylab = "Population frequency (%)", border = NA)
axis(side = 1, at = seq(1,ncol(data))*1.2-0.5, labels = colnames(data), las = 2)
legend(x = ncol(data)*1.2+0.5, y = 100, legend = row.names(data), fill = my_colors, bty = "n", border = NA)
dev.off()

write.table(data, "/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellTypeFrequencies_B_cell_byACRType_new_colors.txt", sep="\t", quote=F)


# Make plots for BCR chain usage
# IGL IGK IGH frequncy barplot
data_list=NULL
top_ten_list=NULL
expt_list=unique(M@meta.data[["orig.ident"]])
for(expt in 1:length(expt_list)){
  
  # Set up path and read in Seurat object
  name <- expt_list[expt]
  path <- paste0("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/", name)
  M_temp <- subset(x = M, subset = orig.ident == name)
  
  # Pull in BCR data files
  clonotypes <- read.csv(paste0("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/", name, "/clonotypes.csv"), sep=',', header=T)
  filtered_contig <- read.csv(paste0("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/", name, "/filtered_contig_annotations.csv"), sep=',', header=T)
  merged_BCR=left_join(filtered_contig, clonotypes, by=c("raw_clonotype_id" = "clonotype_id"))
  
  # Add full BCR data to Seurat object
  M_temp@misc[["BCR_full_table"]]<-merged_BCR
  
  # Add BCR data to each cell
  AB_table=as.data.frame(matrix(nrow=length(unique(merged_BCR$barcode)),ncol=4))
  colnames(AB_table)=c("sample" , "IGL", "IGK", "IGH")
  row.names(AB_table)=unique(merged_BCR$barcode)
  for(i in unique(merged_BCR$barcode)){
    temp=merged_BCR[which(merged_BCR$barcode==i),] #subset table to just unique cell barcode i
    temp=temp[,c(1,3,6,23,28)] #reduce to columns we care about, for ease of visualization when coding
    IGL=temp[which(temp$chain=="IGL"),] #pull out column(s) that are IGL
    IGK=temp[which(temp$chain=="IGK"),] #pull out column(s) that are IGK
    IGH=temp[which(temp$chain=="IGH"),] #pull out column(s) that are IGH
    AB_table[i,1]=name
    AB_table[i,2]=paste(IGL$cdr3, collapse=",")
    AB_table[i,3]=paste(IGK$cdr3, collapse=",")
    AB_table[i,4]=paste(IGH$cdr3, collapse=",")
  }
  AB_table=as.data.frame(cbind(barcodes=row.names(AB_table), AB_table))
  data_to_add=as.data.frame(cbind(barcodes_full=M_temp@assays[["RNA"]]@data@Dimnames[[2]], barcodes=M_temp@assays[["RNA"]]@data@Dimnames[[2]], samples=M_temp@meta.data$orig.ident))
  data_to_add[,2]=gsub("_.*", "", data_to_add[,2])
  data_to_add=left_join(data_to_add, AB_table, by=c("barcodes"="barcodes", "samples"="sample"))
  data_to_add[which(data_to_add[,4] == ""),4]<-NA #deal with blank spaces from the above paste
  data_to_add[which(data_to_add[,5] == ""),5]<-NA
  data_to_add[which(data_to_add[,6] == ""),6]<-NA
  M_temp@meta.data[["IGL"]]<-data_to_add$IGL
  M_temp@meta.data[["IGK"]]<-data_to_add$IGK
  M_temp@meta.data[["IGH"]]<-data_to_add$IGH
  
  # Save data to make histogram of IGL, IGK, IGH presence
  hist_data=data_to_add[,4:6]
  row.names(hist_data)=data_to_add$barcodes_full
  hist_data[which(!is.na(hist_data$IGL)),1]<-1 #binarize
  hist_data[which(!is.na(hist_data$IGK)),2]<-1
  hist_data[which(!is.na(hist_data$IGH)),3]<-1
  hist_data[which(is.na(hist_data$IGL)),1]<-0
  hist_data[which(is.na(hist_data$IGK)),2]<-0
  hist_data[which(is.na(hist_data$IGH)),3]<-0
  hist_data$IGL=as.numeric(hist_data$IGL)
  hist_data$IGK=as.numeric(hist_data$IGK)
  hist_data$IGH=as.numeric(hist_data$IGH)
  
  
  IGL_IGK_IGH_Freq=as.data.frame(table(hist_data))$Freq
  if(sum(colSums(hist_data))>0){
    name_list=NULL
    vals=c("IGL", "IGK", "IGH")
    for(rrow in 1:nrow(as.data.frame(table(hist_data)))){
      a=as.data.frame(table(hist_data))[rrow,]
      temp=as.numeric(lapply(a[1:3], as.character))
      if(sum(temp)==0){
        name_list=c(name_list, "None")
      }else{
        temp=paste(vals[which(temp != 0)], collapse="+")
        name_list=c(name_list, temp)
      }
    }
    names(IGL_IGK_IGH_Freq)=name_list
    M_temp@misc[["IGL_IGK_IGH_Freq"]]<-IGL_IGK_IGH_Freq
    data_list[[expt]]<-IGL_IGK_IGH_Freq
  }else{
    data_list[[expt]]<-IGL_IGK_IGH_Freq
    names(data_list[[expt]])<-"None"
  }
  
  top_ten=clonotypes[order(clonotypes$frequency, decreasing=T),"cdr3s_aa"][1:10]
  top_ten_list[[expt]]<-top_ten
  
  # Color UMAP by clone size
  clonesize_table=unique(merged_BCR[,c(1,32)])
  clonesize_table=left_join(data_to_add, clonesize_table, by=c("barcodes" = "barcode"))
  M_temp@meta.data[["clonesize_B"]]<-clonesize_table$frequency
  clones_table=unique(merged_BCR[,c(1,34)])
  clones_table=left_join(data_to_add, clones_table, by=c("barcodes" = "barcode"))
  M_temp@meta.data[["clones_B"]]<-clones_table$cdr3s_aa
  clones_top_ten<-M_temp@meta.data[["clones_B"]]
  clones_top_ten[which(!clones_top_ten %in% top_ten)]<-NA
  M_temp@meta.data[["clones_top_ten_B"]]<-clones_top_ten
  
}

names(data_list)=expt_list
data_matrix=do.call(rbind, lapply(data_list, function(x) x[match(names(data_list[[1]]), names(x))]))
data_matrix[which(is.na(data_matrix))]<-0

# Without the "None" group
pdf(file = "/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/barplot_BCRonly.pdf", width = 10, height = 8.5)
par(mar=c(10, 6, 4,4))
barplot(height = t(data_matrix)[2:8,], 
        beside = TRUE, 
        col = hue_pal()(7), 
        las=2,
        ylab="Barcodes",
        ylim=c(0,150))
legend("topright", legend=colnames(data_matrix)[2:8], fill = hue_pal()(7))
dev.off()

#K/L usage only (from cell ranger analysis)
data_matrix_IGLIGK=cbind(IGK=rowSums(data_matrix[,grep("IGK", colnames(data_matrix))]), IGL=rowSums(data_matrix[,grep("IGL", colnames(data_matrix))]))
pdf(file = "/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/barplot_BCRonly_KL.pdf", width = 10, height = 8.5)
par(mar=c(10, 6, 4,4))
barplot(height = t(data_matrix_IGLIGK), 
        beside = TRUE, 
        col = hue_pal()(2), 
        las=2,
        ylab="Barcodes",
        ylim=c(0,150))
legend("topright", legend=colnames(data_matrix_IGLIGK), fill = hue_pal()(2))
dev.off()

write.table(t(data_matrix_IGLIGK), "/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/barplot_BCRonly_KL.txt", sep="\t", quote=F)

# Heavy chain usage (from Krish analysis)
data=read.table("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/all_bcr_mut_updated_namecorrection.csv", sep=",", header=T)
data[,3]=paste(data[,1], data[,2], sep="-")
data[,3]=gsub("-$", "", data[,3])
data=data[,c(3,5)]
data=data[data$source %in% expt_list,]
data_matrix=data.frame(unclass(table(data)))
colnames(data_matrix)[1]<-"Unknown"

# Without the "Unknown" group
pdf(file = "/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/barplot_BCRonly_Heavy.pdf", width = 10, height = 8.5)
par(mar=c(10, 6, 4,4))
barplot(height = t(data_matrix)[2:15,], 
        beside = TRUE, 
        col = hue_pal()(14), 
        las=2,
        ylab="Barcodes",
        ylim=c(0,1000))
legend("topright", legend=colnames(data_matrix)[2:15], fill = hue_pal()(14))
dev.off()

# Heavy chain usage (from Krish analysis) — reduced to IGC, IGD, IGG, and IGM level only
data=read.table("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/all_bcr_mut_updated_namecorrection.csv", sep=",", header=T)
data[,3]=paste(data[,1], data[,2], sep="-")
data[,3]=gsub("-$", "", data[,3])
data=data[,c(3,5)]
data=data[data$source %in% expt_list,]
data[,2]=substr(gsub("HC", "", data[,2]), start = 1, stop = 3)
data[,2]=gsub("IG", "Ig", data[,2])
data_matrix=data.frame(unclass(table(data)))
colnames(data_matrix)[1]<-"Unknown"

# Without the "Unknown" group
pdf(file = "/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/barplot_BCRonly_Heavy_Redu.pdf", width = 10, height = 8.5)
par(mar=c(10, 6, 4,4))
barplot(height = t(data_matrix)[2:5,], 
        beside = TRUE, 
        col = hue_pal()(4), 
        las=2,
        ylab="Barcodes",
        ylim=c(0,1000))
legend("topright", legend=colnames(data_matrix)[2:5], fill = hue_pal()(4))
dev.off()

write.table(t(data_matrix)[2:5,], "/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/barplot_BCRonly_Heavy_Redu.txt", sep="\t", quote=F)


# Normalize data using log for visualization purposes
M_B <- NormalizeData(M_B, assay="RNA")
all.genes <- rownames(M_B)
M_B <- ScaleData(M_B, features = all.genes, assay="RNA")

myGenes=c("FCER2", "CD1C", "HLA-DRB1", "HLA-DRA", "CD72", "CD5", "CD9")

pdf(file = "Bcell_Manual_RNA_Stacked_Violin_Names.pdf", width = 6, height = 6)
par(mar=c(2, 2, 2, 2))
VlnPlot(M_B, assay="RNA", group.by="cluster_names_new", stack=T, features=myGenes, flip=T, sort=T) + NoLegend()
dev.off()

# Reorder levels for dotplot
M_B@meta.data[["cluster_names_new"]]<-factor(M_B@meta.data[["cluster_names_new"]], levels=c("Naïve B", "Transitional B"))

pdf(file = "Bcell_Manual_RNA_DotPlot_Names_log.pdf", width = 8, height = 4)
par(mar=c(4, 4, 4, 4))
DotPlot(object=M_B, assay="RNA", features = myGenes, dot.scale=10, scale=F) + labs(y =NULL, x=NULL) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text.y = element_text(angle = 90, vjust = 1, hjust=1))
dev.off()

p1<-DotPlot(object=M_B, assay="RNA", features = myGenes, dot.scale=10, scale=F) + labs(y =NULL, x=NULL) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text.y = element_text(angle = 90, vjust = 1, hjust=1))
write.table(p1[["data"]], "/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/B_Manual_RNA_DotPlot_Names_log.txt", sep="\t", quote=F)


#myGenes2=c("CD19", "CD24", "CD38", "CD27", "PTPRC", "CD44", "A4GALT", "MME", "IL10", "IGHM", "IGHD", "IGHA1", "IGHG1", "MS4A1")
myGenes2=c("MS4A1", "CD19", "PTPRC", "IGHM", "IGHD", "IGHA1", "IGHG1", "CD24", "CD38", "CD27", "FCER2", "CD1C", "CD72", "CD5", "CD9", "A4GALT")
pdf(file = "Bcell_Manual2_RNA_DotPlot_Names_log.pdf", width = 8, height = 4)
par(mar=c(4, 4, 4, 4))
DotPlot(object=M_B, assay="RNA", features = myGenes2, dot.scale=10, scale=F) + labs(y =NULL, x=NULL)+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5))
dev.off()

