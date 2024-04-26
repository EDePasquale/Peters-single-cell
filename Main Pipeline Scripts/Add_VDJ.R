##############################
#                            #
#        Add VDJ Data        #
# Peters - Single Cell Liver #
#     Erica DePasquale       #
#        03 Oct 2023         #
#                            #
##############################

# Load libraries
library(Seurat)
library(dplyr)
library(plyr)
library(scales)
library(pheatmap)
library(ggplot2)

# Read in new Seurat object
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23")
M<-readRDS("Seurat_Liver_30_subcluster_names_meta.rds")
expt_list=unique(M@meta.data$orig.ident)

#################
#               #
#      TCR      #
#               #
#################

# Make empty data.frames and lists to use in the loop
AB_table_full=as.data.frame(matrix(nrow=1, ncol=4))
colnames(AB_table_full)=c("barcodes", "sample", "TRA", "TRB")

merged_TCR_full=as.data.frame(matrix(nrow=1, ncol=38))
colnames(merged_TCR_full)=c("sample", "barcode",               "is_cell",               "contig_id",             "high_confidence",      
                            "length",                "chain",                 "v_gene",                "d_gene",               
                            "j_gene",                "c_gene",                "full_length",           "productive",           
                            "fwr1",                  "fwr1_nt",               "cdr1",                  "cdr1_nt",              
                            "fwr2",                  "fwr2_nt",               "cdr2",                  "cdr2_nt",              
                            "fwr3",                  "fwr3_nt",               "cdr3",                  "cdr3_nt",              
                            "fwr4",                  "fwr4_nt",               "reads",                 "umis",                 
                            "raw_clonotype_id",      "raw_consensus_id",      "exact_subclonotype_id", "frequency",            
                            "proportion",            "cdr3s_aa",              "cdr3s_nt",              "inkt_evidence",        
                            "mait_evidence")
data_list=NULL
top_ten_list=NULL

# Loop through each sample to add data to the shared data frame and make individual plots
for(expt in 1:length(expt_list)){
  
  ###############
  # FULL OBJECT #
  ###############
  
  # Set up path and read in Seurat object
  name <- expt_list[expt]
  path <- paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/", name)
  
  # Pull in TCR data files
  clonotypes <- read.csv(paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/", name, "/clonotypes.csv"), sep=',', header=T)
  filtered_contig <- read.csv(paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/", name, "/filtered_contig_annotations.csv"), sep=',', header=T)
  merged_TCR=left_join(filtered_contig, clonotypes, by=c("raw_clonotype_id" = "clonotype_id"))
  merged_TCR_full=rbind(merged_TCR_full, cbind(sample=rep(name, nrow(merged_TCR)),merged_TCR))
  
  # Add TCR data to each cell
  AB_table=as.data.frame(matrix(nrow=length(unique(merged_TCR$barcode)),ncol=3))
  colnames(AB_table)=c("sample", "TRA", "TRB")
  row.names(AB_table)=unique(merged_TCR$barcode)
  for(i in unique(merged_TCR$barcode)){
    temp=merged_TCR[which(merged_TCR$barcode==i),] #subset table to just unique cell barcode i
    temp=temp[,c(1,3,6,23,28)] #reduce to columns we care about, for ease of visualization when coding
    TRA=temp[which(temp$chain=="TRA"),] #pull out column(s) that are TRA
    TRB=temp[which(temp$chain=="TRB"),] #pull out column(s) that are TRB
    AB_table[i,1]=name
    AB_table[i,2]=paste(TRA$cdr3, collapse=",")
    AB_table[i,3]=paste(TRB$cdr3, collapse=",")
  }
  AB_table=as.data.frame(cbind(barcodes=row.names(AB_table), AB_table))
  AB_table_full=rbind(AB_table_full, AB_table)
  
  #####################
  # INDIVIDUAL OBJECT #
  #####################
  M_sub=subset(x = M, subset = orig.ident == name)
  data_to_add=as.data.frame(cbind(barcodes_full=M_sub@assays[["RNA"]]@data@Dimnames[[2]], barcodes=M_sub@assays[["RNA"]]@data@Dimnames[[2]]))
  data_to_add[,2]=gsub("_.*", "", data_to_add[,2])
  data_to_add=left_join(data_to_add, AB_table)
  data_to_add[which(data_to_add[,3] == ""),3]<-NA #deal with blank spaces from the above paste
  data_to_add[which(data_to_add[,4] == ""),4]<-NA
  data_to_add[which(data_to_add[,5] == ""),5]<-NA
  M_sub@meta.data[["TRA"]]<-data_to_add$TRA
  M_sub@meta.data[["TRB"]]<-data_to_add$TRB
  
  # Save data to make histogram of TRA and TRB presence
  hist_data=data_to_add[,4:5]
  row.names(hist_data)=data_to_add$barcodes_full
  hist_data[which(!is.na(hist_data$TRA)),1]<-1 #binarize
  hist_data[which(!is.na(hist_data$TRB)),2]<-1
  hist_data[which(is.na(hist_data$TRA)),1]<-0
  hist_data[which(is.na(hist_data$TRB)),2]<-0
  hist_data$TRA=as.numeric(hist_data$TRA)
  hist_data$TRB=as.numeric(hist_data$TRB)
  
  TRA_TRB_Freq=as.data.frame(table(hist_data))$Freq
  names(TRA_TRB_Freq)=c("None", "TRA Only", "TRB Only", "TRA and TRB")
  M_sub@misc[["TRA_TRB_Freq"]]<-TRA_TRB_Freq
  data_list[[expt]]<-TRA_TRB_Freq
  
  top_ten=clonotypes[order(clonotypes$frequency, decreasing=T),"cdr3s_aa"][1:10]
  top_ten_list[[expt]]<-top_ten
  
  #Barcode Frequency plot
  pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/",  name, "/", name, "_barcode_frequency.pdf"), width = 11, height = 8.5)
  par(mar=c(6,6, 4,4))
    barplot(sort(clonotypes$frequency, decreasing=T), ylab="Barcode Frequency", xlab="Clonotypes", col="darkblue", main=name)
    legend("topright", top_ten)
  dev.off()

  #Barcode Frequency plot, scaled
  pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/",  name, "/", name, "_barcode_frequency_scaled.pdf"), width = 11, height = 8.5)
  par(mar=c(6,6, 4,4))
    barplot(sort(clonotypes$frequency, decreasing=T), ylab="Barcode Frequency", xlab="Clonotypes", col="darkblue", ylim=c(0,35), main=name)
    legend("topright", top_ten)
  dev.off()
  
  # Color UMAP by clone size
  clonesize_table=unique(merged_TCR[,c(1,32)])
  clonesize_table=left_join(data_to_add, clonesize_table, by=c("barcodes" = "barcode"))
  M_sub@meta.data[["clonesize"]]<-clonesize_table$frequency
  clones_table=unique(merged_TCR[,c(1,34)])
  clones_table=left_join(data_to_add, clones_table, by=c("barcodes" = "barcode"))
  M_sub@meta.data[["clones"]]<-clones_table$cdr3s_aa
  clones_top_ten<-M_sub@meta.data[["clones"]]
  clones_top_ten[which(!clones_top_ten %in% top_ten)]<-NA
  M_sub@meta.data[["clones_top_ten"]]<-clones_top_ten
  
  pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/",  name, "/", name, "_UMAP_clonesize.pdf"), width = 8, height = 8)
  par(mar=c(2, 2, 2, 2))
    print(FeaturePlot(object = M_sub, reduction = "umap", pt.size=0.5, features="clonesize", cols=c("lightgrey", "red")))
  dev.off()

  # Color UMAP by top 10 clones
  pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/",  name, "/", name, "_UMAP_clones.pdf"), width = 8, height = 11)
  par(mar=c(2, 2, 2, 2))
    print(DimPlot(object = M_sub, reduction = "umap", pt.size=0.5, group.by="clones_top_ten", order = rev(top_ten)) +
    theme(legend.position="bottom") +
    guides(color=guide_legend(ncol=1)))
  dev.off()

  saveRDS(M_sub, file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/QC/", name, "/", name, "_filtered_matrices/Seurat.rds"))

}

# Modify AB_table and merged_TCR table to remove duplicate barcodes
AB_table_full2=AB_table_full[!duplicated(AB_table_full),] #remove duplicate rows
AB_table_full3=AB_table_full2[ !duplicated(AB_table_full2$barcodes), ] #take the first row where there are cell name conflicts (only 29 cells)

merged_TCR_full2=merged_TCR_full[!duplicated(merged_TCR_full),] #remove duplicate rows
merged_TCR_full3=merged_TCR_full2[ !duplicated(merged_TCR_full2$barcode), ] #take the first row where there are cell name conflicts (only 29 cells)

# Add AB_table data to Seurat object
data_to_add=as.data.frame(cbind(barcodes_full=M@assays[["SCT"]]@data@Dimnames[[2]], barcodes=M@assays[["SCT"]]@data@Dimnames[[2]], samples=M@meta.data$orig.ident))
data_to_add[,2]=gsub("_.*", "", data_to_add[,2])
data_to_add=left_join(data_to_add, AB_table_full3, by=c("barcodes"="barcodes", "samples"="sample"))
data_to_add[which(data_to_add[,4] == ""),4]<-NA #deal with blank spaces from the above paste
data_to_add[which(data_to_add[,5] == ""),5]<-NA
M@meta.data[["TRA"]]<-data_to_add$TRA
M@meta.data[["TRB"]]<-data_to_add$TRB

# Add clone size data to Seurat object
clonesize_table=unique(merged_TCR_full3[,c(1,2,33)])
clonesize_table=left_join(data_to_add, clonesize_table, by=c("barcodes" = "barcode", "samples"="sample"))
M@meta.data[["clonesize"]]<-clonesize_table$frequency

# Make UMAP and TSNE plots
pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/UMAP_clonesize.pdf"), width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
  print(FeaturePlot(object = M, reduction = "umap", pt.size=0.3, features="clonesize", cols=c("lightgrey", "red")))
dev.off()

pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/TSNE_clonesize.pdf"), width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
  print(FeaturePlot(object = M, reduction = "tsne", pt.size=0.3, features="clonesize", cols=c("lightgrey", "red")))
dev.off()

# Make plots for TRA and TRB distribution
names(data_list)=expt_list
data_matrix=do.call("rbind", data_list)

# With the "None" group
pdf(file = "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/barplot.pdf", width = 10, height = 8.5)
par(mar=c(10, 6, 4,4))
barplot(height = t(data_matrix), 
        beside = TRUE, 
        col = c("grey", hue_pal()(3)), 
        las=2,
        ylab="Barcodes")
legend("topright", legend=colnames(data_matrix), fill = c("grey", hue_pal()(3)))
dev.off()

# Without the "None" group
pdf(file = "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/barplot_TCRonly.pdf", width = 10, height = 8.5)
par(mar=c(10, 6, 4,4))
barplot(height = t(data_matrix)[2:4,], 
        beside = TRUE, 
        col = hue_pal()(3), 
        las=2,
        ylab="Barcodes",
        ylim=c(0,300))
legend("topright", legend=colnames(data_matrix)[2:4], fill = hue_pal()(3))
dev.off()

# Make plots for shared clones
names(top_ten_list)=expt_list
top_ten_matrix=do.call("rbind", top_ten_list)
unique_clones=unique(as.character(top_ten_matrix))
top_ten_matrix_wide=as.data.frame(matrix(nrow=30, ncol=length(unique_clones)))
row.names(top_ten_matrix_wide)=expt_list
colnames(top_ten_matrix_wide)=unique_clones
for(i in 1:length(expt_list)){
  top_ten_matrix_wide[i,top_ten_matrix[i,]]<-1
}
top_ten_matrix_wide[is.na(top_ten_matrix_wide)] <- 0

# Reorder to group samples
top_ten_matrix_wide=top_ten_matrix_wide[c(1,2,3,4,5,6,7,23,8,9,10,11,12,13,14,24,15,16,17,25,26,18,27,28,29,30,19,20,21,22),]

pheatmap(top_ten_matrix_wide,
         cluster_rows = F,
         color=c("gray92", "red"),
         gaps_row=c(5,6,8,11,13,16,21,22,24,26,27,28,29),
         width=40,
         height = 11.5,
         filename = "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/heatmap_top_ten_clones.pdf")

# Save new RDS file
saveRDS(M, file = "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Seurat_Liver_30_subcluster_names_meta_TCR.rds")


#################
#               #
#      BCR      #
#               #
#################

# Make empty data.frames and lists to use in the loop
AB_table_full=as.data.frame(matrix(nrow=1, ncol=5))
colnames(AB_table_full)=c("barcodes", "sample", "IGL", "IGK", "IGH")

merged_BCR_full=as.data.frame(matrix(nrow=1, ncol=36))
colnames(merged_BCR_full)=c("sample", "barcode",               "is_cell",               "contig_id",             "high_confidence",       "length",                "chain",                
                            "v_gene",                "d_gene",                "j_gene",                "c_gene",                "full_length",           "productive",           
                            "fwr1",                  "fwr1_nt",               "cdr1",                  "cdr1_nt",               "fwr2",                  "fwr2_nt",              
                            "cdr2",                  "cdr2_nt",               "fwr3",                  "fwr3_nt",               "cdr3",                  "cdr3_nt",              
                            "fwr4",                  "fwr4_nt",               "reads",                 "umis",                  "raw_clonotype_id",      "raw_consensus_id",     
                            "exact_subclonotype_id", "frequency",             "proportion",            "cdr3s_aa",              "cdr3s_nt") 
data_list=NULL
top_ten_list=NULL

# Loop through each sample to add data to the shared data frame and make individual plots
for(expt in 1:length(expt_list)){
  
  ###############
  # FULL OBJECT #
  ###############
  
  # Set up path and read in Seurat object
  name <- expt_list[expt]
  path <- paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/", name)
  
  # Pull in BCR data files
  clonotypes <- read.csv(paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/", name, "/clonotypes.csv"), sep=',', header=T)
  filtered_contig <- read.csv(paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/", name, "/filtered_contig_annotations.csv"), sep=',', header=T)
  merged_BCR=left_join(filtered_contig, clonotypes, by=c("raw_clonotype_id" = "clonotype_id"))
  merged_BCR_full=rbind(merged_BCR_full, cbind(sample=rep(name, nrow(merged_BCR)),merged_BCR))
  
  # Add BCR data to each cell
  AB_table=as.data.frame(matrix(nrow=length(unique(merged_BCR$barcode)),ncol=4))
  colnames(AB_table)=c("sample", "IGL", "IGK", "IGH")
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
  AB_table_full=rbind(AB_table_full, AB_table)
  
  #####################
  # INDIVIDUAL OBJECT #
  #####################
  
  M_sub=subset(x = M, subset = orig.ident == name)
  data_to_add=as.data.frame(cbind(barcodes_full=M_sub@assays[["RNA"]]@data@Dimnames[[2]], barcodes=M_sub@assays[["RNA"]]@data@Dimnames[[2]], samples=M_sub@meta.data$orig.ident))
  data_to_add[,2]=gsub("_.*", "", data_to_add[,2])
  data_to_add=left_join(data_to_add, AB_table, by=c("barcodes"="barcodes", "samples"="sample"))
  data_to_add[which(data_to_add[,4] == ""),4]<-NA #deal with blank spaces from the above paste
  data_to_add[which(data_to_add[,5] == ""),5]<-NA
  data_to_add[which(data_to_add[,6] == ""),6]<-NA
  M_sub@meta.data[["IGL"]]<-data_to_add$IGL
  M_sub@meta.data[["IGK"]]<-data_to_add$IGK
  M_sub@meta.data[["IGH"]]<-data_to_add$IGH
  
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
    M_sub@misc[["IGL_IGK_IGH_Freq"]]<-IGL_IGK_IGH_Freq
    data_list[[expt]]<-IGL_IGK_IGH_Freq
  }else{
    data_list[[expt]]<-IGL_IGK_IGH_Freq
    names(data_list[[expt]])<-"None"
  }
  
  top_ten=clonotypes[order(clonotypes$frequency, decreasing=T),"cdr3s_aa"][1:10]
  top_ten_list[[expt]]<-top_ten
  
  #Barcode Frequency plot
  pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/",  name, "/", name, "_barcode_frequency.pdf"), width = 11, height = 8.5)
  par(mar=c(6,6, 4,4))
    barplot(sort(clonotypes$frequency, decreasing=T), ylab="Barcode Frequency", xlab="Clonotypes", col="darkblue", main=name)
    legend("topright", top_ten)
  dev.off()
  
  #Barcode Frequency plot, scaled
  pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/",  name, "/", name, "_barcode_frequency_scaled.pdf"), width = 11, height = 8.5)
  par(mar=c(6,6, 4,4))
    barplot(sort(clonotypes$frequency, decreasing=T), ylab="Barcode Frequency", xlab="Clonotypes", col="darkblue", ylim=c(0,120), main=name)
    legend("topright", top_ten)
  dev.off()
  
  # Color UMAP by clone size
  clonesize_table=unique(merged_BCR[,c(1,32)])
  clonesize_table=left_join(data_to_add, clonesize_table, by=c("barcodes" = "barcode"))
  M_sub@meta.data[["clonesize_B"]]<-clonesize_table$frequency
  clones_table=unique(merged_BCR[,c(1,34)])
  clones_table=left_join(data_to_add, clones_table, by=c("barcodes" = "barcode"))
  M_sub@meta.data[["clones_B"]]<-clones_table$cdr3s_aa
  clones_top_ten<-M_sub@meta.data[["clones_B"]]
  clones_top_ten[which(!clones_top_ten %in% top_ten)]<-NA
  M_sub@meta.data[["clones_top_ten_B"]]<-clones_top_ten
  
  if(sum(colSums(hist_data))>0){
    pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/",  name, "/", name, "_UMAP_clonesize.pdf"), width = 8, height = 8)
    par(mar=c(2, 2, 2, 2))
      print(FeaturePlot(object = M_sub, reduction = "umap", pt.size=0.5, features="clonesize_B", cols=c("lightgrey", "red")))
    dev.off()
  }
  
  saveRDS(M_sub, file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/QC/", name, "/", name, "_filtered_matrices/Seurat.rds"))
  
}

# Modify AB_table and merged_TCR table to remove duplicate barcodes
AB_table_full2=AB_table_full[!duplicated(AB_table_full),] #remove duplicate rows
AB_table_full3=AB_table_full2[ !duplicated(AB_table_full2$barcodes), ] #take the first row where there are cell name conflicts (only 7 cells)

merged_BCR_full2=merged_BCR_full[!duplicated(merged_BCR_full),] #remove duplicate rows
merged_BCR_full3=merged_BCR_full2[ !duplicated(merged_BCR_full2$barcode), ] #this loss doesn't matter because only clone size is needed and that is shared by IGH, IGK, and IGL

# Add AB_table data to Seurat object
data_to_add=as.data.frame(cbind(barcodes_full=M@assays[["RNA"]]@data@Dimnames[[2]], barcodes=M@assays[["RNA"]]@data@Dimnames[[2]], samples=M@meta.data$orig.ident))
data_to_add[,2]=gsub("_.*", "", data_to_add[,2])
data_to_add=left_join(data_to_add, AB_table_full3, by=c("barcodes"="barcodes", "samples"="sample"))
data_to_add[which(data_to_add[,4] == ""),4]<-NA #deal with blank spaces from the above paste
data_to_add[which(data_to_add[,5] == ""),5]<-NA
data_to_add[which(data_to_add[,6] == ""),6]<-NA
M@meta.data[["IGL"]]<-data_to_add$IGL
M@meta.data[["IGK"]]<-data_to_add$IGK
M@meta.data[["IGH"]]<-data_to_add$IGH

# Add clone size data to Seurat object
clonesize_table=unique(merged_BCR_full3[,c(1,2,33)])
clonesize_table=left_join(data_to_add, clonesize_table, by=c("barcodes" = "barcode", "samples"="sample")) #kept only 960/2736. lose nearly all of the largest clones because the barcode doesn't exist in our data
M@meta.data[["clonesize_B"]]<-clonesize_table$frequency

# Make barplots for clone size
temp=as.data.frame(table(clonesize_table$frequency))
pdf(file = "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/clonesize_barplot_full.pdf", width = 6, height = 4.25)
par(mar=c(10, 6, 4,4))
  print(ggplot(temp, aes(x=Var1, y=Freq)) + 
        geom_bar(stat="identity", fill="lightblue",color="black", 
                 position=position_dodge()) +theme_bw() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") + 
        labs(y = "Frequency", x = "Clone Size"))
dev.off()

pdf(file = "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/clonesize_barplot_redu.pdf", width = 6, height = 4.25)
par(mar=c(10, 6, 4,4))
  print(ggplot(temp, aes(x=Var1, y=Freq)) + 
        geom_bar(stat="identity", fill="lightblue",color="black", 
                 position=position_dodge()) +theme_bw() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") + 
        labs(y = "Frequency", x = "Clone Size") + ylim(0, 40)) 
dev.off()

# Make UMAP and TSNE plots
pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/UMAP_B_clonesize.pdf"), width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
  print(FeaturePlot(object = M, reduction = "umap", pt.size=0.3, features="clonesize_B", cols=c("lightgrey", "red")))
dev.off()

pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/TSNE_B_clonesize.pdf"), width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
  print(FeaturePlot(object = M, reduction = "tsne", pt.size=0.3, features="clonesize_B", cols=c("lightgrey", "red")))
dev.off()

# Save new RDS file
saveRDS(M, file = "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Seurat_Liver_30_subcluster_names_meta_TCR_BCR.rds")



##################
#                #
#      BCR       #
# Krish Pipeline #
#                #
##################

# Pull in BCR data and make minor modifications
data=read.table("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/all_bcr_mut_updated_namecorrection.csv", sep=",", header=T)
data=data[,-c(3,6)]
data[which(data$sample==""),2]<-NA
bar <- apply(cbind(data$subject, data$sample), 1, 
             function(x) paste(x[!is.na(x)], collapse = "-"))
data=cbind(sample_ID=bar, data)
data=cbind(cell=gsub("_.*","", data$read), data)
data=data[-c(which(data$lineage=="")),]

cell_names=M@assays[["SCT"]]@data@Dimnames[[2]]
cell_names=gsub("_.*", "", cell_names)
cell_names_sample=paste0(cell_names, "_", M@meta.data$orig.ident)
sample_names=unique(M@meta.data[["orig.ident"]])
data2=data[data$sample_ID %in% sample_names,]

# Calculate clone size (based off "lineage" column)
CS_long=data.frame(table(data2$sample_ID, data2$lineage))
CS_wide=data.frame(unclass(table(data2$sample_ID, data2$lineage)))

# Reduce Krish's data to only those cells in the Seurat object (real cells that passed QC)
data3=cbind(cell_sample=paste0(data2$cell, "_", data2$sample_ID), data2)
data_redu=data3[data3$cell_sample %in% cell_names_sample,]

# Attach clone size to data
data_redu=left_join(data_redu, CS_long, by=c("lineage"="Var2", "sample_ID"="Var1"))

data_redu2=data_redu[which(data_redu$Freq>1),]
data_redu2=data_redu2[order(data_redu2$Freq, decreasing=T),]
barplot(data_redu2$Freq)

mycolor=hue_pal()(length(unique(data_redu2$sample_ID))) # 30 samples

temp=mapvalues(data_redu2$sample_ID,
               from=unique(data_redu2$sample_ID),
               to=mycolor)

pdf(file = "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/Clonesize_Barplot.pdf", width = 10, height = 8)
par(mar=c(2, 2, 2, 2))
barplot(data_redu2$Freq,
        col=as.character(mapvalues(data_redu2$sample_ID,
                                   from=unique(data_redu2$sample_ID),
                                   to=mycolor)),
        ylim=c(0,80),
        ylab="Frequency",
        main="Clonotype Size Colored by Sample")
legend("topright", legend=unique(data_redu2$sample_ID), fill=mycolor)
dev.off()

# Make barplot for clone size > 2
CS_long_redu=CS_long[which(CS_long$Freq>2),]
CS_long_redu=CS_long_redu[order(CS_long_redu$Freq, decreasing=T),]
mycolor=hue_pal()(length(unique(CS_long_redu$Var1))) # 22 unique samples
temp=mapvalues(CS_long_redu$Var1,
               from=unique(CS_long_redu$Var1),
               to=mycolor)
pdf(file = "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/BCRs/Clonesize_Barplot_Expanded.pdf", width = 10, height = 8)
par(mar=c(2, 2, 2, 2))
barplot(CS_long_redu$Freq, 
        col=as.character(mapvalues(CS_long_redu$Var1,
                                   from=unique(CS_long_redu$Var1),
                                   to=mycolor)),
        ylim=c(0,120),
        ylab="Frequency",
        main="Clonotype Size Colored by Sample (3 or Larger)")
legend("topright", legend=unique(CS_long_redu$Var1), fill=mycolor)
dev.off()

# Add Krish clone size to Seurat object
temp=data.frame(cbind(cell=cell_names, sample=M@meta.data[["orig.ident"]]))
temp=left_join(temp, data_redu, by=c("cell"="cell", "sample"="sample_ID"))
M@meta.data[["clonesize_B_Krish"]]<-temp$Freq

# Look to see how well clone size agrees
cor.test(M@meta.data[["clonesize_B"]], M@meta.data[["clonesize_B_Krish"]]) #0.9602526

# Make UMAP and TSNE plots
pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/UMAP_B_Krish_clonesize.pdf"), width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
print(FeaturePlot(object = M, reduction = "umap", pt.size=0.3, features="clonesize_B_Krish", cols=c("lightgrey", "red")))
dev.off()

pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/TSNE_B_Krish_clonesize.pdf"), width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
print(FeaturePlot(object = M, reduction = "tsne", pt.size=0.3, features="clonesize_B_Krish", cols=c("lightgrey", "red")))
dev.off()

M@meta.data[["mutation_level"]]<-temp$mutation_level

# Make UMAP and TSNE plots
pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/UMAP_B_Krish_mutation.pdf"), width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
print(FeaturePlot(object = M, reduction = "umap", pt.size=0.3, features="mutation_level", cols=c("lightgrey", "blue")))
dev.off()

pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/TSNE_B_Krish_mutation.pdf"), width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
print(FeaturePlot(object = M, reduction = "tsne", pt.size=0.3, features="mutation_level", cols=c("lightgrey", "blue")))
dev.off()

# Save new RDS file
saveRDS(M, file = "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")
