# Load libraries
library(immunarch)
library(dplyr)
library(RColorBrewer)

# https://immunarch.com/articles/web_only/v8_tracking.html 

# Set working directory for saving plots
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/")

# Set up colors
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#########################################
#                                       #
# <><><><> SHARED AND EXPANDED <><><><> #
#                                       #
#########################################

##############################
#                            #
# Patient 1 (i.e., ALP-0003) #
#                            #
##############################

expt_list=c("ALP-0003-BX1",
            "ALP-0003-BX2",
            "ALP-0003-BX3",
            "ALP-0003-BX4",
            "ALP-0003-BX5")

# Find the shared, expanded clonotypes and set colors accordingly
shared=NULL
for(i in 1:length(expt_list)){
  path=paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/", expt_list[i], "/clonotypes.csv")
  Y=read.table(path, sep=",", header=T)
  assign(paste0(expt_list[i],"_clonotypes"),Y)
  Z=Y[which(Y$frequency > 2),]
  shared=c(shared, Z$cdr3s_aa)
}
shared=unique(shared)
mycols=sample(col_vector, length(shared))
# [1] "#B3CDE3" "#80B1D3" "#1F78B4" "#CCCCCC" "#FF7F00" "#E5C494" "#FFFFB3" "#E78AC3" "#FDDAEC"
# [10] "#F781BF" "#984EA3" "#FFED6F" "#D95F02" "#E5D8BD" "#F4CAE4" "#E41A1C" "#FB9A99" "#B2DF8A"
# [19] "#A65628"
mycols=c("#B3CDE3", "#80B1D3", "#1F78B4", "#CCCCCC", "#FF7F00", "#E5C494", "#FFFFB3", "#E78AC3", "#FDDAEC", 
         "#F781BF", "#984EA3", "#FFED6F", "#D95F02", "#E5D8BD", "#F4CAE4", "#E41A1C", "#FB9A99", "#B2DF8A", 
         "#A65628")

# Create data frames to save the clonotype proportions (used in barplot) and frequency (saved as table)
results=as.data.frame(matrix(nrow=length(shared), ncol=length(expt_list)+1))
colnames(results)=c("CDR3.aa", expt_list)
results[,1]=shared

results2=as.data.frame(matrix(nrow=length(shared), ncol=length(expt_list)+1))
colnames(results2)=c("CDR3.aa", expt_list)
results2[,1]=shared

for(i in 1:length(expt_list)){
  Y=get(paste0(expt_list[i],"_clonotypes"))[,2:4]
  results[,i+1]=left_join(results, Y[,2:3], by=c("CDR3.aa"="cdr3s_aa"))[,7]
  results2[,i+1]=left_join(results2, Y[,c(1,3)], by=c("CDR3.aa"="cdr3s_aa"))[,7]
}

results[is.na(results)] <- 0
results2[is.na(results2)] <- 0
write.table(results2, file="ALP-0003_shared_expanded_clones_freq.txt", row.names=F, quote=F, sep="\t")

# Change class to work with visualization function and plot
class(results)<-c("immunr_dynamics", "data.table",  "data.frame")
A=vis(results, .plot = "smooth")

pdf(file = "ALP-0003_shared_expanded_clones_bar.pdf", width = 12, height = 8)
par(mar=c(2, 2, 2, 2))
  A + scale_fill_manual(values=mycols)
dev.off()

#global stats table
results3=as.data.frame(matrix(nrow=2, ncol=length(expt_list)))
colnames(results3)=c(expt_list)
row.names(results3)=c("Total Clonotypes", "Expanded Clonotypes")
for(i in 1:length(expt_list)){
  Y=get(paste0(expt_list[i],"_clonotypes"))
  results3[1,i]=nrow(Y)
  results3[2,i]=length(which(Y$frequency > 2))
}
write.table(results3, file="ALP-0003_shared_expanded_clones_freq_part2.txt", row.names=F, quote=F, sep="\t")


################################
#                              #
# Patient 17 (i.e., ALP-00036) # without BX1 (removed 4.18.23)
#                              #
################################

expt_list=c("ALP-00036",
            "ALP-00036-BX2",
            "ALP-00036-BX3",
            "ALP-00036-BX4",
            "ALP-00036-BX5")

# Find the shared, expanded clonotypes and set colors accordingly
shared=NULL
for(i in 1:length(expt_list)){
  path=paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/", expt_list[i], "/clonotypes.csv")
  Y=read.table(path, sep=",", header=T)
  assign(paste0(expt_list[i],"_clonotypes"),Y)
  Z=Y[which(Y$frequency > 2),]
  shared=c(shared, Z$cdr3s_aa)
}
shared=unique(shared)
mycols=sample(col_vector, length(shared))
# [1] "#666666" "#FDBF6F" "#CAB2D6" "#984EA3" "#FDB462" "#CBD5E8" "#FF7F00" "#FED9A6" "#80B1D3"
# [10] "#FFED6F" "#FFFFCC" "#E31A1C" "#1F78B4" "#CCEBC5" "#E5D8BD" "#FFFF99" "#A6D854" "#F0027F"
mycols=c("#666666", "#FDBF6F", "#CAB2D6", "#984EA3", "#FDB462", "#CBD5E8", "#FF7F00", "#FED9A6", "#80B1D3",
         "#FFED6F", "#FFFFCC", "#E31A1C", "#1F78B4", "#CCEBC5", "#E5D8BD", "#FFFF99", "#A6D854", "#F0027F")

# Create data frames to save the clonotype proportions (used in barplot) and frequency (saved as table)
results=as.data.frame(matrix(nrow=length(shared), ncol=length(expt_list)+1))
colnames(results)=c("CDR3.aa", expt_list)
results[,1]=shared

results2=as.data.frame(matrix(nrow=length(shared), ncol=length(expt_list)+1))
colnames(results2)=c("CDR3.aa", expt_list)
results2[,1]=shared

for(i in 1:length(expt_list)){
  Y=get(paste0(expt_list[i],"_clonotypes"))[,2:4]
  results[,i+1]=left_join(results, Y[,2:3], by=c("CDR3.aa"="cdr3s_aa"))[,7]
  results2[,i+1]=left_join(results2, Y[,c(1,3)], by=c("CDR3.aa"="cdr3s_aa"))[,7]
}

results[is.na(results)] <- 0
results2[is.na(results2)] <- 0
write.table(results2, file="ALP-00036_shared_expanded_clones_freq.txt", row.names=F, quote=F, sep="\t")

# Change class to work with visualization function and plot
class(results)<-c("immunr_dynamics", "data.table",  "data.frame")
B=vis(results, .plot = "smooth")

pdf(file = "ALP-00036_shared_expanded_clones_bar.pdf", width = 12, height = 8)
par(mar=c(2, 2, 2, 2))
  B + scale_fill_manual(values=mycols)
dev.off()

#global stats table
results3=as.data.frame(matrix(nrow=2, ncol=length(expt_list)))
colnames(results3)=c(expt_list)
row.names(results3)=c("Total Clonotypes", "Expanded Clonotypes")
for(i in 1:length(expt_list)){
  Y=get(paste0(expt_list[i],"_clonotypes"))
  results3[1,i]=nrow(Y)
  results3[2,i]=length(which(Y$frequency > 2))
}
write.table(results3, file="ALP-00036_shared_expanded_clones_freq_part2.txt", row.names=F, quote=F, sep="\t")


###############################
#                             #
#         ALP-00017           #
#                             #
###############################

expt_list=c("ALP-00017",
            "ALP-00017-BX2")

# Find the shared, expanded clonotypes and set colors accordingly
shared=NULL
for(i in 1:length(expt_list)){
  path=paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/", expt_list[i], "/clonotypes.csv")
  Y=read.table(path, sep=",", header=T)
  assign(paste0(expt_list[i],"_clonotypes"),Y)
  Z=Y[which(Y$frequency > 2),]
  shared=c(shared, Z$cdr3s_aa)
}
shared=unique(shared)
mycols=sample(col_vector, length(shared))
# [1] "#1B9E77" "#E6F5C9" "#FC8D62" "#386CB0" "#FFFF33" "#FBB4AE" "#33A02C" "#7FC97F" "#CCCCCC"
# [10] "#FFFF99" "#E6AB02" "#999999"
mycols=c("#1B9E77", "#E6F5C9", "#FC8D62", "#386CB0", "#FFFF33", "#FBB4AE", "#33A02C", "#7FC97F", "#CCCCCC",
         "#FFFF99", "#E6AB02", "#999999")

# Create data frames to save the clonotype proportions (used in barplot) and frequency (saved as table)
results=as.data.frame(matrix(nrow=length(shared), ncol=length(expt_list)+1))
colnames(results)=c("CDR3.aa", expt_list)
results[,1]=shared

results2=as.data.frame(matrix(nrow=length(shared), ncol=length(expt_list)+1))
colnames(results2)=c("CDR3.aa", expt_list)
results2[,1]=shared

for(i in 1:length(expt_list)){
  Y=get(paste0(expt_list[i],"_clonotypes"))[,2:4]
  results[,i+1]=left_join(results, Y[,2:3], by=c("CDR3.aa"="cdr3s_aa"))[,4]
  results2[,i+1]=left_join(results2, Y[,c(1,3)], by=c("CDR3.aa"="cdr3s_aa"))[,4]
}

results[is.na(results)] <- 0
results2[is.na(results2)] <- 0
write.table(results2, file="ALP-00017_shared_expanded_clones_freq.txt", row.names=F, quote=F, sep="\t")

# Change class to work with visualization function and plot
class(results)<-c("immunr_dynamics", "data.table",  "data.frame")
C=vis(results, .plot = "smooth")

pdf(file = "ALP-00017_shared_expanded_clones_bar.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
  C + scale_fill_manual(values=mycols)
dev.off()

#global stats table
results3=as.data.frame(matrix(nrow=2, ncol=length(expt_list)))
colnames(results3)=c(expt_list)
row.names(results3)=c("Total Clonotypes", "Expanded Clonotypes")
for(i in 1:length(expt_list)){
  Y=get(paste0(expt_list[i],"_clonotypes"))
  results3[1,i]=nrow(Y)
  results3[2,i]=length(which(Y$frequency > 2))
}
write.table(results3, file="ALP-00017_shared_expanded_clones_freq_part2.txt", row.names=F, quote=F, sep="\t")


###############################
#                             #
#         ALP-00023           #
#                             #
###############################

expt_list=c("ALP-00023-BX1",
            "ALP-00023-BX2",
            "ALP-00023-BX3")

# Find the shared, expanded clonotypes and set colors accordingly
shared=NULL
for(i in 1:length(expt_list)){
  path=paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/", expt_list[i], "/clonotypes.csv")
  Y=read.table(path, sep=",", header=T)
  assign(paste0(expt_list[i],"_clonotypes"),Y)
  Z=Y[which(Y$frequency > 2),]
  shared=c(shared, Z$cdr3s_aa)
}
shared=unique(shared)
mycols=sample(col_vector, length(shared))
# [1] "#CCEBC5" "#FBB4AE" "#4DAF4A" "#377EB8" "#FDCDAC" "#BEAED4" "#FFED6F" "#BF5B17" "#80B1D3"
# [10] "#CCCCCC" "#F781BF" "#D95F02" "#FFFF99" "#CBD5E8" "#FFF2AE"
mycols=c("#CCEBC5", "#FBB4AE", "#4DAF4A", "#377EB8", "#FDCDAC", "#BEAED4", "#FFED6F", "#BF5B17", "#80B1D3",
         "#CCCCCC", "#F781BF", "#D95F02", "#FFFF99", "#CBD5E8", "#FFF2AE")

# Create data frames to save the clonotype proportions (used in barplot) and frequency (saved as table)
results=as.data.frame(matrix(nrow=length(shared), ncol=length(expt_list)+1))
colnames(results)=c("CDR3.aa", expt_list)
results[,1]=shared

results2=as.data.frame(matrix(nrow=length(shared), ncol=length(expt_list)+1))
colnames(results2)=c("CDR3.aa", expt_list)
results2[,1]=shared

for(i in 1:length(expt_list)){
  Y=get(paste0(expt_list[i],"_clonotypes"))[,2:4]
  results[,i+1]=left_join(results, Y[,2:3], by=c("CDR3.aa"="cdr3s_aa"))[,5]
  results2[,i+1]=left_join(results2, Y[,c(1,3)], by=c("CDR3.aa"="cdr3s_aa"))[,5]
}

results[is.na(results)] <- 0
results2[is.na(results2)] <- 0
write.table(results2, file="ALP-00023_shared_expanded_clones_freq.txt", row.names=F, quote=F, sep="\t")

# Change class to work with visualization function and plot
class(results)<-c("immunr_dynamics", "data.table",  "data.frame")
D=vis(results, .plot = "smooth")

pdf(file = "ALP-00023_shared_expanded_clones_bar.pdf", width = 10, height = 8)
par(mar=c(2, 2, 2, 2))
  D + scale_fill_manual(values=mycols)
dev.off()

#global stats table
results3=as.data.frame(matrix(nrow=2, ncol=length(expt_list)))
colnames(results3)=c(expt_list)
row.names(results3)=c("Total Clonotypes", "Expanded Clonotypes")
for(i in 1:length(expt_list)){
  Y=get(paste0(expt_list[i],"_clonotypes"))
  results3[1,i]=nrow(Y)
  results3[2,i]=length(which(Y$frequency > 2))
}
write.table(results3, file="ALP-00023_shared_expanded_clones_freq_part2.txt", row.names=F, quote=F, sep="\t")


###############################
#                             #
#         ALP-00033           # # Not enough expanded clones (2) to keep this figure in the paper
#                             #
###############################

expt_list=c("ALP-00033-BX2",
            "ALP-00033-BX3",
            "ALP-00033-BX4")

# Find the shared, expanded clonotypes and set colors accordingly
shared=NULL
for(i in 1:length(expt_list)){
  path=paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/", expt_list[i], "/clonotypes.csv")
  Y=read.table(path, sep=",", header=T)
  assign(paste0(expt_list[i],"_clonotypes"),Y)
  Z=Y[which(Y$frequency > 2),]
  shared=c(shared, Z$cdr3s_aa)
}
shared=unique(shared)
mycols=sample(col_vector, length(shared))
# [1] "#CCEBC5" "#FBB4AE" "#4DAF4A" "#377EB8" "#FDCDAC" "#BEAED4" "#FFED6F" "#BF5B17" "#80B1D3"
# [10] "#CCCCCC" "#F781BF" "#D95F02" "#FFFF99" "#CBD5E8" "#FFF2AE"
# mycols=c("#CCEBC5", "#FBB4AE", "#4DAF4A", "#377EB8", "#FDCDAC", "#BEAED4", "#FFED6F", "#BF5B17", "#80B1D3",
#          "#CCCCCC", "#F781BF", "#D95F02", "#FFFF99", "#CBD5E8", "#FFF2AE")

# Create data frames to save the clonotype proportions (used in barplot) and frequency (saved as table)
results=as.data.frame(matrix(nrow=length(shared), ncol=length(expt_list)+1))
colnames(results)=c("CDR3.aa", expt_list)
results[,1]=shared

results2=as.data.frame(matrix(nrow=length(shared), ncol=length(expt_list)+1))
colnames(results2)=c("CDR3.aa", expt_list)
results2[,1]=shared

for(i in 1:length(expt_list)){
  Y=get(paste0(expt_list[i],"_clonotypes"))[,2:4]
  results[,i+1]=left_join(results, Y[,2:3], by=c("CDR3.aa"="cdr3s_aa"))[,5]
  results2[,i+1]=left_join(results2, Y[,c(1,3)], by=c("CDR3.aa"="cdr3s_aa"))[,5]
}

results[is.na(results)] <- 0
results2[is.na(results2)] <- 0
write.table(results2, file="ALP-00033_shared_expanded_clones_freq.txt", row.names=F, quote=F, sep="\t")

# Change class to work with visualization function and plot
class(results)<-c("immunr_dynamics", "data.table",  "data.frame")
D=vis(results, .plot = "smooth")

pdf(file = "ALP-00033_shared_expanded_clones_bar.pdf", width = 10, height = 8)
par(mar=c(2, 2, 2, 2))
D + scale_fill_manual(values=mycols)
dev.off()

#global stats table
results3=as.data.frame(matrix(nrow=2, ncol=length(expt_list)))
colnames(results3)=c(expt_list)
row.names(results3)=c("Total Clonotypes", "Expanded Clonotypes")
for(i in 1:length(expt_list)){
  Y=get(paste0(expt_list[i],"_clonotypes"))
  results3[1,i]=nrow(Y)
  results3[2,i]=length(which(Y$frequency > 2))
}
write.table(results3, file="ALP-00033_shared_expanded_clones_freq_part2.txt", row.names=F, quote=F, sep="\t")

