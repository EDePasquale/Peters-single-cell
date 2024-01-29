##############################
#                            #
#       QC and Filtering     #
# Peters - Single Cell Liver #
#     Erica DePasquale       #
#        30 Aug 2023         #
#                            #
##############################

# Scripts modified from Peter van Galen with permission

# Options
options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)

# Load libraries
library("Matrix")
library("Matrix.utils")
library(data.table)
library(gdata)
library(Rtsne)
library(umap)
library(dplyr)
library(DropletUtils)

# Define Seurat object list
exp_list=c("ALP0003/10X-Peters-ALP-0003-BX1-20210430-5v1-1hg/ALP-0003-BX1",
           "ALP0003/10X-Peters-ALP-0003-BX2-20210430-5v1-1hg/ALP-0003-BX2",
           "ALP0003/10X-Peters-ALP-003-BX3-20210428-5v1-1hg/ALP-0003-BX3",
           "ALP0003/10X-Peters-ALP-00037-20191122-5v1hg/ALP-00037", #ALP-0003-BX4
           "ALP0003/10X-Peters-ALP-0003-BX5-20210901-5v1-1hg/ALP-0003-BX5",
           "ALP00011/10X-Peters-ALP-00011-BX2-20211108-5v1-1hg/ALP-00011-BX2",
           "ALP00017/10X-Peters-ALP-00017-20210514-5v1-1hg/ALP-00017",
           "ALP00023/10X-Peters-ALP-00023-BX1-20210901-5v1-1hg/ALP-00023-BX1",
           "ALP00023/10X-Peters-ALP-00023-BX2-20210819-5v1-1hg/ALP-00023-BX2",
           "ALP00023/10X-Peters-ALP-00023-BX3-20210819-5v1-1hg/ALP-00023-BX3",
           "ALP00030/10X-Peters-ALP-00030-BX1-20210820-5v1-1hg/ALP-00030-BX1",
           "ALP00030/10X-Peters-ALP-00030-BX4-20210820-5v1-1hg/ALP-00030-BX4",
           "ALP00033/10X-Peters-ALP-00033-BX2-20211105-5v1-1hg/ALP-00033-BX2",
           "ALP00033/10X-Peters-ALP-00033-BX3-20211105-5v1-1hg/ALP-00033-BX3",
           "ALP00036/10X-Peters-ALP-00036-20210728-5v1-1hg/ALP-00036",
           "ALP00036/10X-Peters-ALP-00036-BX2-20220224-5v1-1-hg/ALP-00036-BX2",
           "ALP00036/10X-Peters-ALP-00036-BX3-20220224-5v1-1-hg/ALP-00036-BX3",
           "ALP00041/10X-Peters-ALP-00041-TXP-20211101-5v1-1hg/ALP-00041-TXP",
           "ALP00048/10X-Peters-ALP-00048-20210514-5v1-1hg/ALP-00048",
           "ALP00066/10X-Peters-ALP-00066-20210407-5v1-1hg/ALP-00066",
           "ALP00069/10X-Peters-ALP-00069-20191121-5v1hg/ALP-00069",
           "ALP00078/10X-Peters-ALP-00078-20210322-5v1-1hg/ALP-00078",
           "ALP00017/10X-Peters-ALP-00017-BX2-20220407-5v1-1-hg/ALP-00017-BX2",
           "ALP00033/10X-Peters-ALP-00033-BX4-20220406-5v1-1-hg/ALP-00033-BX4",
           "ALP00036/10X-Peters-ALP-00036-BX4-20220228-5v1-1-hg/ALP-00036-BX4",
           "ALP00036/10X-Peters-ALP-00036-BX5-20220406-5v1-1-hg/ALP-00036-BX5",
           "ALP00042/10X-Peters-ALP-00042-BX2-20220506-5v1-1-hg/ALP-00042-BX2",
           "ALP00042/10X-Peters-ALP-00042-TXP-20220506-5v1-1-hg/ALP-00042-TXP",
           "ALP00045/10X-Peters-ALP-00045-BX1-20220505-5v1-1-hg/ALP-00045-BX1",
           "ALP00045/10X-Peters-ALP-00045-TXP-20220505-5v1-1-hg/ALP-00045-TXP")

# Create QC output folder
dir.create("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/QC_08.30.23/")

for(i in 1:length(exp_list)){
  
  print(paste0("Running experiment ", exp_list[i]))
  print(paste0("...loading files ", exp_list[i]))
  
  # Arguments
  if(basename(exp_list[i])=="ALP-00037"){name="ALP-0003-BX4"}else{name <- basename(exp_list[i])} # this takes care of the weirdly named sample
  in.dir="/data/LiverTransplantSingleCell/data"
  out.dir=paste0("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/QC_08.30.23/", name, "/")
  dir.create(out.dir)
  setwd(out.dir) #folder that contains outs
  
  genome <- "/data/GI-Informatics/DePasquale/References/refdata-gex-GRCh38-2020-A"
  min.umis <- 500 # example = 2000 # QC filter umis detected
  min.genes <- 250 # example = 1000 # QC filter genes detected
  max.mito <- 0.3 # example = 0.2 # Filter cells that exceed this fraction of mito alignments (0.2 - 1)
  sign <- T # example = F # Plot signatures?
  
  # Load auxiliary files
  colors=c("#111111", "#AAAAAA", "#B10DC9", "#FF4136", "#FFDC00", "#2ECC40", "#0074D9", "#FF851B",
           "#85144b", "#39CCCC", "#01FF70", "#3D9970", "#7FDBFF", "#F012BE", "#001f3f", "#DDDDDD")
  if (sign == T) { signatures <- read.xls("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/220420_signatures_top150.xlsx", header = T, sheet = 1) }
  source("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/200116_FunctionsGeneral.R")
  geneInfo.dt <- fread("/data/GI-Informatics/DePasquale/References/gene_info_GRCh38_PvG200204.txt")
  
  # Set up
  sampleSheet=as.data.frame(cbind(library_id=name, molecule_h5=" ", Donor=paste(unlist(strsplit(name, "-"))[1:2], collapse="-")))
  sampleSheet$WellCol <- colors[1:nrow(sampleSheet)]
  sampleSheet$DonorCol <- colors[as.numeric(as.factor(sampleSheet$Donor))+6]
  
  
  #=================================================
  # Load count matrix from CellRanger
  #=================================================
  
  samples=NULL
  sample_matrices=list(rep(NA, nrow(sampleSheet)))
  for(x in 1:nrow(sampleSheet)){
    print(paste0("Loading matrix from ", sampleSheet[x,1], "/outs/filtered_feature_bc_matrix"))
    barcode.path = paste0(in.dir, "/", exp_list[i], "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
    features.path = paste0(in.dir, "/", exp_list[i], "/outs/filtered_feature_bc_matrix/features.tsv.gz")
    matrix.path = paste0(in.dir, "/", exp_list[i], "/outs/filtered_feature_bc_matrix/matrix.mtx.gz")
    CM2=readMM(file = matrix.path)
    feature.names = read.delim(features.path, header = F)
    barcode.names = read.delim(barcode.path, header = F)
    colnames(CM2) = barcode.names$V1
    rownames(CM2) = cutf(feature.names$V2, d = "\\.", f = 1)
    CM1 <- aggregate.Matrix(CM2, groupings = rownames(CM2), fun = "sum")
    sample_matrices[x]=CM1
    rm(CM1, CM2, barcode.path, features.path, matrix.path, feature.names, barcode.names)
    samples=c(samples, sampleSheet[x,1])
  }
  names(sample_matrices)=paste0(samples, "_CM1")
  
  
  #=================================================
  # Make statistics file & QC plots
  #=================================================
  
  print(paste0("...QC plots ", exp_list[i]))
  
  # Make vectors with genes of interest
  if (genome == "/data/GI-Informatics/DePasquale/References/refdata-gex-GRCh38-2020-A") {
    Ychr.ch <- geneInfo.dt$gene_name[geneInfo.dt$chr == "Y"]
    Ychr.ch <- intersect(Ychr.ch, rownames(sample_matrices[[x]]))
    chrMgenes.ch <- geneInfo.dt$gene_name[geneInfo.dt$chr == "MT"]
  } else if (genome == "/data/GI-Informatics/DePasquale/References/refdata-gex-GRCm38-2020-A") {
    Ychr.ch <- geneInfo.dt$gene_name[geneInfo.dt$chr == "chrY"]
    Ychr.ch <- intersect(Ychr.ch, rownames(sample_matrices[[x]]))
    chrMgenes.ch <- geneInfo.dt$gene_name[geneInfo.dt$chr == "chrM"]
  }
  
  # Make file for plotting
  stats.dt=list(rep(NA, nrow(sampleSheet)))
  for(x in 1:length(sample_matrices)){
    stats1.dt <- data.table(cell = colnames(sample_matrices[[x]]),
                            umis = colSums(sample_matrices[[x]]),
                            genes = colSums(sample_matrices[[x]] != 0),
                            chrYfraction = colSums(sample_matrices[[x]][Ychr.ch[Ychr.ch %in% row.names(sample_matrices[[x]])],]) / colSums(sample_matrices[[x]]),
                            chrMfraction = colSums(sample_matrices[[x]][chrMgenes.ch[chrMgenes.ch %in% row.names(sample_matrices[[x]])],]) / colSums(sample_matrices[[x]]),
                            chrMcount = colSums(sample_matrices[[x]][chrMgenes.ch[chrMgenes.ch %in% row.names(sample_matrices[[x]])],]),
                            wellcol = sampleSheet$WellCol[[x]],
                            donorcol = sampleSheet$DonorCol[[x]])
    setorder(stats1.dt, -umis)
    stats.dt[[x]]=stats1.dt
    rm(stats1.dt)
  }
  names(stats.dt)=paste0(samples, "_stats")
  
  # Define cells that pass umi and gene thresholds
  values=list(rep(NA, nrow(sampleSheet)))
  for(x in 1:nrow(sampleSheet)){
    umiGenePass.log <- stats.dt[[x]][,"genes"] >= min.genes & stats.dt[[x]][,"umis"] >= min.umis
    mitoPass.log <- stats.dt[[x]][,"chrMfraction"] <= max.mito
    mycol <- ifelse(umiGenePass.log, yes = ifelse(mitoPass.log, yes = "black", no = "red"), no = "grey")
    cexsize <- round(max(c(0.3, 0.5-nrow(stats.dt[[x]])/40000)),2)
    values[[x]]=list(umiGenePass.log=umiGenePass.log, mitoPass.log=mitoPass.log, mycol=mycol, cexsize=cexsize)
  }
  names(values)=paste0(samples, "_values")
  
  for(x in 1:nrow(sampleSheet)){
    
    print(paste0("...QC ", x, "/", nrow(sampleSheet), " ", exp_list[i]))
    # Generate QC plots
    
    # Plot
    pdf(file = paste0(out.dir, "/", samples[x], "_QC.pdf"), width = 6, height = 6)
    par(oma=c(0,0,2,0))
    par(mfrow=c(1,1),mar=c(4,4,4,4))
    
    # cells vs umis
    plot(stats.dt[[x]]$umis, pch = 16, cex = values[[x]]$cexsize, log = "xy", col = values[[x]]$mycol, ylim = c(1, max(stats.dt[[x]]$umis)),
         xlab = "Cells", ylab = "umis", main=samples[x])
    abline(h = min.umis, col = "black", lty = 2)
    title("Cell vs UMI", outer=TRUE)
    
    # cells vs genes
    plot(stats.dt[[x]]$genes, pch = 16, cex = values[[x]]$cexsize, log = "xy", col = values[[x]]$mycol, ylim = c(1, max(stats.dt[[x]]$genes)),
         xlab = "Cells", ylab = "genes", main=samples[x])
    abline(h = min.genes, col = "black", lty = 2)
    title("Cell vs Gene", outer=TRUE)
    
    # Fraction of reads that align to chrY (color by CM well)
    plot(stats.dt[[x]]$chrYfraction, pch = 16, cex = values[[x]]$cexsize, col = stats.dt[[x]]$donorcol, 
         ylab = "Fraction chrY umis", xlab = "Cells", log = "x", main=samples[x])
    text(x=1, y=max(stats.dt[[x]]$chrYfraction)-max(stats.dt[[x]]$chrYfraction)*0.05, labels = sampleSheet$library_id[[x]], col = sampleSheet$DonorCol[[x]], pos = 4)
    title("Fraction of reads that align to chrY", outer=TRUE)
    
    # Fraction of reads that align to chrM
    plot(stats.dt[[x]]$chrMfraction,  pch = 16, cex = values[[x]]$cexsize, col = values[[x]]$mycol, ylab = "Fraction chrM umis", xlab = "Cells", log = "x", ylim = c(0, 1), main=samples[x])
    abline(h = max.mito, col = "black", lty = 2)
    title("Fraction of reads that align to chrM", outer=TRUE)
    
    # Total reads that align to chrM
    plot(stats.dt[[x]]$chrMcount, pch = 16, cex = values[[x]]$cexsize, col = values[[x]]$mycol, ylab = "chrM umi count", xlab = "Cells", log = "x", main=samples[x])
    title("Total reads that align to chrM", outer=TRUE)
    
    # Plot: genes vs UMIs
    plot(stats.dt[[x]][,c("umis", "genes")], pch = 16, cex = values[[x]]$cexsize, col = values[[x]]$mycol, main=samples[x])
    abline(v = min.umis, col = "black", lty = 2)
    abline(h = min.genes, col = "black", lty = 2)
    text(x = max(stats.dt[[x]]$umis)*0.8, y = max(stats.dt[[x]]$genes)*0.6, labels = sum(values[[x]]$mycol == "black"), col = "black")
    text(x = max(stats.dt[[x]]$umis)*0.8, y = max(stats.dt[[x]]$genes)*0.55, labels = sum(values[[x]]$mycol == "green"), col = "green")
    text(x = max(stats.dt[[x]]$umis)*0.8, y = max(stats.dt[[x]]$genes)*0.5, labels = sum(values[[x]]$mycol == "red"), col = "red")
    text(x = max(stats.dt[[x]]$umis)*0.8, y = max(stats.dt[[x]]$genes)*0.45, labels = sum(values[[x]]$mycol == "grey"), col = "grey")
    title("Gene vs UMI", outer=TRUE)
    dev.off()
    
    save.image(file = paste0(out.dir, "/", name,".RData", sep=""))
  }
  
  
  #=================================================
  # Filters
  #=================================================
  
  print(paste0("...filters ", exp_list[i]))
  sample_matrices_filtered=list(rep(NA, nrow(sampleSheet)))
  for(x in 1:nrow(sampleSheet)){
    print(x)
    # Filter cells by UMI and number of genes
    CM1=sample_matrices[[x]]
    message(paste0("\n", ncol(CM1), " cells and ", nrow(CM1), " genes in count matrix."))
    CM_XYM.dgm <- CM1[,stats.dt[[x]]$cell[values[[x]][["umiGenePass.log"]] == T]]
    message(paste0(ncol(CM_XYM.dgm), " cells and ", nrow(CM_XYM.dgm), " genes remaining after filtering for ", min.umis, " UMIs and ", min.genes, " genes."))
    
    # Filter cells based on chrM alignments
    CM.dgm <- CM_XYM.dgm[,stats.dt[[x]]$cell[values[[x]][["umiGenePass.log"]] == T & values[[x]][["mitoPass.log"]] == T]]
    CM.dgm <- CM.dgm[! rownames(CM.dgm) %in% chrMgenes.ch,]
    message(paste0(ncol(CM.dgm), " cells and ", nrow(CM.dgm), " genes remaining after filtering for cells with <", max.mito, " chrM alignments and chrM removal."))
    
    # Filter genes that map to X or Y
    if (genome == "/data/GI-Informatics/DePasquale/References/refdata-gex-GRCh38-2020-A") {
      XYchr.ch <- geneInfo.dt$gene_name[geneInfo.dt$chr == "Y" | geneInfo.dt$chr == "X"]
    } else if (genome == "/data/GI-Informatics/DePasquale/References/refdata-gex-GRCm38-2020-A") {
      XYchr.ch <- geneInfo.dt$gene_name[geneInfo.dt$chr == "chrY" | geneInfo.dt$chr == "chrX"]
    }
    CM.dgm <- CM.dgm[! rownames(CM.dgm) %in% XYchr.ch,]
    message(paste0(ncol(CM.dgm), " cells and ", nrow(CM.dgm), " genes remaining after chrXY removal."))
    
    # Subset stats table
    stats.dt[[x]] <- stats.dt[[x]][stats.dt[[x]]$cell %in% colnames(CM.dgm),]
    
    # Save CM.dgm for future use
    sample_matrices_filtered[x] <- CM.dgm
    
  }
  names(sample_matrices_filtered)=paste0(samples, "_filtered")
  
  #=================================================
  # Save sparse matrices (MM format)
  #=================================================
  
  for(x in 1:length(unique(sampleSheet$library_id))){
    
    # Combine expression matrices for (genes and order match between datasets)
    temp_i=which(sampleSheet$library_id == unique(sampleSheet$library_id)[x])
    temp_m=NULL
    for(j in temp_i){
      temp_m2=sample_matrices_filtered[[j]]
      colnames(temp_m2)=paste0(colnames(sample_matrices_filtered[[j]]), "_", j)
      temp_m=cbind(temp_m, temp_m2)
    }
    
    # Write out expression matrix (10x format)
    DropletUtils:::write10xCounts(paste0("./", unique(sampleSheet$library_id)[x], "_filtered_matrices"), temp_m, version="3")
  }
  
  
  #=================================================
  # Generate Dimensionality Plots
  #=================================================
  
  for(x in 1:nrow(sampleSheet)){
    
    print(paste0("...variable genes ", x, "/", nrow(sampleSheet), " ", exp_list[i]))
    
    # Normalize count matrix
    CM.norm <- t(t(as.matrix(sample_matrices_filtered[[x]]))/colSums(as.matrix(sample_matrices_filtered[[x]]))) * 10000   # Normalize count matrix to 10,000 transcripts per cell
    
    # The following is basically identical to the function variableGenes that is used for Seq-Well
    CM.norm.min <- 0.01
    sd.cutoff <- 1
    
    # Calculate mean and coefficient of variation per gene, remove genes that are not expressed at all
    CM.norm.mean <- rowMeans(CM.norm)
    CM.norm.cv <- apply(CM.norm, 1, sd) / CM.norm.mean
    CM.norm.mean <- sort(CM.norm.mean[!is.na(CM.norm.cv)])
    CM.norm.cv <- CM.norm.cv[names(CM.norm.mean)]
    
    # Plot
    message(paste("\nGenerating", paste0(out.dir, "/", samples[x], "_VarGenes.pdf")))
    pdf(file = paste0(out.dir, "/", samples[x], "_VarGenes.pdf"), width = 6, height = 6)
    par(mar=c(4,4,4,4))
    plot(CM.norm.mean, CM.norm.cv, pch = ".", log = "xy", xlab = "Mean expression", ylab = "Coefficient of variation")
    abline(v = CM.norm.min, col = "blue")
    
    # Linear fit
    CM.norm.lm <- lm(y~x, subset = CM.norm.mean>=CM.norm.min, data.frame(x = log10(CM.norm.mean), y = log10(CM.norm.cv)))
    CM.norm.cv.predict <- predict(CM.norm.lm, data.frame(x = log10(CM.norm.mean)))
    lines(CM.norm.mean, 10^CM.norm.cv.predict, col = "red")
    
    # This is from Volker's function, I do not have a complete understanding
    CM.norm.cv.residue <- log10(CM.norm.cv) - CM.norm.cv.predict
    plot(CM.norm.mean, CM.norm.cv.residue, log="x", pch=".", xlab="Mean", ylab="Residue")
    abline(v=CM.norm.min, col="green")
    abline(h=0, col="red")
    abline(h=sd(CM.norm.cv.residue[CM.norm.mean>=CM.norm.min])*c(1, 2, 3), col="red", lty=2)
    axis(4, at = round(sd(CM.norm.cv.residue[CM.norm.mean>=CM.norm.min])*c(1, 2, 3), 3), las = 2)
    
    plot(seq(-0.25+0.005, 1-0.005, 0.01), as.vector(table(cut(CM.norm.cv.residue[CM.norm.mean>=CM.norm.min], seq(-0.25, 1, 0.01)))), xlab="Residue", ylab="Frequency", type="b")
    abline(v = 0, col = "red")
    abline(v = sd(CM.norm.cv.residue[CM.norm.mean >= CM.norm.min])*c(1, 2, 3), col="red", lty = 2)
    axis(3, at = round(sd(CM.norm.cv.residue[CM.norm.mean>=CM.norm.min])*c(1, 2, 3), 3), las = 2)
    
    plot(ecdf(CM.norm.cv.residue[CM.norm.mean>=CM.norm.min]), main="ECDF")
    qqnorm(CM.norm.cv.residue[CM.norm.mean>=CM.norm.min])
    
    # This is the key object
    vargenes.num <- sort((CM.norm.cv.residue)[CM.norm.mean>=CM.norm.min], decreasing = T)
    attr(vargenes.num, which = "sd") <- sd(CM.norm.cv.residue[CM.norm.mean>=CM.norm.min])
    VarGenes1.ch <- names(which(vargenes.num > sd.cutoff*attr(vargenes.num, "sd")))
    
    # Remove cell cycle genes to generate VarGenes.ch
    VarGenes.ch <- setdiff(VarGenes1.ch, c("ASPM", "CENPE", "CENPF", "DLGAP5", "MKI67", "NUSAP1", "PCLAF", "STMN1", "TOP2A", "TUBB"))
    
    plot(CM.norm.mean, CM.norm.cv, log="xy", pch=".", col=ifelse(names(CM.norm.cv) %in% names(which(vargenes.num > sd.cutoff*attr(vargenes.num, which = "sd"))), "red", "grey"), xlab="Mean", ylab="CV")
    abline(v=CM.norm.min, col="blue")
    lines(CM.norm.mean, 10^CM.norm.cv.predict, col="red")
    text(x = 50, y = 20, labels = paste(length(VarGenes.ch), "genes"), col = "red")
    dev.off()
    
    
    #=================================================
    # tSNE
    #=================================================
    
    print(paste0("...tSNE ", x, "/", nrow(sampleSheet), " ", exp_list[i]))
    # Adjustments to count matrix
    CM.norm.log <- log(CM.norm+1)                         # Normalize & log
    CM.norm.log.var <- CM.norm.log[VarGenes.ch,]          # Normalize & log and only keep variable genes 
    
    # Based on 2_tSNE.PvG190810.R and 1d_runTSNE.vh20190325.R
    message("\nCalculating tSNE coordinates")
    
    # Make tSNE with default values & add values to stats.dt
    set.seed(42)
    tryCatch(
      expr = {
        CM.tsne <- Rtsne(as.matrix(t(CM.norm.log.var)), max_iter = 1000, theta = 0.5, perplexity = 30, verbose = TRUE)
      },
      error = function(e){ 
        tryCatch(
          expr={
            CM.tsne <- Rtsne(as.matrix(t(CM.norm.log.var)), max_iter = 1000, theta = 0.5, perplexity = 5, verbose = TRUE)
          },
          error = function(e){
            CM.tsne <- Rtsne(as.matrix(t(CM.norm.log.var)), max_iter = 1000, theta = 0.5, perplexity = 1, verbose = TRUE)
          }
        )
      }
    )
    stopifnot( all(colnames(CM.norm.log) == stats.dt[[x]]$cell) )# just checking...
    stats.dt$tSNEx <- CM.tsne$Y[,1]
    stats.dt$tSNEy <- CM.tsne$Y[,2]
    
    message(paste("\nGenerating", paste0(out.dir, "/", samples[x], "_tSNE.pdf")))
    pdf(file = paste0(out.dir, "/", samples[x], "_tSNE.pdf"), width = 6, height = 6)
    par(mar=c(4, 4, 4, 4))
    
    # Color by UMI
    plotTSNE(cbind(stats.dt$tSNEx, stats.dt$tSNEy), pch = 16, cex = cexsize, col = colCustom(stats.dt[[x]]$umis, c(2000, 20000)),
             main = paste0("UMIs (", min(stats.dt[[x]]$umis), "-", max(stats.dt[[x]]$umis), ", mean ", round(mean(stats.dt[[x]]$umis), 0), ")"))
    axis(side = 1, at = 0, label = paste(nrow(stats.dt[[x]]), "cells"), tick = F)
    
    # Color by Donor
    plotTSNE(cbind(stats.dt$tSNEx, stats.dt$tSNEy), cex = cexsize, col = stats.dt[[x]]$donorcol, main = "Donor")
    sapply(1:length(unique(sampleSheet$Donor)), function(z) axis(side = 1, at = ifelse(length(unique(sampleSheet$Donor)) == 1,
                                                                                       yes = 0, no = min(stats.dt$tSNEx) + (max(stats.dt$tSNEx)-min(stats.dt$tSNEx))/(length(unique(sampleSheet$Donor))-1)*(z-1)), 
                                                                 tick = F, col.axis = unique(sampleSheet$DonorCol)[z], labels = unique(sampleSheet$Donor)[z]))
    
    # Color by CM well
    stats.dt[[x]]$wellcol=gsub(".*-","",stats.dt[[x]][["cell"]])
    stats.dt[[x]]$wellcol=gsub("1", sampleSheet$WellCol[1], stats.dt[[x]]$wellcol)
    stats.dt[[x]]$wellcol=gsub("2", sampleSheet$WellCol[2], stats.dt[[x]]$wellcol)
    stats.dt[[x]]$wellcol=gsub("3", sampleSheet$WellCol[3], stats.dt[[x]]$wellcol)
    plotTSNE(cbind(stats.dt$tSNEx, stats.dt$tSNEy), cex = cexsize, col = stats.dt[[x]]$wellcol, main = "Library")
    plot(x = rep(1, nrow(sampleSheet)), y = nrow(sampleSheet):1, pch = 21, bg = sampleSheet$WellCol, cex = 3, xlim = c(0,15), axes = F, ylab = "", xlab = "")
    text(x = rep(1.5, nrow(sampleSheet)), y = nrow(sampleSheet):1, labels = sampleSheet$library_id, pos = 4)
    
    # Color by signature score
    if (sign == T) {
      message("\nCalculating signature scores")
      # Signatures
      signatures <- lapply(as.list(signatures), function(y) y[!is.na(y)])
      signatures <- lapply(signatures, intersect, rownames(CM.norm.log))
      # Average gene expression for scoreSignature
      CM.mean <- rowMeans(CM.norm.log)
      signScore <- lapply(names(signatures), function(g) {
        message(g)
        scoreSignature(CM = as.matrix(CM.norm.log), signatures = signatures[[g]], CM.mean = CM.mean, verbose = T)
      })
      names(signScore) <- names(signatures)
      
      # Plot
      for (n in names(signScore)) {
        mycol <- colItay(signScore[[n]])
        db=cbind(stats.dt$tSNEx, stats.dt$tSNEy)
        plotTSNE(db, pch = 16, cex = cexsize, col = mycol, main = n)
      }
    }
    dev.off()
    
    
    #=================================================
    # UMAP
    #=================================================
    
    print(paste0("...UMAP ", x, "/", nrow(sampleSheet), " ", exp_list[i]))
    # Some standard parameters to plot YMAP coordinates
    message("plotUMAP()")
    plotUMAP <- function(x, pch = 16, ...) {
      par(mfrow=c(1,1),mar=c(4,4,4,4))
      plot(x, xlab = NA, ylab = NA, tck = F, yaxt = "n", xaxt = "n", pch = pch, ...)
    }
    
    message("\nCalculating UMAP coordinates")
    
    # Make tSNE with default values & add values to stats.dt
    set.seed(42)
    tryCatch(
      expr={
        CM.umap <- umap(as.matrix(t(CM.norm.log.var)))
      },
      error = function(e){
        print("NO")
      }
    )
    
    stopifnot( all(colnames(CM.norm.log) == stats.dt$cell) )# just checking...
    stats.dt$umapx <- CM.umap$layout[,1]
    stats.dt$umapy <- CM.umap$layout[,2]
    
    message(paste("\nGenerating", paste0(out.dir, "/", samples[x], "_UMAP.pdf")))
    pdf(file = paste0(out.dir, "/", samples[x], "_UMAP.pdf"), width = 6, height = 6)
    par(mar=c(4, 4, 4, 4))
    
    # Color by UMI
    plotUMAP(cbind(stats.dt$umapx, stats.dt$umapy), pch = 16, cex = cexsize, col = colCustom(stats.dt[[x]]$umis, c(2000, 20000)),
             main = paste0("UMIs (", min(stats.dt[[x]]$umis), "-", max(stats.dt[[x]]$umis), ", mean ", round(mean(stats.dt[[x]]$umis), 0), ")"))
    axis(side = 1, at = 0, label = paste(nrow(stats.dt[[x]]), "cells"), tick = F)
    
    # Color by Donor
    plotUMAP(cbind(stats.dt$umapx, stats.dt$umapy), cex = cexsize, col = stats.dt[[x]]$donorcol, main = "Donor")
    sapply(1:length(unique(sampleSheet$Donor)), function(z) axis(side = 1, at = ifelse(length(unique(sampleSheet$Donor)) == 1,
                                                                                       yes = 0, no = min(stats.dt$umapx) + (max(stats.dt$umapx)-min(stats.dt$umapx))/(length(unique(sampleSheet$Donor))-1)*(z-1)), 
                                                                 tick = F, col.axis = unique(sampleSheet$DonorCol)[z], labels = unique(sampleSheet$Donor)[z]))
    
    # Color by CM well
    stats.dt[[x]]$wellcol=gsub(".*-","",stats.dt[[x]][["cell"]])
    stats.dt[[x]]$wellcol=gsub("1", sampleSheet$WellCol[1], stats.dt[[x]]$wellcol)
    stats.dt[[x]]$wellcol=gsub("2", sampleSheet$WellCol[2], stats.dt[[x]]$wellcol)
    stats.dt[[x]]$wellcol=gsub("3", sampleSheet$WellCol[3], stats.dt[[x]]$wellcol)
    plotUMAP(cbind(stats.dt$umapx, stats.dt$umapy), cex = cexsize, col = stats.dt[[x]]$wellcol, main = "Library")
    plot(x = rep(1, nrow(sampleSheet)), y = nrow(sampleSheet):1, pch = 21, bg = sampleSheet$WellCol, cex = 3, xlim = c(0,15), axes = F, ylab = "", xlab = "")
    text(x = rep(1.5, nrow(sampleSheet)), y = nrow(sampleSheet):1, labels = sampleSheet$library_id, pos = 4)
    
    # Color by signature score
    if (sign == T) {
      # Plot
      for (n in names(signScore)) {
        mycol <- colItay(signScore[[n]])
        plotUMAP(cbind(stats.dt$umapx, stats.dt$umapy), pch = 16, cex = cexsize, col = mycol, main = n)
      }
    }
    dev.off()
  }
  save.image(file = paste0(out.dir, "/", name,"_QC.RData", sep=""))
}

