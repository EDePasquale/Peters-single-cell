##############################
#                            #
#     Make Seurat Objects    #
# Peters - Single Cell Liver #
#     Erica DePasquale       #
#        30 Aug 2023         #
#                            #
##############################

# Largely copied from the Seurat 3k PBMC vignette

# Load libraries
library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
library(base)
library(methods)
library(utils)
library(stats)
library(gdata)
library(graphics)
library(grDevices)

# Set list of samples
expt_list=c("ALP-0003-BX1",
            "ALP-0003-BX2",
            "ALP-0003-BX3",
            "ALP-0003-BX4",
            "ALP-0003-BX5",
            "ALP-00011-BX2",
            "ALP-00017",
            "ALP-00023-BX1",
            "ALP-00023-BX2",
            "ALP-00023-BX3",
            "ALP-00030-BX1",
            "ALP-00030-BX4",
            "ALP-00033-BX2",
            "ALP-00033-BX3",
            "ALP-00036",
            "ALP-00036-BX2",
            "ALP-00036-BX3",
            "ALP-00041-TXP",
            "ALP-00048",
            "ALP-00066",
            "ALP-00069",
            "ALP-00078",
            "ALP-00017-BX2",
            "ALP-00033-BX4",
            "ALP-00036-BX4",
            "ALP-00036-BX5",
            "ALP-00042-BX2",
            "ALP-00042-TXP",
            "ALP-00045-BX1",
            "ALP-00045-TXP")

for (expt in 1:length(expt_list)){
  
  # Set up path
  name <- expt_list[expt]
  path <- paste0("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/QC_08.30.23/", name, "/", name, "_filtered_matrices")
  
  # Load the dataset
  setwd(path)
  M.data <- Read10X(data.dir = path)
  
  #Examine the memory savings between regular and sparse matrices
  dense.size <- object.size(as.matrix(M.data))
  dense.size
  sparse.size <- object.size(M.data,sparseMatrixClass='Matrix')
  sparse.size
  dense.size/sparse.size
  
  # Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
  # Initialize the Seurat object with the raw (non-normalized data).
  M <- CreateSeuratObject(counts = M.data, project = name, min.cells = 3, min.features = 200)
  M
  
  #nGene and nUMI are automatically calculated for every object by Seurat. For non-UMI data, nUMI represents the sum of the non-normalized values within a cell
  # We calculate the percentage of mitochondrial genes here and store it in percent.mito using the AddMetaData. The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
  M[["percent.mt"]] <- PercentageFeatureSet(object = M, pattern = "^MT-")
  pdf("QCstats.pdf", width = 12, height = 15)
  VlnPlot(object = M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  dev.off()
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  
  plot1 <- FeatureScatter(object = M, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(object = M, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  pdf("Features-plot.pdf", width = 15, height = 15)
    CombinePlots(plots = list(plot1, plot2))
  dev.off()
  
  #M <- subset(x = M, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt <40)
  ###Normalizing the data
  ###After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the gene expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
  M <- NormalizeData(object = M, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Identification of highly variable features (feature selection)
  # We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). We and others have found that focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.
  M <- FindVariableFeatures(object = M, selection.method = "vst", nfeatures = 2000)
  
  # Identify the 10 most highly variable genes
  top10 <- head(x = VariableFeatures(object = M), 10)
  write.table(top10,file="top10-high-variable-genes.txt",sep="\t",col.names= NA)
  
  # plot variable features with and without labels
  pdf("Variable-Feature-Plot.pdf", width = 13, height = 15)
  plot1 <- VariableFeaturePlot(object = M)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    CombinePlots(plots = list(plot1, plot2))
  dev.off()
  
  ##Scaling the data and removing unwanted sources of variation
  all.genes <- rownames(x = M)
  M <- ScaleData(object = M, features = all.genes)
  
  ## Perform linear dimensional reduction
  M <- RunPCA(object = M, features = VariableFeatures(object = M))
  M
  
  # Examine  and visualize PCA results a few different ways
  pdf("PCA.pdf", width = 8, height = 8)
    VizDimLoadings(object = M, dims = 1:2, reduction = "pca")
  dev.off()
  
  pdf("PCA-plot.pdf", width = 10, height = 8)
    DimPlot(object = M, reduction = "pca")
  dev.off()
  
  pdf("PC-Heatmap.pdf", width = 8, height = 8)
    DimHeatmap(object = M, dims = 1, cells = 500, balanced = TRUE)
    DimHeatmap(object = M, dims = 1:15, cells = 500, balanced = TRUE)
  dev.off()
  
  ####Determine statistically significant principal components
  M <- JackStraw(object = M, num.replicate = 100)
  M <- ScoreJackStraw(object = M, dims = 1:20)
  png("JackStraw.png", width = 8, height = 8, units = 'in', res = 600)
  pdf("JackStraw.pdf", width = 8, height = 8)
    JackStrawPlot(object = M, dims = 1:20)
  dev.off()
  
  png("PCE-bowl-plot.png", width = 8, height = 8, units = 'in', res = 600)
  pdf("PCE-bowl-plot.pdf", width = 8, height = 8)
    ElbowPlot(object = M)
  dev.off()
  
  ## Cluster the cells
  #FindNeighbors function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs).
  M <- FindNeighbors(object = M, dims = 1:20)
  M <- FindClusters(object = M, resolution = 0.5)
  head(x = Idents(object = M), 5)
  
  ####Run non-linear dimensional reduction (UMAP)####
  pdf("UMAP-plot.pdf", width = 10, height = 8)
  M <- RunUMAP(object = M, dims = 1:20)
    DimPlot(object = M, reduction = "umap",label= TRUE)
  dev.off()
  
  
  ##############################
  #                            #
  # Score Cells for Signatures #
  #                            #
  ##############################
  
  source("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/200116_FunctionsGeneral.R")
  
  signatures <- read.xls("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/220420_signatures_top150.xlsx", header = T, sheet = 1)
  cexsize <- round(max(c(0.3, 0.5-nrow(M@assays[["RNA"]]@data)/40000)),2)
  
  # Color by signature score
  message("\nCalculating signature scores")
  # Signatures
  signatures <- lapply(as.list(signatures), function(y) y[!is.na(y)])
  signatures <- lapply(signatures, intersect, rownames(M@assays[["RNA"]]@data))
  # Average gene expression for scoreSignature
  CM.mean <- rowMeans(M@assays[["RNA"]]@data)
  signScore <- lapply(names(signatures), function(g) {
    message(g)
    scoreSignature(CM = as.matrix(M@assays[["RNA"]]@data), signatures = signatures[[g]], CM.mean = CM.mean, verbose = T)
  })
  names(signScore) <- names(signatures)
  
  # Plot all signatures - UMAP
  pdf(file = paste0(name, "_UMAP_Seurat.pdf"), width = 6, height = 6)
  par(mar=c(4, 4, 4, 4))
  for (n in names(signScore)) {
    mycol <- colItay(signScore[[n]])
    plot(M@reductions[["umap"]]@cell.embeddings, pch = 16, cex = cexsize, col = mycol, main = n)
  }
  dev.off()
  
  
  ##################################
  #                                #
  # Save New Plots and Save Object #
  #                                #
  ##################################
  
  # Plot Seurat UMAP
  pdf(file = paste0(name, "_UMAP_Seurat_Groups.pdf"), width = 6, height = 6)
  par(mar=c(2, 2, 2, 2))
  DimPlot(object = M, reduction = "umap",label= TRUE)
  dev.off()
  
  saveRDS(M, file = "Seurat.rds")
}
