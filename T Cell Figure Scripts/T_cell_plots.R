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
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/")
M <- readRDS("Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")

# Subset to T-cells only
Idents(object = M) <- "cluster_names_new"
Clusters_Keep=c("CD8 Effector Memory T", "CD4 naive T", "CD8 TRM Activated", "Gamma Delta T", "MAIT", "CD8 Effector T", "Cycling T")
M_T=subset(x = M, idents = Clusters_Keep)

# Make TSNE plot
M_T <- RunPCA(M_T, npcs = 30, verbose = FALSE)
M_T <- RunUMAP(M_T, reduction = "pca", dims = 1:30)
M_T <- RunTSNE(M_T, reduction = "pca", dims = 1:30)

# Pull in new colors
myColors=read.table("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Color_Guide.txt", sep = "\t", header=T, comment.char="*")
myColors[4,2]<-"Naïve B"
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

# pal <- viridis(n = 10, option = "C", direction = 1)
# 
# # Plot by Clone size
# pdf(file = "TSNE_Tcell_Clonesize_Viridis.pdf", width = 8, height = 8)
# par(mar=c(2, 2, 2, 2))
# FeaturePlot(M_T, reduction="tsne", features="clonesize", raster=F, pt.size=1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PuOr")))
# dev.off()
# 
# # Full clusters clone size
# pdf(file = "TSNE_Clonesize_Viridis.pdf", width = 8, height = 8)
# par(mar=c(2, 2, 2, 2))
# FeaturePlot(M, reduction="tsne", features ="clonesize", raster=F, cols=c("lightgrey", "red"))
# dev.off()

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

# Make plot
pdf("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellTypeFrequencies_T_cell_byACRType_new_colors.pdf", width = 9, height = 9)
par(mar = c(8,4,4,16), xpd = T)
  barplot(as.matrix(data[nrow(data):1,]), col = rev(my_colors), xaxt = "n", ylab = "Population frequency (%)", border = NA)
  axis(side = 1, at = seq(1,ncol(data))*1.2-0.5, labels = colnames(data), las = 2)
  legend(x = ncol(data)*1.2+0.5, y = 100, legend = row.names(data), fill = my_colors, bty = "n", border = NA)
dev.off()

# Normalize data using log for visualization purposes
M_T <- NormalizeData(M_T, assay="RNA")
all.genes <- rownames(M_T)
M_T <- ScaleData(M_T, features = all.genes, assay="RNA")

#myGenes=unique(c("CD8A", "CD4", "CD3G", "TRAV1-2", "PTPRC", "CD27", "CCR7", "XCL1", "TCF7"))
myGenes2=unique(c("CD3G", "CD4", "CCR7", "CD8A", "TRAV1-2", "CD27", "CD44", "CD69", "XCL1", "TRDV1"))
#myGenes3=unique(c("CD3G", "CD4", "CCR7", "CD8A", "TRAV1-2", "CD27", "CD44", "CD69", "XCL1", "TRDV1", "GZMK", "GZMB", "PRF1", "CX3CR1", "CXCR6", "TBX21"))

# pdf(file = "Tcell_Manual_RNA_Stacked_Violin_Names.pdf", width = 11, height = 6)
# par(mar=c(2, 2, 2, 2))
#   VlnPlot(M_T, assay="RNA", group.by="cluster_names_new", stack=T, features=myGenes, flip=T, sort=T) + NoLegend()
# dev.off()

# New Code for Dot Plot (to display the values as non-scaled to best represent the actual expression values)
########
DotPlot <- function(
    object,
    assay = NULL,
    features,
    cols = c("lightgrey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = 'radius',
    scale.min = NA,
    scale.max = NA) {
  assay <- assay %||% SeuratObject:::DefaultAssay(object = object)
  SeuratObject:::DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning(
        "Some feature groups are unnamed.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = SeuratObject:::CellsByIdentities(object = object, idents = idents))
  
  data.features <- SeuratObject:::FetchData(object = object, vars = features, cells = cells, slot="data")
  data.features$id <- if (is.null(x = group.by)) {
    SeuratObject:::Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = Seurat:::PercentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(
      what = rbind,
      args = lapply(X = data.plot, FUN = unlist)
    )
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning(
      "Only one identity present, the expression values will be not scaled",
      call. = FALSE,
      immediate. = TRUE
    )
  } else if (ngroup < 5 & scale) {
    warning(
      "Scaling data with a low number of groups may produce misleading results",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      if (scale) {
        data.use <- scale(x = data.use)
        data.use <- Seurat:::MinMax(data = data.use, min = col.min, max = col.max)
      } else {
        data.use <- log1p(x = data.use)
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(
      X = as.character(x = data.plot$id),
      FUN = gsub,
      FUN.VALUE = character(length = 1L),
      pattern =  paste0(
        '^((',
        paste(sort(x = levels(x = object), decreasing = TRUE), collapse = '|'),
        ')_)'
      ),
      replacement = '',
      USE.NAMES = FALSE
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        return(colorRampPalette(colors = c('grey', color))(20)[value])
      },
      color = cols[splits.use],
      value = avg.exp.scaled
    )
  }
  color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp') #changed
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(
      x = 'Features',
      y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
    ) +
    theme_cowplot()
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(
      facets = ~feature.groups,
      scales = "free_x",
      space = "free_x",
      switch = "y"
    ) + theme(
      panel.spacing = unit(x = 1, units = "lines"),
      strip.background = element_blank()
    )
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  } else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
  }
  return(plot)
}

########

pdf(file = "Tcell_Manual_RNA_DotPlot_Names.pdf", width = 10, height = 5)
par(mar=c(4, 4, 4, 4))
  DotPlot(object=M_T, assay="RNA", features = myGenes2, dot.scale=10, scale=F) + labs(y ="Cluster", x=NULL) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

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
  clonotypes <- read.csv(paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/", expt, "/clonotypes.csv"), sep=',', header=T)
  filtered_contig <- read.csv(paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/TCRs/", expt, "/filtered_contig_annotations.csv"), sep=',', header=T)
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
clones_matrix_wide_sums=clones_matrix_wide_sums[,order(-clones_matrix_wide_sums[nrow(clones_matrix_wide_sums),])]

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
         filename = "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/heatmap_all_shared_clones_updated.pdf")

pheatmap(t(clones_matrix_wide_redu),
         cluster_cols = F,
         color=c("gray92", "red"),
         gaps_col=c(5,6,8,11,13,16,21,22,24,26,27,28,29),
         width=11.5,
         height = 20,
         angle_col=45,
         filename = "/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/heatmap_all_shared_clones_long_updated.pdf")















#################
# Additional Violin Plots for Anna

# dir.create("New_Violin_Plots")
# setwd("New_Violin_Plots")
# 
# M_T <- NormalizeData(M_T, assay="RNA")
# all.genes <- rownames(M_T)
# M_T <- ScaleData(M_T, features = all.genes, assay="RNA")
# 
# # T cell subset
# myGenes=unique(c("PDCD1", "CXCR6", "IFNG", "CD27", "CD70", "TNFRSF9", "EOMES", "TBX21", "TCF7", "XCL1", "GZMB", "GZMK", "PRF1", "CX3CR1", "TRAV1-2", "TIGIT", "NCAM1", "PRDM1"))
# 
# pdf(file = "Tcell_subset_SCT_Stacked_Violin_Names.pdf", width = 9, height = 9)
# par(mar=c(2, 2, 2, 2))
#   VlnPlot(M_T, assay="SCT", group.by="cluster_names_new", stack=T, features=myGenes, flip=T, sort=T) + NoLegend()
# dev.off()
# 
# pdf(file = "Tcell_subset_RNA_Stacked_Violin_Names.pdf", width = 9, height = 9)
# par(mar=c(2, 2, 2, 2))
#   VlnPlot(M_T, assay="RNA", group.by="cluster_names_new", stack=T, features=myGenes, flip=T, sort=T) + NoLegend()
# dev.off()
# 
# # Activated
# myGenes=unique(c("PFR1", "GZMB", "IFNG", "HLA-DRA", "CX3CR1", "TBX21", "MKI67", "PCNA"))
# 
# pdf(file = "Tcell_activated_SCT_Stacked_Violin_Names.pdf", width = 9, height = 9)
# par(mar=c(2, 2, 2, 2))
#   VlnPlot(M_T, assay="SCT", group.by="cluster_names_new", stack=T, features=myGenes, flip=T, sort=T) + NoLegend()
# dev.off()
# 
# pdf(file = "Tcell_activated_RNA_Stacked_Violin_Names.pdf", width = 9, height = 9)
# par(mar=c(2, 2, 2, 2))
#   VlnPlot(M_T, assay="RNA", group.by="cluster_names_new", stack=T, features=myGenes, flip=T, sort=T) + NoLegend()
# dev.off()
# 
# # Resident/Circulation
# myGenes=unique(c("ZNF683", "PRDM1", "CD69", "ITGAE", "CXCR6", "S1PR1", "S1PR5", "SELL"))
# 
# pdf(file = "Tcell_residentCirculation_SCT_Stacked_Violin_Names.pdf", width = 9, height = 9)
# par(mar=c(2, 2, 2, 2))
#   VlnPlot(M_T, assay="SCT", group.by="cluster_names_new", stack=T, features=myGenes, flip=T, sort=T) + NoLegend()
# dev.off()
# 
# pdf(file = "Tcell_residentCirculation_RNA_Stacked_Violin_Names.pdf", width = 9, height = 9)
# par(mar=c(2, 2, 2, 2))
#   VlnPlot(M_T, assay="RNA", group.by="cluster_names_new", stack=T, features=myGenes, flip=T, sort=T) + NoLegend()
# dev.off()
# 
# # Exhaustion
# myGenes=unique(c("TOX", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "NR4A1", "NKG7"))
# 
# pdf(file = "Tcell_exhaustion_SCT_Stacked_Violin_Names.pdf", width = 9, height = 9)
# par(mar=c(2, 2, 2, 2))
#   VlnPlot(M_T, assay="SCT", group.by="cluster_names_new", stack=T, features=myGenes, flip=T, sort=T) + NoLegend()
# dev.off()
# 
# pdf(file = "Tcell_exhaustion_RNA_Stacked_Violin_Names.pdf", width = 9, height = 9)
# par(mar=c(2, 2, 2, 2))
#   VlnPlot(M_T, assay="RNA", group.by="cluster_names_new", stack=T, features=myGenes, flip=T, sort=T) + NoLegend()
# dev.off()
# 
# 
# ##############################
# # Even more plots for Anna
# 
# dir.create("New_Plots_pt2")
#setwd("New_Plots_pt2")
# 
# M_T <- NormalizeData(M_T, assay="RNA")
# all.genes <- rownames(M_T)
# M_T <- ScaleData(M_T, features = all.genes, assay="RNA")
# 
# Add new classifications to metadata
# normal=c("ALP-00036", "ALP-00041-TXP", "ALP-00042-TXP", "ALP-00045-TXP", "ALP-00066", "ALP-00078")
# notACR=c("ALP-00042-BX2", "ALP-00045-BX1")
# indexRejection=c("ALP-0003-BX1", "ALP-00017", "ALP-00023-BX1", "ALP-00023-BX2", "ALP-00036-BX2")
# followUp=c("ALP-0003-BX3", "ALP-0003-BX4", "ALP-0003-BX5", "ALP-00017-BX2", "ALP-00023-BX3", "ALP-00036-BX3", "ALP-00036-BX4", "ALP-00036-BX5")
# comb=c(normal, notACR, indexRejection, followUp)

# normal=c("ALP-00036", "ALP-00041-TXP", "ALP-00042-TXP", "ALP-00045-TXP", "ALP-00066", "ALP-00078")
# notACR=c("ALP-00042-BX2", "ALP-00045-BX1")
# indexRejection=c("ALP-0003-BX1", "ALP-00017", "ALP-00023-BX1", "ALP-00030-BX1", "ALP-00048", "ALP-00069")
# recurrent=c("ALP-0003-BX2", "ALP-00011-BX2", "ALP-00017-BX2", "ALP-00023-BX2", "ALP-00023-BX3", "ALP-00033-BX2", "ALP-00033-BX3", "ALP-00033-BX4", "ALP-00036-BX2", "ALP-00036-BX3", "ALP-00036-BX4", "ALP-00036-BX5")
# resolved=c("ALP-0003-BX3", "ALP-0003-BX4", "ALP-0003-BX5", "ALP-00030-BX4")

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

# M_T@meta.data$Classification <- factor(M_T@meta.data$Classification, levels=c("Normal", "Not ACR", "Index Rejection", "Follow Up"))
# 
# # Activated
# myGenes_activated=unique(c("PFR1", "GZMB", "IFNG", "HLA-DRA", "CX3CR1", "TBX21", "MKI67", "PCNA"))
# # Resident/Circulation
# myGenes_resident=unique(c("ZNF683", "PRDM1", "CD69", "ITGAE", "CXCR6", "S1PR1", "S1PR5", "SELL"))
# # Exhaustion
# myGenes_exhausted=unique(c("TOX", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "NR4A1", "NKG7"))
# 
# 
# ################
# #              #
# #     ALL      #
# #              #
# ################
# 
# pdf(file = "Tcell_all_activated_RNA_Violin.pdf", width = 6, height = 9)
# par(mar=c(2, 2, 2, 2))
#   VlnPlot(M_T, assay="RNA", group.by="Classification", stack=T, features=myGenes_activated, flip=T, sort=F) + NoLegend()
# dev.off()
# 
# pdf(file = "Tcell_all_resident_RNA_Violin.pdf", width = 6, height = 9)
# par(mar=c(2, 2, 2, 2))
#   VlnPlot(M_T, assay="RNA", group.by="Classification", stack=T, features=myGenes_resident, flip=T, sort=F) + NoLegend()
# dev.off()
# 
# pdf(file = "Tcell_all_exhausted_RNA_Violin.pdf", width = 6, height = 9)
# par(mar=c(2, 2, 2, 2))
#   VlnPlot(M_T, assay="RNA", group.by="Classification", stack=T, features=myGenes_exhausted, flip=T, sort=F) + NoLegend()
# dev.off()
# 
# 
# #######################
# #                     #
# #     SIDE BY SIDE    # 
# #                     #
# #######################
# 
expanded=rep("Not Expanded", nrow(M_T@meta.data))
expanded[which(M_T$clonesize>2)]<-"Expanded"
M_T[["expanded"]] <-expanded
table(M_T@meta.data$expanded)

pdf(file = "Tcell_sidebyside_activated_RNA_Violin.pdf", width = 8, height = 9)
par(mar=c(2, 2, 2, 2))
  VlnPlot(M_T, assay="RNA", group.by="Classification", split.by="expanded", stack=T, features=myGenes_activated, flip=T, sort=F)
dev.off()

pdf(file = "Tcell_sidebyside_resident_RNA_Violin.pdf", width = 8, height = 9)
par(mar=c(2, 2, 2, 2))
  VlnPlot(M_T, assay="RNA", group.by="Classification", split.by="expanded", stack=T, features=myGenes_resident, flip=T, sort=F)
dev.off()

pdf(file = "Tcell_sidebyside_exhausted_RNA_Violin.pdf", width = 8, height = 9)
par(mar=c(2, 2, 2, 2))
  VlnPlot(M_T, assay="RNA", group.by="Classification", split.by="expanded", stack=T, features=myGenes_exhausted, flip=T, sort=F)
dev.off()

myGenes_new=c("GZMB", "HLA-DRA", "LAG3",  "TIGIT", "NKG7", "CD69", "CXCR6")
pdf(file = "Tcell_sidebyside_new_RNA_Violin_newClassifications_justSerial.pdf", width = 8, height = 9)
par(mar=c(2, 2, 2, 2))
  VlnPlot(M_T, assay="RNA", group.by="Classification", split.by="expanded", stack=T, features=myGenes_new, flip=T, sort=F)
dev.off()

