##############################
#                            #
#     Chord Diagrams for     #
#        CellphoneDB         #
# Peters - Single Cell Liver #
#     Erica DePasquale       #
#        9 Oct 2023          #
#                            #
##############################

# Libraries
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(chorddiag)  #devtools::install_github("mattflor/chorddiag")
library(scales)

myScopes=c("All", "Expanded", "NotExpanded")
myClusters=c("CD1C+ Kupffer", "C1Q+ Kupffer", "LST1+ Kupffer", "PTPRC+ Kupffer", "IFI27+ Kupffer", "Monocyte-Derived Macrophage")

for(clus in myClusters){
  
  for(sco in myScopes){
    
    # Load dataset
    my_data <- read.table(paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellPhoneDB/",  clus, "/", sco, "_samples/Peters_DEG_results/relevant_interactions.txt"), sep="\t", header=F)
    my_data[1,]=gsub("\\|", "_", my_data[1,])
    
    my_data2 <- my_data[,-c(1:11)]
    names_list=apply(my_data2[1,], 1, strsplit, "_")
    my_data3=as.data.frame(do.call(rbind, names_list[[1]])) #not sure why this would be rbind and not cbind...
    colnames(my_data3)<-c("rowname", "key")
    a=as.matrix(my_data2[2:nrow(my_data2),])
    a<- matrix(as.numeric(a),    # Convert to numeric matrix
               ncol = ncol(a))
    my_data3=cbind(my_data3, value=colSums(a))
    
    # parameters
    circos.clear()
    circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
    par(mar = rep(0, 4))
    
    # color palette
    myColors=read.table("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Color_Guide.txt", sep = "\t", header=T, comment.char="*")
    myColors[4,2]<-"Naïve B"
    myColors=myColors[which(myColors$Cluster_Name %in% c("CD8 Effector Memory T", "CD8 Effector T", "MAIT", "CD8 TRM Activated", "Cycling T", clus)),]
    myColors=myColors[order(myColors$Cluster_Name, decreasing=F),]
    mycolor=myColors$Cluster_Color
    
    pdf(file = paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellPhoneDB/",  clus, "/", sco, "_samples/Peters_DEG_results/Chord_allInteractions.pdf"), width = 8, height = 8)
    par(mar=c(2, 2, 2, 2))
    # Base plot
    chordDiagram(
      x=my_data3,
      grid.col = mycolor,
      transparency = 0.25,
      directional = 1,
      direction.type = c("arrows", "diffHeight"), 
      diffHeight  = -0.04,
      annotationTrack = "grid", 
      annotationTrackHeight = c(0.05, 0.1),
      link.arr.type = "big.arrow", 
      link.sort = TRUE, 
      link.largest.ontop = TRUE)
    
    # Add text and axis
    circos.trackPlotRegion(
      track.index = 1, 
      bg.border = NA, 
      panel.fun = function(x, y) {
        
        xlim = get.cell.meta.data("xlim")
        sector.index = get.cell.meta.data("sector.index")
        
        # Add names to the sector. 
        circos.text(
          x = mean(xlim), 
          y = 3.2, 
          labels = sector.index, 
          facing = "bending", 
          cex = 0.8
        )
        
        # Add graduation on axis
        circos.axis(
          h = "top", 
          labels.niceFacing = FALSE)
      }
    )
    dev.off()
  }
}


