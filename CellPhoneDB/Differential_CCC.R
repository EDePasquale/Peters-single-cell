##############################
#                            #
#      Differential CCC      #
#        CellphoneDB         #
# Peters - Single Cell Liver #
#     Erica DePasquale       #
#         9 Oct 2023         #
#                            #
##############################

# Load libraries
library(dplyr)

myClusters=c("CD1C+ Kupffer", "C1Q+ Kupffer", "LST1+ Kupffer", "PTPRC+ Kupffer", "IFI27+ Kupffer", "Monocyte-Derived Macrophage")

for(clus in myClusters){
  
  # Load datasets
  Exp <- read.table(paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellPhoneDB/",  clus, "/Expanded_samples/Peters_DEG_results/relevant_interactions.txt"), sep="\t", header=T)
  Exp[1,]=gsub("\\|", "_", Exp[1,])
  NotExp <- read.table(paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellPhoneDB/",  clus, "/NotExpanded_samples/Peters_DEG_results/relevant_interactions.txt"), sep="\t", header=T)
  NotExp[1,]=gsub("\\|", "_", NotExp[1,])
  
  Shared=intersect(Exp$id_cp_interaction, NotExp$id_cp_interaction)
  UniqueExp=setdiff(Exp$id_cp_interaction, NotExp$id_cp_interaction)
  UniqueNotExp=setdiff(NotExp$id_cp_interaction, Exp$id_cp_interaction)
  
  ##########
  # Shared #
  ##########
  
  Exp_Shared=Exp[Exp$id_cp_interaction %in% Shared,]
  row.names(Exp_Shared)<-Shared
  NotExp_Shared=NotExp[NotExp$id_cp_interaction %in% Shared,]
  row.names(NotExp_Shared)<-Shared
  
  my_results=as.data.frame(matrix(ncol=14, nrow=1))
  colnames(my_results)<-c("Cluster_Pair", "id_cp_interaction", "interacting_pair", "partner_a", "partner_b", "gene_a", "gene_b", 
                          "secreted", "receptor_a", "receptor_b", "annotation_strategy", "is_integrin", "Expanded", "NotExpanded")
  for(i in 13:ncol(Exp_Shared)){
    Ind=which(Exp_Shared[,i] != NotExp_Shared[,i])
    my_results=rbind(my_results, cbind(Cluster_Pair=rep(colnames(Exp_Shared)[i], length(Ind)),
                                       Exp_Shared[Ind,1:11],
                                       Expanded=Exp_Shared[Ind,i],
                                       NotExpanded=NotExp_Shared[Ind,i]))
  }
  my_results<-my_results[-1,]
  
  ##########
  # Unique #
  ##########
  
  # Unique interactions to Expanded. Set all expanded to 1 and not expanded to 0
  X=Exp[Exp$id_cp_interaction %in% UniqueExp,]
  for(i in 1:nrow(X)){
    Ind=which(X[i,]==1)
    my_results=rbind(my_results, cbind(Cluster_Pair=colnames(X)[Ind],
          do.call("rbind", replicate(length(Ind), X[i, 1:11], simplify = FALSE)),
          Expanded=rep(1, length(Ind)),
          NotExpanded=rep(0, length(Ind))))
  }
  
  # Unique interactions to NotExpanded. Set all not expanded to 1 and expanded to 0
  Y=NotExp[NotExp$id_cp_interaction %in% UniqueNotExp,]
  for(i in 1:nrow(Y)){
    Ind=which(Y[i,]==1)
    my_results=rbind(my_results, cbind(Cluster_Pair=colnames(Y)[Ind],
                                       do.call("rbind", replicate(length(Ind), Y[i, 1:11], simplify = FALSE)),
                                       Expanded=rep(0, length(Ind)),
                                       NotExpanded=rep(1, length(Ind))))
  }
  
  
  write.table(my_results, paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellPhoneDB_", clus, "_Differential_Interactions.txt"), sep="\t", quote=F)
}