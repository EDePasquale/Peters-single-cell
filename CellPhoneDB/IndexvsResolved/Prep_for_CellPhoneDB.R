##############################
#                            #
#     Prepare files for      #
#        CellphoneDB         #
# Peters - Single Cell Liver #
#     Erica DePasquale       #
#        9 Oct 2023          #
#                            #
##############################

library(Seurat)
library(SeuratObject)
library(Matrix)

# From: https://github.com/ventolab/CellphoneDB/blob/master/notebooks/0_prepare_your_data_from_Seurat.ipynb
dir.create("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellPhoneDB_IndexvsResolved")

myScopes=c("Expanded", "NotExpanded")
myScopes2=c("IndexRejection", "Resolved")
myClusters=c("CD1C+ Kupffer", "C1Q+ Kupffer", "LST1+ Kupffer", "PTPRC+ Kupffer", "IFI27+ Kupffer", "Monocyte-Derived Macrophage")

for(clus in myClusters){
  
  for(sco2 in myScopes2){
    
    for(sco in myScopes){
      
      # 0. Make working directory
      dir.create(paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellPhoneDB_IndexvsResolved/", clus))
      dir.create(paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellPhoneDB_IndexvsResolved/", clus, "/", sco, "_", sco2, "_samples"))
      setwd(paste0("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellPhoneDB_IndexvsResolved/", clus, "/", sco, "_", sco2, "_samples"))
      dir.create("Peters_5PrimeTCRBCR_counts_mtx")
      
      # 1. Load seurat object
      so = readRDS("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")
      
      # Add new classifications to metadata
      normal=c("ALP-00036", "ALP-00041-TXP", "ALP-00042-TXP", "ALP-00045-TXP", "ALP-00066", "ALP-00078")
      notACR=c("ALP-00042-BX2", "ALP-00045-BX1")
      indexRejection=c("ALP-0003-BX1", "ALP-00017", "ALP-00023-BX1", "ALP-00023-BX2", "ALP-00036-BX2")
      recurrent=c("ALP-0003-BX2",  "ALP-00023-BX3", "ALP-00036-BX3", "ALP-00036-BX4", "ALP-00036-BX5")
      resolved=c("ALP-0003-BX3", "ALP-0003-BX4", "ALP-0003-BX5", "ALP-00017-BX2")
      
      comb=c(normal, notACR, indexRejection, recurrent, resolved)
      
      Idents(object = so) <- "orig.ident"
      so=subset(x = so, idents = comb)
      
      temp=so@meta.data$orig.ident
      temp[temp %in% normal]<-"Pre-TXP"
      temp[temp %in% notACR]<-"No ACR"
      temp[temp %in% indexRejection]<-"Index Rejection"
      temp[temp %in% recurrent]<-"Recurrent"
      temp[temp %in% resolved]<-"Resolved"
      so@meta.data[["Classification"]]<-temp
      
      Idents(object = so) <- "Classification"
      Idents(object = so) <- factor(x = Idents(so), levels = c("Pre-TXP", "No ACR", "Index Rejection", "Recurrent", "Resolved"))
      
      Idents(so) = so$cluster_names_new
      
      so<-subset(so, idents=c("CD8 Effector Memory T", "CD8 Effector T", "MAIT", "CD8 TRM Activated", "Cycling T", clus)) # CD8 only
      expanded=rep("no", nrow(so@meta.data))
      expanded[which(so$clonesize>2)]<-"yes"
      expanded[which(so$cluster_names_new==clus)]<-"kupffer"
      so[["expanded"]] <-expanded
      table(so@meta.data$expanded)
      
      # 1B. Reduce to exp or not exp
      if(sco=="Expanded"){
        Idents(object = so) <- "expanded"
        so<-subset(so, idents=c("yes", "kupffer")) # FOR SUBSETTING
      }else if(sco=="NotExpanded"){
        Idents(object = so) <- "expanded"
        so<-subset(so, idents=c("no", "kupffer")) # FOR SUBSETTING
      }
      Idents(object = so) <- "cluster_names_new"
      table(so@meta.data$expanded)
      
      # 1C. Reduce index rejection or resolved
      if(sco2=="IndexRejection"){
        Idents(object = so) <- "Classification"
        so<-subset(so, idents=c("Index Rejection")) # FOR SUBSETTING
      }else if(sco2=="Resolved"){
        Idents(object = so) <- "Classification"
        so<-subset(so, idents=c("Resolved")) # FOR SUBSETTING
      }
      Idents(object = so) <- "cluster_names_new"
      table(so@meta.data$Classification)
      
      # 2. Write gene expression in mtx format
      writeMM(so@assays$RNA@counts, file = 'Peters_5PrimeTCRBCR_counts_mtx/matrix.mtx')
      write(x = rownames(so@assays$RNA@counts), file = "Peters_5PrimeTCRBCR_counts_mtx/features.tsv")
      write(x = colnames(so@assays$RNA@counts), file = "Peters_5PrimeTCRBCR_counts_mtx/barcodes.tsv")
      
      # 3. Generate your meta
      #In this example, our input is an anndata containing the cluster/celltype information in metadat$'cell_type'
      #The object also has metadat$'lineage' information wich will be used below for a hierarchical DEGs approach.
      so@meta.data$Cell = rownames(so@meta.data)
      df = so@meta.data[, c('Cell', 'cluster_names_new')]
      write.table(df, file ='Peters_5PrimeTCRBCR_meta.tsv', sep = '\t', quote = F, row.names = F)
      
      # 4. Compute DEGs (optional)
      ## OPTION 1 - compute DEGs for all cell types
      DEGs <- FindAllMarkers(so,
                             test.use = 'LR',
                             verbose = T,
                             only.pos = T,
                             random.seed = 1,
                             logfc.threshold = 0.2,
                             min.pct = 0.1,
                             return.thresh = 0.05)
      
      fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_log2FC > 0.1)
      fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')] 
      write.table(fDEGs, file ='Peters_5PrimeTCRBCR_DEGs.tsv', sep = '\t', quote = F, row.names = F)
      
    }
    
  }
  
}

# 5. Run cellphoneDB (not in R)
# cellphonedb method degs_analysis  \
#   Peters_5PrimeTCRBCR_meta.tsv  \
#   Peters_5PrimeTCRBCR_counts_mtx  \
#   Peters_5PrimeTCRBCR_DEGs.tsv  \
#   --microenvs Peters_5PrimeTCRBCR_microenviroments.tsv  \ #optional
#   --counts-data hgnc_symbol  \
#   --database database/database/cellphonedb_user_2021-06-29-11_41.db \
#   --threshold 0.1
