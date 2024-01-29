library(Seurat)

############
# EXPANDED #
############

# Read in object
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/")
M <- readRDS("Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")

# Subset to CD8
M<-subset(M, idents=c("CD8 Effector Memory T", "CD8 Effector T", "MAIT", "CD8 TRM Activated", "Cycling T")) # CD8 only

# Subset to expanded
expanded=rep("no", nrow(M@meta.data))
expanded[which(M$clonesize>2)]<-"yes"
M[["expanded"]] <-expanded
table(M@meta.data$expanded)
# no  yes 
# 2949  566
M<-subset(M, subset = expanded == "yes")


###############################
# Number of cells with 2 TRAs #
###############################

myfunc<-function(myString){ # function to parse the TRA or TRB value (split by comma) and tell me how many TRA or TRB it contains
  temp=unlist(strsplit(myString, split=","))
  if(is.na(temp[1])){return(NA)}else{length(temp)}
}
NumTRA=sapply(M@meta.data[["TRA"]], myfunc)
table(NumTRA, useNA = "always")
#     1    2  <NA> 
#   441   58   67 <- 58 cells of 566 (10.25%) with 2 TRAs

table(NumTRA, M@meta.data$ACRType, useNA = "always") # By ACR type
#         Donor   Late ACR Not ACR    Resolved Late ACR
#   1     22      334      39         46                
#   2     0       49       0          9 <- Number of cells in each ACR type that have 2 TRAs               
#  <NA>   8       51       5          3                 

# How many of these have a TRB?
NumTRB=sapply(M@meta.data[["TRB"]], myfunc)
X=cbind(NumTRA, NumTRB)
Y=X[which(NumTRA>1),]

########################################
# Number of cells with 1 TRA and 1 TRB #
########################################

NumTRB=sapply(M@meta.data[["TRB"]], myfunc)
table(NumTRB, useNA = "always")
#     1    2  <NA> 
#   539    6   21

results=as.data.frame(cbind(NumTRA, NumTRB, Has1OfEach=rep(0, length(NumTRA))))
ind=which(results[,1]=="1" & results[,2]=="1")
results$Has1OfEach[ind]<-1

table(results$Has1OfEach, useNA = "always")
#   0    1      NA> 
#   152  414    0 <- 414 cells of 566 (73.14%) with 1 TRA and 1 TRB

table(results$Has1OfEach, M@meta.data$ACRType, useNA = "always") # By ACR type
#           Donor   Late ACR   Not ACR   Resolved Late ACR
#   0       10      118        5         19    
#   1       20      316        39        39    <- Number of cells in each ACR type that have 1 TRA and 1 TRB
#   <NA>    0       0          0         0     


#############
#  ALL CD8  #
#############

# Read in object
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/")
M <- readRDS("Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")

# Subset to CD8
M<-subset(M, idents=c("CD8 Effector Memory T", "CD8 Effector T", "MAIT", "CD8 TRM Activated", "Cycling T")) # CD8 only

# Subset to expanded
expanded=rep("no", nrow(M@meta.data))
expanded[which(M$clonesize>2)]<-"yes"
M[["expanded"]] <-expanded
table(M@meta.data$expanded)
# no  yes 
# 2949  566
# M<-subset(M, subset = expanded == "yes")


###############################
# Number of cells with 2 TRAs #
###############################

myfunc<-function(myString){ # function to parse the TRA or TRB value (split by comma) and tell me how many TRA or TRB it contains
  temp=unlist(strsplit(myString, split=","))
  if(is.na(temp[1])){return(NA)}else{length(temp)}
}
NumTRA=sapply(M@meta.data[["TRA"]], myfunc)
table(NumTRA, useNA = "always")
# 1     2   <NA> 
# 1829  185 1501 <- 185 cells of 3515 (5.26%) with 2 TRAs

table(NumTRA, M@meta.data$ACRType, useNA = "always") # By ACR type
#         Donor   Late ACR  Not ACR   Resolved Late ACR
#   1     130     1335      112       252   
#   2     10      147       3         25  <- Number of cells in each ACR type that have 2 TRAs  
# <NA>    321     798       136       246 

# How many of these have a TRB?
NumTRB=sapply(M@meta.data[["TRB"]], myfunc)
X=cbind(NumTRA, NumTRB)
Y=X[which(NumTRA>1),]

########################################
# Number of cells with 1 TRA and 1 TRB #
########################################

NumTRB=sapply(M@meta.data[["TRB"]], myfunc)
table(NumTRB, useNA = "always")
#   1       2   <NA> 
#   2198    29  1288 

results=as.data.frame(cbind(NumTRA, NumTRB, Has1OfEach=rep(0, length(NumTRA))))
ind=which(results[,1]=="1" & results[,2]=="1")
results$Has1OfEach[ind]<-1

table(results$Has1OfEach, useNA = "always")
#   0    1      NA> 
#   1774 1741    0  <- 1741 cells of 3515 (49.53%) with 1 TRA and 1 TRB

table(results$Has1OfEach, M@meta.data$ACRType, useNA = "always") # By ACR type
#   Donor   Late ACR  Not ACR   Resolved  Late ACR 
#   0       338       999       143       294   
#   1       123       1281      108       229   <- Number of cells in each ACR type that have 1 TRA and 1 TRB
#   <NA>    0         0         0         0    



################
#  ALL TCELLS  #
################

# Read in object
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/")
M <- readRDS("Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")

# Reduce to T cells
M<-subset(M, idents=c("CD8 Effector Memory T", "CD8 Effector T", "MAIT", "CD8 TRM Activated", "Cycling T", "CD4 naive T", "Gamma Delta T")) 


###############################
# Number of cells with 2 TRAs #
###############################

myfunc<-function(myString){ # function to parse the TRA or TRB value (split by comma) and tell me how many TRA or TRB it contains
  temp=unlist(strsplit(myString, split=","))
  if(is.na(temp[1])){return(NA)}else{length(temp)}
}
NumTRA=sapply(M@meta.data[["TRA"]], myfunc)
table(NumTRA, useNA = "always")
#      1    2 <NA> 
#   2189  226 1932  <- 226 cells of 4247 (5.20%) with 2 TRAs

table(NumTRA, M@meta.data$ACRType, useNA = "always") # By ACR type
#       Donor    Late ACR   Not ACR   Resolved Late ACR
# 1      160     1589       141       299    
# 2       13     174        7         32    <- Number of cells in each ACR type that have 2 TRAs  
# <NA>   358     1112       164       298    

# How many of these have a TRB?
NumTRB=sapply(M@meta.data[["TRB"]], myfunc)
X=cbind(NumTRA, NumTRB)
Y=X[which(NumTRA>1),]

########################################
# Number of cells with 1 TRA and 1 TRB #
########################################

NumTRB=sapply(M@meta.data[["TRB"]], myfunc)
table(NumTRB, useNA = "always")
#   1       2   <NA> 
#   2631    43  1673  

results=as.data.frame(cbind(NumTRA, NumTRB, Has1OfEach=rep(0, length(NumTRA))))
ind=which(results[,1]=="1" & results[,2]=="1")
results$Has1OfEach[ind]<-1

table(results$Has1OfEach, useNA = "always")
#   0    1      NA> 
#   2273 2074   0  <- 2074 cells of 4347 (47.71%) with 1 TRA and 1 TRB

table(results$Has1OfEach, M@meta.data$ACRType, useNA = "always") # By ACR type
#        Donor   Late ACR  Not ACR   Resolved  Late ACR 
# 0      381     1358      177       357  
# 1      150     1517      135       272  <- Number of cells in each ACR type that have 1 TRA and 1 TRB
# <NA>   0       0         0         0    

#####################################################
# Of cells with a TRA or TRB, how many have just 1? #
#####################################################

toRemove=intersect(which(is.na(X[,1])), which(is.na(X[,2])))
Z=X[-toRemove,]

sum(is.na(Z[,2])) #76 TRA only, 76/2750 = 2.76%
sum(is.na(Z[,1])) #335 TRB only, 335/2750 = 12.18%

#have both?
toKeep=intersect(which(!is.na(X[,1])), which(!is.na(X[,2])))
A=X[toKeep,] #2339/2750 = 85.05%

