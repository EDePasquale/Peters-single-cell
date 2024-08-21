# Libraries
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
# BiocManager::install("ComplexHeatmap")
#devtools::install_github("saeyslab/nichenetr")
#devtools::install_github("saeyslab/multinichenetr")
library(igraph)
library(nichenetr)
library(multinichenetr)
library(Seurat)
library(SeuratObject)
library(plyr)
library(circlize)
library(writexl)

# Options
options(timeout = 500)

# Constants
organism = "human"
sample_id = "orig.ident"
group_id = "expansion_diagnosis"
celltype_id = "ident"
batches = NA
covariates = NA
min_cells = 3 # was 10
min_sample_prop = 0.1 # was 50%, made it 10%, i.e., really only needs to be expressed in one sample
fraction_cutoff = 0.05 # was 5%
logFC_threshold = 0.50 # was 50%
p_val_threshold = 0.05
p_val_adj = FALSE 
top_n_target = 250
n.cores = 2
scenario = "regular"
ligand_activity_down = FALSE

# Sets
myScopes=c("Expanded", "NotExpanded")
myScopes2=c("IndexRejection", "Resolved")
myClusters=c("CD1C+ Kupffer", "C1Q+ Kupffer", "LST1+ Kupffer", "PTPRC+ Kupffer", "IFI27+ Kupffer", "Monocyte-Derived Macrophage", "CD8 Effector T")
myClusters_nice=c("CD1CKupffer", "C1QKupffer", "LST1Kupffer", "PTPRCKupffer", "IFI27Kupffer", "MonocyteDerivedMacrophage", "CD8EffectorT")

# Define networks
if(organism == "human"){
  
  # lr_network_all = 
  #   readRDS(url(
  #     "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
  # )) %>% 
  lr_network_all = readRDS("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/MultiNicheNet/lr_network_human_allInfo_30112033.rds"
  ) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  # ligand_target_matrix = readRDS(url(
  #   "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
  # ))
  ligand_target_matrix = readRDS("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/MultiNicheNet/ligand_target_matrix_nsga2r_final.rds")
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
} else if(organism == "mouse"){
  
  lr_network_all = readRDS(url(
    "https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds"
  )) %>%
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
  ))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
}

# Test SCE object
# sce = readRDS(url(
#   "https://zenodo.org/record/8010790/files/sce_subset_breastcancer.rds"
# ))
# sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()
#as.SingleCellExperiment(x, assay = NULL, ...)

# Read in our data
so = readRDS("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/Seurat_Liver_30_subcluster_names_meta_TCR_BCR_Krish.rds")
Idents(so) = so$cluster_names_new
DefaultAssay(so) <- "RNA"

# Fix cell type and sample names
so@meta.data$cluster_names_new=mapvalues(so@meta.data$cluster_names_new, from=myClusters, to=myClusters_nice)
so@meta.data$orig.ident<-make.names(so@meta.data$orig.ident)

# Add new classifications to metadata
normal=c("ALP.00036", "ALP.00041.TXP", "ALP.00042.TXP", "ALP.00045.TXP", "ALP.00066", "ALP.00078")
notACR=c("ALP.00042.BX2", "ALP.00045.BX1")
indexRejection=c("ALP.0003.BX1", "ALP.00017", "ALP.00023.BX1", "ALP.00023.BX2", "ALP.00036.BX2")
recurrent=c("ALP.0003.BX2",  "ALP.00023.BX3", "ALP.00036.BX3", "ALP.00036.BX4", "ALP.00036.BX5")
resolved=c("ALP.0003.BX3", "ALP.0003.BX4", "ALP.0003.BX5", "ALP.00017.BX2")
comb=c(normal, notACR, indexRejection, recurrent, resolved)

Idents(object = so) <- "orig.ident"
so=subset(x = so, idents = comb) # include only samples for which we have classifications

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

# Reduce to CD8T and Kupffer
so<-subset(so, idents=c(myClusters_nice))
expanded=rep("no", nrow(so@meta.data))
expanded[which(so$clonesize>2)]<-"yes"
expanded[which(so$cluster_names_new %in% setdiff(myClusters_nice, "CD8EffectorT"))]<-"kupffer"
so[["expanded"]] <-expanded

# Create object list
Seurat_list=list()

# Reduce to relevant data sets and save in list
for(sco in myScopes){
  for(sco2 in myScopes2){
    # 1A. Reduce to exp or not exp
    if(sco=="Expanded"){
      Idents(object = so) <- "expanded"
      so_temp<-subset(so, idents=c("yes", "kupffer")) # FOR SUBSETTING
    }else if(sco=="NotExpanded"){
      Idents(object = so) <- "expanded"
      so_temp<-subset(so, idents=c("no", "kupffer")) # FOR SUBSETTING
    }
    
    # 1B. Reduce index rejection or resolved
    if(sco2=="IndexRejection"){
      Idents(object = so_temp) <- "Classification"
      so_temp<-subset(so_temp, idents=c("Index Rejection")) # FOR SUBSETTING
    }else if(sco2=="Resolved"){
      Idents(object = so_temp) <- "Classification"
      so_temp<-subset(so_temp, idents=c("Resolved")) # FOR SUBSETTING
    }
    Idents(object = so_temp) <- "cluster_names_new"
    
    # 1C. Add meta data
    so_temp[["expansion_diagnosis"]] <- rep(paste0(sco, "_", sco2), nrow(so_temp@meta.data))
    so_temp[["orig.ident"]] <- sapply(so_temp[["orig.ident"]], paste0, ".", sco, ".", sco2)
    
    # 1D. Save object
    Seurat_list[[paste0(sco, "_", sco2)]]<-so_temp
  }
}

# Now I need to merge these into a single object and set the group level meta data. This will mean that the same cells (kupffer) will
# be present in two different samples, the expanded and non-expanded sample. I think this will be fine, but I'm still questioning it.
# I also think I will do this merging at the Seurat level and then convert the merged object with relevant meta data slots into the SCE.

# Merge objects
Seurat_combined <- merge(x = Seurat_list[[1]], y = c(Seurat_list[[2]], Seurat_list[[3]], Seurat_list[[4]]), add.cell.ids = c("EI", "ER", "NeI","NeR"))

# Convert to SCE
SCE_combined=as.SingleCellExperiment(Seurat_combined)

# Define senders and receivers
senders_oi = SummarizedExperiment::colData(SCE_combined)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(SCE_combined)[,celltype_id] %>% unique()
SCE_combined = SCE_combined[, SummarizedExperiment::colData(SCE_combined)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]

# Rename groups for simplicity of following tutorial
SummarizedExperiment::colData(SCE_combined)[,group_id][SummarizedExperiment::colData(SCE_combined)[,group_id] == "Expanded_IndexRejection"] = "G1.T1" # G1 = expanded, T1 = IR
SummarizedExperiment::colData(SCE_combined)[,group_id][SummarizedExperiment::colData(SCE_combined)[,group_id] == "NotExpanded_IndexRejection"] = "G2.T1" # G2 = not expanded, T1 = IR
SummarizedExperiment::colData(SCE_combined)[,group_id][SummarizedExperiment::colData(SCE_combined)[,group_id] == "Expanded_Resolved"] = "G1.T2" # G1 = expanded, T2 = Resolved
SummarizedExperiment::colData(SCE_combined)[,group_id][SummarizedExperiment::colData(SCE_combined)[,group_id] == "NotExpanded_Resolved"] = "G2.T2" # G2 = not expanded, T2 = Resolved



####################################
# <><><><> Simple test #1 <><><><> #
####################################

## 1) G1.T2-G1.T1: In expanded cells, what is the difference between IR and Resolved?

# Set contrasts and table
contrasts_oi = c("'G1.T2-G1.T1','G1.T1-G1.T2'") 

contrast_tbl = tibble(
  contrast = c("G1.T2-G1.T1", "G1.T1-G1.T2"), 
  group = c("G1.T2","G1.T1")
)

#*Cell-type filtering: determine which cell types are sufficiently present*
abundance_info = get_abundance_info(
  sce = SCE_combined, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  batches = batches
)
abundance_info$abund_plot_sample
# There are many kupffer subtype/sample combinations that don't have a sufficent number of cells to do analysis (min = 10)

#*Gene filtering: determine which genes are sufficiently expressed in each present cell type*
frq_list = get_frac_exprs(
  sce = SCE_combined, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)
# Note: make sure the DefaultAssay is set to "RNA" at the beginning of the process

# Now only keep genes that are expressed by at least one cell type:
genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% pull(gene) %>% unique() 
SCE_combined = SCE_combined[genes_oi, ]

# *Pseudobulk expression calculation: determine and normalize per-sample pseudobulk expression levels for each expressed gene in each present cell type*
abundance_expression_info = process_abundance_expression_info(
  sce = SCE_combined, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)
#Warning message:
#In DGEList.default(pb@assays@data[[celltype_oi]]) :
# At least one library size is zero

# *Differential expression (DE) analysis: determine which genes are differentially expressed*
DE_info_group1 = get_DE_info(
  sce = SCE_combined, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)
#[1] "DE analysis is done:"
# [1] "included cell types are:"
# [1] "LST1Kupffer"
# [1] "excluded cell types are:"
# [1] "CD1CKupffer"               "C1QKupffer"                "MonocyteDerivedMacrophage" "CD8EffectorT"             
# [5] "PTPRCKupffer"              "IFI27Kupffer"             
# TODO: this is obviously not okay!
# I'm going to try cutting min_cells in half (10 -> 5)
# [1] "DE analysis is done:"
# [1] "included cell types are:"
# [1] "MonocyteDerivedMacrophage" "LST1Kupffer"               "CD8EffectorT"              "IFI27Kupffer"             
# [1] "excluded cell types are:"
# [1] "CD1CKupffer"  "C1QKupffer"   "PTPRCKupffer"
# This is much better, but I still hate losing some. I will try min 3.
# [1] "DE analysis is done:"
# [1] "included cell types are:"
# [1] "CD1CKupffer"               "C1QKupffer"                "MonocyteDerivedMacrophage" "LST1Kupffer"              
# [5] "CD8EffectorT"              "IFI27Kupffer"             
# [1] "excluded cell types are:"
# [1] "PTPRCKupffer"
# Even better, and I don't think I really want to go below this

# Check DE output information in table with logFC and p-values for each gene-celltype-contrast:
DE_info_group1$celltype_de$de_output_tidy %>% head()
DE_info_group1$hist_pvals

# The p-value distributions are okay but not great, so we will try something else

empirical_pval = FALSE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info_group1$celltype_de$de_output_tidy)
  celltype_de_group1 = DE_info_emp$de_output_tidy_emp %>% dplyr::select(-p_val, -p_adj) %>% dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de_group1 = DE_info_group1$celltype_de$de_output_tidy
}
#DE_info_emp$hist_pvals_emp

# I actually think the first set were better, so I will go back to that

# *Combine DE information for ligand-senders and receptors-receivers*
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de_group1,
  receiver_de = celltype_de_group1,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

sender_receiver_de %>% head(20)

#*Ligand activity prediction: use the DE analysis output to predict the activity of ligands in receiver cell types and infer their potential target genes*
# We will first inspect the geneset_oi-vs-background ratios for the default tresholds:
geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de_group1, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment

ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de_group1,
    receivers_oi = intersect(receivers_oi, celltype_de_group1$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = TRUE, 
    n.cores = n.cores
  )
))

#*Prioritization: rank cell-cell communication patterns through multi-criteria prioritization*
sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(SCE_combined) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down
))

#*Compile the MultiNicheNet output object*
multinichenet_output_group1_t2vst1 = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de_group1,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = tibble()
) 
multinichenet_output_group1_t2vst1 = make_lite_output(multinichenet_output_group1_t2vst1)

#*Summarizing ChordDiagram circos plots*
# We will look here at the top 50 predictions across all contrasts, senders, and receivers of interest.
prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output_group1_t2vst1$prioritization_tables, 
  top_n = 50, 
  rank_per_group = FALSE
)
prioritized_tbl_oi = 
  multinichenet_output_group1_t2vst1$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

dir.create("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/MultiNicheNet/simple_results")
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/MultiNicheNet/simple_results")

# Save results table
write_xlsx(prioritized_tbl_oi, "InExp_IRvsRes_prioritized_tbl_oi.xlsx")

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

######## modified make_circos_groups_comparison
make_circos_group_comparison = function(prioritized_tbl_oi, colors_sender, colors_receiver){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("circlize")
  
  prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup() # if grouped: things will be messed up downstream
  
  # Link each cell type to a color
  grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
  grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
  
  # Make the plot for each condition
  groups_oi = prioritized_tbl_oi$group %>% unique()
  all_plots = groups_oi %>% lapply(function(group_oi){
    
    # Make the plot for condition of interest - title of the plot
    title = group_oi
    circos_links_oi = prioritized_tbl_oi %>% dplyr::filter(group == group_oi)
    
    # deal with duplicated sector names
    # dplyr::rename the ligands so we can have the same ligand in multiple senders (and receptors in multiple receivers)
    # only do it with duplicated ones!
    # this is in iterative process because changing one ligand to try to solve a mistake here, can actually induce another problem
    circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)
    
    df = circos_links
    
    
    ligand.uni = unique(df$ligand)
    for (i in 1:length(ligand.uni)) {
      df.i = df[df$ligand == ligand.uni[i], ]
      sender.uni = unique(df.i$sender)
      for (j in 1:length(sender.uni)) {
        df.i.j = df.i[df.i$sender == sender.uni[j], ]
        df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
        df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
      }
    }
    receptor.uni = unique(df$receptor)
    for (i in 1:length(receptor.uni)) {
      df.i = df[df$receptor == receptor.uni[i], ]
      receiver.uni = unique(df.i$receiver)
      for (j in 1:length(receiver.uni)) {
        df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
        df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
        df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
      }
    }
    
    intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
    
    while(length(intersecting_ligands_receptors) > 0){
      df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
      df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
      df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
      df = dplyr::bind_rows(df_unique, df_duplicated)
      intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
    }
    
    circos_links = df
    
    # Link ligands/Receptors to the colors of senders/receivers
    circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
    links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
    ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
    grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
    receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
    grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
    grid_col =c(grid_ligand_color,grid_receptor_color)
    
    # give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
    transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency
    
    # Define order of the ligands and receptors and the gaps
    ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
      ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
    }) %>% unlist()
    
    receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
      receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
    }) %>% unlist()
    
    order = c(ligand_order,receptor_order)
    
    width_same_cell_same_ligand_type = 0.275
    width_different_cell = 3
    width_ligand_receptor = 9
    width_same_cell_same_receptor_type = 0.275
    
    # sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    #   sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
    #   gap = width_different_cell
    #   return(c(sector,gap))
    # }) %>% unlist()
    # sender_gaps = sender_gaps[-length(sender_gaps)]
    # 
    # receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    #   sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
    #   gap = width_different_cell
    #   return(c(sector,gap))
    # }) %>% unlist()
    # receiver_gaps = receiver_gaps[-length(receiver_gaps)]
    # 
    # gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
    # 
    # if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
    #   warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
    # }
    
    links_circle$weight[links_circle$weight == 0] = 0.01
    circos.clear()
    # circos.par(gap.degree = gaps)
    circos.par(gap.degree=2)
    chordDiagram(links_circle,
                 directional = 1,
                 order=order,
                 link.sort = TRUE,
                 link.decreasing = TRUE,
                 grid.col = grid_col,
                 transparency = transparency,
                 diffHeight = 0.0075,
                 direction.type = c("diffHeight", "arrows"),
                 link.visible = links_circle$weight > 0.01,
                 annotationTrack = "grid",
                 preAllocateTracks = list(track.height = 0.175),
                 grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
                 reduce = 0,
                 scale = TRUE)
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
    }, bg.border = NA) #
    
    title(title)
    p_circos = recordPlot()
    return(p_circos)
    
  })
  names(all_plots) = groups_oi
  
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  # grid_col_all = c(colors_receiver, colors_sender)
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Receiver")
  ComplexHeatmap::draw(legend, just = c("left", "bottom"))
  
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Sender")
  ComplexHeatmap::draw(legend, just = c("left", "top"))
  
  p_legend = grDevices::recordPlot()
  
  all_plots$legend = p_legend
  
  return(all_plots)
}
#######

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
#TODO Error in rep(width_same_cell_same_ligand_type, times = (circos_links %>%  : invalid 'times' argument
# removed the lines relating to variable gaps (as per https://github.com/saeyslab/nichenetr/issues/5)

pdf(file = "InExp_IRvsRes_Expanded_Index_Rejection.pdf", width = 7, height = 7)
par(mar=c(2, 2, 2, 2))
  circos_list[1] 
dev.off()

pdf(file = "InExp_IRvsRes_Expanded_Resolved.pdf", width = 7, height = 7)
par(mar=c(2, 2, 2, 2))
  circos_list[2]
dev.off()

pdf(file = "InExp_IRvsRes_Legend.pdf", width = 7, height = 7)
par(mar=c(2, 2, 2, 2))
  circos_list[3]
dev.off()

# nearly all interactions are G1.T1 and one in G1.T2, meaning the expanded IR vs expanded Resolved, as expected

# Bubble plots
group_oi = "G1.T1"

prioritized_tbl_oi_50 = get_top_n_lr_pairs(
  multinichenet_output_group1_t2vst1$prioritization_tables, 
  top_n = 50, rank_per_group = FALSE) %>% filter(group == group_oi)

prioritized_tbl_oi_50_omnipath = prioritized_tbl_oi_50 %>% 
  inner_join(lr_network_all)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output_group1_t2vst1$prioritization_tables, 
  prioritized_tbl_oi_50_omnipath)

pdf(file = "InExp_IRvsRes_BubblePlot_G1.T1.pdf", width = 18, height = 14)
par(mar=c(2, 2, 2, 2))
  plot_oi
dev.off()

# Bubble plots
group_oi = "G1.T2"

prioritized_tbl_oi_50 = get_top_n_lr_pairs(
  multinichenet_output_group1_t2vst1$prioritization_tables, 
  top_n = 50, rank_per_group = FALSE) %>% filter(group == group_oi)

prioritized_tbl_oi_50_omnipath = prioritized_tbl_oi_50 %>% 
  inner_join(lr_network_all)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output_group1_t2vst1$prioritization_tables, 
  prioritized_tbl_oi_50_omnipath)

pdf(file = "InExp_IRvsRes_BubblePlot_G1.T2.pdf", width = 18, height = 14)
par(mar=c(2, 2, 2, 2))
plot_oi
dev.off()


####################################
# <><><><> Simple test #2 <><><><> #
####################################

## 2) G2.T2-G2.T1: In non-expanded cells, what is the difference between IR and Resolved?

# Set contrasts and table
contrasts_oi = c("'G2.T2-G2.T1','G2.T1-G2.T2'") 

contrast_tbl = tibble(
  contrast = c("G2.T2-G2.T1", "G2.T1-G2.T2"), 
  group = c("G2.T2","G2.T1")
)

#*Cell-type filtering: determine which cell types are sufficiently present*
abundance_info = get_abundance_info(
  sce = SCE_combined, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  batches = batches
)
abundance_info$abund_plot_sample
# There are many kupffer subtype/sample combinations that don't have a sufficent number of cells to do analysis (min = 10)

#*Gene filtering: determine which genes are sufficiently expressed in each present cell type*
frq_list = get_frac_exprs(
  sce = SCE_combined, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)
# Note: make sure the DefaultAssay is set to "RNA" at the beginning of the process

# Now only keep genes that are expressed by at least one cell type:
genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% pull(gene) %>% unique() 
SCE_combined = SCE_combined[genes_oi, ]

# *Pseudobulk expression calculation: determine and normalize per-sample pseudobulk expression levels for each expressed gene in each present cell type*
abundance_expression_info = process_abundance_expression_info(
  sce = SCE_combined, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)
#Warning message:
#In DGEList.default(pb@assays@data[[celltype_oi]]) :
# At least one library size is zero
#TODO: not sure what to make of this warning yet

# *Differential expression (DE) analysis: determine which genes are differentially expressed*
DE_info_group2 = get_DE_info(
  sce = SCE_combined, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)
# [1] "DE analysis is done:"
# [1] "included cell types are:"
# [1] "CD1CKupffer"               "C1QKupffer"                "MonocyteDerivedMacrophage" "LST1Kupffer"              
# [5] "CD8EffectorT"              "IFI27Kupffer"             
# [1] "excluded cell types are:"
# [1] "PTPRCKupffer"

# Check DE output information in table with logFC and p-values for each gene-celltype-contrast:
DE_info_group2$celltype_de$de_output_tidy %>% head()
DE_info_group2$hist_pvals

# This is awful, but only for the non-T cells. This make sense based on the comparison I'm doing, because for other
# cell types the cells are exactly the same, so the p-value for the DE will be 1. Hopefully this doesn't cause
# downstream issues

# empirical_pval = TRUE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info_group2$celltype_de$de_output_tidy)
  celltype_de_group2 = DE_info_emp$de_output_tidy_emp %>% dplyr::select(-p_val, -p_adj) %>% dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de_group2 = DE_info_group2$celltype_de$de_output_tidy
}
#DE_info_emp$hist_pvals_emp

# *Combine DE information for ligand-senders and receptors-receivers*
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de_group2,
  receiver_de = celltype_de_group2,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

sender_receiver_de %>% head(20)

#*Ligand activity prediction: use the DE analysis output to predict the activity of ligands in receiver cell types and infer their potential target genes*
# We will first inspect the geneset_oi-vs-background ratios for the default tresholds:
geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de_group2, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment

# lots of non TRUE/TRUE, but only for non-T types, as expected

ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de_group2,
    receivers_oi = intersect(receivers_oi, celltype_de_group2$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = TRUE, 
    n.cores = n.cores
  )
))

#*Prioritization: rank cell-cell communication patterns through multi-criteria prioritization*
sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(SCE_combined) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down
))

#*Compile the MultiNicheNet output object*
multinichenet_output_group2_t2vst1 = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de_group2,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = tibble()
) 
multinichenet_output_group2_t2vst1 = make_lite_output(multinichenet_output_group2_t2vst1)

#*Summarizing ChordDiagram circos plots*
# We will look here at the top 50 predictions across all contrasts, senders, and receivers of interest.
prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output_group2_t2vst1$prioritization_tables, 
  top_n = 50, 
  rank_per_group = FALSE
)
prioritized_tbl_oi = 
  multinichenet_output_group2_t2vst1$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

# Save results table
write_xlsx(prioritized_tbl_oi, "InNotExp_IRvsRes_prioritized_tbl_oi.xlsx")

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

######## modified make_circos_groups_comparison
make_circos_group_comparison = function(prioritized_tbl_oi, colors_sender, colors_receiver){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("circlize")
  
  prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup() # if grouped: things will be messed up downstream
  
  # Link each cell type to a color
  grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
  grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
  
  # Make the plot for each condition
  groups_oi = prioritized_tbl_oi$group %>% unique()
  all_plots = groups_oi %>% lapply(function(group_oi){
    
    # Make the plot for condition of interest - title of the plot
    title = group_oi
    circos_links_oi = prioritized_tbl_oi %>% dplyr::filter(group == group_oi)
    
    # deal with duplicated sector names
    # dplyr::rename the ligands so we can have the same ligand in multiple senders (and receptors in multiple receivers)
    # only do it with duplicated ones!
    # this is in iterative process because changing one ligand to try to solve a mistake here, can actually induce another problem
    circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)
    
    df = circos_links
    
    
    ligand.uni = unique(df$ligand)
    for (i in 1:length(ligand.uni)) {
      df.i = df[df$ligand == ligand.uni[i], ]
      sender.uni = unique(df.i$sender)
      for (j in 1:length(sender.uni)) {
        df.i.j = df.i[df.i$sender == sender.uni[j], ]
        df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
        df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
      }
    }
    receptor.uni = unique(df$receptor)
    for (i in 1:length(receptor.uni)) {
      df.i = df[df$receptor == receptor.uni[i], ]
      receiver.uni = unique(df.i$receiver)
      for (j in 1:length(receiver.uni)) {
        df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
        df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
        df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
      }
    }
    
    intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
    
    while(length(intersecting_ligands_receptors) > 0){
      df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
      df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
      df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
      df = dplyr::bind_rows(df_unique, df_duplicated)
      intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
    }
    
    circos_links = df
    
    # Link ligands/Receptors to the colors of senders/receivers
    circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
    links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
    ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
    grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
    receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
    grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
    grid_col =c(grid_ligand_color,grid_receptor_color)
    
    # give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
    transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency
    
    # Define order of the ligands and receptors and the gaps
    ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
      ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
    }) %>% unlist()
    
    receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
      receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
    }) %>% unlist()
    
    order = c(ligand_order,receptor_order)
    
    width_same_cell_same_ligand_type = 0.275
    width_different_cell = 3
    width_ligand_receptor = 9
    width_same_cell_same_receptor_type = 0.275
    
    # sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    #   sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
    #   gap = width_different_cell
    #   return(c(sector,gap))
    # }) %>% unlist()
    # sender_gaps = sender_gaps[-length(sender_gaps)]
    # 
    # receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    #   sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
    #   gap = width_different_cell
    #   return(c(sector,gap))
    # }) %>% unlist()
    # receiver_gaps = receiver_gaps[-length(receiver_gaps)]
    # 
    # gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
    # 
    # if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
    #   warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
    # }
    
    links_circle$weight[links_circle$weight == 0] = 0.01
    circos.clear()
    # circos.par(gap.degree = gaps)
    circos.par(gap.degree=2)
    chordDiagram(links_circle,
                 directional = 1,
                 order=order,
                 link.sort = TRUE,
                 link.decreasing = TRUE,
                 grid.col = grid_col,
                 transparency = transparency,
                 diffHeight = 0.0075,
                 direction.type = c("diffHeight", "arrows"),
                 link.visible = links_circle$weight > 0.01,
                 annotationTrack = "grid",
                 preAllocateTracks = list(track.height = 0.175),
                 grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
                 reduce = 0,
                 scale = TRUE)
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
    }, bg.border = NA) #
    
    title(title)
    p_circos = recordPlot()
    return(p_circos)
    
  })
  names(all_plots) = groups_oi
  
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  # grid_col_all = c(colors_receiver, colors_sender)
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Receiver")
  ComplexHeatmap::draw(legend, just = c("left", "bottom"))
  
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Sender")
  ComplexHeatmap::draw(legend, just = c("left", "top"))
  
  p_legend = grDevices::recordPlot()
  
  all_plots$legend = p_legend
  
  return(all_plots)
}
#######

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
#TODO Error in rep(width_same_cell_same_ligand_type, times = (circos_links %>%  : invalid 'times' argument
# removed the lines relating to variable gaps (as per https://github.com/saeyslab/nichenetr/issues/5)

pdf(file = "InNotExp_IRvsRes_Expanded_Index_Rejection.pdf", width = 7, height = 7)
par(mar=c(2, 2, 2, 2))
  circos_list[1] 
dev.off()

pdf(file = "InNotExp_IRvsRes_Expanded_Resolved.pdf", width = 7, height = 7)
par(mar=c(2, 2, 2, 2))
  circos_list[2]
dev.off()

pdf(file = "InNotExp_IRvsRes_Legend.pdf", width = 7, height = 7)
par(mar=c(2, 2, 2, 2))
  circos_list[3]
dev.off()

# nearly all interactions are G2.T2 with some in G1.T2, meaning more interactions in non-expanded T than T, which might be due to cell numbers

# Bubble plots #1
group_oi = "G2.T1"

prioritized_tbl_oi_50 = get_top_n_lr_pairs(
  multinichenet_output_group2_t2vst1$prioritization_tables, 
  top_n = 50, rank_per_group = FALSE) %>% filter(group == group_oi)

prioritized_tbl_oi_50_omnipath = prioritized_tbl_oi_50 %>% 
  inner_join(lr_network_all)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output_group2_t2vst1$prioritization_tables, 
  prioritized_tbl_oi_50_omnipath)

pdf(file = "InNotExp_IRvsRes_BubblePlot_G2.T1.pdf", width = 18, height = 14)
par(mar=c(2, 2, 2, 2))
  plot_oi
dev.off()

# Bubble plots #1
group_oi = "G2.T2"

prioritized_tbl_oi_50 = get_top_n_lr_pairs(
  multinichenet_output_group2_t2vst1$prioritization_tables, 
  top_n = 50, rank_per_group = FALSE) %>% filter(group == group_oi)

prioritized_tbl_oi_50_omnipath = prioritized_tbl_oi_50 %>% 
  inner_join(lr_network_all)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output_group2_t2vst1$prioritization_tables, 
  prioritized_tbl_oi_50_omnipath)

pdf(file = "InNotExp_IRvsRes_BubblePlot_G2.T2.pdf", width = 18, height = 14)
par(mar=c(2, 2, 2, 2))
plot_oi
dev.off()



####################################
# <><><><> Simple test #3 <><><><> #
####################################

## 1) G1.T2-G2.T2: 

# Set contrasts and table
contrasts_oi = c("'G1.T2-G2.T2','G2.T2-G1.T2'") 

contrast_tbl = tibble(
  contrast = c("G1.T2-G2.T2", "G2.T2-G1.T2"), 
  group = c("G1.T2","G2.T2")
)

#*Cell-type filtering: determine which cell types are sufficiently present*
abundance_info = get_abundance_info(
  sce = SCE_combined, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  batches = batches
)
abundance_info$abund_plot_sample
# There are many kupffer subtype/sample combinations that don't have a sufficent number of cells to do analysis (min = 10)

#*Gene filtering: determine which genes are sufficiently expressed in each present cell type*
frq_list = get_frac_exprs(
  sce = SCE_combined, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)
# Note: make sure the DefaultAssay is set to "RNA" at the beginning of the process

# Now only keep genes that are expressed by at least one cell type:
genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% pull(gene) %>% unique() 
SCE_combined = SCE_combined[genes_oi, ]

# *Pseudobulk expression calculation: determine and normalize per-sample pseudobulk expression levels for each expressed gene in each present cell type*
abundance_expression_info = process_abundance_expression_info(
  sce = SCE_combined, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)
#Warning message:
#In DGEList.default(pb@assays@data[[celltype_oi]]) :
# At least one library size is zero

# *Differential expression (DE) analysis: determine which genes are differentially expressed*
DE_info_group1 = get_DE_info(
  sce = SCE_combined, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)
#[1] "DE analysis is done:"
# [1] "included cell types are:"
# [1] "LST1Kupffer"
# [1] "excluded cell types are:"
# [1] "CD1CKupffer"               "C1QKupffer"                "MonocyteDerivedMacrophage" "CD8EffectorT"             
# [5] "PTPRCKupffer"              "IFI27Kupffer"             
# TODO: this is obviously not okay!
# I'm going to try cutting min_cells in half (10 -> 5)
# [1] "DE analysis is done:"
# [1] "included cell types are:"
# [1] "MonocyteDerivedMacrophage" "LST1Kupffer"               "CD8EffectorT"              "IFI27Kupffer"             
# [1] "excluded cell types are:"
# [1] "CD1CKupffer"  "C1QKupffer"   "PTPRCKupffer"
# This is much better, but I still hate losing some. I will try min 3.
# [1] "DE analysis is done:"
# [1] "included cell types are:"
# [1] "CD1CKupffer"               "C1QKupffer"                "MonocyteDerivedMacrophage" "LST1Kupffer"              
# [5] "CD8EffectorT"              "IFI27Kupffer"             
# [1] "excluded cell types are:"
# [1] "PTPRCKupffer"
# Even better, and I don't think I really want to go below this

# Check DE output information in table with logFC and p-values for each gene-celltype-contrast:
DE_info_group1$celltype_de$de_output_tidy %>% head()
DE_info_group1$hist_pvals

# The p-value distributions are okay but not great, so we will try something else

empirical_pval = FALSE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info_group1$celltype_de$de_output_tidy)
  celltype_de_group1 = DE_info_emp$de_output_tidy_emp %>% dplyr::select(-p_val, -p_adj) %>% dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de_group1 = DE_info_group1$celltype_de$de_output_tidy
}
#DE_info_emp$hist_pvals_emp

# I actually think the first set were better, so I will go back to that

# *Combine DE information for ligand-senders and receptors-receivers*
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de_group1,
  receiver_de = celltype_de_group1,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

sender_receiver_de %>% head(20)

#*Ligand activity prediction: use the DE analysis output to predict the activity of ligands in receiver cell types and infer their potential target genes*
# We will first inspect the geneset_oi-vs-background ratios for the default tresholds:
geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de_group1, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment

ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de_group1,
    receivers_oi = intersect(receivers_oi, celltype_de_group1$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = TRUE, 
    n.cores = n.cores
  )
))

#*Prioritization: rank cell-cell communication patterns through multi-criteria prioritization*
sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(SCE_combined) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down
))

#*Compile the MultiNicheNet output object*
multinichenet_output_group1_t2vst1 = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de_group1,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = tibble()
) 
multinichenet_output_group1_t2vst1 = make_lite_output(multinichenet_output_group1_t2vst1)

#*Summarizing ChordDiagram circos plots*
# We will look here at the top 50 predictions across all contrasts, senders, and receivers of interest.
prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output_group1_t2vst1$prioritization_tables, 
  top_n = 50, 
  rank_per_group = FALSE
)
prioritized_tbl_oi = 
  multinichenet_output_group1_t2vst1$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

dir.create("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/MultiNicheNet/simple_results")
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/MultiNicheNet/simple_results")

# Save results table
write_xlsx(prioritized_tbl_oi, "InRes_ExpvsNotExp_prioritized_tbl_oi.xlsx")

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

######## modified make_circos_groups_comparison
make_circos_group_comparison = function(prioritized_tbl_oi, colors_sender, colors_receiver){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("circlize")
  
  prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup() # if grouped: things will be messed up downstream
  
  # Link each cell type to a color
  grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
  grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
  
  # Make the plot for each condition
  groups_oi = prioritized_tbl_oi$group %>% unique()
  all_plots = groups_oi %>% lapply(function(group_oi){
    
    # Make the plot for condition of interest - title of the plot
    title = group_oi
    circos_links_oi = prioritized_tbl_oi %>% dplyr::filter(group == group_oi)
    
    # deal with duplicated sector names
    # dplyr::rename the ligands so we can have the same ligand in multiple senders (and receptors in multiple receivers)
    # only do it with duplicated ones!
    # this is in iterative process because changing one ligand to try to solve a mistake here, can actually induce another problem
    circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)
    
    df = circos_links
    
    
    ligand.uni = unique(df$ligand)
    for (i in 1:length(ligand.uni)) {
      df.i = df[df$ligand == ligand.uni[i], ]
      sender.uni = unique(df.i$sender)
      for (j in 1:length(sender.uni)) {
        df.i.j = df.i[df.i$sender == sender.uni[j], ]
        df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
        df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
      }
    }
    receptor.uni = unique(df$receptor)
    for (i in 1:length(receptor.uni)) {
      df.i = df[df$receptor == receptor.uni[i], ]
      receiver.uni = unique(df.i$receiver)
      for (j in 1:length(receiver.uni)) {
        df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
        df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
        df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
      }
    }
    
    intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
    
    while(length(intersecting_ligands_receptors) > 0){
      df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
      df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
      df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
      df = dplyr::bind_rows(df_unique, df_duplicated)
      intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
    }
    
    circos_links = df
    
    # Link ligands/Receptors to the colors of senders/receivers
    circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
    links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
    ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
    grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
    receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
    grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
    grid_col =c(grid_ligand_color,grid_receptor_color)
    
    # give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
    transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency
    
    # Define order of the ligands and receptors and the gaps
    ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
      ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
    }) %>% unlist()
    
    receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
      receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
    }) %>% unlist()
    
    order = c(ligand_order,receptor_order)
    
    width_same_cell_same_ligand_type = 0.275
    width_different_cell = 3
    width_ligand_receptor = 9
    width_same_cell_same_receptor_type = 0.275
    
    # sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    #   sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
    #   gap = width_different_cell
    #   return(c(sector,gap))
    # }) %>% unlist()
    # sender_gaps = sender_gaps[-length(sender_gaps)]
    # 
    # receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    #   sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
    #   gap = width_different_cell
    #   return(c(sector,gap))
    # }) %>% unlist()
    # receiver_gaps = receiver_gaps[-length(receiver_gaps)]
    # 
    # gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
    # 
    # if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
    #   warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
    # }
    
    links_circle$weight[links_circle$weight == 0] = 0.01
    circos.clear()
    # circos.par(gap.degree = gaps)
    circos.par(gap.degree=2)
    chordDiagram(links_circle,
                 directional = 1,
                 order=order,
                 link.sort = TRUE,
                 link.decreasing = TRUE,
                 grid.col = grid_col,
                 transparency = transparency,
                 diffHeight = 0.0075,
                 direction.type = c("diffHeight", "arrows"),
                 link.visible = links_circle$weight > 0.01,
                 annotationTrack = "grid",
                 preAllocateTracks = list(track.height = 0.175),
                 grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
                 reduce = 0,
                 scale = TRUE)
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
    }, bg.border = NA) #
    
    title(title)
    p_circos = recordPlot()
    return(p_circos)
    
  })
  names(all_plots) = groups_oi
  
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  # grid_col_all = c(colors_receiver, colors_sender)
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Receiver")
  ComplexHeatmap::draw(legend, just = c("left", "bottom"))
  
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Sender")
  ComplexHeatmap::draw(legend, just = c("left", "top"))
  
  p_legend = grDevices::recordPlot()
  
  all_plots$legend = p_legend
  
  return(all_plots)
}
#######

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
#TODO Error in rep(width_same_cell_same_ligand_type, times = (circos_links %>%  : invalid 'times' argument
# removed the lines relating to variable gaps (as per https://github.com/saeyslab/nichenetr/issues/5)

pdf(file = "InRes_ExpvsNotExp_NotExpanded_Resolved.pdf", width = 7, height = 7)
par(mar=c(2, 2, 2, 2))
circos_list[1] 
dev.off()

pdf(file = "InRes_ExpvsNotExp_Expanded_Resolved.pdf", width = 7, height = 7)
par(mar=c(2, 2, 2, 2))
circos_list[2]
dev.off()

pdf(file = "InRes_ExpvsNotExp_Legend.pdf", width = 7, height = 7)
par(mar=c(2, 2, 2, 2))
circos_list[3]
dev.off()

# nearly all interactions are G1.T1 and one in G1.T2, meaning the expanded IR vs expanded Resolved, as expected

# Bubble plots
group_oi = "G1.T2"

prioritized_tbl_oi_50 = get_top_n_lr_pairs(
  multinichenet_output_group1_t2vst1$prioritization_tables, 
  top_n = 50, rank_per_group = FALSE) %>% filter(group == group_oi)

prioritized_tbl_oi_50_omnipath = prioritized_tbl_oi_50 %>% 
  inner_join(lr_network_all)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output_group1_t2vst1$prioritization_tables, 
  prioritized_tbl_oi_50_omnipath)

pdf(file = "InRes_ExpvsNotExp_BubblePlot_G1.T2.pdf", width = 18, height = 14)
par(mar=c(2, 2, 2, 2))
  plot_oi
dev.off()

# Bubble plots
group_oi = "G2.T2"

prioritized_tbl_oi_50 = get_top_n_lr_pairs(
  multinichenet_output_group1_t2vst1$prioritization_tables, 
  top_n = 50, rank_per_group = FALSE) %>% filter(group == group_oi)

prioritized_tbl_oi_50_omnipath = prioritized_tbl_oi_50 %>% 
  inner_join(lr_network_all)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output_group1_t2vst1$prioritization_tables, 
  prioritized_tbl_oi_50_omnipath)

pdf(file = "InRes_ExpvsNotExp_BubblePlot_G2.T2.pdf", width = 18, height = 14)
par(mar=c(2, 2, 2, 2))
plot_oi
dev.off()




####################################
# <><><><> Simple test #4 <><><><> #
####################################

## 1) G1.T1-G2.T1: 

# Set contrasts and table
contrasts_oi = c("'G1.T1-G2.T1','G2.T1-G1.T1'") 

contrast_tbl = tibble(
  contrast = c("G1.T1-G2.T1", "G2.T1-G1.T1"), 
  group = c("G1.T1","G2.T1")
)

#*Cell-type filtering: determine which cell types are sufficiently present*
abundance_info = get_abundance_info(
  sce = SCE_combined, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  batches = batches
)
abundance_info$abund_plot_sample
# There are many kupffer subtype/sample combinations that don't have a sufficent number of cells to do analysis (min = 10)

#*Gene filtering: determine which genes are sufficiently expressed in each present cell type*
frq_list = get_frac_exprs(
  sce = SCE_combined, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)
# Note: make sure the DefaultAssay is set to "RNA" at the beginning of the process

# Now only keep genes that are expressed by at least one cell type:
genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% pull(gene) %>% unique() 
SCE_combined = SCE_combined[genes_oi, ]

# *Pseudobulk expression calculation: determine and normalize per-sample pseudobulk expression levels for each expressed gene in each present cell type*
abundance_expression_info = process_abundance_expression_info(
  sce = SCE_combined, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)
#Warning message:
#In DGEList.default(pb@assays@data[[celltype_oi]]) :
# At least one library size is zero

# *Differential expression (DE) analysis: determine which genes are differentially expressed*
DE_info_group1 = get_DE_info(
  sce = SCE_combined, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)
#[1] "DE analysis is done:"
# [1] "included cell types are:"
# [1] "LST1Kupffer"
# [1] "excluded cell types are:"
# [1] "CD1CKupffer"               "C1QKupffer"                "MonocyteDerivedMacrophage" "CD8EffectorT"             
# [5] "PTPRCKupffer"              "IFI27Kupffer"             
# TODO: this is obviously not okay!
# I'm going to try cutting min_cells in half (10 -> 5)
# [1] "DE analysis is done:"
# [1] "included cell types are:"
# [1] "MonocyteDerivedMacrophage" "LST1Kupffer"               "CD8EffectorT"              "IFI27Kupffer"             
# [1] "excluded cell types are:"
# [1] "CD1CKupffer"  "C1QKupffer"   "PTPRCKupffer"
# This is much better, but I still hate losing some. I will try min 3.
# [1] "DE analysis is done:"
# [1] "included cell types are:"
# [1] "CD1CKupffer"               "C1QKupffer"                "MonocyteDerivedMacrophage" "LST1Kupffer"              
# [5] "CD8EffectorT"              "IFI27Kupffer"             
# [1] "excluded cell types are:"
# [1] "PTPRCKupffer"
# Even better, and I don't think I really want to go below this

# Check DE output information in table with logFC and p-values for each gene-celltype-contrast:
DE_info_group1$celltype_de$de_output_tidy %>% head()
DE_info_group1$hist_pvals

# The p-value distributions are okay but not great, so we will try something else

empirical_pval = FALSE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info_group1$celltype_de$de_output_tidy)
  celltype_de_group1 = DE_info_emp$de_output_tidy_emp %>% dplyr::select(-p_val, -p_adj) %>% dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de_group1 = DE_info_group1$celltype_de$de_output_tidy
}
#DE_info_emp$hist_pvals_emp

# I actually think the first set were better, so I will go back to that

# *Combine DE information for ligand-senders and receptors-receivers*
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de_group1,
  receiver_de = celltype_de_group1,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

sender_receiver_de %>% head(20)

#*Ligand activity prediction: use the DE analysis output to predict the activity of ligands in receiver cell types and infer their potential target genes*
# We will first inspect the geneset_oi-vs-background ratios for the default tresholds:
geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de_group1, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment

ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de_group1,
    receivers_oi = intersect(receivers_oi, celltype_de_group1$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = TRUE, 
    n.cores = n.cores
  )
))

#*Prioritization: rank cell-cell communication patterns through multi-criteria prioritization*
sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(SCE_combined) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down
))

#*Compile the MultiNicheNet output object*
multinichenet_output_group1_t2vst1 = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de_group1,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = tibble()
) 
multinichenet_output_group1_t2vst1 = make_lite_output(multinichenet_output_group1_t2vst1)

#*Summarizing ChordDiagram circos plots*
# We will look here at the top 50 predictions across all contrasts, senders, and receivers of interest.
prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output_group1_t2vst1$prioritization_tables, 
  top_n = 50, 
  rank_per_group = FALSE
)
prioritized_tbl_oi = 
  multinichenet_output_group1_t2vst1$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

dir.create("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/MultiNicheNet/simple_results")
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/MultiNicheNet/simple_results")

# Save results table
write_xlsx(prioritized_tbl_oi, "InIR_ExpvsNotExp_prioritized_tbl_oi.xlsx")

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

######## modified make_circos_groups_comparison
make_circos_group_comparison = function(prioritized_tbl_oi, colors_sender, colors_receiver){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("circlize")
  
  prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup() # if grouped: things will be messed up downstream
  
  # Link each cell type to a color
  grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
  grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
  
  # Make the plot for each condition
  groups_oi = prioritized_tbl_oi$group %>% unique()
  all_plots = groups_oi %>% lapply(function(group_oi){
    
    # Make the plot for condition of interest - title of the plot
    title = group_oi
    circos_links_oi = prioritized_tbl_oi %>% dplyr::filter(group == group_oi)
    
    # deal with duplicated sector names
    # dplyr::rename the ligands so we can have the same ligand in multiple senders (and receptors in multiple receivers)
    # only do it with duplicated ones!
    # this is in iterative process because changing one ligand to try to solve a mistake here, can actually induce another problem
    circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)
    
    df = circos_links
    
    
    ligand.uni = unique(df$ligand)
    for (i in 1:length(ligand.uni)) {
      df.i = df[df$ligand == ligand.uni[i], ]
      sender.uni = unique(df.i$sender)
      for (j in 1:length(sender.uni)) {
        df.i.j = df.i[df.i$sender == sender.uni[j], ]
        df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
        df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
      }
    }
    receptor.uni = unique(df$receptor)
    for (i in 1:length(receptor.uni)) {
      df.i = df[df$receptor == receptor.uni[i], ]
      receiver.uni = unique(df.i$receiver)
      for (j in 1:length(receiver.uni)) {
        df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
        df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
        df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
      }
    }
    
    intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
    
    while(length(intersecting_ligands_receptors) > 0){
      df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
      df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
      df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
      df = dplyr::bind_rows(df_unique, df_duplicated)
      intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
    }
    
    circos_links = df
    
    # Link ligands/Receptors to the colors of senders/receivers
    circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
    links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
    ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
    grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
    receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
    grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
    grid_col =c(grid_ligand_color,grid_receptor_color)
    
    # give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
    transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency
    
    # Define order of the ligands and receptors and the gaps
    ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
      ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
    }) %>% unlist()
    
    receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
      receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
    }) %>% unlist()
    
    order = c(ligand_order,receptor_order)
    
    width_same_cell_same_ligand_type = 0.275
    width_different_cell = 3
    width_ligand_receptor = 9
    width_same_cell_same_receptor_type = 0.275
    
    # sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    #   sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
    #   gap = width_different_cell
    #   return(c(sector,gap))
    # }) %>% unlist()
    # sender_gaps = sender_gaps[-length(sender_gaps)]
    # 
    # receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    #   sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
    #   gap = width_different_cell
    #   return(c(sector,gap))
    # }) %>% unlist()
    # receiver_gaps = receiver_gaps[-length(receiver_gaps)]
    # 
    # gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
    # 
    # if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
    #   warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
    # }
    
    links_circle$weight[links_circle$weight == 0] = 0.01
    circos.clear()
    # circos.par(gap.degree = gaps)
    circos.par(gap.degree=2)
    chordDiagram(links_circle,
                 directional = 1,
                 order=order,
                 link.sort = TRUE,
                 link.decreasing = TRUE,
                 grid.col = grid_col,
                 transparency = transparency,
                 diffHeight = 0.0075,
                 direction.type = c("diffHeight", "arrows"),
                 link.visible = links_circle$weight > 0.01,
                 annotationTrack = "grid",
                 preAllocateTracks = list(track.height = 0.175),
                 grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
                 reduce = 0,
                 scale = TRUE)
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
    }, bg.border = NA) #
    
    title(title)
    p_circos = recordPlot()
    return(p_circos)
    
  })
  names(all_plots) = groups_oi
  
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  # grid_col_all = c(colors_receiver, colors_sender)
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Receiver")
  ComplexHeatmap::draw(legend, just = c("left", "bottom"))
  
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Sender")
  ComplexHeatmap::draw(legend, just = c("left", "top"))
  
  p_legend = grDevices::recordPlot()
  
  all_plots$legend = p_legend
  
  return(all_plots)
}
#######

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
#TODO Error in rep(width_same_cell_same_ligand_type, times = (circos_links %>%  : invalid 'times' argument
# removed the lines relating to variable gaps (as per https://github.com/saeyslab/nichenetr/issues/5)

pdf(file = "InIR_ExpvsNotExp_NotExpanded_IR.pdf", width = 7, height = 7)
par(mar=c(2, 2, 2, 2))
circos_list[1] 
dev.off()

pdf(file = "InIR_ExpvsNotExp_Expanded_IR.pdf", width = 7, height = 7)
par(mar=c(2, 2, 2, 2))
circos_list[2]
dev.off()

pdf(file = "InIR_ExpvsNotExp_Legend.pdf", width = 7, height = 7)
par(mar=c(2, 2, 2, 2))
circos_list[3]
dev.off()

# nearly all interactions are G1.T1 and one in G1.T2, meaning the expanded IR vs expanded Resolved, as expected

# Bubble plots
group_oi = "G1.T1"

prioritized_tbl_oi_50 = get_top_n_lr_pairs(
  multinichenet_output_group1_t2vst1$prioritization_tables, 
  top_n = 50, rank_per_group = FALSE) %>% filter(group == group_oi)

prioritized_tbl_oi_50_omnipath = prioritized_tbl_oi_50 %>% 
  inner_join(lr_network_all)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output_group1_t2vst1$prioritization_tables, 
  prioritized_tbl_oi_50_omnipath)

pdf(file = "InIR_ExpvsNotExp_BubblePlot_G1.T1.pdf", width = 18, height = 13)
par(mar=c(2, 2, 2, 2))
plot_oi
dev.off()

# Bubble plots
group_oi = "G2.T1"

prioritized_tbl_oi_50 = get_top_n_lr_pairs(
  multinichenet_output_group1_t2vst1$prioritization_tables, 
  top_n = 50, rank_per_group = FALSE) %>% filter(group == group_oi)

prioritized_tbl_oi_50_omnipath = prioritized_tbl_oi_50 %>% 
  inner_join(lr_network_all)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output_group1_t2vst1$prioritization_tables, 
  prioritized_tbl_oi_50_omnipath)

pdf(file = "InIR_ExpvsNotExp_BubblePlot_G2.T1.pdf", width = 18, height = 13)
par(mar=c(2, 2, 2, 2))
plot_oi
dev.off()
