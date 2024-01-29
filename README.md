# Peters lab publication scripts

## Paper Citation
[to be determined]

## Project Description
In this project, the Peters lab performed 5' single-cell RNA-sequencing paired with TCR and BCR detection using the 10X Genomics Chromium Single Cell Immune Profiling kit on liver samples collected from pediatric patients with liver rejection following transplantation. For each patient there is at least 1 sample taken from a biopsy time point following transplantation, but there are frequently multiple time points gathered as well as sequecing of the healthy donor prior to transplantation. The goal of this project is to identify T-cell clones that may be driving rejection and to explore their gene expression differences across responses to treatment and interactions with other immune cells. Secondary goals include exploration of the BCR data and identification of novel cell states during rejection.

## Processing Steps
### Main Pipeline Scripts
The following scripts were used in order to perform the following tasks:  
1. CellRanger for GEX and TCR/BCR data (cellranger.sh, cellranger_TCR.sh, cellranger_BCR.sh)
2. QC and filtering (QC_and_Filtering.sh and .R)
3. Make Seurat objects (Make_Seurat_Objects.sh and .R)
4. SCT and integration (SCT_and_Integration.sh and .R)
5. Subcluster (Subcluster.R)
6. Add cluster names (Add_Cluster_Names.R)
7. Add meta data (Add_Meta_Data.R)
8. Prop plots (Prop_Plots.R)
9. Add VDJ (Add_VDJ.R)

### T cell Figure Scripts and Non T Cell Figures Scripts
We next investigated the T cells in the dataset and used the following scripts to generate the T-cell related figures:  

10. T and Non-T cell figure scripts  
      a. (T_cell_plots.R, Shared_clones_bar.R, TRA_TRB_Stats.R, Shared_clones_bar.R)  
      b. (B_cell_plots.R, NK_cell_plots.R, Kupffer_cell_plots.R)

### Stats For Cell Type Proportions
The next set of scripts were used to perform all of the differential cell type proportion statistics associated with each stacked bar plot for the various cell types:  

11. Stats for cell type proportions  
      a. Figure 1 (Stats_for_prop_all.R)  
      b. T-cells (Stats_for_prop_Tcells.R)  
      c. B-cells (Stats_for_prop_Bcells.R)  
      d. Kupffer cells (Stats_for_prop_Kcells.R)  
      e. NK cells (Stats_for_prop_NKcells.R)

### CellPhoneDB
Finally, we used CellPhoneDB to investigate the interactions between T cells and Kupffer cells:

12. CellPhoneDB  
      a. Prep for running CellPhoneDB, i.e., getting all of the files into the correct format and performing the differential expression (Prep_for_CellPhoneDB.R)  
      b. Run CellPhoneDB for each Kupffer cell or Monocyte-Derived Macrophage population and all CD8+ T-cell populations. This was done for all CD8+ T-cells, only those expanded, and only those not expanded (e.g., Peters_DEG_cellphonedb_PTPRC_NotExp.sh)  
      c. Post CellPhoneDB script to generate figures (Post_CellPhoneDB.R)  
      d. Script to reduce the relevant interactions tables to only those interactions that are different between expanded and not expanded (Differential_CCC.R)

#### Note: All additional figure panels were generated outside R (e.g., in Prism) and the code is not included in this repository.
