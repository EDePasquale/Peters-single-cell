# Defining the T cell transcriptional landscape in pediatric liver transplant rejection at single cell resolution
![GraphicalAbstract](https://github.com/user-attachments/assets/a974bb77-ac7b-426e-a38e-1b4a73236f42)

## Paper Citation
BioRxiv: https://www.biorxiv.org/content/10.1101/2024.02.26.582173v1.full.pdf

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
We next investigated the T cells in the dataset and used the following scripts to generate the T cell related figures. We also did the same with non-T cells:  

10. T and Non-T cell figure scripts  
      a. T cell figure scripts (T_cell_plots.R, Shared_clones_bar.R, TRA_TRB_Stats.R)  
      b. Non-T cell figure scripts (B_cell_plots.R, NK_cell_plots.R, Kupffer_cell_plots.R)

### Splicing Scripts
This set of scripts was used to explore the CD45RA/RO usage in the CD4+ clusters to separate Naive T from TCM:

11. Splicing Scripts
      a. Generate cluster-barcode associations for splitting by cell type (Generate_Clusters_Splicing.R) 
      b. Merge the BAM files from all samples (BAM_merge.sh) 
      c. Split the merged BAM file by cell type (BAMtoSubBAM.sh and extract-reads-pysam.py) 

### Stats For Cell Type Proportions
The next set of scripts were used to perform all of the differential cell type proportion statistics associated with each stacked bar plot for the various cell types:  

12. Stats for cell type proportions  
      a. Figure 1 (Stats_for_prop_all.R)  
      b. T cells (Stats_for_prop_Tcells.R)  
      c. B cells (Stats_for_prop_Bcells.R)  
      d. Kupffer cells (Stats_for_prop_Kcells.R)  
      e. NK cells (Stats_for_prop_NKcells.R)

### MultiNicheNet
This script was used to run MultiNicheNet to investigate the interactions between T cells and Kupffer/Macrophage cells:

13. NicheNet  
      a. Run MultiNicheNet (MultiNicheNet_4simple.R)  

### SoupOrCell
This set of script runs SoupOrCell to identify donor-derived and recipient-derived cells in multiple biopsy samples:

14. SoupOrCell  
      a. Combine BAM files from multiple biopsies, including the matching donor sample (mergebams_Patient4.sh, mergebams_Patient5.sh, mergebams_Patient6.sh)
      b. Run SoupOrCell and generate figures (SouporCell_12Jun2024_Patient4.sh, SouporCell_12Jun2024_Patient5.sh, SouporCell_12Jun2024_Patient6.sh)
      
#### Note: All additional figure panels were generated outside R (e.g., in Prism) and the code is not included in this repository.
