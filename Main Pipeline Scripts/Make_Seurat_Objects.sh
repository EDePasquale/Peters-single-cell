#BSUB -W 48:00
#BSUB -o Make_Seurat_Objects.out
#BSUB -J Make_Seurat_Objects
#BSUB -M 300000

module load R/4.1.1

Rscript Make_Seurat_Objects.R