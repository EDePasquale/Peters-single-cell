#BSUB -W 48:00
#BSUB -o QC_and_Filtering.out
#BSUB -J QC_and_Filtering
#BSUB -M 300000

module load R/4.1.1

Rscript QC_and_Filtering.R