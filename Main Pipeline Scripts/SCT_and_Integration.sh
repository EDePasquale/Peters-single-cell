#BSUB -W 48:00
#BSUB -o SCT_and_Integration.out
#BSUB -J SCT_and_Integration
#BSUB -M 300000

module load R/4.1.1

Rscript SCT_and_Integration.R