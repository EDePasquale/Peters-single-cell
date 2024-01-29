#BSUB -W 24:00
#BSUB -o PT1-BX1.out
#BSUB -J PT1-BX1
#BSUB -M 4000

module load cellranger/6.0.0

cellranger count --id=PT1-BX1 --fastqs=fastqs --transcriptome=/data/GI-Informatics/DePasquale/References/refdata-gex-GRCh38-2020-A --jobmode=lsf --maxjobs=100
