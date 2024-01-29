#BSUB -W 24:00
#BSUB -o PT1-BX1-TCR.out
#BSUB -J PT1-BX1-TCR
#BSUB -M 4000

/data/hildemanlab/programs/cellranger-6.0.0/cellranger vdj --id=PT1-BX1-TCR --chain TR --fastqs=fastqs --reference=/data/hildemanlab/programs/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 --jobmode=lsf --maxjobs=100
