#BSUB -W 24:00
#BSUB -o /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/SoupOrCell_7Jun2024/Patient6/results/SouporCell_merged.out
#BSUB -J SouporCell_merged_Patient6
#BSUB -M 100000
#BSUB -n 4
#BSUB -R "span[hosts=1]"

module load anaconda3
module load souporcell/1.0.0
source activate souporcell
module load bedtools/2.30.0

/usr/local/souporcell/1.0.0/souporcell_pipeline.py \
    -i /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/SoupOrCell_7Jun2024/Patient6/out_sorted_bam.bam \
    -b /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/SoupOrCell_7Jun2024/Patient6/out_barcodes.tsv \
    -f /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/SouporCell_test/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa \
    -t 4 \
    -o /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/SoupOrCell_7Jun2024/Patient6/results \
    -k 2 \
    --common_variants /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/SoupOrCell_7Jun2024/filtered_2p_1kgenomes_unchr.vcf

#In data directory
#cp barcodes.tsv.gz barcodes2.tsv.gz 
#gunzip barcodes2.tsv.gz 
#mv barcodes2.tsv barcodes.tsv

#In home directory
#bsub < /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Scripts_Used/Peters-scLiver/SouporCell_12Jun2024_Patient6.sh