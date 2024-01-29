#BSUB -W 24:00
#BSUB -o /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Scripts_Used/Peters-single-cell/CellPhoneDB/Peters_CD1C_All.out
#BSUB -J Peters_CD1C_All
#BSUB -M 100000

module load python3
source ~/cpdb/bin/activate
module load R/4.1.1

cellphonedb method degs_analysis  \
    /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellPhoneDB/CD1C+\ Kupffer/All_samples/Peters_5PrimeTCRBCR_meta.tsv  \
    /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellPhoneDB/CD1C+\ Kupffer/All_samples/Peters_5PrimeTCRBCR_counts_mtx  \
    /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellPhoneDB/CD1C+\ Kupffer/All_samples/Peters_5PrimeTCRBCR_DEGs.tsv  \
    --output-path /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/CellPhoneDB/CD1C+\ Kupffer/All_samples/Peters_DEG_results\
    --verbose \
    --counts-data hgnc_symbol

#bsub < /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Scripts_Used/Peters-single-cell/CellPhoneDB/Peters_DEG_cellphonedb_CD1C_All.sh