#BSUB -W 24:00
#BSUB -o /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Scripts_Used/Peters-scLiver/mergebams_patient5.out
#BSUB -J mergebams_patient5
#BSUB -M 100000
#BSUB -n 4
#BSUB -R "span[hosts=1]"

~/mergebams/target/release/mergebams \
	-i /data/LiverTransplantSingleCell/data/ALP00045/10X-Peters-ALP-00045-TXP-20220505-5v1-1-hg/ALP-00045-TXP/outs/possorted_genome_bam.bam,/data/LiverTransplantSingleCell/data/ALP00045/10X-Peters-ALP-00045-BX1-20220505-5v1-1-hg/ALP-00045-BX1/outs/possorted_genome_bam.bam \
	-l PreTXP_,BX1_ \
	-b /data/LiverTransplantSingleCell/data/ALP00045/10X-Peters-ALP-00045-TXP-20220505-5v1-1-hg/ALP-00045-TXP/outs/filtered_feature_bc_matrix/barcodes.tsv.gz,/data/LiverTransplantSingleCell/data/ALP00045/10X-Peters-ALP-00045-BX1-20220505-5v1-1-hg/ALP-00045-BX1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
	-o /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/SoupOrCell_7Jun2024/Patient5/ \
	-t 4