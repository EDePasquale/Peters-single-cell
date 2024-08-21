#BSUB -W 24:00
#BSUB -o /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Scripts_Used/Peters-scLiver/mergebams_patient4.out
#BSUB -J mergebams_patient4
#BSUB -M 10000
#BSUB -n 4
#BSUB -R "span[hosts=1]"

~/mergebams/target/release/mergebams \
	-i /data/LiverTransplantSingleCell/data/ALP00042/10X-Peters-ALP-00042-TXP-20220506-5v1-1-hg/ALP-00042-TXP/outs/possorted_genome_bam.bam,/data/LiverTransplantSingleCell/data/ALP00042/10X-Peters-ALP-00042-BX2-20220506-5v1-1-hg/ALP-00042-BX2/outs/possorted_genome_bam.bam \
	-l PreTXP_,BX2_ \
	-b /data/LiverTransplantSingleCell/data/ALP00042/10X-Peters-ALP-00042-TXP-20220506-5v1-1-hg/ALP-00042-TXP/outs/filtered_feature_bc_matrix/barcodes.tsv.gz,/data/LiverTransplantSingleCell/data/ALP00042/10X-Peters-ALP-00042-BX2-20220506-5v1-1-hg/ALP-00042-BX2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
	-o /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/SoupOrCell_7Jun2024/Patient4/ \
	-t 4