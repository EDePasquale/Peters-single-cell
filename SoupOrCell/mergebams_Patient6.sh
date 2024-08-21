#BSUB -W 24:00
#BSUB -o mergebams_patient6.out
#BSUB -J mergebams_patient6
#BSUB -M 10000
#BSUB -n 4
#BSUB -R "span[hosts=1]"

~/mergebams/target/release/mergebams \
	-i /data/LiverTransplantSingleCell/data/ALP00036/10X-Peters-ALP-00036-20210728-5v1-1hg/ALP-00036/outs/possorted_genome_bam.bam, \
		/data/LiverTransplantSingleCell/data/ALP00036/10X-Peters-ALP-00036-BX2-20220224-5v1-1-hg/ALP-00036-BX2/outs/possorted_genome_bam.bam, \
		/data/LiverTransplantSingleCell/data/ALP00036/10X-Peters-ALP-00036-BX3-20220224-5v1-1-hg/ALP-00036-BX3/outs/possorted_genome_bam.bam, \
		/data/LiverTransplantSingleCell/data/ALP00036/10X-Peters-ALP-00036-BX4-20220228-5v1-1-hg/ALP-00036-BX4/outs/possorted_genome_bam.bam, \
		/data/LiverTransplantSingleCell/data/ALP00036/10X-Peters-ALP-00036-BX5-20220406-5v1-1-hg/ALP-00036-BX5/outs/possorted_genome_bam.bam \
	-l PreTXP_,BX2_,BX3_,BX4_,BX5_ \
	-b /data/LiverTransplantSingleCell/data/ALP00036/10X-Peters-ALP-00036-20210728-5v1-1hg/ALP-00036/outs/filtered_feature_bc_matrix/barcodes.tsv.gz, \
		/data/LiverTransplantSingleCell/data/ALP00036/10X-Peters-ALP-00036-BX2-20220224-5v1-1-hg/ALP-00036-BX2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz, \
		/data/LiverTransplantSingleCell/data/ALP00036/10X-Peters-ALP-00036-BX3-20220224-5v1-1-hg/ALP-00036-BX3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz, \
		/data/LiverTransplantSingleCell/data/ALP00036/10X-Peters-ALP-00036-BX4-20220228-5v1-1-hg/ALP-00036-BX4/outs/filtered_feature_bc_matrix/barcodes.tsv.gz, \
		/data/LiverTransplantSingleCell/data/ALP00036/10X-Peters-ALP-00036-BX5-20220406-5v1-1-hg/ALP-00036-BX5/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
	-o /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/Seurat_Integration_0.5_SCT_08.30.23/SoupOrCell_7Jun2024/Patient6/ \
	-t 4