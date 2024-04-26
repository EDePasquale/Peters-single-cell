#BSUB -W 48:00
#BSUB -o BAM_merge.out
#BSUB -J BAM_merge
#BSUB -M 10000
#BSUB -n 8
#BSUB -R "span[hosts=1]"

module load samtools/1.9.0

samtools merge --threads 8 /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/splicing/results/merged.bam \
	/data/LiverTransplantSingleCell/data/ALP00011/10X-Peters-ALP-00011-BX2-20211108-5v1-1hg/ALP-00011-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-20210514-5v1-1hg/ALP-00017/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-BX2-20220407-5v1-1-hg/ALP-00017-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00011/10X-Peters-ALP-00011-BX2-20211108-5v1-1hg/ALP-00011-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-20210514-5v1-1hg/ALP-00017/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-BX2-20220407-5v1-1-hg/ALP-00017-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00011/10X-Peters-ALP-00011-BX2-20211108-5v1-1hg/ALP-00011-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-20210514-5v1-1hg/ALP-00017/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-BX2-20220407-5v1-1-hg/ALP-00017-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00011/10X-Peters-ALP-00011-BX2-20211108-5v1-1hg/ALP-00011-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-20210514-5v1-1hg/ALP-00017/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-BX2-20220407-5v1-1-hg/ALP-00017-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00011/10X-Peters-ALP-00011-BX2-20211108-5v1-1hg/ALP-00011-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-20210514-5v1-1hg/ALP-00017/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-BX2-20220407-5v1-1-hg/ALP-00017-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00011/10X-Peters-ALP-00011-BX2-20211108-5v1-1hg/ALP-00011-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-20210514-5v1-1hg/ALP-00017/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-BX2-20220407-5v1-1-hg/ALP-00017-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00011/10X-Peters-ALP-00011-BX2-20211108-5v1-1hg/ALP-00011-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-20210514-5v1-1hg/ALP-00017/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-BX2-20220407-5v1-1-hg/ALP-00017-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00011/10X-Peters-ALP-00011-BX2-20211108-5v1-1hg/ALP-00011-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-20210514-5v1-1hg/ALP-00017/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-BX2-20220407-5v1-1-hg/ALP-00017-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00011/10X-Peters-ALP-00011-BX2-20211108-5v1-1hg/ALP-00011-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-20210514-5v1-1hg/ALP-00017/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00017/10X-Peters-ALP-00017-BX2-20220407-5v1-1-hg/ALP-00017-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00011/10X-Peters-ALP-00011-BX2-20211108-5v1-1hg/ALP-00011-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00036/10X-Peters-ALP-00036-BX2-20220224-5v1-1-hg/ALP-00036-BX2/outs/possorted_genome_bam.bam \
	/data/LiverTransplantSingleCell/data/ALP00036/10X-Peters-ALP-00036-BX3-20220224-5v1-1-hg/ALP-00036-BX3/outs/possorted_genome_bam.bam