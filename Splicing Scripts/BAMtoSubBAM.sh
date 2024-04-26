DIR=$(pwd)

cat <<EOF
#BSUB -W 24:00
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -M 128000
#BSUB -e merged_split_after.err
#BSUB -o merged_split_after.out

cd $DIR

module load python/2.7.5

python /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/splicing/extract-reads-pysam.py \
--i /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/splicing/results/merged.bam \
--o /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/splicing/results/merged_split \
--f /data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/splicing/cluster_files/Clusters_all.txt \

EOF