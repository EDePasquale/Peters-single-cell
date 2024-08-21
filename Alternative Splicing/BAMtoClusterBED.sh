#!/bin/bash

while IFS=$'\t' read SAMPLE BAM; do
	bsub <<EOF
#BSUB -W 10:00
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -M 32000
#BSUB -e results/${SAMPLE}/${SAMPLE}.err
#BSUB -o results/${SAMPLE}/${SAMPLE}.out
#BSUB -J $SAMPLE

module load python/2.7.5

mkdir results/$SAMPLE
cd results/$SAMPLE

python /data/salomonis2/software/AltAnalyze/import_scripts/BAMtoJunctionClusterBarcode.py --i $BAM --species Hs --c ../../cluster_files/Clusters_${SAMPLE}.txt

EOF
done < input_file_locations_all.txt
