#!/bin/bash

for dir in /data/LiverTransplantSingleCell/data/*/*/*/outs/junctions_alt; do
  SAMPLE=$(basename $(dirname $(dirname $dir)))
  for i in $dir/*; do
    RENAMED=$(basename $i)
    cp $i for_altanalyze/$(echo $RENAMED | sed "s/possorted_genome_bam/$SAMPLE/")
  done
done
#./copy_cluster_junctions.sh 