#!/bin/bash
#This script for alignment from trimmed fastq.gz

function docker_run() { docker run --rm=True -u $(id -u):$(id -g) -v $(pwd):/data "$@" ;}
wordir=/home/bioinfo-sanger/Group_work/species/s.pneumo/results_fastqc_gz/trimmed_results/results_trimmed_gz
cd $wordir
for i in $(ls *_1.trimmed.fastq.gz); do
NAME=$(basename $i _1.trimmed.fastq.gz)
echo "$NAME"
j="${NAME}_1.trimmed.fastq.gz"
echo "$j"
k="${NAME}_2.trimmed.fastq.gz"
echo "$k"
docker_run staphb/bwa bwa mem Reference_sequence_GPSC33.fa $j $k > GPSC33bwa.sam
done
