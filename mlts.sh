#!/bin/bash
#File name:mlts.sh
#This script for mlts from trimmed fastq.gz

function docker_run() { docker run --rm=True -u $(id -u):$(id -g) -v $(pwd):/data "$@" ;}
wordir=/home/bioinfo-sanger/Group_work/species/s.pneumo/results_fastqc_gz/trimmed_results/results_trimmed_gz
cd $wordir
docker pull staphb/mlst
docker pull staphb/any2fasta
docker pull ncbi/blast



for i in $(ls *_1.trimmed.fastq.gz); do

NAME=$(basename $i _1.trimmed.fastq.gz)
echo "$NAME"
j="${NAME}_1.trimmed.fastq.gz"
echo "$j"
k="${NAME}_2.trimmed.fastq.gz"
echo "$k"

docker_run staphb/seroba seroba runSerotyping seroba_k71_14082017 ./$j ./$k ${NAME}_serotype_output;

done
