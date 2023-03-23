#!/bin/bash
#Author: Nathalia Portilla
#This script for triming fastq.gz

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

docker_run staphb/trimmomatic trimmomatic PE $j $k ${NAME}_1.trimmed.fastq.gz /dev/null ${NAME}_2.trimmed.fastq.gz /dev/null ILLUMINACLIP:~/Data/Seccion2/Section_two/Adapter_trimming/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36;

done
