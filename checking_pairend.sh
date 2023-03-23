#!/bin/bash
#This script for checking length on Pairends sequences from trimmed fastq.gz
#Author: Nathalia Portilla

function docker_run() { docker run --rm=True -u $(id -u):$(id -g) -v $(pwd):/data "$@" ;}
wordir=/home/bioinfo-sanger/Group_work/species/s.pneumo/results_fastqc_gz/trimmed_results/results_trimmed_gz
cd $wordir

#Create a empty csv to storage results
csv_file=checking_pairends_results.csv
echo "Genome,Forward Length,Reverse Length" > $csv_file

#Get names of genomes with PE on directory
for i in $(ls *_1.trimmed.fastq.gz); do
NAME=$(basename $i _1.trimmed.fastq.gz)
echo "$NAME"
j="${NAME}_1.trimmed.fastq.gz"
echo "$j"
k="${NAME}_2.trimmed.fastq.gz"
echo "$k"

#Count length of each sequences

j1=$(gzip -cdf $j| wc -l)
k1=$(gzip -cdf $k| wc -l)

#Add data to the file
echo "${NAME},$j1,$k1" >> $csv_file

done
