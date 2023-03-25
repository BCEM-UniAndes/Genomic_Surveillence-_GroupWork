#!/bin/bash
#Author: Felipe Sierra
#This script for checking length on Pairends sequences from trimmed fastq.gz
function docker_run() { docker run --rm=True -u $(id -u):$(id -g) -v $(pwd):/data "$@" ;}
wordir=/home/bioinfo-sanger/s.pneumo/results_fastqc_gz/trimmed_results/results_trimmed_gz
cd $wordir
#Create a empty csv to storage results
csv_file=popunk_input.csv
echo "#ID,R1,R2">$csv_file
#Get names of genomes with PE on directory
for i in $(ls *_1.trimmed.fastq.gz); do
NAME=$(basename $i _1.trimmed.fastq.gz)
echo "$NAME"
j=${NAME}_1.trimmed.fastq.gz
echo "$j"
k=${NAME}_2.trimmed.fastq.gz
echo "$k"
#Merge data to the file
echo "${NAME},$j,$k" >> $csv_file;
#convert cvs into  tsv
sed 's/,/\t/g' popunk_input.csv > poppunkdir.tsv
done
