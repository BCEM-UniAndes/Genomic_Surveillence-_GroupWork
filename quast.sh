#!/bin/bash
#This script for doing  de novo-assembly from a list of timmed fastq.gz

function docker_run() { docker run --rm=True -u $(id -u):$(id -g) -v $(pwd):/data "$@" ;}
wordir=/home/bioinfo-sanger/Data/Group_work/s.pneumo/s.pneumo/results_fastqc_gz/trimmed_results/results_trimmed_gz/spades_output_pneumoniae
cd $wordir
mkdir -p quast_output

for i in $(ls *_1.trimmed.fastq.gz); do

NAME=$(basename $i _1.trimmed.fastq.gz)
echo "$NAME"
j="${NAME}_1.trimmed.fastq.gz"
echo "$j"
k="${NAME}_2.trimmed.fastq.gz"
echo "$k"
l="${NAME}_contigs.fasta"

docker_run staphb/quast quast.py ./$l -r Reference_sequence_GPSC46.fa -g PROKKA_03052023.gff -1 ./$j -2 ./$k -o ./quast_output/${NAME};

done
