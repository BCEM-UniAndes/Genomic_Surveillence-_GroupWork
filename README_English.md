# Group Work- Genomic Surveillance

**Author:** Nathalia Portilla & Felipe Sierra

**Github:** https://github.com/portillanath/Genomic_Surveillence-_GroupWork

In this tutorial, we summarize the whole dataset process workflow for the *S.pneumoniae* (20 sequences with their respective forward *_1 and reverse *_2) and *S.agalactiae* (20 sequences with their respective forward *_1 and reverse *_2) species.

# 1.**NGS data quality control**

A typical whole genome sequencing process involves sample preparation, library preparation and sequencing. Errors occurring at each of these steps can negatively impact the quality of the sequence information.We will be using FastQC tool to carry out quality control on reads. First, check fastqc installation with the command `docker pull staphb/fastqc`

Then download the *S.pneumoniae* and *S.agalacteae* zips from Google Drive and transfer them to your working directory, 🌟 as a shortcut if you are using Linux, go to the file explorer and press `ctrl+L` and write: `sftp://youruser@yourIP` , then ingress the password. This is a more accessible way to upload and download files. Once have the archive in your virtual machine `unzip s.agalactea.zip` and `unzip pneumoniae.zip` .You must obtain two folders with fastq.gz for each forward and reverse sequence, then execute fastqc:

```bash
docker_run staphb/fastqc fastqc *.fastq.gz
```

The final product per sequence is an HTML with the analysis (visualizations), and a stats.csv file, exploring the different metrics and quality distributions. It is recommended for good practice to move and store data on folders by their content file separate the fastq.gz and HTML files, also has a consistent nomenclature on data as follows can help

| Raw_Sequence_s.pneumoniae | New_code |
| --- | --- |
| 20925_3#47_1.fastq.gz | pneumoniae1 |
| 20925_3#47_2.fastq.gz | pneumoniae1 |
| 20925_3#79_1.fastq.gz | pneumoniae2 |
| 20925_3#79_2.fastq.gz | pneumoniae2 |
| 22335_7#18_1.fastq.gz | pneumoniae3 |
| 22335_7#18_2.fastq.gz | pneumoniae3 |
| 22335_7#80_1.fastq.gz | pneumoniae4 |
| 22335_7#80_2.fastq.gz | pneumoniae4 |
| 22335_7#82_1.fastq.gz | pneumoniae5 |
| 22335_7#82_2.fastq.gz | pneumoniae5 |
| 22335_7#118_1.fastq.gz | pneumoniae6 |
| 22335_7#118_2.fastq.gz | pneumoniae6 |
| 22335_7#140_1.fastq.gz | pneumoniae7 |
| 22335_7#140_2.fastq.gz | pneumoniae7 |
| 22335_7#148_1.fastq.gz | pneumoniae8 |
| 22335_7#148_2.fastq.gz | pneumoniae8 |
| 22335_7#153_1.fastq.gz | pneumoniae9 |
| 22335_7#153_2.fastq.gz | pneumoniae9 |
| 22335_8#2_1.fastq.gz | pneumoniae10 |
| 22335_8#2_2.fastq.gz | pneumoniae10 |
| 22335_8#12_1.fastq.gz | pneumoniae11 |
| 22335_8#12_2.fastq.gz | pneumoniae11 |
| 22335_8#21_1.fastq.gz | pneumoniae12 |
| 22335_8#21_2.fastq.gz | pneumoniae12 |
| 22335_8#53_1.fastq.gz | pneumoniae13 |
| 22335_8#53_2.fastq.gz | pneumoniae13 |
| 22335_8#116_1.fastq.gz | pneumoniae14 |
| 22335_8#116_2.fastq.gz | pneumoniae14 |
| 28522_3#32_1.fastq.gz | pneumoniae15 |
| 28522_3#32_2.fastq.gz | pneumoniae15 |
| 28522_3#153_1.fastq.gz | pneumoniae16 |
| 28522_3#153_2.fastq.gz | pneumoniae16 |
| 28522_3#198_1.fastq.gz | pneumoniae17 |
| 28522_3#198_2.fastq.gz | pneumoniae17 |
| 28522_3#201_1.fastq.gz | pneumoniae18 |
| 28522_3#201_2.fastq.gz | pneumoniae18 |
| 28522_3#237_1.fastq.gz | pneumoniae19 |
| 28522_3#237_2.fastq.gz | pneumoniae19 |
| 28522_3#241_1.fastq.gz | pneumoniae20 |
| 28522_3#241_2.fastq.gz | pneumoniae20 |

# 2. **Adapter Trimming**

We use **trimmomatic** tool to remove adaptors, to trim low quality reads and to remove short sequences. We need to check of having the respective adapter sequence on the working directory, for our dataset is TruSeq3-PE.fa. 

<aside>
💡 To optimizated coomands, loops repeat instructions, which is very handful to large datasets, we present here one way. Let’s create a new bash script using `nano [trimming.sh](http://trimming.sh)` then exceute with `bash trimming.sh`

</aside>

```bash
#!/bin/bash
#Author: Nathalia Portilla
#This script for triming fastq.gz

function docker_run() { docker run --rm=True -u $(id -u):$(id -g) -v $(pwd):/da>
wordir=/home/bioinfo-sanger/Group_work/species/s.pneumo/results_fastqc_gz/trimm>
cd $wordir

for i in $(ls *_1.trimmed.fastq.gz); do

NAME=$(basename $i _1.trimmed.fastq.gz)
echo "$NAME"
j="${NAME}_1.trimmed.fastq.gz"
echo "$j"
k="${NAME}_2.trimmed.fastq.gz"
echo "$k"

docker_run staphb/trimmomatic trimmomatic PE $j $k ${NAME}_1.trimmed.fastq.gz /dev/null 
${NAME}_2.trimmed.fastq.gz /dev/null ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:6:30 MINLEN:50

done
```

We move all the `mv *_trimmed.fastq.gz results_trimmed` into a new directory. We let here a report of the surviving porcentage sequences after trimming.

| Sample | Both_Survive |
| --- | --- |
| 20925_3#47 | 89.19% |
| 20925_3#79 | 89.03% |
| 22335_7#18 | 91.04% |
| 22335_7#80 | 90.90% |
| 22335_7#82 | 90.64% |
| 22335_7#118 | 90.56% |
| 22335_7#140 | 90.59% |
| 22335_7#148 | 90.55% |
| 22335_7#153 | 91.14% |
| 22335_8#2 | 90.60% |
| 22335_8#12 | 90.86% |
| 22335_8#21 | 90.64% |
| 22335_8#53 | 90.44% |
| 22335_8#116 | 85.43% |
| 28522_3#32 | 91.74% |
| 28522_3#153 | 90.20% |
| 28522_3#198 | 92.92% |
| 28522_3#201 | 92.62% |
| 28522_3#237 | 91.04% |
| 28522_3#241 | 89.08% |

# 3. **FastQC to trimmed sequences**

Run again fastq to check changes on quality for trimmed sequences. As well you can explore htmls files

```bash
docker_run staphb/fastqc fastqc *.trimmed.fastq.gz
```

# 4. Taxonomy Assignation

**Kraken** is a taxonomic sequence classifier to short DNA reads doing an examination of the *k*-mers within a read and querying a database with those *k*-mers, the database contains a mapping of every k-mer in Kraken’s genomic library to the lowest common ancestor (LCA) of all genomes contain that k-mer. First, pull the docker image `docker pull staphb/kraken2` . ⭐It is recommended copy the database inside your working directory  `cp -r ~/Data/Group_work/s.pneumo/k2_standard_08gb_20221209 .` Once this is set up, create a bash script and execute:****

```bash
#!/bin/bash
#Author: Nathalia Portilla
#This script for running Kraken FOR fastq.gz sequences 

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
docker_run staphb/kraken2 kraken2 --use-names --db k2_standard_08gb_20221209 --paired $j $k --report ${NAME}.kraken.txt --output -
done

```

As a good practice for quality, we will calculate the percentage of “x” species, in this case: `grep "Streptococcus pneumoniae" kraken_report.txt` and execute the bash script `bash percentage_calculator.sh`

# 5. De-Novo G**enome Assembly**

**Genome assembly** is the step to determine how reads fit together by looking for overlaps between them. Where a genome is pieced together without a reference sequence to compare is called ***de novo assembly,*** this doesn't produce complete genomes, it is going to generate multiple contiguous sequences ***(contigs)*** with an arbitrary order, but with a close reference, it is possible to improve it.

For this analysis we are going to run **SPADES**, please pull the docker image `docker pull staphb/spades` To run into multiple strains use:

```bash
#!/bin/bash
#Author: Nathalia Portilla
#This script for doing  de novo-assembly from a list of trimmed fastq.gz
function docker_run() { docker run --rm=True -u $(id -u):$(id -g) -v $(pwd):/data "$@" ;}
wordir=/home/bioinfo-sanger/Data/Group_work/s.pneumo/s.pneumo/results_fastqc_gz/trimmed_results/results_trimmed_gz
cd $wordir
for i in $(ls *_1.trimmed.fastq.gz); do
NAME=$(basename $i _1.trimmed.fastq.gz)
echo "$NAME"
j="${NAME}_1.trimmed.fastq.gz"
echo "$j"
k="${NAME}_2.trimmed.fastq.gz"
echo "$k"
docker_run staphb/spades spades.py -1 $j -2 $k -t 16 --careful --cov-cutoff auto -o ${NAME};
done
```

# 6. Quast

Quast is used for quality assembly with or without a reference genomes. In order to run with multiple assemblies, generate a new bash script `nano [quast.sh](http://quast.sh)` and execute. ⭐**Note** as a recommendation move in your working directory which contains your fastq.gz pairend files and contigs.fasta from assembly, the reference genomes for your species (*S.pneumoniae* Reference_sequence_GPSC46.fa and for *S.agalacteae* Reference_sequence_GBS_cc17.fasta), the PROKKA_03052023.gff both avaible on Data/Section_three

```bash
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
```

# 7.Prokka

Prokka annotated a draft genome sequence. Annotates features like protein coding regions, tRNA, rARN genes. First use protein coding regions found with prodigal and then the function is predicted  by protein databases. First pull the docker image  `docker pull staphb/prokka` make a new `nano [prokka.sh](http://prokka.sh)` and check that all your contgis.fasta are inside the work directory

```bash
#!/bin/bash
#This script for doing  de novo-assembly from a list of timmed fastq.gz

function docker_run() { docker run --rm=True -u $(id -u):$(id -g) -v $(pwd):/data "$@" ;}
wordir=/home/bioinfo-sanger/Data/Group_work/s.pneumo/s.pneumo/results_fastqc_gz/trimmed_results/results_trimmed_gz
cd $wordir

for i in $(ls *_contigs.fasta); do
docker_run staphb/prokka prokka ./$i --outdir prokka_out$i;
done
```

# 8. **Serotyping**

**8.1 *S.pneumoniae* serotyping pipeline**

Let’s distinguish that each species has a different serotype pipeline. First, we run the analysis for the *S.pneunomiae* dataset. For this we are gonna use **SeroBA** tool, make sure that is installed `docker pull staphb/seroba` and storage the k-mer database **seroba_k71_14082017 serobamak** ([https://drive.google.com/drive/folders/1Jy3UDUOqBiIvwqtc6FMA2PAuxs6kjnni](https://drive.google.com/drive/folders/1Jy3UDUOqBiIvwqtc6FMA2PAuxs6kjnni)) on the same folder than your trimmed sequences (working directory)

Create a new bash script [`nano Serotype.sh`](http://Serotype.sh) and execute with `bash Serotype.sh`

```bash
#!/bin/bash
#File name:Serotype.sh
#This script for serotyping S.pneumoniae from trimmed fastq.gz

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

docker_run staphb/seroba seroba runSerotyping seroba_k71_14082017 ./$j ./$k ${NAME}_serotype_output;

done
```

We move all the outputs into a new folder  `mv *_serotype_output serotype_results`, then compile all data inside the new directory with the command `docker_run staphb/seroba seroba summary ./` , must get a tsv file.

**8.2 *S.agalactaeae* serotyping** 

For running srtst2 you will use the raw data (not the trimmed) 

```bash
#!/bin/bash
#This script for serotyping multiple genome reads of s.agalactaeae

function docker_run() { docker run --rm=True -u $(id -u):$(id -g) -v $(pwd):/data "$@" ;}
wordir=/home/bioinfo-sanger/Data/Group_work/s.pneumo/s.agalactiae/lanes_2.txt/
cd $wordir
mkdir -p serotyping_output

for i in $(ls *_1.fastq.gz); do

NAME=$(basename $i _1.fastq.gz)
echo "$NAME"
j="${NAME}_1.fastq.gz"
echo "$j"
k="${NAME}_2.fastq.gz"
echo "$k"

docker_run staphb/srst2 srst2 --input_pe $j $k --output ./serotyping_out/${NAME}_output --log --gene_db analysis/clean_data/GBS-SBG.fasta;

done
```

# 9.**MLTS**

**9.1 MLTS tool (Recommended, run this)** ✅

Download the following images using the command: 

`docker pull staphb/mlst`

`docker pull staphb/any2fasta`

`docker pull ncbi/blast`

Then run the tool using the contigs of all the genomes on your workdir and obtein a csv file

```bash
##for running with s.pneumonae
docker_run staphb/mlst mlst --legacy --scheme spneumoniae *contigs.fasta > mlst_results.csv

##for running with s.agalactiae
docker_run staphb/mlst mlst --legacy --scheme sagalactiae *contigs.fasta > mlst_results.csv
```

**9.2 SRST2 tool**

Check the docker image using `docker pull staphb/srst2`  and a fasta sequence database of all alleles and a tab delim file to match ST profiles as a combination of alleles, download from `docker_run staphb/srst2 [getmlst.py](http://getmlst.py/) --species "Streptococcus pneumonia"`

Must get an alleles_fasta, mlst_data_download_Streptococcus_pneumoniae_None.log, profiles_csv and Streptococcus_pneumoniae.fasta files on your working directory

Then create a bash script `nano [srst2.sh](http://srst2.sh)` and execute `bash srst2.sh`

```bash

#!/bin/bash
#This script for doing MLTS using srst2 tool from trimmed fastq.gz

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
docker_run staphb/srst2 srst2 --input_pe $j $k --output srst2output --log --mlst_db Streptococcus_pneumoniae.fasta --mlst_defin>
done
```

# 10. AMR

ABRicate tool carries out mass screening of contigs for antimicrobial resistance or virulence genes using multiple databases. It only supports contigs, It only detects acquired resistance genes, **NOT** point mutations, It uses a DNA sequence database, not protein and It needs BLAST+ >= 2.7 and any2fasta to be installed. 

Make sure of having the following docker images, by command:

`docker pull staphb/abricate`

`docker pull ncbi/blast`

`docker pull staphb/any2fasta`

Using `docker_run staphb/abricate abricate --list` you can explore the different databases loaded. Now, for running the tool use the command: `docker_run staphb/abricate abricate --db resfinder --quiet *contigs.fasta > amr_results.tab`

# 11. Mapping

We are going to use bwa (Burrows-Wheeler Aligner) and samtools, please download it from the repository using the commands:

`docker pull staphb/bwa`

`docker pull staphb/samtools`

**11.1 *S.pneumoniae* pipeline**

a) **Get the references:** A single reference genome is designated to represent a species for comparative analysis. A complete reference genome should be of high-quality annotation and meets the highest level of experimental support for structural and functional annotation. For *S.pneumoniae* we are going to use sequences GPSC46 and GPSC33, we will use the alignment results to determine which lineage the reads belong to. Make sure that these are in your work directory. 

Reference_sequence_GPSC46.fa

Reference_sequence_GPSC33.fa

b**) Check pair-end sequences**

💡It’s important to count the number of lines in each of the read files and **check they have the same number,** if your two paired files are unbalanced, something has gone wrong (e.g., filtering/trimming went wrong and corrupted the output, or maybe files from different samples are being used)**.** For this, create a new script `nano hecking_pairend.sh`  and execute with `bash checking_pairend.sh`

```bash
#!/bin/bash
#Author: Nathalia Portilla
#This script for checking length on Pairends sequences from trimmed fastq.gz
function docker_run() { docker run --rm=True -u $(id -u):$(id -g) -v $(pwd):/data "$@" ;}
wordir=/home/bioinfo-sanger/Group_work/species/s.pneumo/results_fastqc_gz/trimmed_results/results_trimmed
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
```

c) ****Creating an alignment:**** 

Create an index of the reference sequence `docker_run staphb/bwa bwa index Reference_sequence_GPSC46.fa` and `docker_run staphb/bwa bwa index Reference_sequence_GPSC33.fa`

Then create a script for each linage, you can make a copy of the same script and just replace Reference_sequence_GPSC33.fa and the name output.

```bash
#!/bin/bash
#This script for alignment from trimmed fastq.gz
function docker_run() { docker run --rm=True -u $(id -u):$(id -g) -v $(>
wordir=/home/bioinfo-sanger/Group_work/species/s.pneumo/results_fastqc_gz/trimmed_results/results_trimmed_gz
cd $wordir
for i in $(ls *_1.trimmed.fastq.gz); do
NAME=$(basename $i _1.trimmed.fastq.gz)
echo "$NAME"
j="${NAME}_1.trimmed.fastq.gz"
echo "$j"
k="${NAME}_2.trimmed.fastq.gz"
echo "$k"
docker_run staphb/bwa bwa mem Reference_sequence_GPSC46.fa $j $k > GPSC46bwa.sam
done
```

You are going to get a SAM file, and inspect it `head GPSC46bwa.sam`  file content (reference name, alignment location, CIGAR string, the read sequence itself, quality etc..). . It is good practice to convert your SAM files to BAM files (Binary Alignment Map) using `docker_run staphb/samtools samtools sort GPSC46bwa.sam -o GPSC46bwa.bam`

The next step is to index the BAM file; indexing, which relies on sorted data, enables faster searches downstream. Type this command: `docker_run staphb/samtools samtools index GPSC46bwa.bam`  use **`ls -alh`** to check the size, everything must be >0. 

d) ****Assessing the alignment:**** The 2nd column in the SAM file contains the flag for the read alignment. If the flag includes the number 4 flag in its makeup then the read is unmapped, if it doesn't include the number 4 then it is mapped.

- For knowing  read alignment mapped `docker_run staphb/samtools samtools view -c -F4 GPSC46bwa.bam`For knowing unmapped `docker_run staphb/samtools samtools view -c -f4 GPSC46bwa.bam`

<aside>
🛑 BUT! these are not reads, a read could be mapped equally well to multiple positions, 0ne will be called the primary alignment, and others secondary alignment and in another hand, a read could be split into two parts, So to get the true number of mapped reads you need to count the alignment that do not have flags 4 (unmapped), 256 (not primary), and 2048 (supplementary) = 4+256+2048 = 2308

</aside>

- For knowing reads mapped `docker_run staphb/samtools samtools view -c -F2308 GPSC46bwa.bam`
- Statics  `docker_run staphb/samtools samtools stats GPSC46bwa.bam > GPSC46bwa_bamstats.txt`

# 12. Variant Calling

map sequence reads to a reference sequence and then identify the variants (both single nucleotide polymorphisms (SNPs) and small insertions and deletions (indels)) in the sequence data. First, load `docker pull staphb/snippydocker pull staphb/snippy`

### Add the reference Genome

```bash
#!/bin/bash
#This script for alignment from trimmed fastq.gz
function docker_run() { docker run --rm=True -u $(id -u):$(id -g) -v $(pwd):/data "$@" ;}
wordir=/home/bioinfo-sanger/Group_work/species/s.pneumo/results_fastqc_gz/trimmed_results/results_tr>
cd $wordir
for i in $(ls *_1.trimmed.fastq.gz); do
NAME=$(basename $i _1.trimmed.fastq.gz)
echo "$NAME"
j="${NAME}_1.trimmed.fastq.gz"
echo "$j"
k="${NAME}_2.trimmed.fastq.gz"
echo "$k"
docker_run staphb/snippy snippy --outdir ${NAME}_GPSC46_folder --R1 $j --R2 $k --ref Reference_seque>
done
```

Use this resource for helping you interpretate  VCF files [https://www.internationalgenome.org/wiki/Analysis/vcf4.0](https://www.internationalgenome.org/wiki/Analysis/vcf4.0)

- For doing a **summary of variants data** use: `head -n5 [snps.tab](http://snps.tab)` inside your directory per sequence
- You can also generate a pseudogenome doing remplacing the reference databases with the snips identified, you can explore that on the file  `head snps.consensus.subs.fa`

****creating multiple sequence alignment****

We need to create a Tab separated input file (which contains a list of the forward and reverse reads in the following format: ID, names of R1 reads and names of R2 reads). 

```bash
#!/bin/bash
#Author: Nathalia Portilla
#This script for generate a tab delimiter list
function docker_run() { docker run --rm=True -u $(id -u):$(id -g) -v $(pwd):/data "$@" ;}
wordir=/home/bioinfo-sanger/Group_work/species/s.pneumo/results_fastqc_gz/trimmed_results/results_tr>
cd $wordir
#Create a empty tab to storage results
tab_file=tabdelim.csv
echo "#name,R1,R2" > $tab_file
#Get names of genomes with PE on directory
for i in $(ls *_1.trimmed.fastq.gz); do
NAME=$(basename $i _1.trimmed.fastq.gz)
echo "$NAME"
j="${NAME}_1.trimmed.fastq.gz"
echo "$j"
k="${NAME}_2.trimmed.fastq.gz"
echo "$k"
#Add data to the file
echo "${NAME},$j,$k" >> $tab_file
done
#Convert that into a tab delim txt
sed 's/,/\t/g' tabdelim.csv >  tabdelim.tab
```

Now generate a output script (make sure of having your respective reference genome on your directory)  `docker_run staphb/snippy snippy-multi tabdelim-tab --ref Reference_sequence_GPSC46.fa --cpus 2 > [runme.sh](http://runme.sh)` And then, generate the alignment `docker_run staphb/snippy sh ./runme.sh` Use `ls -l core.*` to list the alignment files.

![Untitled](Group%20Work-%20Genomic%20Surveillance%20fd65548faba249b4b127be28e4ac5f8f/Untitled.png)

The file contained the multiple sequence alignment have the name `core.full.aln` but this mask sequences with low confidence. To rescue it use `docker_run staphb/snippy snippy-clean_full_aln core.full.aln > clean.full.aln`

# 13. Phylogenetics

Phylogenetics is quite important, it is the study of evolutionary relationships among biological entities – often species, individuals or genes (which may be referred to as taxa).There are many tools for making phylogenetic analysis. Here we describe many that you can use and how to use them:

- **phylogenetic tree using FastTree**
    1. ****Make a SNP-only alignment using snp-sites****
    
    **We use SNPs in order to make the phylogenetic assignation in a faster way from the multialigned genomes** 
    
    *First, remove all the invariant sites and create a SNP-only multiple sequence alignment. We will use output from snippy runs described in the previous page. Run the command:*
    
    `docker_run staphb/snp-sites snp-sites -o clean.full.SNPs.aln clean.full.aln`
    
    **We can see how many invariant sites were removed (and what proportion of A, T, G, C they were) using:**
    
    `docker_run staphb/snp-sites snp-sites -C clean.full.aln`
    
    1. Now we use FasTree
    
    `docker_run staphb/fasttree FastTree -nt -gtr clean.full.SNPs.aln > clean.full.SNPs.aln.tree`
    
- Recombination ussing Gubbins

- Clustering by kmers with PopPUNK
    
    1.Create a file which lists your sample names and paths to their sequence data.
    
    ```bash
    #!/bin/bash
    #Author: Felipe Sierra
    #This script for checking length on Pairends sequences from trimmed fastq.gz
    function docker_run() { docker run --rm=True -u $(id -u):$(id -g) -v $(pwd):/data "$@" ;}
    wordir=/home/bioinfo-sanger/Group_work/species/s.agalactiae/lanes_2.txt/analysis/clean_data
    cd $wordir
    
    #Create a empty csv to storage results
    tsv_file=popunk_input.tsv
    rm popunk_input.tsv
    #Get names of genomes with PE on directory
    for i in $(ls *_1.trimmed.fastq.gz); do
    	NAME=$(basename $i _1.trimmed.fastq.gz)
    	echo "$NAME"
    	j=$wordir/${NAME}_1.trimmed.fastq.gz
    	echo $j
    	k=$wordir/${NAME}_2.trimmed.fastq.gz
    	echo "$k"
    	#Count length of each sequences
    	# j1=$($j| "${NAME}_1.trimmed.fastq.gz")
    	# k1=$($k| "${NAME}_2.trimmed.fastq.gz")
    	#Add data to the file
    	echo -e ${NAME}\\t$j\\t$k >> popunk_input.tsv;
    
    done
    ```
    
    - [ ]  Presentación miércoles
    - [ ]  Presentación Jueves
    - [ ]  Subir archivos al github
    - [ ]  Traducir el notion a español
    - [ ]  filogenias