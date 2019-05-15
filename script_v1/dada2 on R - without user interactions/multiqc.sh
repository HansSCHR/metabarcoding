#!/bin/bash
#$ -q normal.q
#$ -N multiqc
#$ -M hans.schrieke@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=6G
#$ -V
#$ -cwd

#SCHRIEKE Hans 
#Multiqc on R 

module load system/conda/5.1.0

source activate dada2

multiqc --version

cd /homedir/schrieke/Fastq
mkdir fastqc 

fastqc /homedir/schrieke/Fastq/*fastq.gz

mv /homedir/schrieke/Fastq/*.zip /homedir/schrieke/Fastq/fastqc
mv /homedir/schrieke/Fastq/*.html /homedir/schrieke/Fastq/fastqc

cd /homedir/schrieke/Fastq/fastqc
multiqc .
