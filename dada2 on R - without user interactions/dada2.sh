#!/bin/bash
#$ -q normal.q
#$ -N dada2
#$ -M hans.schrieke@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=6G
#$ -V
#$ -cwd


#SCHRIEKE Hans 
#DADA2 on R 


module load system/conda/5.1.0

source activate dada2

cd /homedir/schrieke/Fastq/
mkdir tables

Rscript dada2.R 

mv /homedir/schrieke/Fastq/stats.csv /homedir/schrieke/Fastq/tables
mv /homedir/schrieke/Fastq/taxa.csv /homedir/schrieke/Fastq/tables
mv /homedir/schrieke/Fastq/seqtab.csv /homedir/schrieke/Fastq/tables

cd /homedir/schrieke/Fastq/tables

#Rscript phyloseq.R
