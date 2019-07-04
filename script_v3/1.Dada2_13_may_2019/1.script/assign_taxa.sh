#!/bin/bash
#$ -q bigmem.q
#$ -N assign_taxa
#$ -M hans.schrieke@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=6G
#$ -V
#$ -cwd


#SCHRIEKE Hans 
#Dada2 on R 
# assignTaxonomy(minBoot=80)

module load system/conda/5.1.0

source activate dada2-check

Rscript assign_taxa.R 
