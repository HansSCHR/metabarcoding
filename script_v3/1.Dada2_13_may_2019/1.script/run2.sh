#!/bin/bash
#$ -q bigmem.q
#$ -N run2
#$ -M hans.schrieke@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=6G
#$ -V
#$ -cwd


#SCHRIEKE Hans 
#Dada2 on R 
# trunc(245,240)
# maxEE(2,4)
# assignTaxonomy(minBoot=80)

module load system/conda/5.1.0

source activate dada2-check

Rscript run2.R 
