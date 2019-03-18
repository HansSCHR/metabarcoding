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
Rscript dada2.R 
