#!/bin/bash
#$ -q long.q
#$ -N rmarkdown
#$ -M hans.schrieke@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=10G
#$ -V
#$ -cwd

. /etc/profile.d/modules.sh
module purge
module load system/conda/4.4.10
source activate dada2-check
conda list

Rscript metabarcoding.Rmd 