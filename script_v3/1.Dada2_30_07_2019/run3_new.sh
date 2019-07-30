#!/bin/bash
#$ -q long.q
#$ -N run3_new
#$ -M hans.schrieke@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=10G
#$ -V
#$ -cwd


#SCHRIEKE Hans 
#Dada2 on R 
# trunc(245,240)
# maxEE(2,4)
# assignTaxonomy(minBoot=80)

. /etc/profile.d/modules.sh
module purge
module load system/conda/4.4.10
source activate dada2-check
conda list

Rscript run3_new.R 
