#!/bin/bash

#$ -q normal.q
#$ -N assign_taxo
#$ -M hans.schrieke@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=6G
#$ -V
#$ -cwd

# LOAD MODULE ENVIRONMENT
. /etc/profile.d/modules.sh

module purge
module load system/conda/5.1.0
source activate qiime2-2019.4
#conda list

TMP_PATH=$"/homedir/schrieke/stage/oligotyping/assignement/qiime/training-feature-classifiers"

echo "TMP_PATH : " $TMP_PATH

# Import data to qiime2 artefact

qiime tools import \--type 'FeatureData[Sequence]' \--input-path S18_S1_L001_MERGED_uppercase.fasta \--output-path merged.qza

echo "import data to qiime2 artefact : DONE"

# Assign Taxonomy to ASV

qiime feature-classifier classify-sklearn \
  --i-classifier /homedir/schrieke/stage/oligotyping/assignement/qiime/training-feature-classifiers/silva-132-99-nb-classifier.qza \
  --i-reads /homedir/schrieke/stage/oligotyping/assignement/qiime/training-feature-classifiers/merged.qza  \
  --o-classification /homedir/schrieke/stage/oligotyping/assignement/qiime/training-feature-classifiers/taxonomy.qza \

echo "assign taxonomy : DONE"

qiime metadata tabulate \
  --m-input-file /homedir/schrieke/stage/oligotyping/assignement/qiime/training-feature-classifiers/taxonomy.qza  \
  --o-visualization /homedir/schrieke/stage/oligotyping/assignement/qiime/training-feature-classifiers/taxonomy.qzv

echo "create visualization file : DONE"

source deactivate 
