#!/bin/bash

#$ -q long.q
#$ -N final_assign_taxo
#$ -M hans.schrieke@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=6G
#$ -V
#$ -cwd

# LOAD MODULE ENVIRONMENT
. /etc/profile.d/modules.sh

module purge
module load system/conda/4.4.10
source activate qiime2-2019.4
#conda list

TMP_PATH=$"/homedir/schrieke/stage/oligotyping/assignement/qiime/training-feature-classifiers/sequence_merged_format_good"

echo "TMP_PATH : " $TMP_PATH
cd $TMP_PATH




# Import data to qiime2 artefact

for fasta in *.fasta
do
  fullfilename=$(basename $fasta)
  sample=$(echo $fullfilename | cut -d'_' -f1)
  qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path $fasta \
    --output-path merged_$sample.qza
  #echo $sample " is succesfully imported"
done

echo "import data to qiime2 artefact : DONE"




# Assign Taxonomy to ASV

#mkdir taxonomy

for qza in *.qza
do
  fullfilename=$(basename $qza)
  sample=$(echo $fullfilename | cut -d'_' -f2 | cut -d'.' -f1)
  qiime feature-classifier classify-sklearn \
    --i-classifier /homedir/schrieke/stage/oligotyping/assignement/qiime/training-feature-classifiers/silva-132-99-nb-classifier.qza \
    --i-reads /homedir/schrieke/stage/oligotyping/assignement/qiime/training-feature-classifiers/sequence_merged_format_good/$qza  \
    --o-classification /homedir/schrieke/stage/oligotyping/assignement/qiime/training-feature-classifiers/sequence_merged_format_good/taxonomy/taxonomy_$sample.qza \
  #echo $sample " is succesfully assigned"
done

echo "assign taxonomy : DONE"




# Export fasta from qza

cd $TMP_PATH/taxonomy

mkdir sequence_output

for qza in *.qza
do 
  fullfilename=$(basename $qza)
  sample=$(echo $fullfilename | cut -d'_' -f2 | cut -d'.' -f1)
    qiime tools export \
        --input-path /homedir/schrieke/stage/oligotyping/assignement/qiime/training-feature-classifiers/sequence_merged_format_good/taxonomy/$qza \
        --output-path /homedir/schrieke/stage/oligotyping/assignement/qiime/training-feature-classifiers/sequence_merged_format_good/taxonomy/sequence_output/taxonomy_$sample
    cd /homedir/schrieke/stage/oligotyping/assignement/qiime/training-feature-classifiers/sequence_merged_format_good/taxonomy/sequence_output/taxonomy_$sample
    mv *.tsv taxonomy_$sample.tsv
    cd $TMP_PATH/taxonomy
    echo $qza " is now exported"
done

echo "fasta export : DONE" 





# Create visualization files 

#cd $TMP_PATH/taxonomy

for qza in *.qza
do
  fullfilename=$(basename $qza)
  sample=$(echo $fullfilename | cut -d'_' -f2 | cut -d'.' -f1)
  qiime metadata tabulate \
    --m-input-file /homedir/schrieke/stage/oligotyping/assignement/qiime/training-feature-classifiers/sequence_merged_format_good/taxonomy/$qza \
    --o-visualization /homedir/schrieke/stage/oligotyping/assignement/qiime/training-feature-classifiers/sequence_merged_format_good/taxonomy/taxonomy_$sample.qzv
  #echo $sample " is now available to watch"
done

echo "create visualization files : DONE"






source deactivate
