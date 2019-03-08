#!/bin/bash
#$ -q short.q
#$ -N multiqc
#$ -M hans.schrieke@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=6G
#$ -V
#$ -cwd



# Load required modules 

module load system/conda/5.1.0



# Install qiime2 

wget https://data.qiime2.org/distro/core/qiime2-2018.11-py35-linux-conda.yml
conda env create -n qiime --file qiime2-2018.11-py35-linux-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2018.11-py35-linux-conda.yml



# Activate qiime2 environment

conda activate qiime 



# Create a work folder "dada2"

mkdir dada2


# Import data to qiime2 artefact (Casava format)

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path data \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza


# Quality Score with qiime2 (demultiplexed sequences)

qiime demux summarize \
  --i-data demux-paired-end.qza \
  --o-visualization demux.qzv

# Adapters removal 
qiime cutadapt trim-paired \
        --i-demultiplexed-sequences demux-paired-end.qza \
        --p-cores 1 \
        --p-front-f AATGATACGGACCACCGAGATCTACAC \
        --p-front-r CAAGCAGAAGACGGCATACGAGAT \
        --p-error-rate 0.1 \
        --o-trimmed-sequences demux-paired-end-trimmed.qza \

# PCR primers removal 




# DADA2

#qiime dada2 denoise-single \
#  --i-demultiplexed-seqs demux-paired-end.qza \
#  --p-trim-left 0 \
#  --p-trunc-len 240 \
#  --o-representative-sequences rep-seqs-dada2.qza \
#  --o-table table-dada2.qza \
#  --o-denoising-stats stats-dada2.qza

qiime dada2 denoise-paired --i-demultiplexed-seqs demux-paired-end.qza \
                           --p-trunc-len-f 0 \
                           --p-trunc-len-r 240 \ #Quality Score decreases from 240pb for the reverse reads
                           --p-max-ee 2.0 \ #default value : all the reads with number of exepcted errors higher than 2.0 will be discarded
                           --p-trunc-q 10 \ #reads are truncated at the first instance of a quality score less than or equal to 30
                           --p-n-reads-learn 1000000 \ #default value : it's the number of read to use during the training of error model 
                           --p-n-threads ${NSLOTS} \ #it uses all the available cores
                           --p-chimera-method consensus\ #default value : chimeras are detected in samples individually, and sequences found chimeric in a sufficient fraction of samples are removed
                           --o-representative-sequences rep-seq-dada2.qza \
                           --o-table table-dada2.qza \
                           --o-denoising-stats stats-dada2.qza \
                           --verbose



# Build ASV phylogeny with FastTree

mkdir phylogeny

qiime phylogeny align-to-tree-mafft-fasttree \
  --p-n-threads 0\
  --i-sequences rep-seqs-dada2.qza \
  --o-alignment phylogeny/aligned-rep-seqs.qza \
  --o-masked-alignment phylogeny/masked-aligned-rep-seqs.qza \
  --o-tree phylogeny/unrooted-tree.qza \
  --o-rooted-tree phylogeny/rooted-tree.qza \
  --verbose

# Assign Taxonomy to ASV

mkdir taxonomy

qiime feature-classifier classify-sklearn \
  --i-classifier /taxonomy/silva-132-99-nb-classifier.qza \
  --i-reads rep-seqs-dada2.qza  \
  --o-classification taxonomy/taxonomy.qza \
  --p-n-jobs default \
  --verbose

qiime metadata tabulate \
  --m-input-file taxonomy/taxonomy.qza  \
  --o-visualization taxonomy/taxonomy.qzv




# Export ASV table in biom format with taxonomy so we can import it easely in phyloseq/R
# see https://forum.qiime2.org/t/exporting-and-modifying-biom-tables-e-g-adding-taxonomy-annotations/3630
