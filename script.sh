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

# Dada2 
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux-paired-end.qza \
  --p-trim-left 0 \
  --p-trunc-len 240 \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-table table-dada2.qza \
  --o-denoising-stats stats-dada2.qza



