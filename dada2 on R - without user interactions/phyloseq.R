# Phyloseq 
# SCHRIEKE Hans

# Set work directory 
path <- "/gs7k1/home/schrieke/Fastq/tables"
setwd(path)

# Data import : OTU table, taxa table and data frame (metadata)
OTU <- read.csv("seqtab.csv", sep=";", dec=",")
taxa <- read.csv("taxa.csv", sep=";", dec=",")
metada <- read.csv("metadata.csv", sep=";", dec=",")

OTU_table <- as.matrix(OTU)
taxa_table <- as.matrix(taxa)

# Phyloseq installation
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')

# Phyloseq and ggplot2 loading 
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw()) #set ggplot2 graphic theme 

ps <- phyloseq(otu_table(OTU_table.nochim, taxa_are_rows=FALSE), 
               sample_data(metada), 
               tax_table(taxa_table))
