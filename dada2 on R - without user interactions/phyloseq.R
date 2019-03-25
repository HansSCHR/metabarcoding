# Phyloseq 
# SCHRIEKE Hans

# Set work directory 
path <- "/gs7k1/home/schrieke/Fastq/tables"
setwd(path)


# Data import : OTU table, taxa table and data frame (metadata)
otu <- read.csv("seqtabnochim.csv", sep=";", dec=",")
taxa <- read.csv("taxa.csv", sep=";", dec=",")
metadata <- read.csv("metadata_run1.csv", sep=",", row.names = 1)

metadata[metadata=="N"] <- NA # negative controls 

otu <- as.matrix(otu)
#print("otu table (matrix) : ")
#head(otu)

taxa <- as.matrix(taxa)
#print("taxa table (matrix) : ")
#head(taxa)

metadata <- as.data.frame(metadata)
#print ("metadata (data frame) : ")
#head(metadata)

# Phyloseq installation
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')

# Phyloseq and ggplot2 loading 
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw()) #set ggplot2 graphic theme 

# Phyloseq object creation
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(taxa)
SAM = sample_data(metadata)

ps <- phyloseq(OTU, TAX, SAM) 


# Alpha diversity plot 
dir.create("phyloseq_plot")
path2 <- "/gs7k1/home/schrieke/Fastq/tables/phyloseq_plot"
setwd(path2)

pdf("5-alpha_diversity.pdf")
plot_richness(ps, measures=c("Shannon", "Simpson"), color="LOCATION")
dev.off()


# Bray-Curtis NMDS
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

pdf("6-NMDS.pdf")
plot_ordination(ps.prop, ord.nmds.bray, color="LOCATION", title="Bray NMDS")
dev.off()

setwd(path) 

# Save the session
install.packages("session", repos = "http://cran.us.r-project.org")
library(session)
save.session("phyloseq.Rda")

