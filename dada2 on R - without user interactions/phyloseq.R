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
metadata <- na.omit(metadata) # remove N/A lines from metadata

otu <- as.matrix(otu)
taxa <- as.matrix(taxa)
metadata <- as.data.frame(metadata)


# Normalization (min-max and log)
otu_norm <- (otu -min(otu))/(max(otu)-min(otu))
otu_log <- log(otu +1)
otu_norm <- as.matrix(otu_norm)
otu_log <- as.matrix(otu_log)

# Phyloseq installation
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')


# Phyloseq and ggplot2 loading 
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw()) #set ggplot2 graphic theme 


# Phyloseq object creation
OTU = otu_table(otu, taxa_are_rows = TRUE)
OTU2 = otu_table(otu_log, taxa_are_rows = TRUE)
TAX = tax_table(taxa)
SAM = sample_data(metadata)

ps <- phyloseq(OTU, TAX, SAM) 
ps2 <- phyloseq(OTU2, TAX, SAM)


# Alpha diversity plot 
dir.create("phyloseq_plot")
path2 <- "/gs7k1/home/schrieke/Fastq/tables/phyloseq_plot"
setwd(path2)

dir.create("alpha-diversity")
path3 <- "D:/data/phyloseq_plot/alpha-diversity"
setwd(path3)

pdf("alpha_diversity.pdf")
plot_richness(ps, measures=c("Shannon", "Simpson"), color="LOCATION")
dev.off()

setwd(path2)

# Bray-Curtis NMDS with otu_norm
ps.prop <- transform_sample_counts(ps, function(otu_norm) otu_norm/sum(otu_norm))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

dir.create("NMDS-norm")
path3 <- "D:/data/phyloseq_plot/NMDS-norm"
setwd(path3)
                                   
pdf("NMDS-location-norm.pdf")
plot_ordination(ps.prop, ord.nmds.bray, color="LOCATION", title="Bray NMDS (location)", label="SAMPLE")
dev.off()

pdf("NMDS-date-norm.pdf")
plot_ordination(ps.prop, ord.nmds.bray, color="DATE", title="Bray NMDS (date)", label="SAMPLE")
dev.off()

pdf("NMDS-organism-norm.pdf")
plot_ordination(ps.prop, ord.nmds.bray, color="ORGANISM", title="Bray NMDS (organism)", label="SAMPLE")
dev.off()

setwd(path2)
                                   
# Bray-Curtis NMDS with otu_log
ps2.prop <- transform_sample_counts(ps2, function(otu_log) otu_log/sum(otu_log))
ord.nmds.bray <- ordinate(ps2.prop2, method="NMDS", distance="bray")

dir.create("NMDS-log")
path3 <- "D:/data/phyloseq_plot/NMDS-log"
setwd(path3)
                                    
pdf("NMDS-location-log.pdf")
plot_ordination(ps2.prop, ord.nmds.bray, color="LOCATION", title="Bray NMDS (location)", label="SAMPLE")
dev.off()

pdf("NMDS-date-log.pdf")
plot_ordination(ps2.prop, ord.nmds.bray, color="DATE", title="Bray NMDS (date)", label="SAMPLE")
dev.off()

pdf("NMDS-organism-log.pdf")
plot_ordination(ps2.prop, ord.nmds.bray, color="ORGANISM", title="Bray NMDS (organism)", label="SAMPLE")
dev.off()

setwd(path) 

                                    
# Adonis 
library(vegan)
library(dist)
dist.jac <- vegdist(otu_norm, method="jaccard", binary=TRUE)
adonis(formula = dist.jac, ps, data = metadata, permutation = 9999)
                                    
                                    
# Save the session
install.packages("session", repos = "http://cran.us.r-project.org")
library(session)
save.session("phyloseq.Rda")

