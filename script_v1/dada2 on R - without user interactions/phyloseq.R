# Phyloseq 
# SCHRIEKE Hans

# SET WORK DIRECTORY 
path <- "/gs7k1/home/schrieke/Fastq/tables"
setwd(path)


# DATA IMPORT : OTU table, taxa table and data frame (metadata)
otu <- read.csv("seqtabnochim.csv", sep=";", dec=",")
taxa <- read.csv("taxa.csv", sep=";", dec=",")
metadata <- read.csv("metadata_run1.csv", sep=",", row.names = 1)

metadata[metadata=="N"] <- NA # negative controls
metadata <- na.omit(metadata) # remove N/A lines from metadata

otu <- as.matrix(otu)
taxa <- as.matrix(taxa)
metadata <- as.data.frame(metadata)


# Normalization (min-max / log / DESeq2)
otu_norm <- (otu -min(otu))/(max(otu)-min(otu))
otu_norm <- as.matrix(otu_norm)

otu_log <- log(otu +1)
otu_log <- as.matrix(otu_log)


# Normalization (%)

otu_percent <- (otu*100/colSums(otu))
as.matrix(otu_percent)


# DEseq2
rownames(metadata)
colnames(otu)

otu2 <- otu[,c(1:96)]
colnames(otu2)


metadata2 <- metadata[order(metadata$SAMPLE),]
metadata2 <- metadata2[,c(2:9)]

sort(metadata2)
ncol(otu2)
nrow(metadata)

otu2 <- sort(otu2)
otu2<- otu2+1 #allow to remove the zero --> DESeq doesn't work with zero 

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=otu2, colData=metadata2, design=~DATE+LOCATION+ORGANISM)
dds <- DESeq(dds)
ts <- counts(dds, normalized=TRUE)

setwd(path2)
dir.create("DESeq")
path3 <- "D:/data/phyloseq_plot/DESeq"
setwd(path3)

pdf("dispersion-plot-DESeq.pdf")
plotDispEsts(dds, main="Dispersion plot")
dev.off



# PHYLOSEQ INSTALLATION
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')


# PHYLOSEQ AND GGPLOT 2 LOADING
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw()) #set ggplot2 graphic theme 


# PHYLOSEQ OBJECT CREATION
OTU0 = otu_table(otu, taxa_are_rows =TRUE) # without normalization
OTU = otu_table(otu, taxa_are_rows = TRUE)
OTU2 = otu_table(otu_log, taxa_are_rows = TRUE)
OTU3 = otu_table(otu_percent, taxa_are_rows = TRUE)

TAX = tax_table(taxa)
SAM = sample_data(metadata)

ps0 <- phyloseq(OTU0, TAX, SAM) # without normalization
ps <- phyloseq(OTU, TAX, SAM) 
ps2 <- phyloseq(OTU2, TAX, SAM)
ps3 <- phyloseq(OTU3, TAX, SAM)


# ALPHA DIVERSITY PLOT 
dir.create("phyloseq_plot")
path2 <- "/gs7k1/home/schrieke/Fastq/tables/phyloseq_plot"
setwd(path2)

dir.create("alpha-diversity")
path3 <- "/gs7k1/home/schrieke/Fastq/tables/phyloseq_plot/alpha-diversity"
setwd(path3)

pdf("alpha_diversity.pdf")
plot_richness(ps, measures=c("Shannon", "Simpson"), color="LOCATION")
dev.off()

setwd(path2)


# BRAY CURTIS NMDS with OTU_NORM
ps.prop <- transform_sample_counts(ps, function(otu_norm) otu_norm/sum(otu_norm))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

dir.create("NMDS-norm")
path3 <- "/gs7k1/home/schrieke/Fastq/tables/phyloseq_plot/NMDS-norm"
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
                      
                                   
# BRAY CURTIS NMDS with OTU_LOG
ps2.prop <- transform_sample_counts(ps2, function(otu_log) otu_log/sum(otu_log))
ord.nmds.bray <- ordinate(ps2.prop, method="NMDS", distance="bray")

dir.create("NMDS-log")
path3 <- "/gs7k1/home/schrieke/Fastq/tables/phyloseq_plot/NMDS-log"
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

setwd(path2) 

                                    
# BRAY CURTIS NMDS with OTU_PERCENT
ps3.prop <- transform_sample_counts(ps3, function(otu_percent) otu_log/sum(otu_percent))
ord.nmds.bray <- ordinate(ps3.prop, method="NMDS", distance="bray")

dir.create("NMDS-percent")
path3 <- "/gs7k1/home/schrieke/Fastq/tables/phyloseq_plot/NMDS-percent"
setwd(path3)
                                    
pdf("NMDS-location-percent.pdf")
plot_ordination(ps2.prop, ord.nmds.bray, color="LOCATION", title="Bray NMDS (location)", label="SAMPLE")
dev.off()

pdf("NMDS-date-percent.pdf")
plot_ordination(ps2.prop, ord.nmds.bray, color="DATE", title="Bray NMDS (date)", label="SAMPLE")
dev.off()

pdf("NMDS-organism-percent.pdf")
plot_ordination(ps2.prop, ord.nmds.bray, color="ORGANISM", title="Bray NMDS (organism)", label="SAMPLE")
dev.off()

setwd(path) 
        
                                    
# ADONIS 
library(vegan)
library(dist)
dist.jac <- vegdist(otu_norm, method="jaccard", binary=TRUE)
dist.jac <- as.matrix(dist.jac)
metadata <- as.data.frame(metadata)
as.recursive(ps)
is.atomic(ps)
adonis(formula = dist.jac ~LOCATION, data=metadata2, permutation = 9999)
                                    

# RAREFACTION CURVE

source('functions.R') # import amp_rarecurve, plot_composition, ggrare

library(dplyr)
library(vegan)
library(magrittr)

ps.E <- subset_samples(ps, ORGANISM == "E") # full organism 
ps.I <- subset_samples(ps, ORGANISM == "I") # intestine
ps.GS <- subset_samples(ps, ORGANISM == "GS") # salivary gland
ps.O <- subset_samples(ps, ORGANISM == "O") # ovaries
ps.P <- subset_samples(ps, ORGANISM == "P") # organs pull
ps.T <- subset_samples(ps, LOCATION == "T") # blanks
ps.wT <- subset_samples(ps, LOCATION != "T") # without blanks



amp_rarecurve(subset_samples(ps, LOCATION =="N"),
              step=100,
              label = T,
              legend.position = "bottomright",
              legend = T)



readsumsdf = data.frame(nreads = sort(taxa_sums(ps.wT), TRUE), sorted = 1:ntaxa(ps.wT), 
                                      type = "OTUs")

readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps.wT), TRUE), sorted = 1:nsamples(ps.wT), 
                                          type = "Samples"))

ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity") + ggtitle("Total number of reads before Preprocessing (Sans Témoins)") + scale_y_log10() +
  facet_wrap(~type, ncol = 1, scales = "free") #+ scale_y_log10()

graph2ppt(file=paste0("Preanalysis",".ppt"),append=T,width=9,aspectr=sqrt(2))




# PLOT COMPOSITION 

library(scales)
library(reshape2)
library(ggplot2)

subset_samples(ps, LOCATION == "N") %>%
  plot_composition("Kingdom", "Bacteria", "Species", numberOfTaxa =
                     1000, fill = "Species")

subset_samples(ps, LOCATION == "N") %>%
  plot_richness(ps, measures=c("Observed","Shannon","ACE"))


                                    
# GGRARE

ggrare(ps, step = 7, label = NULL, color = NULL,
       plot = TRUE, parallel = FALSE, se = TRUE)
                                    
                                    
                                    
# SAVE SESSION
install.packages("session", repos = "http://cran.us.r-project.org")
library(session)
save.session("phyloseq.Rda")

