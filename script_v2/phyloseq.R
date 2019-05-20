# SCHRIEKE Hans 
# phyloseq
# input : decontamed otu_table, asv_table, metadata 
# output : rarefaction curve, distribution plot, richness plot, NMDS, PCoA, Adonis



#--------------------------------------------------------------------------------------------#
#--------------------------------------PATH SETTING------------------------------------------#
#--------------------------------------------------------------------------------------------#

path = "D:/stage/data/runs_new"
setwd(path)

dir.create("plot") # folder for plot
path2 <- "D:/stage/data/runs_new/plot"




#--------------------------------------------------------------------------------------------#
#-------------------------------------LOAD PACKAGES------------------------------------------#
#--------------------------------------------------------------------------------------------#

library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library("DESeq2")
library("microbiome")
library("hrbrthemes")
#library("dplyr")

scripts <- c("graphical_methods.R",
             "tree_methods.R",
             "plot_merged_trees.R",
             "specificity_methods.R",
             "ternary_plot.R",
             "richness.R",
             "edgePCA.R",
             "copy_number_correction.R",
             "import_frogs.R",
             "prevalence.R",
             "compute_niche.R",
             "NCM_fit.R")
urls <-
  paste0("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/",
         scripts)

for (url in urls) {
  source(url)
}

url2 <- "https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/graphical_methods.R"
source(url2)





#--------------------------------------------------------------------------------------------#
#--------------------------------------IMPORT DATA-------------------------------------------#
#--------------------------------------------------------------------------------------------#


count_tab <- read.table("otu_decontam.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")

tax_tab <- as.matrix(read.table("tax_decontam.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

metadata <- read.table("metadata_decontam.tsv", header=T, row.names=1,
                              check.names=F, sep="\t")

dim(count_tab)
dim(metadata)
dim(tax_tab)

setwd(path2)

OTU <- otu_table(count_tab, taxa_are_rows=T)
TAX <- tax_table(tax_tab)
SAM <- sample_data(metadata)
ps <- phyloseq(OTU, TAX, SAM)







#--------------------------------------------------------------------------------------------#
#-----------------------------------RAREFACTION CURVE----------------------------------------#
#--------------------------------------------------------------------------------------------#

# method 1
#rarecurve(t(count_tab), step=100, lwd=2, col=metadata$color_run, ylab="ASVs", label=F)
#abline(v=(min(rowSums(t(count_tab)))))


# method 2
pdf("rarecurve.pdf")
#jpeg("rarecurve.jpg")
ggrare(ps, step = 2000, label = NULL, color = "run", plot = TRUE, parallel = FALSE, se = TRUE)
dev.off()



### NOTE ###
# Sample effort is good for all the samples 






#--------------------------------------------------------------------------------------------#
#--------------------------------------DISTRIBUTION------------------------------------------#
#--------------------------------------------------------------------------------------------#

ps.F <- subset_samples(ps, dna_from == "Full_body") # full organism 
ps.I <- subset_samples(ps, dna_from == "Intestine") # intestine
ps.SG <- subset_samples(ps, dna_from == "Salivary_gland") # salivary gland
ps.O <- subset_samples(ps, dna_from == "Ovary") # ovaries
ps.P <- subset_samples(ps, dna_from == "Organs_pull") # organs pull
ps.B <- subset_samples(ps, dna_from == "Blank") # blanks
ps.wB <- subset_samples(ps, dna_from != "Blank") # without blanks



readsumsdf = data.frame(nreads = sort(taxa_sums(ps.wB), TRUE), sorted = 1:ntaxa(ps.wB), 
                        type = "OTUs")

readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps.wB), TRUE), sorted = 1:nsamples(ps.wB), 
                                          type = "Samples"))



pdf("distribution.pdf")
#jpeg("distribution.jpg")
ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity") + ggtitle("Total number of reads before Preprocessing (Sans Témoins)") + scale_y_log10() +
  facet_wrap(~type, ncol = 1, scales = "free") #+ scale_y_log10()
dev.off()






#--------------------------------------------------------------------------------------------#
#---------------------------RICHNESS AND DIVERSITY ESTIMATES---------------------------------#
#--------------------------------------------------------------------------------------------#


# x based on samples
pdf("plot_richness_boxplot.pdf", width=20, height=10)
#jpeg("plot_richness.jpg", width=20, height =10)
plot_richness(ps, x= "run", color = "species", measures=c("Observed", "Chao1", "Shannon")) + 
  geom_boxplot() +
  theme(legend.title = element_blank())
dev.off()

# x based on location
pdf("plot_richness_location_boxplot.pdf", width=20, height=10)
#jpeg("plot_richness_location.jpg", width=20, height =10)
plot_richness(ps, x="location", color="run", measures=c("Observed", "Chao1", "Chao2", "Shannon")) + 
  geom_boxplot() +
  theme(legend.title = element_blank())
dev.off()


# richness estimation
results <- estimate_richness(ps, measures=c("Observed", "Chao1", "Chao2", "Shannon"))
summary(results)



### NOTE ###

# Chao1 is based on abundance -> Sest = Sobs + F2 / 2G 
# Sest = number of species that we want to know 
# Sobs = number of different species observed in a sample
# F = number of singletons (the number of species with only a single occurrence in the sample)
# G = the number of doubletons (the number of species with exactly two occurrences in the sample)
# On these plots, we can observe that "Observed" and "Chao1" are the same
# Maybe there is no singleton or doubleton in the sample

# Shannon is based on presence/absence -> H' = - ?? ((Ni / N) * log2 (Ni / N))
# H' = shannon index
# Ni = number of individual from a i species
# N = total number of individuals 
# H' = 0 when all individuals are from the same species or when each species is represented by a unique individual 
# H' is max when all individuals are evenly ditributed over all species 










#--------------------------------------------------------------------------------------------#
#-------------------------------------NORMALIZATION------------------------------------------#
#--------------------------------------------------------------------------------------------#


# Min/Max
otu_norm <- (count_tab -min(count_tab))/(max(count_tab)-min(count_tab))


# Log
otu_log <- log(count_tab +1)


# %

otu_percent <- cbind(0)
for (i in 1:ncol(count_tab)){
  otu_percent <- cbind(otu_percent,(count_tab[i]/colSums(count_tab[i]))*100)
  #print(otu_col_percent)
  #otu_percent2 <- cbind(otu_col_percent, otu_col_percent)
  #otu_percent2 <- otu_percent2
}

otu_percent <- otu_percent[,-1]

# (count_tab[1]/colSums(count_tab[1]))*100
# otu_percent <- ((count_tab/colSums(count_tab))*100)
# otu_percent <- as.matrix(otu_percent)


# DESeq2
rownames(metadata)
colnames(count_tab)

otu2 <- as.matrix(count_tab)
otu2 <- (otu2+1) # allows to remove the zero --> DESeq doesn't work with zero 

ncol(otu2) # check up
nrow(metadata) # check up

#otu_deseq <- counts(dds, normalized=TRUE) # normalized otu with deseq

count_tab2 <- (count_tab+1)

deseq_counts <- DESeqDataSetFromMatrix(otu2, colData = metadata, design = ~dna_from)
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

otu_deseq <- assay(deseq_counts_vst) # normalized otu_table with deseq







#--------------------------------------------------------------------------------------------#
#----------------------------------TAXONOMIC SUMMARIES---------------------------------------#
#--------------------------------------------------------------------------------------------#




OTU_percent <- otu_table(otu_percent, taxa_are_rows = TRUE)
ps_percent <- phyloseq(OTU_percent, TAX, SAM)


top20 <- names(sort(taxa_sums(ps_percent), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps_percent, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

x11()

pdf("taxa_plot_2.pdf", width=20, height =10)
plot_composition(ps.top20,
                      taxonomic.level = "Family",
                      sample.sort = "sample",
                      x.label = "sample") +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_percent() +
  labs(x = "Samples", y = "Relative abundance (%)",
       title = "Relative abundance data",
       caption = "Caption text.") + 
  theme_ipsum(grid="Y")
dev.off()

print(p)
dev.off()



p.phy <- plot_composition(ps.top20, sample.sort = NULL, otu.sort = NULL,
                          x.label = "sample", plot.type = "barplot", verbose = FALSE)

print(p.phy + scale_fill_brewer(palette = "Paired") + theme_bw())


plot_composition(ps.top20, x.label="Family", plot.type="barplot")

pdf("taxa_plot3.pdf")
plot_composition(ps.top20, "Phylum", "Proteobacteria", "Family", 20, fill = "Family") +
  theme_ipsum(grid="Y") +
  theme(axis.text.x = element_text(color="#993333", 
                                   size=5),
        axis.text.y = element_text(color="#993333", 
                                   size=14))
  #geom_bar(stat="identity", color="black")
dev.off()

pdf("taxa_plot.pdf", width=20, height =10)
#jpeg("taxa_plot.jpg", width=20, height =10)
plot_bar(ps.top20, x="sample", fill="Family")
dev.off()


# plot according to location

pdf("taxa_plot_location.pdf", width=20, height =10)
#jpeg("taxa_plot_location.jpg", width=20, height =10)
plot_bar(ps.top20, x="sample", fill="Family") + facet_wrap(~location, scales="free_x")
dev.off()

# plot according to condition of DNA extraction
pdf("taxa_plot_dnafrom.pdf", width=20, height =10)
#jpeg("taxa_plot_dnafrom.jpg", width=20, height =10)
plot_bar(ps.top20, x="sample", fill="Family") + facet_wrap(~dna_from, scales="free_x")
dev.off()


### NOTE ###
# We can see that the taxonomic composition changes with the location and the condition of DNA extraction
# It's interesting!

# Taxa_plot_location :
# Anaplasmataceae family is present everywhere except in Labo Tetracycline
# Acetobacteraceae is present close everywhere but majoritary in Guadeloupe 
# Enterobacteriaceae is present only in Lavar, Bosc and Labo Tetracycline

# Taxa_plot_dnafrom :
# Anaplasmataceae is present everywhere 
# Ovary and Salivary gland are majoritary composed by Anaplasmataceae
# Full_body and Intestine present the most important diversity 











#--------------------------------------------------------------------------------------------#
#--------------------------------------ORDINATION--------------------------------------------#
#--------------------------------------------------------------------------------------------#

OTU_log <- otu_table(otu_log, taxa_are_rows=T)
ps_log <- phyloseq(OTU_log, SAM)

OTU_percent <- otu_table(otu_percent, taxa_are_rows=T)
ps_percent <- phyloseq(OTU_percent, SAM)


OTU_deseq <- otu_table(otu_deseq, taxa_are_rows=T)
ps_deseq <- phyloseq(OTU_deseq, SAM)


# PCoA log (euclidean)
pdf("PCoA_log_species.pdf")
#jpeg("PCoA_log_species.jpg")
plot_ordination(ps_log, ordinate(ps_log, method ="MDS", distance = "euclidean"), color = "species") +
  geom_point(size = 3) +
  ggtitle("PCoA (log normalization)")
dev.off()


# PCoA deseq (euclidean)
pdf("PCoA_deseq_species.pdf")
#jpeg("PCoA_deseq_species.jpg")
plot_ordination(ps_deseq, ordinate(ps_deseq, method ="MDS", distance = "euclidean"), color = "species") +
  geom_point(size = 3) +
  ggtitle("PCoA (deseq normalization)")
dev.off()


# PCoA percent (bray)
pdf("PCoA_percent_species_corr.pdf")
#jpeg("PCoA_percent_species.jpg")
plot_ordination(ps_percent, ordinate(ps_percent, method ="MDS", distance = "bray"), color = "species") +
  geom_point(size = 3) +
  ggtitle("PCoA (percent normalization)")
dev.off()



### NOTE ###
# We can see an artefact pattern on log and deseq PCoA
# It's a possible effect of 0 in otu_table but log doesn't change the result, so...
# PCoA is maybe not applicable with our data







#--------------------------------------------------------------------------------------------#
#--------------------------------NMMDS (euclidean distance)----------------------------------#
#--------------------------------------------------------------------------------------------#

# phyloseq object creation for each condition
ps.OI <- subset_samples(ps, dna_from != "Full_body")
ps.OI <- subset_samples(ps, dna_from != "Salivary_gland")
ps.OI <- subset_samples(ps, dna_from != "Organs_pull")
ps.OI <- subset_samples(ps, dna_from != "Blank") # ovaries + intestine

#otu_deseq[otu_deseq < 0.0] <- 0.0

# otu_table transformation and ordination
ps.prop.euc <- transform_sample_counts(ps, function(count_tab) count_tab/sum(count_tab))
ord.nmds.euc <- ordinate(ps.prop.euc, method="NMDS", distance="euclidean")

ps.propF.euc <- transform_sample_counts(ps.F, function(count_tab) count_tab/sum(count_tab))
ord.nmds.eucF <- ordinate(ps.propF.euc, method="NMDS", distance="euclidean")

ps.propI.euc <- transform_sample_counts(ps.I, function(count_tab) count_tab/sum(count_tab))
ord.nmds.eucI <- ordinate(ps.propI.euc, method="NMDS", distance="euclidean")

ps.propO.euc <- transform_sample_counts(ps.O, function(count_tab) count_tab/sum(count_tab))
ord.nmds.eucO <- ordinate(ps.propO.euc, method="NMDS", distance="euclidean")

ps.propOI.euc <- transform_sample_counts(ps.OI, function(count_tab) count_tab/sum(count_tab))
ord.nmds.eucOI <- ordinate(ps.propOI.euc, method="NMDS", distance="euclidean")


# NMDS plots (euclidean)

pdf("NMDS_euc_all.pdf")
#jpeg("NMDS_euc_all.jpg")
plot_ordination(ps.prop.euc, ord.nmds.euc, color="species", title="Euclidean NMDS (all)", label="sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_euc_fullbody.pdf")
#jpeg("NMDS_euc_fullbody.jpg")
plot_ordination(ps.propF.euc, ord.nmds.eucF, color="species", title="Euclidean NMDS (full body)", label="sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_euc_intestine.pdf")
#jpeg("NMDS_euc_intestine.jpg")
plot_ordination(ps.propI.euc, ord.nmds.eucI, color="species", shape="species", title="Euclidean NMDS (intestine)", label="sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_euc_ovary.pdf")
#jpeg("NMDS_euc_ovary.jpg")
plot_ordination(ps.propO.euc, ord.nmds.eucO, color="species", title="Euclidean NMDS (ovary)", label="sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_euc_ovary_intestine.pdf")
#jpeg("NMDS_euc_ovary_intestine.jpg")
plot_ordination(ps.propOI.euc, ord.nmds.eucOI, color="species", title="Euclidean NMDS (ovary+intestine)", label="sample")+
  geom_point(size = 3)
dev.off()







#--------------------------------------------------------------------------------------------#
#-----------------------------------NMMDS (bray distance)------------------------------------#
#--------------------------------------------------------------------------------------------#

# otu_table transformation and ordination
ps.prop <- transform_sample_counts(ps, function(count_tab) count_tab/sum(count_tab))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

ps.propF <- transform_sample_counts(ps.F, function(count_tab) count_tab/sum(count_tab))
ord.nmds.brayF <- ordinate(ps.propF, method="NMDS", distance="bray")

ps.propI <- transform_sample_counts(ps.I, function(count_tab) count_tab/sum(count_tab))
ord.nmds.brayI <- ordinate(ps.propI, method="NMDS", distance="bray")

ps.propO <- transform_sample_counts(ps.O, function(count_tab) count_tab/sum(count_tab))
ord.nmds.brayO <- ordinate(ps.propO, method="NMDS", distance="bray")

ps.propOI <- transform_sample_counts(ps.OI, function(count_tab) count_tab/sum(count_tab))
ord.nmds.brayOI <- ordinate(ps.propOI, method="NMDS", distance="bray")


# NMDS plots (bray)

pdf("NMDS_bray_all.pdf")
#jpeg("NMDS_bray_all.jpg")
plot_ordination(ps.prop, ord.nmds.bray, color="species", title="Bray NMDS (all)", label="sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_bray_fullbody.pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(ps.propF, ord.nmds.brayF, color="species", title="Bray NMDS (full body)", label="sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_bray_intestine.pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(ps.propI, ord.nmds.brayI, color="species", shape="species", title="Bray NMDS (intestine)", label="sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_bray_ovary.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(ps.propO, ord.nmds.brayO, color="species", title="Bray NMDS (ovary)", label="sample") +
  geom_point(size = 3)
dev.off()

#pdf("NMDS_bray_ovary_intestine.pdf")
jpeg("NMDS_bray_ovary_intestine.jpg")
plot_ordination(ps.propOI, ord.nmds.brayOI, color="species", title="Bray NMDS (ovary+intestine)", label="sample") +
  geom_point(size = 3)
dev.off()




### NOTE ### 
# Bray-Curtis distance doesn't seem applicable for our data, contrary to euclidean distance
# Why ? I don't really know but I think each data depends on one type of distance for NMDS
# I think we need to choose a distance and show that our results are robust with this distance 







#--------------------------------------------------------------------------------------------#
#-------------------------------------------ADONIS-------------------------------------------#
#--------------------------------------------------------------------------------------------#


adonis(vegdist(t(otu_table(ps_deseq)), method = "euclidean") ~ location,
       data=as(sample_data(ps_deseq), "data.frame"), permutation = 9999)

### NOTE : ###
# Adonis is used to test the difference of abundance between groups

#           Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
#location    5     88105 17621.0  8.2206 0.1624  1e-04 ***
#  Residuals 212    454427  2143.5         0.8376           
#Total     217    542531                 1.0000           
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# R2 = 0.1624 means that location could explain the distance variation 
# p is very significative, that means the observed difference is unlikely to be due to chance




adonis(vegdist(t(otu_table(ps_deseq)), method = "euclidean") ~ dna_from,
       data=as(sample_data(ps_deseq), "data.frame"), permutation = 9999)


adonis(vegdist(t(otu_table(ps_deseq)), method = "euclidean") ~ date,
       data=as(sample_data(ps_deseq), "data.frame"), permutation = 9999)






#--------------------------------------------------------------------------------------------#
#---------------------------------------Save session-----------------------------------------#
#--------------------------------------------------------------------------------------------#

setwd(path)
install.packages("session", repos = "http://cran.us.r-project.org")
library(session); packageVersion("session")
save.session("phyloseq.Rda")


