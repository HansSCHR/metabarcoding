# SCHRIEKE Hans
# DECONTAM
# input : tables from dada2 step
# output : phyloseq objects with decontam tables 



#--------------------------------------------------------------------------------------------#
#--------------------------------------PATH SETTING------------------------------------------#
#--------------------------------------------------------------------------------------------#

path = "D:/stage/data/runs_new2/"
setwd(path)




#--------------------------------------------------------------------------------------------#
#-------------------------------------LOAD PACKAGES------------------------------------------#
#--------------------------------------------------------------------------------------------#

library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
library(DESeq2); packageVersion("DESeq2")

theme_set(theme_gray()) #set ggplot2 graphic theme 





#--------------------------------------------------------------------------------------------#
#---------------------------------------IMPORT DATA------------------------------------------#
#--------------------------------------------------------------------------------------------#

otu <- read.csv("seqtabnochimcor.csv", sep=";", dec=",")
taxa <- read.csv("taxafinal.csv", sep=";", dec=",")
metadata <- read.csv("metadata_runs.csv", sep=",", row.names = 1)



#--------------------------------------------------------------------------------------------#
#--------------------------------EXTRACT GOOD STANDARDS--------------------------------------#
#--------------------------------------------------------------------------------------------#

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- rownames(otu)
asv_headers <- vector(dim(otu)[1], mode="character")
asv_headers2 <- vector(dim(otu)[1], mode="character") 

for (i in 1:dim(otu)[1]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
  asv_headers2[i] <- paste("ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
asv_fasta2 <- cbind(asv_headers2, asv_seqs)
#write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- otu
row.names(asv_tab) <- sub(">", "", asv_headers)

# tax table:
asv_tax <- as.matrix(taxa)
row.names(asv_tax) <- sub(">", "", asv_headers)



#--------------------------------------------------------------------------------------------#
#------------------------------------PHYLOSEQ OBJECT-----------------------------------------#
#--------------------------------------------------------------------------------------------#

OTU = otu_table(asv_tab, taxa_are_rows =TRUE)
TAX = tax_table(asv_tax)
SAM = sample_data(metadata)


ps <- phyloseq(OTU, TAX, SAM)
ps <- subset_samples(ps, Species!="CuT" & Species!="CuP" & Species!="CuN" & Species!="CuG")
ps <- prune_samples(sample_sums(ps) >= 1, ps) # remove counts = 0




#--------------------------------------------------------------------------------------------#
#----------------------------------------DECONTAM--------------------------------------------#
#--------------------------------------------------------------------------------------------#


######################################## PREPROCESS #########################################

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

pdf("LibrarySize.pdf")
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Dna)) + geom_point()
dev.off()

as.numeric(get_variable(ps, "Dna"))
get_variable(ps, "Dna")
sample_data(ps)

sample_data(ps)$Dna <- as.numeric(get_variable(ps, "Dna"))



######################  IDENTIFY CONTAMINANTS - PREVALENCE #################################

sample_data(ps)$is.neg <- sample_data(ps)$Control == "Control sample"

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant) # 80 contaminants

contam_asvs_prev05 <- row.names(contamdf.prev05[contamdf.prev05$contaminant == TRUE, ])
contam <-asv_tax[row.names(asv_tax) %in% contam_asvs_prev05, ]


ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Control == "Control sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Control == "True sample", ps.pa)


# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev05$contaminant)

pdf("Prevalence_Decontam.pdf")
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
dev.off()

# write table of contaminant asv 
write.table(contam, "asv_contam_prev05.tsv",
            sep="\t", quote=F, col.names=NA)





#--------------------------------------------------------------------------------------------#
#---------------------------------------NEW TABLES-------------------------------------------#
#--------------------------------------------------------------------------------------------#

# new otu table 
asv_tab_prev05 <- asv_tab[!row.names(asv_tab) %in% contam_asvs_prev05, ]

# new tax table
asv_tax_prev05 <- asv_tax[!row.names(asv_tax) %in% contam_asvs_prev05, ]

# new fasta  
contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs_prev05))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_prev05 <- asv_fasta[- dont_want]

# write tables
write.table(asv_tab_prev05, "asv_tab_prev05.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_prev05, "asv_tax_prev05.tsv",
            sep="\t", quote=F, col.names=NA)
write(asv_fasta_prev05, "fasta_prev05.fa")


# extract new metadata as dataframe from phyloseq 
meta <- function(x) { # https://rdrr.io/github/microbiome/microbiome/src/R/meta.R
  df <- as(sample_data(x), "data.frame")
  rownames(df) <- sample_names(x)
  df
}
metadata_prev05 <- meta(ps)

# write decontamed metadata
write.table(metadata_prev05, "metadata_decontam.tsv", sep="\t", quote=F, col.names=NA)




#--------------------------------------------------------------------------------------------#
#-----------------------------NEW DECONTAM PHYLOSEQ OBJECT-----------------------------------#
#--------------------------------------------------------------------------------------------#

OTU_decontam <- otu_table(asv_tab_prev05, taxa_are_rows = TRUE)
TAX_decontam <- tax_table(asv_tax_prev05)
SAM_decontam <- sample_data(metadata_prev05)

ps_decontam <- phyloseq(OTU_decontam, TAX_decontam, SAM_decontam)
ps_decontam <- subset_samples(ps_decontam, Control!="Control sample") # remove blanks
ps_decontam <- subset_taxa(ps_decontam, Kingdom == "Bacteria") # select only Bacteria
ps_decontam <- prune_samples(sample_sums(ps_decontam) >= 1, ps_decontam) # remove counts = 0 in otu_table
ps_decontam <- prune_taxa(taxa_sums(ps_decontam) >= 1, ps_decontam)# remove counts = 0 in tax_table




#--------------------------------------------------------------------------------------------#
#-------------------REMOVING REMAINING CHLOROPLAST AND MITOCHONDRIA--------------------------#
#--------------------------------------------------------------------------------------------#
ps_decontam2 <- subset_taxa(ps_decontam, Family!="Mitochondria" & Order!="Chloroplast")
ps_decontam2 <- prune_samples(sample_sums(ps_decontam2) >= 1, ps_decontam2)
ps_decontam2 <- prune_taxa(taxa_sums(ps_decontam2) >= 1, ps_decontam2)



#--------------------------------------------------------------------------------------------#
#-------------------------------------NORMALIZATION------------------------------------------#
#--------------------------------------------------------------------------------------------#

# % 
ps_percent <- transform_sample_counts( ps_decontam2, function(x) x/sum(x)*100 )



#--------------------------------------------------------------------------------------------#
#-----------------------------------EXPLORING WOLBACHIA--------------------------------------#
#--------------------------------------------------------------------------------------------#
ps_wolbachia <- subset_taxa(ps_percent, Genus == "Wolbachia")
ps_wolbachia <- prune_samples(sample_sums(ps_wolbachia) >= 1, ps_wolbachia)
ps_wolbachia <- prune_samples(taxa_sums(ps_wolbachia) >= 1, ps_wolbachia)


ps_proteo <- subset_taxa(ps_percent, Phylum=="Proteobacteria")
ps_proteo <- prune_samples(sample_sums(ps_proteo) >= 1, ps_proteo)
ps_wolbachia <- prune_samples(taxa_sums(ps_proteo) >= 1, ps_proteo)

#--------------------------------------------------------------------------------------------#
#---------------------------------------Save objects-----------------------------------------#
#--------------------------------------------------------------------------------------------#
save(ps, ps_decontam, ps_decontam2, ps_percent, ps_wolbachia, ps_proteo, file = "objects.RData")




#--------------------------------------------------------------------------------------------#
#---------------------------------------Save session-----------------------------------------#
#--------------------------------------------------------------------------------------------#

#install.packages("session", repos = "http://cran.us.r-project.org")
library(session); packageVersion("session")
save.session("decontam2.Rda")
#load("phyloseq.Rda")

