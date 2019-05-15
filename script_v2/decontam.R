# SCHRIEKE Hans 
# decontam 
# input : otu_table and tax_table from dada2, metadata file
# output : tables without contaminants 


#--------------------------------------------------------------------------------------------#
#--------------------------------------PATH SETTING------------------------------------------#
#--------------------------------------------------------------------------------------------#

path = "D:/stage/data/runs_new"
setwd(path)




#--------------------------------------------------------------------------------------------#
#-------------------------------------LOAD PACKAGES------------------------------------------#
#--------------------------------------------------------------------------------------------#

library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")

theme_set(theme_bw()) #set ggplot2 graphic theme 





#--------------------------------------------------------------------------------------------#
#---------------------------------------IMPORT DATA------------------------------------------#
#--------------------------------------------------------------------------------------------#

otu <- read.csv("seqtabnochimcor.csv", sep=";", dec=",")
taxa <- read.csv("taxafinal.csv", sep=";", dec=",")
metadata <- read.csv("metadata_runs.csv", sep=",", row.names = 1)




#--------------------------------------------------------------------------------------------#
#------------------------------REMOVE SAMPLES WITHOUT READS----------------------------------#
#--------------------------------------------------------------------------------------------#


metadata <- as.data.frame(t(metadata))
otu <- otu[, order(colnames(otu))]
metadata <- metadata[, order(colnames(metadata))]


for (i in 1:ncol(otu)){
  if (sum(otu[,i])==0){
    print(names(otu[i]))
  }
}

for (i in 1:ncol(otu)){
  if (sum(otu[,i])==0){
    print(names(otu[i]))
    otu <- otu[,-i]
    print("otu column removed")
    metadata <- metadata[,-i]
    print("metadata column removed")
  }

}
metadata <- t(metadata)

ncol(otu) # check
nrow(metadata) # check

otu <- as.matrix(t(otu))
taxa <- as.matrix(taxa)
metadata <- as.data.frame(metadata)




#--------------------------------------------------------------------------------------------#
#--------------------------------EXTRACT GOOD STANDARDS--------------------------------------#
#--------------------------------------------------------------------------------------------#

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(otu)
asv_headers <- vector(dim(otu)[2], mode="character")

for (i in 1:dim(otu)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
#write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(otu)
row.names(asv_tab) <- sub(">", "", asv_headers)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)





#--------------------------------------------------------------------------------------------#
#------------------------------------PHYLOSEQ OBJECT-----------------------------------------#
#--------------------------------------------------------------------------------------------#

OTU = otu_table(asv_tab, taxa_are_rows =TRUE)
TAX = tax_table(asv_tax)
SAM = sample_data(metadata)


ps <- phyloseq(OTU, TAX, SAM)
ps <- subset_samples(ps,species!="CuT")
ps <- subset_samples(ps, species!="CuP")
ps <- subset_samples(ps, species!="CuN")
ps <- subset_samples(ps, species!="CuG") # cullicoides removed 




#--------------------------------------------------------------------------------------------#
#----------------------------------------DECONTAM--------------------------------------------#
#--------------------------------------------------------------------------------------------#


######################################## PREPROCESS #########################################

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

pdf("LibrarySize.pdf")
ggplot(data=df, aes(x=Index, y=LibrarySize, color=true_or_blank)) + geom_point()
dev.off()

as.numeric(get_variable(ps, "dna_quant"))
get_variable(ps, "dna_quant")
sample_data(ps)

sample_data(ps)$dna_quant <- as.numeric(get_variable(ps, "dna_quant"))




######################  IDENTIFY CONTAMINANTS - FREQUENCY ###################################

contamdf.freq <- isContaminant(ps, method="frequency", conc="dna_quant")
head(contamdf.freq)

table(contamdf.freq$contaminant) # 22 contaminants (27 avant les runs séparés)
head(which(contamdf.freq$contaminant))

contam_asvs_freq <- row.names(contamdf.freq[contamdf.freq$contaminant == TRUE, ])
asv_tax[row.names(asv_tax) %in% contam_asvs_freq, ]


plot_frequency(ps, taxa_names(ps)[c(2,32)], conc="dna_quant") + 
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

set.seed(100)
plot_frequency(ps, taxa_names(ps)[sample(which(contamdf.freq$contaminant),3)], conc="dna_quant") +
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")




######################  IDENTIFY CONTAMINANTS - PREVALENCE #################################

sample_data(ps)$is.neg <- sample_data(ps)$true_or_blank == "Control sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")

table(contamdf.prev$contaminant) # 37 contaminants
head(which(contamdf.prev$contaminant))

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant) # 69 contaminants (68 avant)

contam_asvs_prev05 <- row.names(contamdf.prev05[contamdf.prev05$contaminant == TRUE, ])
contam <-asv_tax[row.names(asv_tax) %in% contam_asvs_prev05, ]
write.table(contam, "contam_prev.tsv",
            sep="\t", quote=F, col.names=NA)


ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$true_or_blank == "Control sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$true_or_blank == "True sample", ps.pa)


# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

pdf("PrevalenceDecontam.pdf")
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
dev.off()





######################  IDENTIFY CONTAMINANTS - COMBINED ###################################

contamdf.combi <- isContaminant(ps, method="combined", conc="dna_quant", neg="is.neg")
table(contamdf.combi$contaminant) # 8 contaminants

contam_asvs_combi <- row.names(contamdf.combi[contamdf.combi$contaminant == TRUE, ])
asv_tax[row.names(asv_tax) %in% contam_asvs_combi, ]


#ps.noncontam <- prune_taxa(!contamdf.combi$contaminant, ps)
#ps.noncontam 

asv_tab <- as(otu_table(ps),"matrix")
dim(asv_tab)





#--------------------------------------------------------------------------------------------#
#---------------------------------------NEW TABLES-------------------------------------------#
#--------------------------------------------------------------------------------------------#

# new otu_table
asv_tab_no_contam_freq <- asv_tab[!row.names(asv_tab) %in% contam_asvs_freq, ]
asv_tab_no_contam_prev05 <- asv_tab[!row.names(asv_tab) %in% contam_asvs_prev05, ]
asv_tab_no_contam_combi <- asv_tab[!row.names(asv_tab) %in% contam_asvs_combi, ]

# new tax_table
asv_tax_no_contam_freq <- asv_tax[!row.names(asv_tax) %in% contam_asvs_freq, ]
asv_tax_no_contam_prev05 <- asv_tax[!row.names(asv_tax) %in% contam_asvs_prev05, ]
asv_tax_no_contam_combi <- asv_tax[!row.names(asv_tax) %in% contam_asvs_combi, ]

# write tables of frequency decontam method
write.table(asv_tab_no_contam_freq, "ASVs_counts-no-contam_freq.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam_freq, "ASVs_taxonomy-no-contam_freq.tsv",
            sep="\t", quote=F, col.names=NA)

# write tables of prevalence (0.5) decontam method
write(asv_fasta_no_contam_prev05, "ASVs-no-contam_prev05.fa")
write.table(asv_tab_no_contam_prev05, "ASVs_counts-no-contam_prev05.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam_prev05, "ASVs_taxonomy-no-contam_prev05.tsv",
            sep="\t", quote=F, col.names=NA)

# write tables of combined decontam method
write(asv_fasta_no_contam_combi, "ASVs-no-contam_combi.fa")
write.table(asv_tab_no_contam_combi, "ASVs_counts-no-contam_combi.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam_combi, "ASVs_taxonomy-no-contam_combi.tsv",
            sep="\t", quote=F, col.names=NA)


# extract new metadata as dataframe from phyloseq 
meta <- function(x) { # https://rdrr.io/github/microbiome/microbiome/src/R/meta.R
  df <- as(sample_data(x), "data.frame")
  rownames(df) <- sample_names(x)
  df
}
metadata_decontam <- meta(ps)

# write decontamed metadata
write.table(metadata_decontam, "metadata_decontam.tsv", sep="\t", quote=F, col.names=NA)





#--------------------------------------------------------------------------------------------#
#--------------------DETECT REMAINING CHLOROPLAST AND MITOCHONDRIA---------------------------#
#--------------------------------------------------------------------------------------------#

# replace NA by 0 in the decontamed otu_table 
asv_tax_no_contam_prev05[is.na(asv_tax_no_contam_prev05)] <- 0

# create lists to storage names of chloroplast and mitochondria ASV 
contam_chloro = c()
contam_mito = c()
count_chloro = 0
count_mito = 0

# detect and storage ASV of chloroplast and mitochondria
for (i in 1:nrow(asv_tax_no_contam_prev05)){
  if (asv_tax_no_contam_prev05[i,"Order"]=="Chloroplast"){
    print(rownames(asv_tax_no_contam_prev05)[i])
    new_element <- rownames(asv_tax_no_contam_prev05)[i]
    contam_chloro <- c(contam_chloro,new_element)
    print(asv_tax_no_contam_prev05[i,4])
    count_chloro = count_chloro+1
    print(count_chloro) # 38 chloroplast
  }
  if (asv_tax_no_contam_prev05[i,"Family"]=="Mitochondria"){
    print(rownames(asv_tax_no_contam_prev05)[i])
    new_element <- rownames(asv_tax_no_contam_prev05)[i]
    contam_mito <- c(contam_mito,new_element)
    print(asv_tax_no_contam_prev05[i,4])
    count_mito = count_mito+1
    print(count_mito) # 2 mitochondria
  }
}

# replace 0 by NA in decontamed otu_table
asv_tax_no_contam_prev05[asv_tax_no_contam_prev05 == 0] <- NA





#--------------------------------------------------------------------------------------------#
#--------------------REMOVE REMAINING CHLOROPLAST AND MITOCHONDRIA---------------------------#
#--------------------------------------------------------------------------------------------#

# merge lists of chloroplast and mitochondria ASV 
contam_mito_chloro <- c(contam_chloro, contam_mito)

# create phyloseq object with decontamed tax and otu tables 
OTU_deseq = otu_table(asv_tab_no_contam_prev05, taxa_are_rows =TRUE)
TAX_deseq = tax_table(asv_tax_no_contam_prev05)
ps_deseq <- phyloseq(OTU_deseq, TAX_deseq, SAM)

# storage good taxa (removing the list of chloroplast and mitochondria ASV)
goodTaxa <- setdiff(taxa_names(ps_deseq), contam_mito_chloro)
length(goodTaxa)

# create new phyloseq object without contaminants ASV
ps_no_contam <- prune_taxa(goodTaxa, ps_deseq) 

# extract otu and tax tables from this new phyloseq object 
otu_decontam <- as(otu_table(ps_no_contam), "matrix")
tax_decontam <- as(tax_table(ps_no_contam), "matrix")

# write the final decontamed otu_table 
write.table(otu_decontam, "otu_decontam.tsv",
            sep="\t", quote=F, col.names=NA)

# write the final decontamed tax_table
write.table(tax_decontam, "tax_decontam.tsv",
            sep="\t", quote=F, col.names=NA)





#--------------------------------------------------------------------------------------------#
#---------------------------------------Save session-----------------------------------------#
#--------------------------------------------------------------------------------------------#

install.packages("session", repos = "http://cran.us.r-project.org")
library(session); packageVersion("session")
save.session("phyloseq.Rda")


