# SCHRIEKE Hans
# DECONTAM
# input : tables from dada2 step
# output : phyloseq objects with decontam tables 



#--------------------------------------------------------------------------------------------#
#--------------------------------------PATH SETTING------------------------------------------#
#--------------------------------------------------------------------------------------------#

path = "D:/stage/data/runs_new2/"
setwd(path)

dir.create("decontam_script")
path_main <- "D:/stage/data/runs_new2/decontam_script/"# folder for plot
setwd(path_main)

dir.create("1.raw_data")
path_raw <- "D:/stage/data/runs_new2/decontam_script/1.raw_data"

dir.create("2.decontam")
path_decontam <- "D:/stage/data/runs_new2/decontam_script/2.decontam"

dir.create("3.decontam2")
path_decontam2 <- "D:/stage/data/runs_new2/decontam/script/3.decontam2"




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
setwd(path)
asv_tab <- read.csv("seqtabnochimcor.csv", sep=";", dec=",")
asv_tax <- read.csv("taxafinal.csv", sep=";", dec=",")
metadata <- read.csv("metadata_runs_update_02_07_19.csv", sep=",", row.names = 1)
load("tree.Rdata")


#--------------------------------------------------------------------------------------------#
#--------------------------------EXTRACT GOOD STANDARDS--------------------------------------#
#--------------------------------------------------------------------------------------------#

setwd(path_raw)

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- rownames(asv_tab)
asv_headers <- vector(dim(asv_tab)[1], mode="character")
asv_headers2 <- vector(dim(asv_tab)[1], mode="character")

for (i in 1:dim(asv_tab)[1]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
  asv_headers2[i] <- paste("ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta_standard <- c(rbind(asv_headers, asv_seqs))
asv_fasta <- cbind(asv_headers2, asv_seqs)
write(asv_fasta_standard, "ASVs_standard.fa")
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab_standard <- as.matrix(asv_tab)
row.names(asv_tab_standard) <- sub(">", "", asv_headers)

# tax table:
asv_tax_standard <- as.matrix(asv_tax)
row.names(asv_tax_standard) <- sub(">", "", asv_headers)

# table asv name - sequence
table_asv_sequence <- cbind(asv_tax_standard[,0], rownames(asv_tax))
colnames(table_asv_sequence) <- c("Sequence")
write(table_asv_sequence, "ASVs_name_sequence.fa")




#--------------------------------------------------------------------------------------------#
#------------------------------------PHYLOSEQ OBJECT-----------------------------------------#
#--------------------------------------------------------------------------------------------#

# Phyloseq object with STANDARD
OTU_std = otu_table(as.matrix(asv_tab_standard), taxa_are_rows =TRUE)
TAX_std = tax_table(as.matrix(asv_tax_standard))
SAM = sample_data(metadata)

ps_std <- phyloseq(OTU_std, TAX_std, SAM)
ps_std1 <- subset_samples(ps_std, Species!="CuT" & Species!="CuP" & Species!="CuN" & Species!="CuG") # remove cullicoides	

asv_tab_std0 <- as(otu_table(ps_std1),"matrix")
write.table(asv_tab_std0, "0.asv_tab_standard0.tsv", sep="\t", quote=F, col.names=NA)

ps_std <- subset_samples(ps_std1, Sample!="NP16" & Sample!="NP17" & Sample!="NP18" & Sample!="NP19" & Sample!="S69" & Sample!="S70" & Sample!="S81" & Sample!="S82" & Sample!="NP20" & Sample!="NP21") # 11 juin 2019
ps_std <- prune_taxa(taxa_sums(ps_std) >= 1, ps_std) # remove asv that are not present in samples 
ps_std <- prune_samples(sample_sums(ps_std) >= 1, ps_std) # remove sample with 0 read


# Phyloseq object WITHOUT STANDARD
OTU = otu_table(asv_tab, taxa_are_rows =TRUE)
TAX = tax_table(as.matrix(asv_tax))
TREE = phy_tree(tree)

ps <- phyloseq(OTU, TAX, SAM, TREE)
ps <- subset_samples(ps, Species!="CuT" & Species!="CuP" & Species!="CuN" & Species!="CuG") # remove cullicoides
ps <- subset_samples(ps, Sample!="NP16" & Sample!="NP17" & Sample!="NP18" & Sample!="NP19" & Sample!="S69" & Sample!="S70" & Sample!="S81" & Sample!="S82" & Sample!="NP20" & Sample!="NP21")
ps <- prune_taxa(taxa_sums(ps) >= 1, ps) # remove asv that are not present in samples
ps <- prune_samples(sample_sums(ps) >= 1, ps) # remove counts = 0

write.table(asv_tab_standard, "1.asv_tab_standard.tsv", sep="\t", quote=F, col.names=NA)

ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2454 taxa and 208 samples ]
# sample_data() Sample Data:       [ 208 samples by 15 sample variables ]
# tax_table()   Taxonomy Table:    [ 2454 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 2454 tips and 2452 internal nodes ]


#--------------------------------------------------------------------------------------------#
#----------------------------------------DECONTAM--------------------------------------------#
#--------------------------------------------------------------------------------------------#

setwd(path_decontam)

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
sample_data(ps_std)$is.neg <- sample_data(ps_std)$Control == "Control sample"

contamdf.prev01 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.1)
contamdf.prev01_std <- isContaminant(ps_std, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev01$contaminant) # 36 contaminants

contam_asvs_prev01 <- row.names(contamdf.prev01[contamdf.prev01$contaminant == TRUE, ])
contam_asvs_prev01_std <- row.names(contamdf.prev01_std[contamdf.prev01_std$contaminant == TRUE, ])

contam_standard <-asv_tax_standard[row.names(asv_tax_standard) %in% contam_asvs_prev01_std, ]
contam <-asv_tax[row.names(asv_tax) %in% contam_asvs_prev01, ]


ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Control == "Control sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Control == "True sample", ps.pa)


# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev01$contaminant)

pdf("Prevalence_Decontam.pdf")
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
dev.off()

# write table of contaminant asv 
write.table(contam, "asv_contam_prev01.tsv",
            sep="\t", quote=F, col.names=NA)





#--------------------------------------------------------------------------------------------#
#---------------------------------------NEW TABLES-------------------------------------------#
#--------------------------------------------------------------------------------------------#

# new otu table 
asv_tab_prev01 <- asv_tab[!row.names(asv_tab) %in% contam_asvs_prev01, ]
asv_tab_prev01_std <- asv_tab_standard[!row.names(asv_tab_standard) %in% contam_asvs_prev01_std, ]

# new tax table
asv_tax_prev01 <- asv_tax[!row.names(asv_tax) %in% contam_asvs_prev01, ]
asv_tax_prev01_std <- asv_tax_standard[!row.names(asv_tax_standard) %in% contam_asvs_prev01_std, ]

# new fasta  
contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs_prev01))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_prev01 <- asv_fasta[- dont_want]

# write tables
write.table(asv_tab_prev01, "asv_tab_prev01.tsv", sep="\t", quote=F, col.names=NA)
write.table(asv_tab_prev01_std, "2.asv_tab_prev01_std.tsv", sep="\t", quote=F, col.names=NA)

write.table(asv_tax_prev01, "asv_tax_prev01.tsv", sep="\t", quote=F, col.names=NA)
write.table(asv_tax_prev01_std, "asv_tax_prev01_std.tsv", sep="\t", quote=F, col.names=NA)

write(asv_fasta_prev01, "fasta_prev01.fa")


# extract new metadata as dataframe from phyloseq 
meta <- function(x) { # https://rdrr.io/github/microbiome/microbiome/src/R/meta.R
  df <- as(sample_data(x), "data.frame")
  rownames(df) <- sample_names(x)
  df
}
metadata_prev01 <- meta(ps)

# write decontamed metadata
write.table(metadata_prev01, "metadata_decontam.tsv", sep="\t", quote=F, col.names=NA)




#--------------------------------------------------------------------------------------------#
#-----------------------------NEW DECONTAM PHYLOSEQ OBJECT-----------------------------------#
#--------------------------------------------------------------------------------------------#

OTU_decontam <- otu_table(asv_tab_prev01_std, taxa_are_rows = TRUE)
TAX_decontam <- tax_table(as.matrix(asv_tax_prev01_std))
SAM_decontam <- sample_data(metadata_prev01)

ps_decontam <- phyloseq(OTU_decontam, TAX_decontam, SAM_decontam)
ps_decontam <- prune_taxa(taxa_sums(ps_decontam) >= 1, ps_decontam)# remove counts = 0 in tax_table
ps_decontam
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2418 taxa and 208 samples ]
# sample_data() Sample Data:       [ 208 samples by 16 sample variables ]
# tax_table()   Taxonomy Table:    [ 2418 taxa by 7 taxonomic ranks ]

ps_decontam <- subset_taxa(ps_decontam, Kingdom == "Bacteria") # select only Bacteria
ps_decontam <- prune_samples(sample_sums(ps_decontam) >= 1, ps_decontam) # remove counts = 0 in otu_table

ps_decontam <- subset_samples(ps_decontam, Control!="Control sample") # remove blanks

#--------------------------------------------------------------------------------------------#
#-------------------REMOVING REMAINING CHLOROPLAST AND MITOCHONDRIA--------------------------#
#--------------------------------------------------------------------------------------------#
ps_decontam2 <- subset_taxa(ps_decontam, Family!="Mitochondria" & Order!="Chloroplast")
ps_decontam2 <- prune_samples(sample_sums(ps_decontam2) >= 1, ps_decontam2)
ps_decontam2 <- prune_taxa(taxa_sums(ps_decontam2) >= 1, ps_decontam2)

asv_decontam2 <- as(otu_table(ps_decontam2),"matrix")

write.table(asv_decontam2, "3.asv_decontam2.tsv",sep="\t", quote=F, col.names=NA)

mitochondria <- subset_taxa(ps_decontam, Family=="Mitochondria") # 1 mitochondria 
chloroplast <- subset_taxa(ps_decontam, Order=="Chloroplast") # 26 chloroplast 

tax_mitochondria <- as(tax_table(mitochondria),"matrix")
tax_chloroplast <- as(tax_table(chloroplast),"matrix")


#--------------------------------------------------------------------------------------------#
#-------------------------------------NORMALIZATION------------------------------------------#
#--------------------------------------------------------------------------------------------#

# % 
ps_percent <- transform_sample_counts( ps_decontam2, function(x) x/sum(x)*100 )
asv_percent <- as(otu_table(ps_percent),"matrix")
write.table(asv_percent, "4.asv_percent.tsv", sep="\t", quote=F, col.names=NA)


# deseq 
ps_deseq <- transform_sample_counts(ps_decontam2, function(x) x+1)

otu_deseq <- as(otu_table(ps_deseq),"matrix")
metadata_deseq <- as(sample_data(ps_deseq),"matrix")

deseq_counts <- DESeqDataSetFromMatrix(otu_deseq, colData = metadata_deseq, design = ~Dna)
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

asv_deseq <- as.matrix(assay(deseq_counts_vst))

write.table(asv_deseq, "5.asv_deseq.tsv", sep="\t", quote=F, col.names=NA)

otu_table(ps_deseq) <- otu_table(asv_deseq, taxa_are_rows = TRUE)

ps_deseq


#--------------------------------------------------------------------------------------------#
#-----------------------------------EXPLORING WOLBACHIA--------------------------------------#
#--------------------------------------------------------------------------------------------#
ps_wolbachia <- subset_taxa(ps_percent, Genus == "Wolbachia")
ps_wolbachia <- prune_taxa(taxa_sums(ps_wolbachia) >= 1, ps_wolbachia)
ps_wolbachia <- prune_samples(sample_sums(ps_wolbachia) >= 1, ps_wolbachia)



ps_proteo <- subset_taxa(ps_percent, Phylum=="Proteobacteria")
ps_wolbachia <- prune_taxa(taxa_sums(ps_proteo) >= 1, ps_proteo)
ps_proteo <- prune_samples(sample_sums(ps_proteo) >= 1, ps_proteo)


#--------------------------------------------------------------------------------------------#
#---------------------------------------Save objects-----------------------------------------#
#--------------------------------------------------------------------------------------------#

setwd(path_main)
save(ps, ps_decontam, ps_decontam2, ps_percent, ps_wolbachia, ps_proteo, ps_deseq, file = "objects.RData")

setwd(path)
save(ps, ps_decontam, ps_decontam2, ps_percent, ps_wolbachia, ps_proteo, ps_deseq, file = "objects.RData")




#--------------------------------------------------------------------------------------------#
#---------------------------------------Save session-----------------------------------------#
#--------------------------------------------------------------------------------------------#

#install.packages("session", repos = "http://cran.us.r-project.org")
library(session); packageVersion("session")
save.session("decontam2.Rda")
#load("phyloseq.Rda")

