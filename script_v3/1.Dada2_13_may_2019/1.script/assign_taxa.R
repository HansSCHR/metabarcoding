# SCHRIEKE Hans
# dada2 - assignTaxonomy
# Input : merged otu table (included 3 runs)


path = "/gs7k1/home/schrieke/stage/dada2/taxa_assignment"
setwd(path)

otu <- read.csv("seqtabnochimcor.csv", sep=";", dec=",")
otu <- t(otu)

library(dada2)

taxa <- assignTaxonomy(otu, file.path(path, "silva_nr_v132_train_set.fa.gz"), multithread=TRUE, minBoot=80)
taxa <- addSpecies(taxa, file.path(path, "silva_species_assignment_v132.fa.gz"))

#taxa.print <- taxa # Removing sequence rownames for display only
#rownames(taxa.print) <- NULL
#head(taxa.print)

write.table(taxa,"taxafinal.csv",sep=";",dec=",")
saveRDS(taxa, file="taxafinal.rds")

print ("Taxonomy done.")

