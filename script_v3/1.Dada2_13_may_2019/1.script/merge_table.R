# SCHRIEKE Hans 
# Merge tables of run1, run2 and run3 from separated dada2 processing 


library(dada2)

setwd("D:/stage/data/runs_new")

# ASV table 
st1 <- readRDS("seqtab1.rds")
st2 <- readRDS("seqtab2.rds")
st3 <- readRDS("seqtab3.rds")

st1 <- t(st1)
st2 <- t(st2)
st3 <- t(st3)

st1 <- as.data.frame(st1)
st2 <- as.data.frame(st2)
st3 <- as.data.frame(st3)

#st1$Undetermined <- st1$Undetermined1

names(st1)[97]
names(st1)[97] <- "Undetermined1"
names(st1)[97]

names(st2)[80]
names(st2)[80] <- "Undetermined2"
names(st2)[80]

names(st3)[78]
names(st3)[78] <- "Undetermined3"
names(st3)[78]

st1 <- t(st1)
st2 <- t(st2)
st3 <- t(st3)

st1 <- as.matrix(st1)
st2 <- as.matrix(st2)
st3 <- as.matrix(st3)


merged_seqtabnochim <- mergeSequenceTables(st1, st2, orderBy = "abundance")
merged_seqtabnochim <- mergeSequenceTables(merged_seqtabnochim, st3, orderBy = "abundance")
merged_seqtabnochim <- t(merged_seqtabnochim)

merged_seqtabnochim <- merged_seqtabnochim[, order(colnames(merged_seqtabnochim))]
write.table(merged_seqtabnochim, "seqtabnochim.csv", sep=";", dec=",")



# Stats 
stat1 <- read.table("stats1.csv", header=T, row.names=1, check.names=F, sep=";", dec=",")
stat2 <- read.table("stats2.csv", header=T, row.names=1, check.names=F, sep=";", dec=",")
stat3 <- read.table("stats3.csv", header=T, row.names=1, check.names=F, sep=";", dec=",")

stats <- rbind(stat1, stat2)
stats <- rbind(stats, stat3)

write.table(stats, "stats.csv", sep=";", dec=",")


# Taxa table 
tx1 <- readRDS("taxa1.rds")
tx2 <- readRDS("taxa2.rds")
tx3 <- readRDS("taxa3.rds")

tx1 <- as.data.frame(tx1)
tx2 <- as.data.frame(tx2)
tx3 <- as.data.frame(tx3)

taxa <- rbind(tx1, tx2)
taxa <- rbind(taxa, tx3)
taxa <- as.data.frame(taxa)

write.table(taxa, "taxa.csv", sep=";", dec=",")


