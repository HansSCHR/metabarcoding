# SCHRIEKE Hans
# dada2 function



library(dada2); packageVersion("dada2")


### FUNCTIONS ### 

plot_folder <- function(path, name, plotype, x, y){
  folder <- sprintf("%s", name)
  dir.create(folder)
  path2 <- getwd()
  path3 <- file.path(path2, folder)
  setwd(path3)


dada2 <- function(path, run){

  if (run == "run1"){
    path <- file.path(path, "run1")
  }
  if (run == "run2"){
    path <- file.path(path, "run2")
  }
  if (run == "run3"){
    path <- file.path(path, "run3")
  }
  if (run == "runs"){
    path <- file.path(path, "runs")
  }
  
  setwd(path)
  print(path)

  # Store forward and reverse reads 
  forward <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
  fwd <<- forward
  reverse <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
  rvs <<- reverse
  # Store sample names 
  sample.names <- sapply(strsplit(basename(forward), "_"), `[`, 1) 

  filtForward <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtReverse <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  
  out <- filterAndTrim(forward, filtForward, reverse, filtReverse, truncLen=c(0,240),
                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=TRUE)
  head(out)

  errForward <- learnErrors(filtForward, multithread=TRUE)
  errFwd <<- errForward
  errReverse <- learnErrors(filtReverse, multithread=TRUE)
  errRvs <<- errReverse
  
  # Dereplication
  derepForward <- derepFastq(filtForward, verbose=TRUE)
  derepReverse <- derepFastq(filtReverse, verbose=TRUE)
  # Name the derep-class objects by the sample names
  names(derepForward) <- sample.names
  names(derepReverse) <- sample.names

  dadaForward <- dada(derepForward, err=errForward, multithread=TRUE)
  dadaReverse <- dada(derepReverse, err=errReverse, multithread=TRUE)

  mergers <- mergePairs(dadaForward, derepForward, dadaReverse, derepReverse, verbose=TRUE)
  head(mergers[[1]])
  print("Reads merged.")
  
  seqtab <- makeSequenceTable(mergers)
  
  # Remove chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  dim(seqtab.nochim)
  sum(seqtab.nochim)/sum(seqtab)
  
  # seqtable csv
  write.table(t(seqtab),  "seqtab.csv", sep=";", dec=",")
  write.table(t(seqtab.nochim),  "seqtabnochim.csv", sep=";", dec=",")
  print("Seqtable done.")

  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaForward, getN), sapply(dadaReverse, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
  
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  head(track)
  
  # stats csv
  write.table(track,"stats.csv",sep=";",dec=",") 
  print("Stats done.")

  taxa <- assignTaxonomy(seqtab.nochim, file.path(path, "silva_nr_v132_train_set.fa.gz"), multithread=TRUE)
  taxa <- addSpecies(taxa, file.path(path, "silva_species_assignment_v132.fa.gz"))
  
  taxa.print <- taxa # Removing sequence rownames for display only
  rownames(taxa.print) <- NULL
  head(taxa.print)
  
  # taxa csv
  write.table(taxa,"taxa.csv",sep=";",dec=",")
  print ("Taxonomy done.")
}




### MAIN ### 

path = "D:/data"
setwd(path)

## RUN 1 ##

# DADA2
dada2(path, "run1")
 
# PLOTS
plot_folder(path, "plot-quality")
pdf("1-plotForward.pdf")
plotQualityProfile(fwd, aggregate = TRUE) #forward reads
dev.off()
pdf("2-plotReverse.pdf")
plotQualityProfile(rvs, aggregate = TRUE) #reverse reads
dev.off()
plot_folder(path, "plot-error")
pdf("3-plotErrorsFwd.pdf")
plotErrors(errFwd, nominalQ=TRUE)
dev.off()
pdf("4-plotErrorsRvs.pdf")
plotErrors(errRvs, nominalQ=TRUE)
dev.off()
  
path = "D:/data"
setwd(path)
  

## RUN 2 ## 
  
# DADA2 #
dada2(path, "run2")
  
# PLOTS
plot_folder(path, "plot-quality")
pdf("1-plotForward.pdf")
plotQualityProfile(fwd, aggregate = TRUE) #forward reads
dev.off()
pdf("2-plotReverse.pdf")
plotQualityProfile(rvs, aggregate = TRUE) #reverse reads
dev.off()
plot_folder(path, "plot-error")
pdf("3-plotErrorsFwd.pdf")
plotErrors(errFwd, nominalQ=TRUE)
dev.off()
pdf("4-plotErrorsRvs.pdf")
plotErrors(errRvs, nominalQ=TRUE)
dev.off()

path = "D:/data"
setwd(path)
  

## RUN 3 ## 

# DADA2 
dada2(path, "run3")

# PLOTS
plot_folder(path, "plot-quality")
pdf("1-plotForward.pdf")
plotQualityProfile(fwd, aggregate = TRUE) #forward reads
dev.off()
pdf("2-plotReverse.pdf")
plotQualityProfile(rvs, aggregate = TRUE) #reverse reads
dev.off()
plot_folder(path, "plot-error")
pdf("3-plotErrorsFwd.pdf")
plotErrors(errFwd, nominalQ=TRUE)
dev.off()
pdf("4-plotErrorsRvs.pdf")
plotErrors(errRvs, nominalQ=TRUE)
dev.off()

path = "D:/data"
setwd(path)


## All RUNS ## 
  
# DADA2 
dada2(path, "runs")

# PLOTS
plot_folder(path, "plot-quality")
pdf("1-plotForward.pdf")
plotQualityProfile(fwd, aggregate = TRUE) #forward reads
dev.off()
pdf("2-plotReverse.pdf")
plotQualityProfile(rvs, aggregate = TRUE) #reverse reads
dev.off()
plot_folder(path, "plot-error")
pdf("3-plotErrorsFwd.pdf")
plotErrors(errFwd, nominalQ=TRUE)
dev.off()
pdf("4-plotErrorsRvs.pdf")
plotErrors(errRvs, nominalQ=TRUE)
dev.off()
  
path = "D:/data"
setwd(path)
  
 
## SAVE SESSION ##   
library(session); packageVersion("session")
save.session("dada2_v2.Rda")
