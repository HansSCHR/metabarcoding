# SCHRIEKE Hans
# dada2
# input : fastq files of the V3-V4 region of the 16S rRNA from mosquito
# output : otu_table, tax_table, stats 


#--------------------------------------------------------------------------------------------#
#---------------------------------------LOAD DADA2-------------------------------------------#
#--------------------------------------------------------------------------------------------#

library(dada2); packageVersion("dada2")




#--------------------------------------------------------------------------------------------#
#----------------------------------------FUNCTIONS-------------------------------------------#
#--------------------------------------------------------------------------------------------#


# CREATE FOLDER (FOR PLOT)

plot_folder <- function(path, name){
  folder <- sprintf("%s", name)
  dir.create(folder)
  path2 <- getwd()
  path3 <- file.path(path2, folder)
  setwd(path3)
  
  # loop for plot doesn't work - I don't know why for now but I will update
  #if (plotype == "quality"){
    #print(x)
    #pdf("1-plotForward.pdf")
    #plotQualityProfile(x, aggregate = TRUE) #forward reads
    #dev.off()
    
    #print(y)
    #pdf("2-plotReverse.pdf")
    #plotQualityProfile(y, aggregate = TRUE) #reverse reads
    #dev.off()
  #}
  
  #if (plotype == "error"){
    #pdf("3-plotErrorsFwd.pdf")
    #plotErrors(x, nominalQ=TRUE)
    #dev.off()
    
    #pdf("4-plotErrorsRvs.pdf")
    #plotErrors(y, nominalQ=TRUE)
    #dev.off()
  #}
  
}





# DADA2 PROCESSING 

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
    path <- file.path(path, "runscopy")
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
  
  library(session); packageVersion("session")
  save.session("dada2_v2_run3.Rda")



  filtForward <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtReverse <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  
  out <- filterAndTrim(forward, filtForward, reverse, filtReverse, truncLen=c(245,235),
                       maxN=0, maxEE=c(3,3), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=TRUE)
  #head(out)

  library(session); packageVersion("session")
  save.session("dada2_v2_run3.Rda")

  errForward <- learnErrors(filtForward, multithread=TRUE)
  errFwd <<- errForward
  errReverse <- learnErrors(filtReverse, multithread=TRUE)
  errRvs <<- errReverse

  library(session); packageVersion("session")
  save.session("dada2_v2_run3.Rda")

  # Dereplication
  derepForward <- derepFastq(filtForward, verbose=TRUE)
  derepReverse <- derepFastq(filtReverse, verbose=TRUE)
  # Name the derep-class objects by the sample names
  names(derepForward) <- sample.names
  names(derepReverse) <- sample.names

  library(session); packageVersion("session")
  save.session("dada2_v2_run3.Rda")

  dadaForward <- dada(derepForward, err=errForward, multithread=TRUE)
  dadaReverse <- dada(derepReverse, err=errReverse, multithread=TRUE)

  library(session); packageVersion("session")
  save.session("dada2_v2_run3.Rda")

  mergers <- mergePairs(dadaForward, derepForward, dadaReverse, derepReverse, verbose=TRUE)
  head(mergers[[1]])
  print("Reads merged.")

  library(session); packageVersion("session")
  save.session("dada2_v2_run3.Rda")


  seqtab <- makeSequenceTable(mergers)
  
  # Remove chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  dim(seqtab.nochim)
  sum(seqtab.nochim)/sum(seqtab)
  
  # csv
  write.table(t(seqtab),  "seqtab.csv", sep=";", dec=",")
  write.table(t(seqtab.nochim),  "seqtabnochim.csv", sep=";", dec=",")
  print("Seqtable done.")

  library(session); packageVersion("session")
  save.session("dada2_v2_run3.Rda")

  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaForward, getN), sapply(dadaReverse, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
  
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  head(track)
  
  # csv
  write.table(track,"stats(maxEE33).csv",sep=";",dec=",") 
  print("Stats done.")

  library(session); packageVersion("session")
  save.session("dada2_v2_run3.Rda")

  taxa <- assignTaxonomy(seqtab.nochim, file.path(path, "silva_nr_v132_train_set.fa.gz"), multithread=TRUE, minBoot=80)
  taxa <- addSpecies(taxa, file.path(path, "silva_species_assignment_v132.fa.gz"))
  
  write.table(taxa,"taxa.csv",sep=";",dec=",")
  print ("Taxonomy done.")
}





#--------------------------------------------------------------------------------------------#
#-------------------------------------------MAIN---------------------------------------------#
#--------------------------------------------------------------------------------------------#



# RUN1 

path = "D:/stage/data"
setwd(path)

dada2(path, "run1")
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



# RUN 2

path = "D:/stage/data"
setwd(path)

dada2(path, "run2")

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



# RUN 3

path = "D:/data"
setwd(path)

dada2(path, "run3")
plot_folder(path, "plot-quality")

pdf("1-plotForward.pdf")
plotQualityProfile(x, aggregate = TRUE) #forward reads
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



# RUNS

path = "D:/stage/data"
setwd(path)

dada2(path, "runs")

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
  



#--------------------------------------------------------------------------------------------#
#---------------------------------------Save session-----------------------------------------#
#--------------------------------------------------------------------------------------------#

install.packages(session)
library(session); packageVersion("session")
save.session("dada2.Rda")

  