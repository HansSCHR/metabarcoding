# SCHRIEKE Hans
# dada2 function
# trunc 245 235 , maxEE 3 4


library(dada2); packageVersion("dada2")



plot_folder <- function(path, name){
  folder <- sprintf("%s", name)
  dir.create(folder)
  path2 <- getwd()
  path3 <- file.path(path2, folder)
  setwd(path3)
  
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


dada2 <- function(path, run){

  if (run == "run1_new"){
    path <- file.path(path, "run1_new")
  }
  if (run == "run2_new"){
    path <- file.path(path, "run2_new")
  }
  if (run == "run3_new"){
    path <- file.path(path, "run3_new")
  }
  if (run == "run4_new"){
    path <- file.path(path, "run4_new")
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
  
  out <- filterAndTrim(forward, filtForward, reverse, filtReverse, truncLen=c(245,235),
                       maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=TRUE)
  head(out)
  print("filter and trim ok")


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
  
  # csv
  write.table(t(seqtab),  "seqtab3.csv", sep=";", dec=",")
  write.table(t(seqtab.nochim),  "seqtabnochim3.csv", sep=";", dec=",")
  saveRDS(seqtab.nochim, file="seqtab3.rds")
  seqtabnochim3 <<- seqtab.nochim
  print("Seqtable done.")



  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaForward, getN), sapply(dadaReverse, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
  
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  head(track)
  
  # csv
  write.table(track,"stats3.csv",sep=";",dec=",") 
  saveRDS(track, file="stats3.rds")
  stats3 <<- track
  print("Stats done.")



  taxa <- assignTaxonomy(seqtab.nochim, file.path(path, "silva_nr_v132_train_set.fa.gz"), multithread=TRUE, minBoot=80)
  taxa <- addSpecies(taxa, file.path(path, "silva_species_assignment_v132.fa.gz"))
  
  #taxa.print <- taxa # Removing sequence rownames for display only
  #rownames(taxa.print) <- NULL
  #head(taxa.print)
  write.table(taxa,"taxa3.csv",sep=";",dec=",")
  saveRDS(taxa, file="taxa3.rds")
  taxa3 <<- taxa
  print ("Taxonomy done.")
  
  # Phylogenetic tree (with decipher and phangorn packages)
  #library(DECIPHER)
  #library(phangorn)
  
  #seqs <- getSequences(seqtab.nochim)
  #names(seqs) <- seqs
  #alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
  
  #phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
  #dm <- dist.ml(phang.align)
  #treeNJ <- NJ(dm)
  #fit = pml(treeNJ, data=phang.align)
  
  #fitGTR <- update(fit, k=4, inv=0.2)
  #fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      #rearrangement = "stochastic", control = pml.control(trace = 0))
  #detach("package:phangorn", unload=TRUE)
  
  #write.table(fitGTR$tree, "phylotree.csv", sep=";", dec=",")
  #print("Job done")

}




#path = "/gs7k1/home/schrieke/stage/dada2"
#setwd(path)

#dada2(path, "run1")
#plot_folder(path, "plot-quality")

#pdf("1-plotForward.pdf")
#plotQualityProfile(fwd, aggregate = TRUE) #forward reads
#dev.off()

#pdf("2-plotReverse.pdf")
#plotQualityProfile(rvs, aggregate = TRUE) #reverse reads
#dev.off()

#plot_folder(path, "plot-error")

#pdf("3-plotErrorsFwd.pdf")
#plotErrors(errFwd, nominalQ=TRUE)
#dev.off()

#pdf("4-plotErrorsRvs.pdf")
#plotErrors(errRvs, nominalQ=TRUE)
#dev.off()





#path = "/gs7k1/home/schrieke/stage/dada2"
#setwd(path)

#dada2(path, "run2v4")

#plot_folder(path, "plot-quality")

#pdf("1-plotForward.pdf")
#plotQualityProfile(fwd, aggregate = TRUE) #forward reads
#dev.off()

#pdf("2-plotReverse.pdf")
#plotQualityProfile(rvs, aggregate = TRUE) #reverse reads
#dev.off()

#plot_folder(path, "plot-error")

#pdf("3-plotErrorsFwd.pdf")
#plotErrors(errFwd, nominalQ=TRUE)
#dev.off()

#pdf("4-plotErrorsRvs.pdf")
#plotErrors(errRvs, nominalQ=TRUE)
#dev.off()



#path = "/gs7k1/home/schrieke/DADA2"
#setwd(path)

#dada2(path, "run3")
#plot_folder(path, "plot-quality")

#pdf("1-plotForward.pdf")
#plotQualityProfile(fwd, aggregate = TRUE) #forward reads
#dev.off()

#pdf("2-plotReverse.pdf")
#plotQualityProfile(rvs, aggregate = TRUE) #reverse reads
#dev.off()

#plot_folder(path, "plot-error")

#pdf("3-plotErrorsFwd.pdf")
#plotErrors(errFwd, nominalQ=TRUE)
#dev.off()

#pdf("4-plotErrorsRvs.pdf")
#plotErrors(errRvs, nominalQ=TRUE)
#dev.off()


path = "/gs7k1/home/schrieke/stage/dada2"
setwd(path)

dada2(path, "run3_new")

plot_folder(path, "plot-quality3")

pdf("1-plotForward3.pdf")
plotQualityProfile(fwd, aggregate = TRUE) #forward reads
dev.off()

pdf("2-plotReverse3.pdf")
plotQualityProfile(rvs, aggregate = TRUE) #reverse reads
dev.off() 

plot_folder(path, "plot-error3")

pdf("3-plotErrorsFwd3.pdf")
plotErrors(errFwd, nominalQ=TRUE)
dev.off()

pdf("4-plotErrorsRvs3.pdf")
plotErrors(errRvs, nominalQ=TRUE)
dev.off()

setwd(path)
  
install.packages("session", repos="http://cran.r-project.org")
library(session); packageVersion("session")
save.session("run3.Rda")

  
