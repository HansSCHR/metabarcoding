# Set your working directory 
path <- "/gs7k1/home/schrieke/stage/dada2/run3"
setwd(path)

# Load packages 
library(dada2); packageVersion("dada2")

# Store forward and reverse reads 
forward <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
reverse <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Store sample names 
sample.names <- sapply(strsplit(basename(forward), "_"), `[`, 1) 

# Prepare folds and names of filter reads
filtForward <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtReverse <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Filter and trim process
out <- filterAndTrim(forward, filtForward, reverse, filtReverse, truncLen=c(245,235),
                     maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
  
  
  # Check output
  head(out)

# Learn errors
errForward <- learnErrors(filtForward, multithread=TRUE)
errReverse <- learnErrors(filtReverse, multithread=TRUE)

# Dereplication
derepForward <- derepFastq(filtForward, verbose=TRUE)
derepReverse <- derepFastq(filtReverse, verbose=TRUE)

#Name the derep-class objects by the sample names
names(derepForward) <- sample.names
names(derepReverse) <- sample.names

# Sample inference
dadaForward <- dada(derepForward, err=errForward, multithread=TRUE)
dadaReverse <- dada(derepReverse, err=errReverse, multithread=TRUE)

# Merge
mergers <- mergePairs(dadaForward, derepForward, dadaReverse, derepReverse, verbose=TRUE)

# Check merge
head(mergers[[1]])

# Construct table 
seqtab <- makeSequenceTable(mergers)

# Construct table without chimera
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Percentage of removing chimera
sum(seqtab.nochim)/sum(seqtab)*100

# Set directory to output
dir.create("output")
path_output <- "/gs7k1/home/schrieke/stage/dada2/run3/output"
setwd(path_output)

# Write tables to disk 
write.table(t(seqtab),  "seqtab3.csv", sep=";", dec=",")
write.table(t(seqtab.nochim),  "seqtabnochim3.csv", sep=";", dec=",")
saveRDS(seqtab.nochim, file="seqtab3.rds")

# Set work directory
setwd(path)

# Sum function
getN <- function(x) sum(getUniques(x))

# Create table with sum of dada forward reads, dada reverse reads, merged reads and count of reads
track <- cbind(out, sapply(dadaForward, getN), sapply(dadaReverse, getN), 
               sapply(mergers, getN),rowSums(seqtab.nochim))

# Rename columns and rows
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

# Check track 
head(track)

# Set directory to output
setwd(path_output)

# Write track on disk 
write.table(track,"stats3.csv",sep=";",dec=",") 
saveRDS(track, file="stats3.rds")

# Set work directory 
setwd(path)

