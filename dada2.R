# DADA2 


# Prerequies : 
# R 3.5.0 or more 

# Set the working directory 
#setwd("C:/Users/U117-F435/Desktop/data")
setwd("C:/Users/SCHRIEKE Hans/Desktop/DADA2/data")


# Install dada2 package and load it 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.8")



# Check and load dada2 version
library(dada2); packageVersion("dada2")



# Assign data path 
#path <- /homedir/schrieke/data 
#path <- "C:/Users/U117-F435/Desktop/data"
path <- "C:/Users/SCHRIEKE Hans/Desktop/DADA2/data"
files <- list.files(path)


# Unzip fastq.gz 
install.packages("R.utils")
library("R.utils")
for (i in 1:length(files)) {gunzip(files[i])}
  

# Store forward and reverse reads 
forward <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
reverse <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))


# Store sample names 
sample.names <- sapply(strsplit(basename(forward), "_"), `[`, 1) 


# Reads quality profil to know where truncate
plotQualityProfile(forward[1:2]) #forward reads
plotQualityProfile(reverse[1:2]) #reverse reads


# Filter and trim
filtForward <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtReverse <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(forward, filtForward, reverse, filtReverse, truncLen=c(0,240),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
head(out)

# Learn the error rates
errForward <- learnErrors(filtForward, multithread=TRUE)
errReverse <- learnErrors(filtReverse, multithread=TRUE)

plotErrors(errForward, nominalQ=TRUE)

# Dereplication
derepForward <- derepFastq(filtForward, verbose=TRUE)
derepReverse <- derepFastq(filtReverse, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepForward) <- sample.names
names(derepReverse) <- sample.names


# Sample Inference
dadaForward <- dada(derepForward, err=errForward, multithread=TRUE)
dadaReverse <- dada(derepReverse, err=errReverse, multithread=TRUE)


# Merge paired reads
mergers <- mergePairs(dadaForward, derepForward, dadaReverse, derepReverse, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


# Sequence table 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaForward, getN), sapply(dadaReverse, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


# Assign taxonomy
# taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/U117-F435/Desktop/silva_nr_v132_train_set.fa.gz", multithread=FALSE)
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/SCHRIEKE Hans/Desktop/DADA2/silva_nr_v132_train_set.fa.gz", multithread=FALSE)
taxa <- addSpecies(taxa, "C:/Users/SCHRIEKE Hans/Desktop/DADA2/silva_species_assignment_v132.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
print(taxa)


# Accuracy 
