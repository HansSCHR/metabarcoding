library("DECIPHER")
library("dada2")
library("phangorn")

path = "D:/stage/data/runs_new2"
setwd(path)

otu <- read.csv("seqtabnochimcor.csv", sep=";", dec=",")
otu <- as.matrix(t(otu))

seqs <- getSequences(otu)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)



phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)



treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)



tree <- fitGTR$tree

library("session")
save.session("tree.Rda")
save(tree, file = "tree.RData")
