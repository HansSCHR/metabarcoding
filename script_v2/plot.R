#--------------------------------------------------------------------------------------------#
#--------------------------------------PATH SETTING------------------------------------------#
#--------------------------------------------------------------------------------------------#

path = "D:/stage/data/runs_new"
setwd(path)

dir.create("plotv3") # folder for plot
path2 <- "D:/stage/data/runs_new/plotv3"




#--------------------------------------------------------------------------------------------#
#-------------------------------------LOAD PACKAGES------------------------------------------#
#--------------------------------------------------------------------------------------------#

library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library("DESeq2")
library("plotly")

#library("microbiome")
library("hrbrthemes")
#library("dplyr")

scripts <- c("graphical_methods.R",
             "tree_methods.R",
             "plot_merged_trees.R",
             "specificity_methods.R",
             "ternary_plot.R",
             "richness.R",
             "edgePCA.R",
             "copy_number_correction.R",
             "import_frogs.R",
             "prevalence.R",
             "compute_niche.R",
             "NCM_fit.R")
urls <-
  paste0("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/",
         scripts)

for (url in urls) {
  source(url)
}


th <- theme_set(theme_bw()) #set ggplot2 graphic theme 

load("objects.Rdata")



#--------------------------------------------------------------------------------------------#
#---------------------------------CHECK NP10 NP11 NP12---------------------------------------#
#--------------------------------------------------------------------------------------------#

otu_ps <- as(otu_table(ps), "matrix")
otu_decontam <- as(otu_table(ps_decontam_bacteria), "matrix")
otu_decontam2 <- as(otu_table(ps_decontam2_bacteria), "matrix")


pull_NP10 <- subset_samples(ps_decontam2_bacteria, Sample=="NP10")
pull_NP11 <- subset_samples(ps_decontam2_bacteria, Sample=="NP11")
pull_NP12 <- subset_samples(ps_decontam2_bacteria, Sample=="NP12")

otu_pullNP10 <- as(otu_table(pull_NP10), "matrix")
otu_pullNP11 <- as(otu_table(pull_NP11), "matrix")
otu_pullNP12 <- as(otu_table(pull_NP12), "matrix")

otu_pullNP10 <- as.data.frame(otu_pullNP10)
otu_pullNP11 <- as.data.frame(otu_pullNP11)
otu_pullNP12 <- as.data.frame(otu_pullNP12)

tax_pullNP10 <- as(tax_table(pull_NP10), "matrix")
tax_pullNP11 <- as(tax_table(pull_NP11), "matrix")
tax_pullNP12 <- as(tax_table(pull_NP12), "matrix")

#rownames(otu_pull)[1,1]

#print(names(otu_pull[2]))

#colSums(otu_pull)
asv_NP10 <- c()
asv_NP11 <- c()
asv_NP12 <- c()

for (i in 1:nrow(otu_pullNP10)){
  #print(i)
  if (otu_pullNP10[i,1]!=0){
    print(rownames(otu_pullNP10)[i])
    print(otu_pullNP10[i,1])
    new <- rownames(otu_pullNP10)[i]
    asv_NP10 <- c(asv_NP10,new) #85 ASV pour 1 022 615 count
  }
}

colSums(otu_pullNP10)
asv_NP10 <- as.matrix(asv_NP10)
dim(asv_NP10)

for (i in 1:nrow(otu_pullNP11)){
  #print(i)
  if (otu_pullNP11[i,1]!=0){
    print(rownames(otu_pullNP11)[i])
    print(otu_pullNP11[i,1])
    new <- rownames(otu_pullNP11)[i]
    asv_NP11 <- c(asv_NP11,new) #95 ASV pour 1 208 708 count
  }
}

colSums(otu_pullNP11)
asv_NP11 <- as.matrix(asv_NP11)
dim(asv_NP11)

for (i in 1:nrow(otu_pullNP12)){
  #print(i)
  if (otu_pullNP12[i,1]!=0){
    print(rownames(otu_pullNP12)[i])
    print(otu_pullNP12[i,1])
    new <- rownames(otu_pullNP12)[i]
    asv_NP12 <- c(asv_NP12,new) #50 ASV pour 1748 count
  }
}

colSums(otu_pullNP12)
asv_NP12 <- as.matrix(asv_NP12)
dim(asv_NP12)

sum(otu_decontam2)
dim(as(tax_table(ps_decontam2_bacteria),"matrix")) #??? 2894 ASV pour 26 479 026 counts


# TOTAL ASV pour NP10, NP11, NP12 : 230 --> 8% des ASV 
# TOTAL COUNT NP10, NP11, NP12 : 2 233 071 --> 8% des counts 







# rownames(tax_pull)[1]
# print(asv_NP10)
# tax_pull[2,]

# tax_NP10 <- c()
for (i in 1:nrow(asv_NP10)){
  print(tax_pullNP10)[i,]
  # new <- tax_pull[i,]
  # tax_NP10 <- c(tax_NP10,new)
}

for (i in 1:nrow(asv_NP11)){
  print(tax_pullNP11)[i,]
}

for (i in 1:nrow(asv_NP12)){
  print(tax_pullNP12)[i,]
}




#--------------------------------------------------------------------------------------------#
#-----------------------------------RAREFACTION CURVE----------------------------------------#
#--------------------------------------------------------------------------------------------#

setwd(path2)

# I 
# p <- ggrare(ps,
#             step = 500,
#             color = "true_or_blank",
#             plot = T,
#             parallel = F,
#             se = T)
# 
# p <- p + 
#   facet_wrap(~ dna_from ) + 
#   geom_vline(xintercept = min(sample_sums(ps)), 
#              color = "gray60") +
#   ggtitle("Raw data")
# 
# pdf("rarecurve_ps.pdf")
# plot(p)
# dev.off()
# 
# 
# p2 <- ggrare(ps_decontam_bacteria,
#             step = 500,
#             color = "true_or_blank",
#             plot = T,
#             parallel = F,
#             se = T)
# 
# p2 <- p2 + 
#   facet_wrap(~ dna_from ) + 
#   geom_vline(xintercept = min(sample_sums(ps)), 
#              color = "gray60") +
#   ggtitle("After decontam")
# 
# pdf("rarecurve_ps_decontam.pdf")
# plot(p2)
# dev.off()


p3 <- ggrare(ps_decontam2_bacteria,
             step = 500,
             color = "Control",
             plot = T,
             parallel = F,
             se = T)

p3bis <- p3 + 
  facet_wrap(~ Organ ) + 
  geom_vline(xintercept = min(sample_sums(ps)), 
             color = "gray60") +
  ggtitle("After removing mito, chloro, archeae") +
  xlim(0,10000) +
  ylim(0, 170)

p4 <- p3 + 
  facet_wrap(~ Organ) + 
  geom_vline(xintercept = min(sample_sums(ps)), 
             color = "gray60") +
  ggtitle("After removing mito, chloro, archeae") +
  xlim(0,100000) +
  ylim(0, 230)
  

pdf("rarecurve_ps_decontam2.pdf")
plot(p3bis)
dev.off()

pdf("rarecurve2_ps_decontam2.pdf")
plot(p4)
dev.off()



#--------------------------------------------------------------------------------------------#
#--------------------------------------DISTRIBUTION------------------------------------------#
#--------------------------------------------------------------------------------------------#

readsumsdf <- data.frame(nreads = sort(taxa_sums(ps), TRUE),
                         sorted = 1:ntaxa(ps), 
                         type = "OTU")
# readsumsdf %>% head()

ggplot(readsumsdf, 
       aes(x = sorted, y = nreads)) + 
  geom_bar(stat = "identity") + 
  scale_y_log10() 

readsumsdf2 <- data.frame(nreads = sort(sample_sums(ps), TRUE), 
                          sorted = 1:nsamples(ps), 
                          type = "Samples")

readsumsdf3 <- rbind(readsumsdf,readsumsdf2)

# readsumsdf3 %>% head()
# readsumsdf3 %>% tail()

p  <-  ggplot(readsumsdf3, 
              aes(x = sorted, y = nreads)) + 
  geom_bar(stat = "identity")

pdf("distribution.pdf")
p + 
  ggtitle("Total number of reads before Preprocessing") + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
dev.off()





#--------------------------------------------------------------------------------------------#
#-----------------------------------------RICHNESS-------------------------------------------#
#--------------------------------------------------------------------------------------------#

ps_noblank <- subset_samples(ps, Control!="Control sample")
ps_decontam_noblank <- subset_samples(ps_decontam2_bacteria, Control!="Control sample")
ps_decontam2_noblank <- subset_samples(ps_decontam2_bacteria, Control!="Control sample")

data1 <-  filter_taxa(ps_noblank, 
                     function(x) sum(x >= 10) > (1), 
                     prune =  TRUE) 

data2 <-  filter_taxa(ps_decontam_noblank, 
                      function(x) sum(x >= 10) > (1), 
                      prune =  TRUE) 

data3 <-  filter_taxa(ps_decontam2_noblank, 
                      function(x) sum(x >= 10) > (1), 
                      prune =  TRUE) 


p1 <- plot_richness(data1,
                   x="Sample",
                   color="Species",
                   measures=c("Observed","Shannon","ACE", "Chao1"),
                   nrow = 1)+
  ggtitle("Raw data")

# pdf("richness_ps_species.pdf")
# print(p1)
# dev.off()



# p2 <- plot_richness(data1,
#                     x="Sample",
#                     color="Dna",
#                     measures=c("Observed","Shannon","ACE", "Chao1"),
#                     nrow = 1)+
#   ggtitle("Raw data")
# # pdf("richness_ps_dnafrom.pdf")
# # print(p2)
# # dev.off()
# # 
# # p1$data %>% head()
# 
# pdf("richness_ps_species2.pdf")
# ggplot(p1$data,aes(Species,value,colour=Species,shape=Species)) +
#   facet_grid(variable ~ Species, drop=T,scale="free",space="fixed") +
#   geom_boxplot(outlier.colour = NA,alpha=0.8, 
#                position = position_dodge(width=0.9)) + 
#   geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
#   ylab("Diversity index")  + xlab(NULL) + theme_bw() +
#   ggtitle("Raw data (species)")
# dev.off()
# 
# pdf("richness_ps_dnafrom2.pdf")
# ggplot(p1$data,aes(Organ,value,colour=Organ,shape=Species)) +
#   facet_grid(variable ~ Organ, drop=T,scale="free",space="fixed") +
#   geom_boxplot(outlier.colour = NA,alpha=0.8, 
#                position = position_dodge(width=0.9)) + 
#   geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
#   ylab("Diversity index")  + xlab(NULL) + theme_bw() +
#   ggtitle("Raw data (organ)")
# dev.off()
# 
# 
# 
# p1 <- plot_richness(data2, 
#                     x="Sample", 
#                     color="Species", 
#                     measures=c("Observed","Shannon","ACE", "Chao1"), 
#                     nrow = 1) +
#   ggtitle("After decontam")
# # pdf("richness_ps_decontam_species.pdf")
# # print(p1)
# # dev.off()
# 
# p2 <- plot_richness(data2, 
#                     x="Sample", 
#                     color="Organ", 
#                     measures=c("Observed","Shannon","ACE", "Chao1"), 
#                     nrow = 1) +
#   ggtitle("After decontam")
# # pdf("richness_ps_decontam_dnafrom.pdf")
# # print(p2)
# # dev.off()
# 
# pdf("richness_ps_decontam_species2.pdf")
# ggplot(p1$data,aes(Species,value,colour=Species,shape=Species)) +
#   facet_grid(variable ~ Species, drop=T,scale="free",space="fixed") +
#   geom_boxplot(outlier.colour = NA,alpha=0.8, 
#                position = position_dodge(width=0.9)) + 
#   geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
#   ylab("Diversity index")  + xlab(NULL) + theme_bw() +
#   ggtitle("After decontam (species)")
# dev.off()
# 
# pdf("richness_ps_decontam_dnafrom2.pdf")
# ggplot(p1$data,aes(Organ,value,colour=Organ,shape=Species)) +
#   facet_grid(variable ~ Organ, drop=T,scale="free",space="fixed") +
#   geom_boxplot(outlier.colour = NA,alpha=0.8, 
#                position = position_dodge(width=0.9)) + 
#   geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
#   ylab("Diversity index")  + xlab(NULL) + theme_bw() +
#   ggtitle("After decontam (dnafrom)")
# dev.off()


p1 <- plot_richness(data3, 
                    x="Sample", 
                    color="Species", 
                    measures=c("Observed","Shannon","ACE", "Chao1"), 
                    nrow = 1) +
  ggtitle("After removing mito, chloro, archeae")
# pdf("richness_ps_decontam2_species.pdf")
# print(p1)
# dev.off()

p2 <- plot_richness(data3, 
                    x="Sample", 
                    color="Organ", 
                    measures=c("Observed","Shannon","ACE", "Chao1"), 
                    nrow = 1) +
  ggtitle("After removing mito, chloro, archeae")
# pdf("richness_ps_decontam2_dnafrom.pdf")
# print(p2)
# dev.off()


pdf("richness_ps_decontam2_species2.pdf")
ggplot(p1$data,aes(Species,value,colour=Species,shape=Species)) +
  facet_grid(variable ~ Species, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.8, 
               position = position_dodge(width=0.9)) + 
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  ylab("Diversity index")  + xlab(NULL) + theme_bw() +
  ggtitle("After removing mito, chloro, archeae (species)")
dev.off()

pdf("richness_ps_decontam2_dnafrom2.pdf")
ggplot(p1$data,aes(Organ,value,colour=Organ,shape=Species)) +
  facet_grid(variable ~ Organ, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.8, 
               position = position_dodge(width=0.9)) + 
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  ylab("Diversity index")  + xlab(NULL) + theme_bw() +
  ggtitle("After removing mito, chloro, archeae (organ)")
dev.off()




#--------------------------------------------------------------------------------------------#
#----------------------------------TAXONOMIC SUMMARIES---------------------------------------#
#--------------------------------------------------------------------------------------------#
ps_deseq <- subset_samples(ps_deseq, Control!="Control sample")
ps_percent <- subset_samples(ps_percent, Control!="Control sample")

p <- plot_composition(ps_deseq,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Phylum", 
                      numberOfTaxa = 20, 
                      fill= "Phylum") +
  facet_wrap(~ Species, scales = "free_x", nrow = 1) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Deseq normalization (species)")

pdf("composition_ps_deseq_species.pdf")
plot(p)
dev.off()

p <- plot_composition(ps_deseq,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Phylum", 
                      numberOfTaxa = 20, 
                      fill= "Phylum") +
  facet_wrap(~ Organ, scales = "free_x", nrow = 1) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Deseq normalization (dnafrom)")

pdf("composition_ps_deseq_organ.pdf")
plot(p)
dev.off()




p <- plot_composition(ps_percent,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Phylum", 
                      numberOfTaxa = 20, 
                      fill= "Phylum") +
  facet_wrap(~ Species, scales = "free_x", nrow = 1) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Percent normalization (species)")

pdf("composition_ps_percent_species.pdf")
plot(p)
dev.off()

p <- plot_composition(ps_percent,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Phylum", 
                      numberOfTaxa = 20, 
                      fill= "Phylum") +
  facet_wrap(~ Organ, scales = "free_x", nrow = 1) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Percent normalization (organ)")

pdf("composition_ps_percent_organ.pdf")
plot(p)
dev.off()



# p <- plot_composition(data3,
#                       taxaRank1 = "Kingdom",
#                       taxaSet1 ="Bacteria",
#                       taxaRank2 = "Phylum", 
#                       numberOfTaxa = 20, 
#                       fill= "Phylum") +
#   facet_wrap(~ species, scales = "free_x", nrow = 1) + 
#   theme(plot.title = element_text(hjust = 0.5)) +
#   ggtitle("After removing mito, chloro, archeae (species)")
# 
# pdf("composition_decontam2_species.pdf")
# plot(p)
# dev.off()
# 
# p <- plot_composition(data3,
#                       taxaRank1 = "Kingdom",
#                       taxaSet1 ="Bacteria",
#                       taxaRank2 = "Phylum", 
#                       numberOfTaxa = 20, 
#                       fill= "Phylum") +
#   facet_wrap(~ dna_from, scales = "free_x", nrow = 1) + 
#   theme(plot.title = element_text(hjust = 0.5)) +
#   ggtitle("After removing mito, chrloro, archeae (dnafrom)")
# 
# pdf("composition_decontam2_dnafrom.pdf")
# plot(p)
# dev.off()






#--------------------------------------------------------------------------------------------#
#------------------------------------------PCOA----------------------------------------------#
#--------------------------------------------------------------------------------------------#


# # PCoA deseq (euclidean)
# pdf("PCoA_deseq_species.pdf")
# #jpeg("PCoA_deseq_species.jpg")
# plot_ordination(ps_deseq, ordinate(ps_deseq, method ="MDS", distance = "euclidean"), color = "species", shape="dna_from") +
#   geom_point(size = 3) +
#   ggtitle("Euclidean PCoA (deseq normalization)")
# dev.off()
# 
# pdf("PCoA_deseq_species2.pdf")
# #jpeg("PCoA_deseq_species.jpg")
# plot_ordination(ps_deseq, ordinate(ps_deseq, method ="MDS", distance = "euclidean"), color = "species", shape="location") +
#   geom_point(size = 3) +
#   ggtitle("Euclidean PCoA (deseq normalization)")
# dev.off()
# 
# pdf("PCoA_deseq_species3.pdf")
# #jpeg("PCoA_deseq_species.jpg")
# plot_ordination(ps_deseq, ordinate(ps_deseq, method ="MDS", distance = "euclidean"), color = "species", shape="date") +
#   geom_point(size = 3) +
#   scale_shape_manual(values=seq(0,15)) +
#   ggtitle("Euclidean PCoA (deseq normalization)")
# dev.off()

# PCoA percent (bray)
pdf("PCoA_percent_species.pdf")
#jpeg("PCoA_percent_species.jpg")
plot_ordination(ps_percent, ordinate(ps_percent, method ="MDS", distance = "bray"), color = "Species", shape="Organ") +
  geom_point(size = 3) +
  ggtitle("Bray PCoA (% normalization)")
dev.off()

pdf("PCoA_percent_species2.pdf")
#jpeg("PCoA_percent_species.jpg")
plot_ordination(ps_percent, ordinate(ps_percent, method ="MDS", distance = "bray"), color = "Species", shape="Location") +
  geom_point(size = 3) +
  ggtitle("Bray PCoA (% normalization)")
dev.off()

pdf("PCoA_percent_species3.pdf")
#jpeg("PCoA_percent_species.jpg")
plot_ordination(ps_percent, ordinate(ps_percent, method ="MDS", distance = "bray"), color = "Species", shape="Date") +
  geom_point(size = 3) +
  th +
  scale_shape_manual(values=seq(0,15)) +
  ggtitle("Bray PCoA (% normalization)")
dev.off()






# #--------------------------------------------------------------------------------------------#
# #--------------------------------NMMDS (euclidean distance)----------------------------------#
# #--------------------------------------------------------------------------------------------#
# 
# # phyloseq object creation for each condition
# ps.F <- subset_samples(ps_deseq, dna_from == "Full_body")
# ps.I <- subset_samples(ps_deseq, dna_from == "Intestine")
# ps.O <- subset_samples(ps_deseq, dna_from == "Ovary")
# ps.OI <- subset_samples(ps_deseq, dna_from != "Full_body")
# ps.OI <- subset_samples(ps_deseq, dna_from != "Salivary_gland")
# ps.OI <- subset_samples(ps_deseq, dna_from != "Organs_pull")
# ps.OI <- subset_samples(ps_deseq, dna_from != "Blank") # ovaries + intestine
# 
# #otu_deseq[otu_deseq < 0.0] <- 0.0
# 
# # otu_table transformation and ordination
# ps.prop.euc <- transform_sample_counts(ps_deseq, function(count_tab) count_tab/sum(count_tab))
# ord.nmds.euc <- ordinate(ps.prop.euc, method="NMDS", distance="euclidean")
# 
# ps.propF.euc <- transform_sample_counts(ps.F, function(count_tab) count_tab/sum(count_tab))
# ord.nmds.eucF <- ordinate(ps.propF.euc, method="NMDS", distance="euclidean")
# 
# ps.propI.euc <- transform_sample_counts(ps.I, function(count_tab) count_tab/sum(count_tab))
# ord.nmds.eucI <- ordinate(ps.propI.euc, method="NMDS", distance="euclidean")
# 
# ps.propO.euc <- transform_sample_counts(ps.O, function(count_tab) count_tab/sum(count_tab))
# ord.nmds.eucO <- ordinate(ps.propO.euc, method="NMDS", distance="euclidean")
# 
# ps.propOI.euc <- transform_sample_counts(ps.OI, function(count_tab) count_tab/sum(count_tab))
# ord.nmds.eucOI <- ordinate(ps.propOI.euc, method="NMDS", distance="euclidean")
# 
# 
# # NMDS plots (euclidean)
# 
# pdf("NMDS_euc_all.pdf")
# #jpeg("NMDS_euc_all.jpg")
# plot_ordination(ps.prop.euc, ord.nmds.euc, color="species", shape="dna_from", title="Euclidean NMDS with all (deseq)", label="sample") +
#   geom_point(size = 3)
# dev.off()
# 
# pdf("NMDS_euc_all2.pdf")
# #jpeg("NMDS_euc_all.jpg")
# plot_ordination(ps.prop.euc, ord.nmds.euc, color="species", shape="location", title="Euclidean NMDS with all (deseq)", label="sample") +
#   geom_point(size = 3)
# dev.off()
# 
# pdf("NMDS_euc_all3.pdf")
# #jpeg("NMDS_euc_all.jpg")
# plot_ordination(ps.prop.euc, ord.nmds.euc, color="species", shape="date", title="Euclidean NMDS with all (deseq)", label="sample") +
#   geom_point(size = 3) +
#   scale_shape_manual(values=seq(0,15))
# dev.off()
# 
# #ord <- metaMDS(ps_deseq, "euclidean")
# 
# pdf("NMDS_euc_fullbody.pdf")
# #jpeg("NMDS_euc_fullbody.jpg")
# plot_ordination(ps.propF.euc, ord.nmds.eucF, color="species", shape ="dna_from", title="Euclidean NMDS with full body (deseq)", label="sample") +
#   geom_point(size = 3)
# dev.off()
# 
# pdf("NMDS_euc_fullbody2.pdf")
# #jpeg("NMDS_euc_fullbody.jpg")
# plot_ordination(ps.propF.euc, ord.nmds.eucF, color="species", shape ="location", title="Euclidean NMDS with full body (deseq)", label="sample") +
#   geom_point(size = 3)
# dev.off()
# 
# pdf("NMDS_euc_fullbody3.pdf")
# #jpeg("NMDS_euc_fullbody.jpg")
# plot_ordination(ps.propF.euc, ord.nmds.eucF, color="species", shape ="date", title="Euclidean NMDS with full body (deseq)", label="sample") +
#   geom_point(size = 3) +
#   scale_shape_manual(values=seq(0,15))
# dev.off()
# 
# pdf("NMDS_euc_intestine.pdf")
# #jpeg("NMDS_euc_intestine.jpg")
# plot_ordination(ps.propI.euc, ord.nmds.eucI, color="species", shape="dna_from", title="Euclidean NMDS with intestine (deseq)", label="sample") +
#   geom_point(size = 3)
# dev.off()
# 
# pdf("NMDS_euc_intestine2.pdf")
# #jpeg("NMDS_euc_intestine.jpg")
# plot_ordination(ps.propI.euc, ord.nmds.eucI, color="species", shape="location", title="Euclidean NMDS with intestine (deseq)", label="sample") +
#   geom_point(size = 3)
# dev.off()
# 
# pdf("NMDS_euc_intestine3.pdf")
# #jpeg("NMDS_euc_intestine.jpg")
# plot_ordination(ps.propI.euc, ord.nmds.eucI, color="species", shape="date", title="Euclidean NMDS with intestine (deseq)", label="sample") +
#   geom_point(size = 3) +
#   scale_shape_manual(values=seq(0,15))
# dev.off()
# 
# pdf("NMDS_euc_ovary.pdf")
# #jpeg("NMDS_euc_ovary.jpg")
# plot_ordination(ps.propO.euc, ord.nmds.eucO, color="species", shape = "dna_from", title="Euclidean NMDS with ovary (deseq)", label="sample") +
#   geom_point(size = 3)
# dev.off()
# 
# pdf("NMDS_euc_ovary2.pdf")
# #jpeg("NMDS_euc_ovary.jpg")
# plot_ordination(ps.propO.euc, ord.nmds.eucO, color="species", shape = "location", title="Euclidean NMDS with ovary (deseq)", label="sample") +
#   geom_point(size = 3)
# dev.off()
# 
# pdf("NMDS_euc_ovary3.pdf")
# #jpeg("NMDS_euc_ovary.jpg")
# plot_ordination(ps.propO.euc, ord.nmds.eucO, color="species", shape = "date", title="Euclidean NMDS with ovary (deseq)", label="sample") +
#   geom_point(size = 3) +
#   scale_shape_manual(values=seq(0,15))
# dev.off()
# 
# pdf("NMDS_euc_ovary_intestine.pdf")
# #jpeg("NMDS_euc_ovary_intestine.jpg")
# plot_ordination(ps.propOI.euc, ord.nmds.eucOI, color="species", shape="dna_from", title="Euclidean NMDS with ovary+intestine (deseq)", label="sample")+
#   geom_point(size = 3)
# dev.off()
# 
# pdf("NMDS_euc_ovary_intestine2.pdf")
# #jpeg("NMDS_euc_ovary_intestine.jpg")
# plot_ordination(ps.propOI.euc, ord.nmds.eucOI, color="species", shape="location", title="Euclidean NMDS with ovary+intestine (deseq)", label="sample")+
#   geom_point(size = 3)
# dev.off()
# 
# pdf("NMDS_euc_ovary_intestine3.pdf")
# #jpeg("NMDS_euc_ovary_intestine.jpg")
# plot_ordination(ps.propOI.euc, ord.nmds.eucOI, color="species", shape="date", title="Euclidean NMDS with ovary+intestine (deseq)", label="sample")+
#   geom_point(size = 3) +
#   scale_shape_manual(values=seq(0,15))
# dev.off()






#--------------------------------------------------------------------------------------------#
#-----------------------------------NMMDS (bray distance)------------------------------------#
#--------------------------------------------------------------------------------------------#

# phyloseq object creation for each condition
ps_percent <- subset_samples(ps_percent, Sample != "S175")
ps.F <- subset_samples(ps_percent, Organ == "Full" | Organ =="Pull")
ps.I <- subset_samples(ps_percent, Organ == "Intestine")
ps.O <- subset_samples(ps_percent, Organ == "Ovary")
ps.OI <- subset_samples(ps_percent, Organ != "Full")
ps.OI <- subset_samples(ps.OI, Organ != "Salivary gland")
ps.OI <- subset_samples(ps.OI, Organ != "Pull")
ps.OI <- subset_samples(ps.OI, Organ != "Blank") # ovaries + intestine
ps.ICE <- subset_samples(ps_percent, Location == "Camping Europe" & Organ =="Intestine")
ps.ICE <- subset_samples(ps.ICE, Date =="30/05/2017" | Date =="28/06/2017")
ps.OCE <- subset_samples(ps_percent, Location == "Camping Europe" & Organ =="Ovary")
ps.OCE <- subset_samples(ps.OCE, Date =="30/05/2017" | Date =="28/06/2017")
ps.I2 <- subset_samples(ps.I, Sample !="NP17" & Sample !="NP20" & Sample !="S81" & Sample != "S82")
ps.O2 <- subset_samples(ps.O, Sample !="NP17" & Sample !="NP20" & Sample !="S81" & Sample != "S82" & Sample !="S103")

# otu_table transformation and ordination
# ps.prop <- transform_sample_counts(ps_percent, function(count_tab) count_tab/sum(count_tab))
# ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

ps.propF <- transform_sample_counts(ps.F, function(count_tab) count_tab/sum(count_tab))
ord.nmds.brayF <- ordinate(ps.propF, method="NMDS", distance="bray")

ps.propI <- transform_sample_counts(ps.I, function(count_tab) count_tab/sum(count_tab))
ord.nmds.brayI <- ordinate(ps.propI, method="NMDS", distance="bray")

ps.propO <- transform_sample_counts(ps.O, function(count_tab) count_tab/sum(count_tab))
ord.nmds.brayO <- ordinate(ps.propO, method="NMDS", distance="bray")

ps.propOI <- transform_sample_counts(ps.OI, function(count_tab) count_tab/sum(count_tab))
ord.nmds.brayOI <- ordinate(ps.propOI, method="NMDS", distance="bray")

ps.propICE <- transform_sample_counts(ps.ICE, function(count_tab) count_tab/sum(count_tab))
ord.nmds.brayICE <- ordinate(ps.propCE, method="NMDS", distance="bray")

ps.propOCE <- transform_sample_counts(ps.OCE, function(count_tab) count_tab/sum(count_tab))
ord.nmds.brayOCE <- ordinate(ps.propOCE, method="NMDS", distance="bray")

ps.propI2 <- transform_sample_counts(ps.I2, function(count_tab) count_tab/sum(count_tab))
ord.nmds.brayI2 <- ordinate(ps.propI2, method="NMDS", distance="bray")

ps.propO2 <- transform_sample_counts(ps.O2, function(count_tab) count_tab/sum(count_tab))
ord.nmds.brayO2 <- ordinate(ps.propO2, method="NMDS", distance="bray")


# NMDS plots (bray)

# pdf("NMDS_bray_all.pdf")
# #jpeg("NMDS_bray_all.jpg")
# plot_ordination(ps.prop, ord.nmds.bray, color="species", shape="dna_from", title="Bray NMDS with all (%)", label="sample") +
#   geom_point(size = 3)
# dev.off()
# 
# pdf("NMDS_bray_all2.pdf")
# #jpeg("NMDS_bray_all.jpg")
# plot_ordination(ps.prop, ord.nmds.bray, color="species", shape="location", title="Bray NMDS with all (%)", label="sample") +
#   geom_point(size = 3)
# dev.off()
# 
# pdf("NMDS_bray_all3.pdf")
# #jpeg("NMDS_bray_all.jpg")
# plot_ordination(ps.prop, ord.nmds.bray, color="species", shape="date", title="Bray NMDS with all (%)", label="sample") +
#   geom_point(size = 3) +
#   scale_shape_manual(values=seq(0,15))
# dev.off()

pdf("NMDS_bray_fullbody.pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(ps.propF, ord.nmds.brayF, color="Species", shape="Organ", title="Bray NMDS with full body (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_bray_fullbody2.pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(ps.propF, ord.nmds.brayF, color="Species", shape="Location", title="Bray NMDS with full body (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_bray_fullbody3.pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(ps.propF, ord.nmds.brayF, color="Species", shape="Date", title="Bray NMDS with full body (%)", label="Sample") +
  geom_point(size = 3) +
  scale_shape_manual(values=seq(0,15))
dev.off()

pdf("NMDS_bray_intestine.pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(ps.propI, ord.nmds.brayI, color="Species", shape="Organ", title="Bray NMDS with intestine (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_bray_intestine2.pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(ps.propI, ord.nmds.brayI, color="Species", shape="Location", title="Bray NMDS with intestine (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_bray_intestine3.pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(ps.propI, ord.nmds.brayI, color="Species", shape="Date", title="Bray NMDS with intestine (%)", label="Sample") +
  geom_point(size = 3) +
  scale_shape_manual(values=seq(0,15))
dev.off()

pdf("NMDS_bray_intestine4(ICE).pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(ps.propICE, ord.nmds.brayICE, color="Date", shape="Species", title="Bray NMDS with intestine (%)", label="Sample") +
  geom_point(size = 1) +
  scale_shape_manual(values=seq(0,15))
dev.off()

pdf("NMDS_bray_intestine4(ICE2).pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(ps.propI2, ord.nmds.brayI2, color="Date", shape="Species", title="Bray NMDS with intestine (%)", label="Sample") +
  geom_point(size = 1) +
  scale_shape_manual(values=seq(0,15))
dev.off()

pdf("NMDS_bray_ovary.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(ps.propO, ord.nmds.brayO, color="Species", shape="Organ", title="Bray NMDS with ovary (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_bray_ovary2.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(ps.propO, ord.nmds.brayO, color="Species", shape="Location", title="Bray NMDS with ovary (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_bray_ovary3.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(ps.propO, ord.nmds.brayO, color="Species", shape="Date", title="Bray NMDS with ovary (%)", label="Sample") +
  geom_point(size = 3) +
  scale_shape_manual(values=seq(0,15))
dev.off()

pdf("NMDS_bray_ovary4(O2).pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(ps.propO2, ord.nmds.brayO2, color="Date", shape="Species", title="Bray NMDS with ovary (%)", label="Sample") +
  geom_point(size = 1) +
  scale_shape_manual(values=seq(0,15))
dev.off()

pdf("NMDS_bray_ovary4(OCE2).pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(ps.propOCE2, ord.nmds.brayOCE2, color="Date", shape="Species", title="Bray NMDS with ovary (%)", label="Sample") +
  geom_point(size = 1) +
  scale_shape_manual(values=seq(0,15))
dev.off()

pdf("NMDS_bray_ovary_intestine.pdf")
#jpeg("NMDS_bray_ovary_intestine.jpg")
plot_ordination(ps.propOI, ord.nmds.brayOI, color="Species", shape="Organ", title="Bray NMDS with ovary+intestine (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_bray_ovary_intestine2.pdf")
#jpeg("NMDS_bray_ovary_intestine.jpg")
plot_ordination(ps.propOI, ord.nmds.brayOI, color="Species", shape="Location", title="Bray NMDS with ovary+intestine (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_bray_ovary_intestine3.pdf")
#jpeg("NMDS_bray_ovary_intestine.jpg")
plot_ordination(ps.propOI, ord.nmds.brayOI, color="Species", shape="Date", title="Bray NMDS with ovary+intestine (%)", label="Sample") +
  geom_point(size = 3) +
  scale_shape_manual(values=seq(0,15))
dev.off()




### NOTE ### 
# Bray-Curtis distance doesn't seem applicable for our data, contrary to euclidean distance
# Why ? I don't really know but I think each data depends on one type of distance for NMDS
# I think we need to choose a distance and show that our results are robust with this distance 






#--------------------------------------------------------------------------------------------#
#-----------------------------------NMMDS (jaccard distance)------------------------------------#
#--------------------------------------------------------------------------------------------#

ps.propF <- transform_sample_counts(ps.F, function(count_tab) count_tab/sum(count_tab))
ord.nmds.jaccF <- ordinate(ps.propF, method="NMDS", distance="jaccard")

ps.propI <- transform_sample_counts(ps.I, function(count_tab) count_tab/sum(count_tab))
ord.nmds.jaccI <- ordinate(ps.propI, method="NMDS", distance="jaccard")

ps.propO <- transform_sample_counts(ps.O, function(count_tab) count_tab/sum(count_tab))
ord.nmds.jaccO <- ordinate(ps.propO, method="NMDS", distance="jaccard")

ps.propOI <- transform_sample_counts(ps.OI, function(count_tab) count_tab/sum(count_tab))
ord.nmds.jaccOI <- ordinate(ps.propOI, method="NMDS", distance="jaccard")


pdf("NMDS_jacc_fullbody.pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(ps.propF, ord.nmds.jaccF, color="Species", shape="Organ", title="Jaccard NMDS with full body (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_jacc_fullbody2.pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(ps.propF, ord.nmds.jaccF, color="Species", shape="Location", title="Jaccard NMDS with full body (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_jacc_fullbody3.pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(ps.propF, ord.nmds.jaccF, color="Species", shape="Date", title="Jaccard NMDS with full body (%)", label="Sample") +
  geom_point(size = 3) +
  scale_shape_manual(values=seq(0,15))
dev.off()

pdf("NMDS_jacc_intestine.pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(ps.propI, ord.nmds.jaccI, color="Species", shape="Organ", title="Jaccard NMDS with intestine (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_jacc_intestine2.pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(ps.propI, ord.nmds.jaccI, color="Species", shape="Location", title="Jaccard NMDS with intestine (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_jacc_intestine3.pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(ps.propI, ord.nmds.jaccI, color="Species", shape="Date", title="Jaccard NMDS with intestine (%)", label="Sample") +
  geom_point(size = 3) +
  scale_shape_manual(values=seq(0,15))
dev.off()

pdf("NMDS_jacc_ovary.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(ps.propO, ord.nmds.jaccO, color="Species", shape="Organ", title="Jaccard NMDS with ovary (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_bray_ovary2.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(ps.propO, ord.nmds.jaccO, color="Species", shape="Location", title="Jaccard NMDS with ovary (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_jacc_ovary3.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(ps.propO, ord.nmds.jaccO, color="Species", shape="Date", title="Jaccard NMDS with ovary (%)", label="Sample") +
  geom_point(size = 3) +
  scale_shape_manual(values=seq(0,15))
dev.off()

pdf("NMDS_jacc_ovary_intestine.pdf")
#jpeg("NMDS_bray_ovary_intestine.jpg")
plot_ordination(ps.propOI, ord.nmds.jaccOI, color="Species", shape="Organ", title="Jaccard NMDS with ovary+intestine (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_jacc_ovary_intestine2.pdf")
#jpeg("NMDS_bray_ovary_intestine.jpg")
plot_ordination(ps.propOI, ord.nmds.jaccOI, color="Species", shape="Location", title="Jaccard NMDS with ovary+intestine (%)", label="Sample") +
  geom_point(size = 3)
dev.off()

pdf("NMDS_jacc_ovary_intestine3.pdf")
#jpeg("NMDS_bray_ovary_intestine.jpg")
plot_ordination(ps.propOI, ord.nmds.jaccOI, color="Species", shape="Date", title="Jaccard NMDS with ovary+intestine (%)", label="Sample") +
  geom_point(size = 3) +
  scale_shape_manual(values=seq(0,15))
dev.off()





#--------------------------------------------------------------------------------------------#
#-------------------------------------------ADONIS-------------------------------------------#
#--------------------------------------------------------------------------------------------#

adonis(vegdist(t(otu_table(ps_percent)), method = "bray") ~ dna_from*location*date,
       data=as(sample_data(ps_percent), "data.frame"), permutation = 9999)

adonis(vegdist(t(otu_table(ps_percent)), method = "euclidean") ~ location,
       data=as(sample_data(ps_percent), "data.frame"), permutation = 9999)

### NOTE : ###
# Adonis is used to test the difference of abundance between groups

#           Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
#location    5     88105 17621.0  8.2206 0.1624  1e-04 ***
#  Residuals 212    454427  2143.5         0.8376           
#Total     217    542531                 1.0000           
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# R2 = 0.1624 means that location could explain the distance variation 
# p is very significative, that means the observed difference is unlikely to be due to chance




adonis(vegdist(t(otu_table(ps_deseq)), method = "euclidean") ~ dna_from,
       data=as(sample_data(ps_deseq), "data.frame"), permutation = 9999)


adonis(vegdist(t(otu_table(ps_deseq)), method = "euclidean") ~ date,
       data=as(sample_data(ps_deseq), "data.frame"), permutation = 9999)



adonis(vegdist(t(otu_table(ps_percent)), method = "bray") ~ dna_from,
       data=as(sample_data(ps_percent), "data.frame"), permutation = 9999)

adonis(vegdist(t(otu_table(ps_percent)), method = "bray") ~ location,
       data=as(sample_data(ps_percent), "data.frame"), permutation = 9999)

adonis(vegdist(t(otu_table(ps_percent)), method = "bray") ~ date,
       data=as(sample_data(ps_percent), "data.frame"), permutation = 9999)





TukeyHSD_Observed <- TukeyHSD(aov(Observed ~ dna_from, data =  rich.plus))
TukeyHSD_Observed_df <- data.frame(TukeyHSD_Observed$dna_from)
TukeyHSD_Observed_df$measure = "Observed"
TukeyHSD_Observed_df$shapiro_test_pval = (shapiro.test(residuals(aov(Observed ~ dna_from, data =  rich.plus))))$p.value
TukeyHSD_Observed_df


library(session); packageVersion("session")
save.session("plot.Rda")



