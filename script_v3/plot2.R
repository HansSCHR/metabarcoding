# SCHRIEKE Hans
# PLOT
# input : phyloseq objects from decontam step
# output : plot and stats 



#--------------------------------------------------------------------------------------------#
#--------------------------------------PATH SETTING------------------------------------------#
#--------------------------------------------------------------------------------------------#

path = "D:/stage/data/runs_new2"
setwd(path)

dir.create("plotv5") # folder for plot
path2 <- "D:/stage/data/runs_new2/plotv5"


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
library("plotly")
library("hrbrthemes")
library("grid")

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


theme_set(theme_bw()) #set ggplot2 graphic theme 



#--------------------------------------------------------------------------------------------#
#--------------------------------------LOAD OBJECTS------------------------------------------#
#--------------------------------------------------------------------------------------------#
load("objects.Rdata")

# ps (raw data)
# ps_decontam (data after decontam, without blanks)
# ps_decontam2 (data after decontam and after removing remaining mitochondria, chloroplast)
# ps_percent (data after percent normalization)
# ps_wolbachia (only wolbachia)
# ps_proteo (only proteobacteria)


setwd(path2)
#--------------------------------------------------------------------------------------------#
#--------------------------------------DISTRIBUTION------------------------------------------#
#--------------------------------------------------------------------------------------------#

readsumsdf <- data.frame(nreads = sort(taxa_sums(ps), TRUE),
                         sorted = 1:ntaxa(ps), 
                         type = "OTU")

ggplot(readsumsdf, 
       aes(x = sorted, y = nreads)) + 
  geom_bar(stat = "identity") + 
  scale_y_log10() 

readsumsdf2 <- data.frame(nreads = sort(sample_sums(ps), TRUE), 
                          sorted = 1:nsamples(ps), 
                          type = "Samples")

readsumsdf3 <- rbind(readsumsdf,readsumsdf2)

p  <-  ggplot(readsumsdf3, 
              aes(x = sorted, y = nreads)) + 
  geom_bar(stat = "identity")+
  theme_gray()

pdf("distribution.pdf")
p + ggtitle("Total number of reads before Preprocessing") + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
dev.off()




#--------------------------------------------------------------------------------------------#
#-----------------------------------RAREFACTION CURVE----------------------------------------#
#--------------------------------------------------------------------------------------------#
p1 <- ggrare(ps_decontam2,
             step = 500,
             #color = "red",
             plot = T,
             parallel = F,
             se = T)

# p2 <- p1 + 
#   facet_wrap(~ Organ) + 
#   geom_vline(xintercept = min(sample_sums(ps)), 
#              color = "gray60") +
#   ggtitle("Rarefaction curve") +
#   xlim(0,100000) +
#   ylim(0, 230) +
#   theme_gray()

p2 <- p1 + 
  facet_wrap(~ Organ) + 
  geom_vline(xintercept = min(sample_sums(ps)), 
             color = "gray60") +
  xlim(0,100000) +
  ylim(0, 230) +
  labs(title = "Suffisant observations have been made for sampling",
       caption = "Rarefaction curve",
       x = "Sample Size", y = "Species Richness")

pdf("rarecurve_decontam2.pdf")
plot(p2)
dev.off()



#--------------------------------------------------------------------------------------------#
#-------------------------------------ALPHA DIVERSITY----------------------------------------#
#--------------------------------------------------------------------------------------------#

data1 <-  filter_taxa(ps_decontam2,
                      function(x) sum(x >= 10) > (1),
                      prune =  TRUE) 


p1 <- plot_richness(data1, 
                    x="Sample", 
                    color="Species", 
                    measures=c("Observed","Shannon", "Chao1"), 
                    nrow = 1) +
  ggtitle("Alpha diversity")
# pdf("richness_ps_decontam2_species.pdf")
# print(p1)
# dev.off()

p2 <- plot_richness(data1, 
                    x="Sample", 
                    color="Organ", 
                    measures=c("Observed","Shannon", "Chao1"), 
                    nrow = 1) +
  ggtitle("Alpha diversity")


pdf("richness_decontam2_species.pdf")
ggplot(p1$data,aes(Species,value,colour=Species,shape=Species)) +
  facet_grid(variable ~ Species, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.8, 
               position = position_dodge(width=0.9)) + 
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  labs(title = "A greater diversity in Guadeloupe",
       caption = "Alpha diversity", y = "Diversity index")
dev.off()

pdf("richness_decontam2_organ.pdf")
ggplot(p1$data,aes(Organ,value,colour=Species,shape=Organ)) +
  facet_grid(variable ~ Organ, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.8, 
               position = position_dodge(width=0.9)) + 
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  labs(title = "A greater diversity in whole organisms",
       caption = "Alpha diversity", y = "Diversity index")
dev.off()




# Culex pipiens - lab vs field
ps_pipiens <- subset_samples(ps_decontam2, Species=="Culex pipiens" & Organ!="Salivary gland")
data4 <-  filter_taxa(ps_pipiens, 
                      function(x) sum(x >= 10) > (1), 
                      prune =  TRUE) 

p1 <- plot_richness(data4, 
                    x="Sample", 
                    color="Field", 
                    measures=c("Observed","Shannon", "Chao1"), 
                    nrow = 1) +
  ggtitle("Alpha diversity")


pdf("richness_pipiens_organs.pdf")
ggplot(p1$data,aes(Species,value,colour=Field)) +
  facet_grid(variable ~ Organ, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.8, 
               position = position_dodge(width=0.9)) + 
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  labs(title = expression(paste("A greater diversity in Ovary and Whole from Lab within ", italic("Culex pipiens"))),
       caption = "Alpha diversity", y = "Diversity index")
dev.off()




# Culex pipiens labo vs Culex quinque field
ps_quinque_fieldvslabo <- subset_samples(ps_decontam2, Species=="Culex quinquefasciatus" & Organ=="Full")

data2 <-  filter_taxa(ps_quinque_fieldvslabo, 
                      function(x) sum(x >= 10) > (1), 
                      prune =  TRUE) 

p1 <- plot_richness(data2, 
                    x="Sample", 
                    color="Field", 
                    measures=c("Observed","Shannon", "Chao1"), 
                    nrow = 1) +
  ggtitle("Alpha diversity")

pdf("richness_quinque_full_labofield.pdf")
ggplot(p1$data,aes(Species,value,colour=Field,shape=Location)) +
  facet_grid(variable ~ Species, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.8, 
               position = position_dodge(width=0.9)) + 
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  labs(title = expression(paste("A greater diversity in Lab within Whole ", italic("Culex quinquefasciatus")," samples")),
       caption = "Alpha diversity", y = "Diversity index")
dev.off()






# # Whole organim - Field *
# ps_decontam2_fullfield <- subset_samples(ps_decontam2, Organ=="Full" & Field=="Field")
# data2 <-  filter_taxa(ps_decontam2_fullfield, 
#                       function(x) sum(x >= 10) > (1), 
#                       prune =  TRUE) 
# p1 <- plot_richness(data2, 
#                     x="Sample", 
#                     color="Species", 
#                     measures=c("Observed","Shannon", "Chao1"), 
#                     nrow = 1) +
#   ggtitle("Alpha diversity")
# # pdf("richness_ps_decontam2_species.pdf")
# # print(p1)
# # dev.off()
# 
# 
# 
# pdf("richness_decontam2_fullfield.pdf")
# ggplot(p1$data,aes(Species,value,colour=Species,shape=Species)) +
#   facet_grid(variable ~ Species, drop=T,scale="free",space="fixed") +
#   geom_boxplot(outlier.colour = NA,alpha=0.8, 
#                position = position_dodge(width=0.9)) + 
#   geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
#   labs(title = "A greater diversity in Guadeloupe within samples from field",
#        caption = "Alpha diversity", y = "Diversity index")
# dev.off()





# # Whole organim - Labo 
# ps_decontam2_labfull <- subset_samples(ps_decontam2, Location=="Lavar" | Location=="Labo Tetracycline" & Organ=="Full")
# data3 <-  filter_taxa(ps_decontam2_labfull, 
#                       function(x) sum(x >= 10) > (1), 
#                       prune =  TRUE) 
# p1 <- plot_richness(data3, 
#                     x="Sample", 
#                     color="Species", 
#                     measures=c("Observed","Shannon", "Chao1"), 
#                     nrow = 1) +
#   ggtitle("Alpha diversity")
# # pdf("richness_ps_decontam2_species.pdf")
# # print(p1)
# # dev.off()
# 
# 
# 
# pdf("richness_decontam2_labfull.pdf")
# ggplot(p1$data,aes(Species,value,colour=Species,shape=Species)) +
#   facet_grid(variable ~ Species, drop=T,scale="free",space="fixed") +
#   geom_boxplot(outlier.colour = NA,alpha=0.8, 
#                position = position_dodge(width=0.9)) + 
#   geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
#   labs(title = "A greater diversity in Guadeloupe within samples from labo",
#        caption = "Alpha diversity", y = "Diversity index")
# dev.off()






# # Ovary, Intestine, Salivary gland - Field
# ps_decontam2_organsfield <- subset_samples(ps_decontam2, Field=="Field" & Organ!="Full")
# data4 <-  filter_taxa(ps_decontam2_organsfield, 
#                       function(x) sum(x >= 10) > (1), 
#                       prune =  TRUE) 
# 
# p1 <- plot_richness(data4, 
#                     x="Sample", 
#                     color="Species", 
#                     measures=c("Observed","Shannon", "Chao1"), 
#                     nrow = 1) +
#   ggtitle("Alpha diversity")
# # pdf("richness_ps_decontam2_species.pdf")
# # print(p1)
# # dev.off()
# 
# p2 <- plot_richness(data4, 
#                     x="Sample", 
#                     color="Organ", 
#                     measures=c("Observed","Shannon", "Chao1"), 
#                     nrow = 1) +
#   ggtitle("Alpha diversity")
# 
# 
# pdf("richness_decontam2_organsfield1.pdf")
# ggplot(p1$data,aes(Species,value,colour=Species,shape=Species)) +
#   facet_grid(variable ~ Species, drop=T,scale="free",space="fixed") +
#   geom_boxplot(outlier.colour = NA,alpha=0.8, 
#                position = position_dodge(width=0.9)) + 
#   geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
#   labs(title = "A greater diversity in Guadeloupe for organs within samples from field",
#        caption = "Alpha diversity", y = "Diversity index")
# dev.off()
# 
# pdf("richness_decontam2_organsfield2.pdf")
# ggplot(p1$data,aes(Species,value,colour=Species,shape=Species)) +
#   facet_grid(variable ~ Organ, drop=T,scale="free",space="fixed") +
#   geom_boxplot(outlier.colour = NA,alpha=0.8, 
#                position = position_dodge(width=0.9)) + 
#   geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
#   labs(title = "A greater diversity for ovary of Culex quinquefasciatus within samples from field",
#        caption = "Alpha diversity", y = "Diversity index")
# dev.off()
# 
# 
# # Ovary, Intestine, Salivary gland - Labo
# ps_decontam2_organslab <- subset_samples(ps_decontam2, Location=="Lavar" | Location=="Labo Tetracycline")
# ps_decontam2_organslab <- subset_samples(ps_decontam2_organslab, Organ!="Full")
# data5 <-  filter_taxa(ps_decontam2_organslab, 
#                       function(x) sum(x >= 10) > (1), 
#                       prune =  TRUE) 
# 
# p1 <- plot_richness(data5, 
#                     x="Sample", 
#                     color="Species", 
#                     measures=c("Observed","Shannon", "Chao1"), 
#                     nrow = 1) +
#   ggtitle("Alpha diversity")
# # pdf("richness_ps_decontam2_species.pdf")
# # print(p1)
# # dev.off()
# 
# p2 <- plot_richness(data5, 
#                     x="Sample", 
#                     color="Organ", 
#                     measures=c("Observed","Shannon", "Chao1"), 
#                     nrow = 1) +
#   ggtitle("Alpha diversity")
# 
# 
# # pdf("richness_decontam2_organslab1.pdf")
# # ggplot(p1$data,aes(Species,value,colour=Species,shape=Species)) +
# #   facet_grid(variable ~ Species, drop=T,scale="free",space="fixed") +
# #   geom_boxplot(outlier.colour = NA,alpha=0.8, 
# #                position = position_dodge(width=0.9)) + 
# #   geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
# #   labs(title = "A greater diversity in Guadeloupe within samples from labo",
# #        caption = "Alpha diversity", y = "Diversity index")
# # dev.off()
# 
# pdf("richness_decontam2_organslab2.pdf")
# ggplot(p1$data,aes(Species,value,colour=Species,shape=Species)) +
#   facet_grid(variable ~ Organ, drop=T,scale="free",space="fixed") +
#   geom_boxplot(outlier.colour = NA,alpha=0.8, 
#                position = position_dodge(width=0.9)) + 
#   geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
#   labs(title = "Same diversity for ovary and intestine from labo Culex pipiens",
#        caption = "Alpha diversity", y = "Diversity index")
# dev.off()

# data4 <-  filter_taxa(ps_wolbachia, 
#                       function(x) sum(x >= 10) > (1), 
#                       prune =  TRUE) 
# 
# p1 <- plot_richness(data4, 
#                     x="Sample", 
#                     color="Field", 
#                     measures=c("Observed","Shannon", "Chao1"), 
#                     nrow = 1) +
#   ggtitle("Alpha diversity")
# 
# 
# pdf("richness_pipiens_organs.pdf")
# ggplot(p1$data,aes(Species,value,colour=Field,)) +
#   facet_grid(variable ~ Field, drop=T,scale="free",space="fixed") +
#   geom_boxplot(outlier.colour = NA,alpha=0.8, 
#                position = position_dodge(width=0.9)) + 
#   geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
#   labs(title = "A greater diversity in lab for Culex quinquefasciatus organs",
#        caption = "Alpha diversity", y = "Diversity index")
# dev.off()


#--------------------------------------------------------------------------------------------#
#----------------------------------TAXONOMIC SUMMARIES---------------------------------------#
#--------------------------------------------------------------------------------------------#

# ps_percent
p <- plot_composition(ps_percent,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Phylum", 
                      numberOfTaxa = 20, 
                      fill= "Phylum") +
  facet_wrap(~ Species, scales = "free_x", nrow = 3) + 
  labs(title = "Proteobacteria is the dominant phylum",
       caption = "Taxonomic composition (20 most abundant phylum)", x="Sample", y = "Abundance") +
  theme_gray()

pdf("composition_percent_species.pdf")
plot(p)
dev.off()

p <- plot_composition(ps_percent,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Phylum", 
                      numberOfTaxa = 20, 
                      fill= "Phylum") +
  facet_wrap(~ Organ, scales = "free_x", nrow = 3) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Proteobacteria is the dominant phylum in organs",
       caption = "Taxonomic composition (20 most abundant phylum)", x="Sample", y = "Abundance") +
  theme_gray()

pdf("composition_percent_organ.pdf")
plot(p)
dev.off()


# ps_proteo 
p <- plot_composition(ps_proteo,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Phylum", 
                      numberOfTaxa = 10, 
                      fill= "Class") +
  facet_wrap(~ Species, scales = "free_x", nrow = 3) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Alphaproteobacteria is the dominant class",
       caption = "Taxonomic composition (class)", x="Sample", y = "Abundance")+
  theme_gray()

pdf("composition_proteobacteria_species.pdf")
plot(p)
dev.off()

p <- plot_composition(ps_proteo,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Phylum", 
                      numberOfTaxa = 10, 
                      fill= "Class") +
  facet_wrap(~ Organ, scales = "free_x", nrow = 3) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Alphaproteobacteria is the dominant class in organs",
       caption = "Taxonomic composition (class)", x="Sample", y = "Abundance")+
  theme_gray()

pdf("composition_proteobacteria_organ.pdf")
plot(p)
dev.off()


# gammaproteo
p <- plot_composition(ps_proteo,
                      taxaRank1 = "Class",
                      taxaSet1 ="Gammaproteobacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 15, 
                      fill= "Genus") +
  facet_wrap(~ Species, scales = "free_x", nrow = 5) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Klebsiella and Erwinia are specific of Culex species in Gammaproteobacteria",
       caption = "Taxonomic composition (15 most abundant class)", x="Sample", y = "Abundance")+
  theme_gray()

pdf("composition_gammaproteobacteria_species.pdf")
plot(p)
dev.off()

p <- plot_composition(ps_proteo,
                      taxaRank1 = "Class",
                      taxaSet1 ="Gammaproteobacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 15, 
                      fill= "Genus") +
  facet_wrap(~ Organ, scales = "free_x", nrow = 5) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Klebsiella and Erwinia are specific of intestine and whole organisms",
       caption = "Taxonomic composition (15 most abundant class)", x="Sample", y = "Abundance") +
  theme_gray()

pdf("composition_gammaproteobacteria_organ.pdf")
plot(p)
dev.off()



# alphaproteo
p <- plot_composition(ps_proteo,
                      taxaRank1 = "Class",
                      taxaSet1 ="Alphaproteobacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 8, 
                      fill= "Genus") +
  facet_wrap(~ Species, scales = "free_x", nrow = 5) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Wolbachia is the dominant class within Alphaproteobacteria",
       caption = "Taxonomic composition (10 most abundant class)", x="Sample", y = "Abundance") +
  theme_gray()

pdf("composition_alphaproteobacteria_species.pdf")
plot(p)
dev.off()



p <- plot_composition(ps_proteo,
                      taxaRank1 = "Class",
                      taxaSet1 ="Alphaproteobacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 10, 
                      fill= "Genus") +
  facet_wrap(~ "Ovary", scales = "free_x", nrow=5) + 
  labs(title = "Wolbachia is the dominant class within Alphaproteobacteria in organism",
       caption = "Taxonomic composition (10 most abundant class)", x="Sample", y = "Abundance") +
  theme_gray()

pdf("composition_alphaproteobacteria_organ.pdf")
plot(p)
dev.off()



#deltaproteo
p <- plot_composition(ps_proteo,
                      taxaRank1 = "Class",
                      taxaSet1 ="Deltaproteobacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 15, 
                      fill= "Genus") +
  facet_wrap(~ Species, scales = "free_x", nrow = 5) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Taxonomic composition - Deltaproteobacteria",
       caption = "Taxonomic composition (8 most abundant class)", x="Sample", y = "Abundance")+
  theme_gray()

pdf("composition_deltaproteobacteria_species.pdf")
plot(p)
dev.off()

p <- plot_composition(ps_proteo,
                      taxaRank1 = "Class",
                      taxaSet1 ="Deltaproteobacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 20, 
                      fill= "Genus") +
  facet_wrap(~ Organ, scales = "free_x", nrow = 5) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Taxonomic composition - Deltaproteobacteria") +
  theme_gray()

pdf("composition_deltaproteobacteria_organ.pdf")
plot(p)
dev.off()

test <- as(tax_table(ps_proteo),"matrix")



#wolbachia
p <- plot_composition(ps_proteo,
                      taxaRank1 = "Genus",
                      taxaSet1 ="Wolbachia",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 15, 
                      fill= "Genus") +
  facet_wrap(~ Species, scales = "free_x", nrow = 5) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Taxonomic composition - Deltaproteobacteria",
       caption = "Taxonomic composition (8 most abundant class)", x="Sample", y = "Abundance")+
  theme_gray()

pdf("composition_deltaproteobacteria_species.pdf")
plot(p)
dev.off()

p <- plot_composition(ps_proteo,
                      taxaRank1 = "Genus",
                      taxaSet1 ="Wolbachia",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 20, 
                      fill= "Genus") +
  facet_wrap(~ Organ, scales = "free_x", nrow = 5) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Taxonomic composition - Deltaproteobacteria") +
  theme_gray()

pdf("composition_deltaproteobacteria_organ.pdf")
plot(p)
dev.off()

test <- as(tax_table(ps_proteo),"matrix")





#--------------------------------------------------------------------------------------------#
#------------------------------------------PCOA----------------------------------------------#
#--------------------------------------------------------------------------------------------#
pdf("PCoA_percent_field1.pdf")
plot_ordination(ps_percent, ordinate(ps_percent, method ="MDS", distance = "bray"), color = "Species", shape="Field") +
  geom_point(size = 4) + 
   labs(title = "Species have specific bacterial community structure",
        caption = "Bray PCoA", x="Axis.1 [37.5%]", y = "Axis.2 [10.2%]") +
  theme_gray()
dev.off()

# with tree
# plot_ordination(ps_percent, ordinate(ps_percent, method ="MDS", distance = "unifrac"), color = "Species", shape="Field") +
#   geom_point(size = 4) +
#   ggtitle("Bray PCoA - Labo vs Field") +
#   theme_gray()

pdf("PCoA_percent_field2.pdf")
plot_ordination(ps_percent, ordinate(ps_percent, method ="MDS", distance = "bray"), color = "Field", shape="Species") +
  geom_point(size = 4) +
  labs(title = "Distinction between laboratory and field samples",
       caption = "Bray PCoA", x="Axis.1 [37.5%]", y = "Axis.2 [10.2%]") +
  theme_gray()
dev.off()







#--------------------------------------------------------------------------------------------#
#-----------------------------------NMMDS (bray distance)------------------------------------#
#--------------------------------------------------------------------------------------------#

# phyloseq objects with the different conditions
ps_percent <- subset_samples(ps_percent, Sample != "S175")


full <- subset_samples(ps_percent, Organ == "Full" | Organ == "Pool")
full_no_aedes <- subset_samples(full, Species!="Aedes aegypti")
full_no_labo <- subset_samples(full_no_aedes, Location!="Labo Tetracycline" & Location!="Lavar")

intestine <- subset_samples(ps_percent, Organ == "Intestine")
intestine_camping <- subset_samples(ps_percent, Location == "Camping Europe" & Organ =="Intestine")
intestine_camping_date <- subset_samples(intestine_camping, Date =="30/05/2017" | Date =="28/06/2017")
intestine_no_lavar <- subset_samples(intestine, Location != "Lavar")
intestine_filter <- subset_samples(intestine_camping_date, Sample !="NP17" & Sample !="NP20" & Sample !="S81" & Sample != "S82")

ovary <- subset_samples(ps_percent, Organ == "Ovary")
ovary_culex_camping <- subset_samples(ovary, Species=="Culex pipiens" & Location == "Camping Europe")
ovary_culex_camping_date <- subset_samples(ovary_culex_camping, Date =="30/05/2017" | Date =="28/06/2017")
ovary_filter <- subset_samples(ovary_culex_camping_date, Sample !="S103")



# NMDS of entire organism
prop.full <- transform_sample_counts(full, function(count_tab) count_tab/sum(count_tab))
bray.full <- ordinate(full, method="NMDS", distance="bray")

prop.full_no_aedes <- transform_sample_counts(full_no_aedes, function(count_tab) count_tab/sum(count_tab))
bray.full_no_aedes <- ordinate(full_no_aedes, method="NMDS", distance="bray")

prop.full_no_labo <- transform_sample_counts(full_no_labo, function(count_tab) count_tab/sum(count_tab))
bray.full_no_labo <- ordinate(full_no_labo, method="NMDS", distance="bray")


pdf("NMDS_bray_full(without full vs pool).pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.full, bray.full, color="Species", title="Bray NMDS with full body", label="Sample") +
  labs(title = "Does species influence bacterial community structure ? ",
       caption = "Bray NMDS", x="NMDS1", y = "NMDS2")+
  geom_point(size = 4) +
  theme_gray()
dev.off()

pdf("NMDS_bray_full(with full vs pool).pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.full, bray.full, color="Organ", shape="Species", title="Bray NMDS with full body - Full vs Pool", label="Sample") +
  labs(title = "Whole and pooled organisms don't show any difference",
       caption = "Bray NMDS", x="NMDS1", y = "NMDS2") +
  geom_point(size = 4) +
  theme_gray()
dev.off()

# pdf("NMDS_bray_full(with aedes).pdf")
# #jpeg("NMDS_bray_fullbody.jpg")
# plot_ordination(prop.full, bray.full, color="Species", shape="Location", title="Bray NMDS with full body - Location", label="Sample") +
#   geom_point(size = 4) +
#   theme_gray()
# dev.off()

pdf("NMDS_bray_full(without aedes).pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.full_no_aedes, bray.full_no_aedes, color="Field", shape="Location", title="Bray NMDS with full body - Location without Aedes aegypti", label="Sample") +
  labs(title = "Do antibiotics influence microbiote ? ",
       caption = "Bray NMDS", x="NMDS1", y = "NMDS2") +
  geom_point(size = 4) +
  theme_gray()
dev.off()

pdf("NMDS_bray_full(without aedes)2.pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.full_no_aedes, bray.full_no_aedes, color="Field", shape="Species", title="Bray NMDS with full body - Location without Aedes aegypti", label="Sample") +
  labs(title = "Does laboratory influence microbiote ?",
       caption = "Bray NMDS", x="NMDS1", y = "NMDS2") +
  geom_point(size = 4) +
  theme_gray()
dev.off()

# pdf("NMDS_bray_full(without aedes and with field).pdf")
# #jpeg("NMDS_bray_fullbody.jpg")
# plot_ordination(prop.full_no_aedes, bray.full_no_aedes, color="Field", shape="Location", title="Bray NMDS with full body - Labo vs Field", label="Sample") +
#   geom_point(size = 4) +
#   theme_gray()
# dev.off()

# pdf("NMDS_bray_full(france vs gwada with aedes and labo).pdf")
# #jpeg("NMDS_bray_fullbody.jpg")
# plot_ordination(prop.full, bray.full, color="Country", shape="Location", title="Bray NMDS with full body - France vs Guadeloupe", label="Sample") +
#   geom_point(size = 4) +
#   scale_shape_manual(values=seq(0,15)) +
#   theme_gray()
# dev.off()

pdf("NMDS_bray_full(france vs gwada without aedes and labo).pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.full_no_labo, bray.full_no_labo, color="Country", shape="Location", title="Bray NMDS with full body - France vs Guadeloupe without Labo and Aedes aegypti", label="Sample") +
  labs(title = "Distinction between France and Guadeloupe",
       caption = "Bray NMDS", x="NMDS1", y = "NMDS2") +
  geom_point(size = 4) +
  scale_shape_manual(values=seq(0,15)) +
  theme_gray()
dev.off()



# NMDS of intestine
prop.intestine <- transform_sample_counts(intestine, function(count_tab) count_tab/sum(count_tab))
bray.intestine <- ordinate(intestine, method="NMDS", distance="bray")

prop.intestine_camping_date <- transform_sample_counts(intestine_camping_date, function(count_tab) count_tab/sum(count_tab))
bray.intestine_camping_date <- ordinate(intestine_camping_date, method="NMDS", distance="bray")

prop.intestine_filter <- transform_sample_counts(intestine_filter, function(count_tab) count_tab/sum(count_tab))
bray.intestine_filter <- ordinate(intestine_filter, method="NMDS", distance="bray")

prop.intestine_no_lavar <- transform_sample_counts(intestine_no_lavar, function(count_tab) count_tab/sum(count_tab))
bray.intestine_no_lavar <- ordinate(intestine_no_lavar, method="NMDS", distance="bray")

pdf("NMDS_bray_intestine(camping europe, date).pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(prop.intestine_camping_date, bray.intestine_camping_date, color="Date", title="Bray NMDS with intestine of Culex pipiens in Camping Europe", label="Sample")+
  labs(title = expression(paste("Does date influence the microbiote of ", italic("Culex pipiens"),"?")),
     caption = "Bray NMDS on intestine at Camping Europe", x="NMDS1", y = "NMDS2") +
  geom_point(size = 4) +
  theme_gray()
dev.off()

#voir quels samples sont écartés et les enlever pour faire ce plot
# pdf("NMDS_bray_intestine.pdf")
# #jpeg("NMDS_bray_intestine.jpg")
# plot_ordination(prop.intestine, bray.intestine, color="Date", shape="Location", title="Bray NMDS with intestine - Date", label="Sample") +
#   geom_point(size = 4) +
#   theme_gray()
# dev.off()
# 
# pdf("NMDS_bray_intestine_no_lavar.pdf")
# #jpeg("NMDS_bray_intestine.jpg")
# plot_ordination(prop.intestine_no_lavar, bray.intestine_no_lavar, color="Date", shape="Location", title="Bray NMDS with intestine - Date without Lavar", label="Sample") +
#   geom_point(size = 4) +
#   theme_gray()
# dev.off()



# NMDS of ovary
prop.ovary <- transform_sample_counts(ovary, function(count_tab) count_tab/sum(count_tab))
bray.ovary <- ordinate(ovary, method="NMDS", distance="bray")

prop.ovary_culex_camping_date <- transform_sample_counts(ovary_culex_camping_date, function(count_tab) count_tab/sum(count_tab))
bray.ovary_culex_camping_date <- ordinate(ovary_culex_camping_date, method="NMDS", distance="bray")

prop.ovary_filter <- transform_sample_counts(ovary_filter, function(count_tab) count_tab/sum(count_tab))
bray.ovary_filter <- ordinate(ovary_filter, method="NMDS", distance="bray")

pdf("NMDS_bray_ovary.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(prop.ovary, bray.ovary, color="Species", shape="Date", title="Bray NMDS with ovary", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()

# S103 = outlier
pdf("NMDS_bray_ovary_camping_culex_date.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(prop.ovary_culex_camping_date, bray.ovary_culex_camping_date, color="Date", title="Bray NMDS with ovary - Culex, Camping Europe, Dates", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()


pdf("NMDS_bray_ovary_filter.pdf") 
#jpeg("NMDS_bray_ovary.jpg") 
plot_ordination(prop.ovary_filter, bray.ovary_filter, color="Date", title="Bray NMDS with ovary - Culex, Camping Europe, Dates", label="Sample") +
  labs(title = expression(paste("Does date influence the microbiote of ", italic("Culex pipiens"),"?")),
       caption = "Bray NMDS on ovary at Camping Europe", x="NMDS1", y = "NMDS2") +
  geom_point(size = 4) +
  theme_gray()
dev.off()







# NMDS of Wolbachia 
ps_wolbachia_ovary_intestine <- subset_samples(ps_wolbachia, Organ=="Intestine" | Organ=="Ovary")
ps_wolbachia_ovary_intestine <- subset_samples(ps_wolbachia, Individuals!="0")
ps_wolbachia_filter <- subset_samples(ps_wolbachia_ovary_intestine, Sample != "S175" & Sample != "S68" & Sample!="S99")
# ps_wolbachia_filter <- subset_samples(ps_wolbachia_ovary_intestine, Sample != "S175" & Sample != "NP38" & Sample!="S99" & Sample!="NP22" & Sample!="NP10" & Sample!="S68" &
#                                         Sample!="NP2" & Sample!="NP11" & Sample!="NP8" & Sample!="NP5" & Sample!="S176")

prop.wolbachia <- transform_sample_counts(ps_wolbachia_ovary_intestine, function(count_tab) count_tab/sum(count_tab))
bray.wolbachia <- ordinate(ps_wolbachia_ovary_intestine, method="NMDS", distance="bray")

prop.wolbachia_filter <- transform_sample_counts(ps_wolbachia_filter, function(count_tab) count_tab/sum(count_tab))
bray.wolbachia_filter <- ordinate(ps_wolbachia_filter, method="NMDS", distance="bray")


pdf("NMDS_bray_wolbachia.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(prop.wolbachia, bray.wolbachia, color="Individuals", shape="Organ", title="Bray NMDS - Wolbachia", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off() # outliers = S68 S175 S99

pdf("NMDS_bray_wolbachia_filter.pdf") # removing outliers 
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(prop.wolbachia_filter, bray.wolbachia_filter, color="Individuals", shape="Organ", title="Bray NMDS - Wolbachia filter", label="Individuals") +
  labs(title = "Does structure of Wolbachia bacterial community depend on individual or organ?  ",
       caption = "Bray NMDS on Wolbachia", x="NMDS1", y = "NMDS2")+
  geom_point(size = 4) +
  theme_gray()
dev.off()

# pdf("NMDS_bray_wolbachia_filter2.pdf")
# #jpeg("NMDS_bray_ovary.jpg")
# plot_ordination(prop.wolbachia_filter, bray.wolbachia_filter, color="Individuals", shape="Location", title="Bray NMDS - Wolbachia filter", label="Sample") +
#   geom_point(size = 4) +
#   theme_gray()
# dev.off()



#--------------------------------------------------------------------------------------------#
#---------------------------------NMMDS (jaccard distance)-----------------------------------#
#--------------------------------------------------------------------------------------------#

# NMDS of whole organism
jacc.full <- ordinate(full, method="NMDS", distance="jaccard")

pdf("NMDS_jacc_full(without full vs pool).pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.full, jacc.full, color="Species", title="Jaccard NMDS with full body", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()

pdf("NMDS_jacc_full(with full vs pool).pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.full, jacc.full, color="Species", shape="Organ", title="Jacc NMDS with full body - Full vs Pool", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()


# NMDS of intestine 
jacc.intestine <- ordinate(intestine, method="NMDS", distance="jaccard")
jacc.intestine_camping_date <- ordinate(intestine_camping_date, method="NMDS", distance="jaccard")
jacc.intestine_filter <- ordinate(intestine_filter, method="NMDS", distance="jaccard")
jacc.intestine_no_lavar <- ordinate(intestine_no_lavar, method="NMDS", distance="jaccard")

pdf("NMDS_jacc_intestine(camping europe, date).pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(prop.intestine_camping_date, jacc.intestine_camping_date, color="Date", title="Jacc NMDS with intestine of Culex pipiens in Camping Europe", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()

#voir quels samples sont écartés et les enlever pour faire ce plot
pdf("NMDS_jacc_intestine.pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(prop.intestine, jacc.intestine, color="Date", shape="Location", title="Jacc NMDS with intestine - Date", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off

pdf("NMDS_jacc_intestine.pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(prop.intestine_no_lavar, jacc.intestine_no_lavar, color="Date", shape="Location", title="Jacc NMDS with intestine - Date without Lavar", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()



#--------------------------------------------------------------------------------------------#
#-------------------------------------------ADONIS-------------------------------------------#
#--------------------------------------------------------------------------------------------#

# ALL DATA
adonis(vegdist(t(otu_table(ps_percent)), method = "bray") ~ Organ*Location*Date*Species,
       data=as(sample_data(ps_percent), "data.frame"), permutation = 9999)

#                 Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Organ            3     6.240 2.08007 13.2550 0.12588 0.0001 ***
# Location         4     9.649 2.41228 15.3719 0.19464 0.0001 ***
# Date             4     2.406 0.60157  3.8335 0.04854 0.0001 ***
# Species          1     1.181 1.18073  7.5241 0.02382 0.0001 ***
# Organ:Location   7     3.819 0.54564  3.4770 0.07705 0.0001 ***
# Organ:Date       5     0.756 0.15120  0.9635 0.01525 0.4851    
# Organ:Species    3     0.571 0.19024  1.2123 0.01151 0.3159    
# Residuals      159    24.951 0.15693         0.50332           
# Total          186    49.574                 1.00000 


ps_whole <- subset_samples(ps_percent, Organ=="Full" & Species!="Aedes aegypti")
adonis(vegdist(t(otu_table(ps_whole)), method = "bray") ~ Location*Date*Species,
       data=as(sample_data(ps_whole), "data.frame"), permutation = 9999)

#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Location   4    9.0662 2.26655 15.0794 0.43306  1e-04 ***
# Date       3    1.9487 0.64957  4.3216 0.09308  1e-04 ***
# Residuals 66    9.9203 0.15031         0.47386           
# Total     73   20.9352                 1.00000           
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


ps_ovary_intestine_ce <- subset_samples(ps_percent, Organ=="Ovary" | Organ=="Intestine" & Location=="Camping Europe")
adonis(vegdist(t(otu_table(ps_ovary_intestine_ce)), method = "bray") ~ Organ*Date*Species,
       data=as(sample_data(ps_ovary_intestine_ce), "data.frame"), permutation = 9999)

#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Organ       1    1.4193 1.41929 13.9943 0.17241 0.0001 ***
#   Date        5    0.3270 0.06539  0.6448 0.03972 0.8757    
# Species     1    0.0058 0.00581  0.0573 0.00071 0.7109    
# Organ:Date  2    0.3950 0.19751  1.9475 0.04798 0.0513 .  
# Residuals  60    6.0851 0.10142         0.73919           
# Total      69    8.2322                 1.00000           
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



ps_intestine_ce <- subset_samples(ps_percent, Organ=="Intestine" & Location=="Camping Europe")
adonis(vegdist(t(otu_table(ps_intestine_ce)), method = "bray") ~ Date,
       data=as(sample_data(ps_intestine_ce), "data.frame"), permutation = 9999)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# Date       2    0.6318 0.31589  0.9882 0.09894 0.4027
# Residuals 18    5.7538 0.31966         0.90106       
# Total     20    6.3856                 1.00000 


ps_ovary_ce <- subset_samples(ps_percent, Organ=="Ovary" & Location=="Camping Europe")
adonis(vegdist(t(otu_table(ps_ovary_ce)), method = "bray") ~ Date,
       data=as(sample_data(ps_ovary_ce), "data.frame"), permutation = 9999)

#           Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)
# Date       2  0.013597 0.0067984 0.51535 0.05146 0.6096
# Residuals 19  0.250645 0.0131918         0.94854       
# Total     21  0.264241                   1.00000 



# WOLBACHIA
adonis(vegdist(t(otu_table(ps_wolbachia)), method = "bray") ~ Organ*Location*Individual,
       data=as(sample_data(ps_wolbachia), "data.frame"), permutation = 9999)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Organ              4    2.7847 0.69618  3.3493 0.12480 0.0121 *
# Location           4    1.5046 0.37615  1.8097 0.06743 0.1383  
# Individual        96   10.4837 0.10921  0.5254 0.46983 0.9927  
# Organ:Location     3    0.2119 0.07064  0.3398 0.00950 0.9380  
# Organ:Individual  32    2.7559 0.08612  0.4143 0.12351 0.9958  
# Residuals         22    4.5728 0.20786         0.20493         
# Total            161   22.3136                 1.00000         
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



adonis(vegdist(t(otu_table(ps_wolbachia_ovary_intestine)), method = "bray") ~Organ*Location*Individual,
       data=as(sample_data(ps_wolbachia_ovary_intestine), "data.frame"), permutation = 9999)


# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# Organ             1    0.9813 0.98131 223.367 0.12341 0.2203
# Location          2    0.0860 0.04301   9.789 0.01082 0.4727
# Individual       41    4.0549 0.09890  22.512 0.50992 0.3205
# Organ:Location    2    0.0695 0.03475   7.910 0.00874 0.5019
# Organ:Individual 32    2.7559 0.08612  19.603 0.34656 0.3447
# Residuals         1    0.0044 0.00439         0.00055       
# Total            79    7.9520                 1.00000 


adonis(vegdist(t(otu_table(ps_wolbachia_ovary_intestine)), method = "bray") ~Organ,
       data=as(sample_data(ps_wolbachia_ovary_intestine), "data.frame"), permutation = 9999)

#           Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)    
# Organ      1    0.9813 0.98131  10.981 0.12341  1e-04 ***
# Residuals 78    6.9706 0.08937         0.87659           
# Total     79    7.9520                 1.00000  


adonis(vegdist(t(otu_table(ps_wolbachia_ovary_intestine)), method = "bray") ~Individual,
       data=as(sample_data(ps_wolbachia_ovary_intestine), "data.frame"), permutation = 9999)

#            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# Individual 43    4.2967 0.099922 0.98411 0.54033 0.4975
# Residuals  36    3.6553 0.101536         0.45967       
# Total      79    7.9520                  1.00000 


adonis(vegdist(t(otu_table(ps_wolbachia_ovary_intestine)), method = "bray") ~Location,
       data=as(sample_data(ps_wolbachia_ovary_intestine), "data.frame"), permutation = 9999)

#            Df SumsOfSqs  MeanSqs F.Model    R2  Pr(>F)
# Location   2    0.0963 0.048154   0.472 0.01211 0.7976
# Residuals 77    7.8556 0.102021         0.98789       
# Total     79    7.9520                  1.00000 








# permanova on whole organism 
ps_whole <- subset_samples(ps_decontam2, Organ=="Full")

adonis(vegdist(t(otu_table(ps_whole)), method = "bray") ~Species*Location,
       data=as(sample_data(ps_whole), "data.frame"), permutation = 9999)

#           Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)    
# Species    2     6.644  3.3220  11.455 0.19675  1e-04 ***
# Location   3     4.504  1.5014   5.177 0.13338  1e-04 ***
# Residuals 78    22.621  0.2900         0.66987           
# Total     83    33.770                 1.00000           
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1





#--------------------------------------------------------------------------------------------#
#---------------------------------------Save session-----------------------------------------#
#--------------------------------------------------------------------------------------------#

#install.packages("session", repos = "http://cran.us.r-project.org")
library(session); packageVersion("session")
save.session("plot2.Rda")
#load("phyloseq.Rda")