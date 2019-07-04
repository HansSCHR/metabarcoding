# SCHRIEKE Hans
# PLOT
# input : phyloseq objects from decontam step
# output : plot and stats 



#--------------------------------------------------------------------------------------------#
#--------------------------------------PATH SETTING------------------------------------------#
#--------------------------------------------------------------------------------------------#

path = "D:/stage/data/runs_new2"
setwd(path)

dir.create("new_02.07.2019") # folder for plot
path2 <- "D:/stage/data/runs_new2/new_02.07.2019"


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

pdf("1-distribution.pdf")
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

pdf("2-rarecurve2_decontam2.pdf")
plot(p2)
dev.off()




#--------------------------------------------------------------------------------------------#
#-------------------------------------ALPHA DIVERSITY----------------------------------------#
#--------------------------------------------------------------------------------------------#

# Wolbachia - 
ps_new <- subset_samples(ps_decontam2, Organ=="Full")
ps_new <- subset_samples(ps_new, Organ!="Pool")

data4 <-  filter_taxa(ps_new, 
                      function(x) sum(x >= 10) > (1), 
                      prune =  TRUE) 

p1 <- plot_richness(data4, 
                    x="Sample", 
                    color="Location", 
                    measures=c("Observed","Shannon", "Chao1"), 
                    nrow = 1) +
  ggtitle("Alpha diversity")

pdf("3-alpha_diversity_full.pdf")
ggplot(p1$data,aes(Organ,value,colour=Location)) +
  facet_grid(variable ~ Species, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.8, 
               position = position_dodge(width=0.9)) + 
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  labs(title = expression(paste("A greater diversity for Aedes and lab samples")),
       caption = "Alpha diversity", y = "Diversity index")
dev.off()



# Bonus : Organs labo vs field
ps_pipiens <- subset_samples(ps_percent, Species=="Culex pipiens" & Organ!="Salivary gland")

data4 <-  filter_taxa(ps_pipiens, 
                      function(x) sum(x >= 10) > (1), 
                      prune =  TRUE) 

p1 <- plot_richness(data4, 
                    x="Sample", 
                    color="Field", 
                    measures=c("Observed","Shannon", "Chao1"), 
                    nrow = 1) +
  ggtitle("Alpha diversity")

pdf("3-bonus_labvsfield.pdf")
ggplot(p1$data,aes(Species,value,colour=Field)) +
  facet_grid(variable ~ Organ, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.8, 
               position = position_dodge(width=0.9)) + 
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  labs(title = expression(paste("A greater diversity in Ovary and Whole from Lab within ", italic("Culex pipiens"))),
       caption = "Alpha diversity", y = "Diversity index")
dev.off()





#--------------------------------------------------------------------------------------------#
#----------------------------------TAXONOMIC SUMMARIES---------------------------------------#
#--------------------------------------------------------------------------------------------#

# Phylum 
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

pdf("4-taxo_phylum_species.pdf")
plot(p)
dev.off()



# Proteobacteria 
ps_proteo_nopool <- subset_samples(ps_proteo, Organ!="Pool")
p <- plot_composition(ps_proteo_nopool,
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

pdf("5-taxo_proteo_species.pdf")
plot(p)
dev.off()

p <- plot_composition(ps_proteo_nopool,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Phylum", 
                      numberOfTaxa = 10, 
                      fill= "Class") +
  facet_wrap(~ Organ + Species, scales = "free_x", ncol=3) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Alphaproteobacteria is the dominant class in organs",
       caption = "Taxonomic composition (class)", x="Sample", y = "Abundance")


pdf("6-taxo_proteo_organ.pdf")
plot(p)
dev.off()


# Alphaproteobacteria 
p <- plot_composition(ps_proteo_nopool,
                      taxaRank1 = "Class",
                      taxaSet1 ="Alphaproteobacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 8, 
                      fill= "Genus") +
  facet_wrap(~ Species, scales = "free_x", nrow = 5) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(strip.text = element_text(face = "italic"))+
  labs(title = "Wolbachia is the dominant class within Alphaproteobacteria",
       caption = "Taxonomic composition (10 most abundant class)", x="Sample", y = "Abundance") +
  theme_gray()

pdf("7-taxo_alphaproteo_species.pdf")
plot(p)
dev.off()




p <- plot_composition(ps_proteo_nopool,
                      taxaRank1 = "Class",
                      taxaSet1 ="Alphaproteobacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 10, 
                      fill= "Genus") +
  facet_wrap(~ Organ + Species, scales = "free_x", nrow=5) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Wolbachia is the dominant genus within Alphaproteobacteria in organism",
       caption = "Taxonomic composition (10 most abundant class)", x="Sample", y = "Abundance")

pdf("8-taxo_alphaproteo_organ.pdf")
plot(p)
dev.off()





#--------------------------------------------------------------------------------------------#
#-----------------------------------NMMDS (bray distance)------------------------------------#
#--------------------------------------------------------------------------------------------#

# Whole

full <- subset_samples(ps_percent, Organ == "Full")

full_culex <- subset_samples(full, Species!="Aedes aegypti")
full_culex_field <- subset_samples(full_culex, Location!="Labo Tetracycline" & Location!="Lavar")

full_pipiens <- subset_samples(full, Species=="Culex pipiens")
full_pipiens_field <- subset_samples(full_pipiens, Field=="Field")

metadata_pipiens_full_field <- as(sample_data(full_pipiens_field),"matrix")


prop.full <- transform_sample_counts(full, function(count_tab) count_tab/sum(count_tab))
bray.full <- ordinate(full, method="NMDS", distance="bray")

prop.full_culex <- transform_sample_counts(full_culex, function(count_tab) count_tab/sum(count_tab))
bray.full_culex <- ordinate(full_culex, method="NMDS", distance="bray")

prop.full_culex_field <- transform_sample_counts(full_culex_field, function(count_tab) count_tab/sum(count_tab))
bray.full_culex_field <- ordinate(full_culex_field, method="NMDS", distance="bray")

prop.pipiens <- transform_sample_counts(full_pipiens, function(count_tab) count_tab/sum(count_tab))
bray.pipiens <- ordinate(full_pipiens, method="NMDS", distance="bray")


pdf("9-NMDS_bray_full.pdf")
plot_ordination(prop.full, bray.full, color="Species", title="Bray NMDS with full body", label="Sample") +
  labs(title = "Does species influence bacterial community structure ? ",
       caption = "Bray NMDS on whole mosquitoes", x="NMDS1", y = "NMDS2")+
  geom_point(size = 4) +
  theme_gray()
dev.off()

pdf("10-NMDS_bray_full_culex.pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.full_culex, bray.full_culex,color="Field", shape="Species", title="Bray NMDS with full body - Location without Aedes aegypti", label="Sample") +
  labs(title = "Do antibiotics influence microbiote ? ",
       caption = "Bray NMDS on whole Culex mosquitoes", x="NMDS1", y = "NMDS2") +
  geom_point(size = 4) +
  theme_gray()
dev.off()
