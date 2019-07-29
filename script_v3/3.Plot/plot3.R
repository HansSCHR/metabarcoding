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

dir.create("new_29.07.2019") # folder for plot
path2 <- "D:/stage/data/runs_new2/new_29.07.2019"


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
library("plyr")

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
sample_data(ps_decontam2)$Organ <- factor(sample_data(ps_decontam2)$Organ, 
                                          levels=c("Whole", "Pool", "Intestine", "Ovary", "Salivary gland"))

sample_data(ps_percent)$Organ <- factor(sample_data(ps_percent)$Organ, 
                                          levels=c("Whole", "Pool", "Intestine", "Ovary", "Salivary gland"))

sample_data(ps_proteo)$Organ <- factor(sample_data(ps_proteo)$Organ, 
                                        levels=c("Whole", "Pool", "Intestine", "Ovary", "Salivary gland"))

sample_data(ps_proteo)$Species <- factor(sample_data(ps_proteo)$Species, 
                                       levels=c("Culex pipiens", "Culex quinquefasciatus", "Aedes aegypti"))

ps_nopool <- subset_samples(ps_decontam2, Organ!="Pool")


p1 <- ggrare(ps_nopool,
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
ps_whole <- subset_samples(ps_decontam2, Organ=="Whole")
ps_whole <- subset_samples(ps_whole, Organ!="Pool")


data4 <-  filter_taxa(ps_whole, 
                      function(x) sum(x >= 10) > (1), 
                      prune =  TRUE) 

p1 <- plot_richness(data4, 
                    x="Sample", 
                    color="Location", 
                    measures=c("Observed","Shannon", "Chao1"), 
                    nrow = 1) +
  ggtitle("Alpha diversity")

pdf("3-alpha_diversity_Whole.pdf")
ggplot(p1$data,aes(Organ,value,colour=Location)) +
  facet_grid(variable ~ Species, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.8, 
               position = position_dodge(width=0.9)) + 
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  labs(title = expression(paste("A greater diversity for Aedes and lab samples")),
       caption = "Alpha diversity", y = "Diversity index")
dev.off()




ps_aedes <- subset_samples(ps_decontam2, Species!="Aedes aegypti" & Organ!="Pool")

data5 <-  filter_taxa(ps_aedes, 
                      function(x) sum(x >= 10) > (1), 
                      prune =  TRUE) 

p2 <- plot_richness(data5, 
                    x="Sample", 
                    color="Location", 
                    measures=c("Observed","Shannon", "Chao1"), 
                    nrow = 1) +
  ggtitle("Alpha diversity")


pdf("3bis.richness_decontam2_organ.pdf")
ggplot(p2$data,aes(Field,value,colour=Location,shape=Field)) +
  facet_grid(variable ~ Organ+Species, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.8,
               position = position_dodge(width=0.9)) +
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  labs(title = "A greater diversity in whole organisms",
       caption = "Alpha diversity", y = "Diversity index")
dev.off()





# Bonus : Organs labo vs field
# ps_pipiens <- subset_samples(ps_percent, Species=="Culex pipiens" & Organ!="Salivary gland")
# 
# data4 <-  filter_taxa(ps_pipiens, 
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
# pdf("3-bonus_labvsfield.pdf")
# ggplot(p1$data,aes(Species,value,colour=Field)) +
#   facet_grid(variable ~ Organ, drop=T,scale="free",space="fixed") +
#   geom_boxplot(outlier.colour = NA,alpha=0.8, 
#                position = position_dodge(width=0.9)) + 
#   geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
#   labs(title = expression(paste("A greater diversity in Ovary and Whole from Lab within ", italic("Culex pipiens"))),
#        caption = "Alpha diversity", y = "Diversity index")
# dev.off()





#--------------------------------------------------------------------------------------------#
#----------------------------------TAXONOMIC SUMMARIES---------------------------------------#
#--------------------------------------------------------------------------------------------#


# physeq2 = filter_taxa(ps_decontam2, function(x) mean(x) > 0.1, TRUE)
# physeq2
# physeq3 = transform_sample_counts(physeq2, function(x) x / sum(x) )
# physeq3
# 
# ps_nopool <- subset_samples(physeq3, Organ!="Pool")
# 
# nopool <- psmelt(ps_nopool)
# 
# glom <- tax_glom(ps_nopool, taxrank = 'Phylum')
# glom # should list # taxa as # phyla
# data <- psmelt(glom) # create dataframe from phyloseq object
# data$Phylum <- as.character(data$Phylum) #convert to character
# 
# #simple way to rename phyla with < 1% abundance
# data$Phylum[data$Abundance < 0.01] <- "< 1% abund."
# 
# medians <- ddply(data, ~Phylum, function(x) c(median=median(x$Abundance)))
# remainder <- medians[medians$median <= 0.01,]$Phylum
# remainder
# 
# data[data$Phylum %in% remainder,]$Phylum <- "Phyla < 1% abund."
# #rename phyla with < 1% relative abundance
# data$Phylum[data$Abundance < 0.01] <- "Phyla < 1% abund."
# 
# 
# p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=Phylum))
# p + geom_bar(aes(), stat="identity", position="stack") +
#   scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
#                                "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
#   theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))




# Phylum 
ps_percent_nopool <- subset_samples(ps_percent, Organ!="Pool")
p <- plot_composition(ps_percent_nopool,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 5, 
                      fill= "Genus") +
  facet_wrap(~ Organ + Species, scales = "free_x", ncol = 3) + 
  labs(title = "Proteobacteria is the dominant phylum",
       caption = "Taxonomic composition (20 most abundant phylum)", x="Sample", y = "Abundance")

pdf("4-taxo_phylum.pdf")
plot(p)
dev.off()



# Proteobacteria 
# 
# ps_proteo_nopool <- subset_samples(ps_proteo, Organ!="Pool")
# p <- plot_composition(ps_proteo_nopool,
#                       taxaRank1 = "Kingdom",
#                       taxaSet1 ="Bacteria",
#                       taxaRank2 = "Phylum", 
#                       numberOfTaxa = 10, 
#                       fill= "Class") +
#   facet_wrap(~ Organ + Species, scales = "free_x", ncol=3) + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   labs(title = "Alphaproteobacteria is the dominant class",
#        caption = "Taxonomic composition (class)", x="Sample", y = "Abundance")
# 
# pdf("5-taxo_proteo.pdf")
# plot(p)
# dev.off()
# 
# 
# 
# # Alphaproteobacteria 
# 
# p <- plot_composition(ps_proteo_nopool,
#                       taxaRank1 = "Class",
#                       taxaSet1 ="Alphaproteobacteria",
#                       taxaRank2 = "Genus", 
#                       numberOfTaxa = 10, 
#                       fill= "Genus") +
#   facet_wrap(~ Organ + Species, scales = "free_x", nrow=5) + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   labs(title = "Wolbachia is the dominant genus within Alphaproteobacteria in organism",
#        caption = "Taxonomic composition (10 most abundant class)", x="Sample", y = "Abundance")
# 
# pdf("6-taxo_alphaproteo.pdf")
# plot(p)
# dev.off()







# 20 genus most abundant 
ps_proteo_nopool <- subset_samples(ps_proteo, Organ!="Pool")
p <- plot_composition(ps_proteo_nopool,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 20, 
                      fill= "Genus") +
  facet_wrap(~ Species+Organ+Location, scales = "free_x", ncol=5) + 
  #facet_grid(Species~Organ) +
  theme(plot.title = element_text(hjust = 0.5), strip.text = element_text(size=8), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) + 
  labs(title = "Wolbachia is the dominant genus",
       caption = "Taxonomic composition (20 most abundant genus)", x="Sample", y = "Abundance")

pdf("7-taxo_20_genus.pdf")
plot(p)
dev.off()


# Without the 20 most abundant genus
ps_proteo_nopool2 <- subset_taxa(ps_proteo_nopool, Genus!="Wolbachia" & Genus!="Acinetobacter" & Genus!="Aeromonas" & Genus!="Asaia" &
                                    Genus!="Enhydrobacter" & Genus!="Erwinia" & Genus!="Haemophilus" & Genus!="Klebsiella" & Genus!="Legionella" &
                                    Genus!="Massilia" & Genus!="Morganella" & Genus!="Providencia" & Genus!="Pseudomonas" & Genus!="Rahnella" &
                                    Genus!="Ralstonia" & Genus!="Serratia" & Genus!="Sphingomonas" & Genus!="Thorsellia" & Genus!="Zymobacter")

p <- plot_composition(ps_proteo_nopool2,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 20, 
                      fill= "Genus") +
  facet_wrap(~ Species+Organ+Location, scales = "free_x", ncol=5) + 
  theme(plot.title = element_text(hjust = 0.5), strip.text = element_text(size=8), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) + 
  labs(title = "Wolbachia is the dominant genus",
       caption = "Taxonomic composition (20 most abundant genus)", x="Sample", y = "Abundance")

pdf("7bis-taxo_wt_20_genus.pdf")
plot(p)
dev.off()







#--------------------------------------------------------------------------------------------#
#-----------------------------------NMMDS (bray distance)------------------------------------#
#--------------------------------------------------------------------------------------------#

# Whole

ps_percent_whole <- subset_samples(ps_percent, Organ == "Whole")

ps_culex_whole <- subset_samples(ps_percent_whole, Species!="Aedes aegypti")
ps_culex_whole_field <- subset_samples(ps_culex_whole, Location!="Labo Tetracycline" & Location!="Lavar")

ps_pipiens_whole <- subset_samples(ps_percent_whole, Species=="Culex pipiens")
ps_pipiens_field <- subset_samples(ps_pipiens_whole, Field=="Field")

ps_proteo_nopool2 <- prune_taxa(taxa_sums(ps_proteo_nopool2) >= 1, ps_proteo_nopool2)
ps_proteo_nopool2 <- prune_samples(sample_sums(ps_proteo_nopool2) >= 1, ps_proteo_nopool2)

#metadata_pipiens_Whole_field <- as(sample_data(Whole_pipiens_field),"matrix")


prop.percent_whole <- transform_sample_counts(ps_percent_whole, function(count_tab) count_tab/sum(count_tab))
bray.percent_whole <- ordinate(ps_percent_whole, method="NMDS", distance="bray")

prop.culex_whole <- transform_sample_counts(ps_culex_whole, function(count_tab) count_tab/sum(count_tab))
bray.culex_whole <- ordinate(ps_culex_whole, method="NMDS", distance="bray")

prop.culex_whole_field <- transform_sample_counts(ps_culex_whole_field, function(count_tab) count_tab/sum(count_tab))
bray.Whole_culex_field <- ordinate(ps_culex_whole_field, method="NMDS", distance="bray")

prop.pipiens_whole <- transform_sample_counts(ps_pipiens_whole, function(count_tab) count_tab/sum(count_tab))
bray.pipiens_whole <- ordinate(ps_pipiens_whole, method="NMDS", distance="bray")
jacc.pipiens_whole <- ordinate(ps_pipiens_whole, method="NMDS", distance="jaccard")

prop.proteo_nopool2 <- transform_sample_counts(ps_proteo_nopool2, function(count_tab) count_tab/sum(count_tab))
bray.proteo_nopool2 <- ordinate(ps_proteo_nopool2, method="NMDS", distance="bray")

pdf("7-NMDS_bray_whole.pdf")
plot_ordination(prop.percent_whole, bray.percent_whole, color="Species", title="Bray NMDS with Whole body", label="Sample") +
  labs(title = "Does species influence bacterial community structure ? ",
       caption = "Bray NMDS on whole mosquitoes", x="NMDS1", y = "NMDS2")+
  geom_point(size = 4) +
  theme_gray()
dev.off()

# pdf("8-NMDS_bray_culex_whole.pdf")
# #jpeg("NMDS_bray_Wholebody.jpg")
# plot_ordination(prop.culex_whole, bray.culex_whole, color="Field", shape="Location", title="Bray NMDS with Whole body - Location without Aedes aegypti", label="Sample") +
#   labs(title = "Do antibiotics influence microbiote ? ",
#        caption = "Bray NMDS on whole Culex mosquitoes", x="NMDS1", y = "NMDS2") +
#   #stat_ellipse(geom = "polygon", level=0.6,alpha = 1/2, aes(fill = Field))+
#   scale_fill_manual(values=c("yellow","green"))+
#   scale_color_manual(values=c("red", "blue"))+
#   geom_point(size = 5) +
#   theme_gray()
# dev.off()


pdf("9-NMDS_bray_culex_whole_ellipse.pdf")
#jpeg("NMDS_bray_Wholebody.jpg")
  plot_ordination(prop.culex_whole, bray.culex_whole, color="Field", shape="Location", title="Bray NMDS with Whole body - Location without Aedes aegypti", label="Sample") +
  labs(title = "Do antibiotics influence microbiote ? ",
       caption = "Bray NMDS on whole Culex mosquitoes", x="NMDS1", y = "NMDS2") +
  stat_ellipse(geom = "polygon", level=0.70,alpha = 1/2, aes(fill = Species))+
  scale_fill_manual(values=c("yellow","green"))+
  scale_color_manual(values=c("red", "blue"))+
  geom_point(size = 5) +
  theme_gray()
dev.off()
  

plot_ordination(prop.proteo_nopool2, bray.proteo_nopool2, color="Organ", title="Bray NMDS with Whole body", label="Sample") +
  labs(title = "Does species influence bacterial community structure ? ",
       caption = "Bray NMDS on whole mosquitoes", x="NMDS1", y = "NMDS2")+
  geom_point(size = 4) +
  theme_gray()
  
plot_ordination(prop.proteo_nopool2, bray.proteo_nopool2, color="Field", shape="Location", title="Bray NMDS with Whole body - Location without Aedes aegypti", label="Sample") +
  labs(title = "Do antibiotics influence microbiote ? ",
       caption = "Bray NMDS on whole Culex mosquitoes", x="NMDS1", y = "NMDS2") +
  #stat_ellipse(geom = "polygon", level=0.70,alpha = 1/2, aes(fill = Species))+
  scale_fill_manual(values=c("yellow","green"))+
  scale_color_manual(values=c("red", "blue"))+
  geom_point(size = 5) +
  theme_gray()


# pdf("10-NMDS_bray_pipiens_whole.pdf")
# #jpeg("NMDS_bray_Wholebody.jpg")
# plot_ordination(prop.pipiens_whole, bray.pipiens_whole, color="Field", shape="Location", title="Bray NMDS with Whole body - Location without Aedes aegypti", label="Sample") +
#   labs(title = expression(paste("Does laboratory influence microbiote of ", italic("Culex pipiens"), "?")),
#        caption = "Bray NMDS", x="NMDS1", y = "NMDS2") +
#   geom_point(size = 4) +
#   theme_gray()
# dev.off()
# 
# 
# 
# pdf("11-NMDS_jaccard_pipiens_whole.pdf")
# #jpeg("NMDS_bray_Wholebody.jpg")
# plot_ordination(prop.pipiens_whole, jacc.pipiens_whole, color="Field", shape="Location", title="Jaccard NMDS with Whole body - Location without Aedes aegypti", label="Sample") +
#   labs(title = expression(paste("Does laboratory influence microbiote of ", italic("Culex pipiens"), "?")),
#        caption = "Jaccard NMDS", x="NMDS1", y = "NMDS2") +
#   geom_point(size = 4) +
#   theme_gray()
# dev.off()


# Organs of Culex pipiens from Cammping Europe 
ps_pipiens_organ <- subset_samples(ps_percent, Organ!="Whole" & Species=="Culex pipiens")
ps_pipiens_organ <- subset_samples(ps_pipiens_organ, Sample!="S175")
ps_pipiens_organ_camping <- subset_samples(ps_pipiens_organ, Location=="Camping Europe")
#test0 <- as(sample_data(ps_whole),"matrix")

prop.pipiens_organ_camping <- transform_sample_counts(ps_pipiens_organ_camping, function(count_tab) count_tab/sum(count_tab))
bray.pipiens_organ_camping <- ordinate(ps_pipiens_organ_camping, method="NMDS", distance="bray")
jacc.pipiens_organ_camping <- ordinate(ps_pipiens_organ_camping, method="NMDS", distance="jaccard")


pdf("11-NMDS_bray_organs_culex_CE.pdf")
plot_ordination(prop.pipiens_organ_camping, bray.pipiens_organ_camping, color="Organ", title="Bray NMDS", label="Sample") +
  labs(title = expression(paste("Organ influences the structure of ", italic("Culex pipiens"), " bacterial community")),
       caption = expression(paste("Bray NMDS on Ovary, Intestine and Salivary Gland of ", italic('Culex pipiens'), " from Camping Europe")), x="NMDS1", y = "NMDS2")+
  geom_point(size = 4) +
  theme_gray()
dev.off()

pdf("12-NMDS_jaccard_organs_culex_CE.pdf")
plot_ordination(prop.pipiens_organ_camping, bray.pipiens_organ_camping, color="Organ", title="Jaccard NMDS", label="Sample") +
  labs(title = expression(paste("Organ influences the structure of ", italic("Culex pipiens"), " bacterial community")),
       caption = expression(paste("Jaccard NMDS on Ovary, Intestine and Salivary Gland of ", italic('Culex pipiens'), " from Camping Europe")), x="NMDS1", y = "NMDS2")+
  geom_point(size = 4) +
  theme_gray()
dev.off()



# Intestine of Culex pipiens from Camping Europe at 2 dates

ps_intestine_camping <- subset_samples(ps_percent, Location == "Camping Europe" & Organ =="Intestine")
ps_intestine_camping_date <- subset_samples(ps_intestine_camping, Date =="30/05/2017" | Date =="28/06/2017")

prop.intestine_camping_date <- transform_sample_counts(ps_intestine_camping_date, function(count_tab) count_tab/sum(count_tab))
bray.intestine_camping_date <- ordinate(ps_intestine_camping_date, method="NMDS", distance="bray")

pdf("13-NMDS_bray_intestine_CE_date.pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(prop.intestine_camping_date, bray.intestine_camping_date, color="Date", title="Bray NMDS with intestine of Culex pipiens in Camping Europe", label="Sample")+
  labs(title = expression(paste("Does time of sampling influence the microbiote of ", italic("Culex pipiens"),"?")),
       caption = "Bray NMDS on intestine at Camping Europe", x="NMDS1", y = "NMDS2") +
  geom_point(size = 4) +
  theme_gray()
dev.off()



# Ovary of Culex pipiens from Camping Europe at 2 dates
ps_percent_ovary <- subset_samples(ps_percent, Organ == "Ovary")
ps_ovary_culex_camping <- subset_samples(ps_percent_ovary, Species=="Culex pipiens" & Location == "Camping Europe")
ps_ovary_culex_camping_date <- subset_samples(ps_ovary_culex_camping, Date =="30/05/2017" | Date =="28/06/2017")
ps_ovary_culex_camping_date <- subset_samples(ps_ovary_culex_camping_date, Sample!="S103")

prop.ovary_culex_camping_date <- transform_sample_counts(ps_ovary_culex_camping_date, function(count_tab) count_tab/sum(count_tab))
bray.ovary_culex_camping_date <- ordinate(ps_ovary_culex_camping_date, method="NMDS", distance="bray")

pdf("14-NMDS_bray_ovary_CE_date.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(prop.ovary_culex_camping_date, bray.ovary_culex_camping_date, color="Date", title="Bray NMDS with ovary - Culex, Camping Europe, Dates", label="Sample") +
  labs(title = expression(paste("Does time of sampling influence the microbiote of ", italic("Culex pipiens"),"?")),
       caption = "Bray NMDS on intestine at Camping Europe", x="NMDS1", y = "NMDS2") +
  geom_point(size = 4) +
  theme_gray()
dev.off()








#--------------------------------------------------------------------------------------------#
#--------------------------------------------HEATMAP-----------------------------------------#
#--------------------------------------------------------------------------------------------#

# Culex pipiens
ps_pipiens <- subset_samples(ps_percent, Species=="Culex pipiens")
ps_pipiens_wolbachia <- subset_taxa(ps_pipiens, Genus=="Wolbachia")
ps_pipiens_wolbachia <- prune_taxa(taxa_sums(ps_pipiens_wolbachia) >= 1, ps_pipiens_wolbachia)
ps_pipiens_wolbachia <- prune_samples(sample_sums(ps_pipiens_wolbachia) >= 1, ps_pipiens_wolbachia)
ps_pipiens_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_pipiens_wolbachia),TRUE)[1:30]), ps_pipiens_wolbachia)

jpeg("15-heatmap_pipiens_wolbachia.jpg", width = 1080, height = 720)
plot_heatmap(ps_pipiens_wolbachia, sample.label="Location", sample.order="Location", low="#000033", high="#00FF00", trans=identity_trans())+
  #scale_fill_gradient(low="#000033", high="#FF3300",breaks=c(100,64,32,0), labels=c("100","64","32","0"))+
  facet_wrap(~ Field + Organ, scales = "free_x", ncol = 3)+
  labs(title = expression(paste("ASV1 is abundant in almost all sequences of ", italic("Culex pipiens"))),
       caption = expression(paste("Heatmap of ", italic('Culex pipiens'), " that contains Wolbachia (log10 transformation of abundance")), x="Field", y = "ASV")
dev.off()




# Culex quinquefasciatus
ps_quinque <- subset_samples(ps_decontam2, Species=="Culex quinquefasciatus")
ps_quinque_wolbachia <- subset_taxa(ps_quinque, Genus=="Wolbachia")
ps_quinque_wolbachia <- prune_taxa(taxa_sums(ps_quinque_wolbachia) >= 1, ps_quinque_wolbachia)
ps_quinque_wolbachia <- prune_samples(sample_sums(ps_quinque_wolbachia) >= 1, ps_quinque_wolbachia)
ps_quinque_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_quinque_wolbachia),TRUE)[1:30]), ps_quinque_wolbachia)
otu_quinque_wolbachia <- as(otu_table(ps_quinque_wolbachia),"matrix")
metadata_quinque_wolbachia <- as(sample_data(ps_quinque_wolbachia),"matrix")
ps_quinque_wolbachia # 30 taxa / 32 samples
sum(as(otu_table(ps_quinque_wolbachia),"matrix")) # 5 406 449

ps_quinque_wolbachia_lab <- subset_samples(ps_quinque_wolbachia, Field=="Lab ")
tax_quinque_wolbachia_lab <- as(tax_table(ps_quinque_wolbachia_lab), "matrix")
ps_quinque_wolbachia_lab <- prune_taxa(taxa_sums(ps_quinque_wolbachia_lab) >= 1, ps_quinque_wolbachia_lab)
ps_quinque_wolbachia_lab <- prune_samples(sample_sums(ps_quinque_wolbachia_lab) >= 1, ps_quinque_wolbachia_lab)
ps_quinque_wolbachia_lab # 4 taxa / 13 samples
sum(as(otu_table(ps_quinque_wolbachia_lab),"matrix")) # 9 265 -> 0,17% of total reads

ps_quinque_wolbachia_field <- subset_samples(ps_quinque_wolbachia, Field=="Field")
ps_quinque_wolbachia_field <- prune_taxa(taxa_sums(ps_quinque_wolbachia_field) >= 1, ps_quinque_wolbachia_field)
ps_quinque_wolbachia_field <- prune_samples(sample_sums(ps_quinque_wolbachia_field) >= 1, ps_quinque_wolbachia_field)
sum(as(otu_table(ps_quinque_wolbachia_field),"matrix")) # 5 397 184 -> 99,83% of total reads
ps_quinque_wolbachia_field # 30 taxa / 19 samples

meta_test <- as(sample_data(ps_quinque_wolbachia), "matrix")
sample_data(ps_quinque_wolbachia)$Organ <- factor(sample_data(ps_quinque_wolbachia)$Organ, levels=c("Intestine", "Ovary", "Salivary gland","Whole", "Pool"))

taxa_test <- as(tax_table(ps_quinque_wolbachia),"matrix")
otu_test <- as(otu_table(ps_quinque_wolbachia),"matrix")

jpeg("16-heatmap_quinque_wolbachia.jpg", width = 1080, height = 720)
plot_heatmap(ps_quinque_wolbachia, sample.label="Location", sample.order="Location", low="#000033", high="#00FF00",trans=identity_trans())+
  #scale_fill_gradient2(low = grey, mid = "#FE1B00", high = "#00FF00", midpoint=625000)+
  #scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue", "black"),values = c(1250000, 1000000, 75000,50000, 25000,0))+
  #scale_fill_gradient(low="gray", high="#FF3300",breaks=c(1250000,1000000,0), labels=c("1250000","750000","0"))+
  facet_wrap(~ Organ, scales = "free_x", ncol = 3)+
  labs(title = expression(paste("ASV1 is abundant in almost all sequences of ", italic("Culex quinquefasciatus"), ", especially in ovary and pool samples")),
       caption = expression(paste("Heatmap of ", italic('Culex quinquefasciatus'), " that contains Wolbachia (log10 transformation of abundance")), x="Field", y = "ASV")
dev.off()

jpeg("16bis-heatmap_quinque_wolbachia.jpg", width = 1080, height = 720)
plot_heatmap(ps_quinque_wolbachia, sample.label="Location", sample.order="Location", low="#000033", high="#00FF00", trans=identity_trans())+
  #scale_fill_gradient(low="#000033", high="#FF3300",breaks=c(1250000,1000000,0), labels=c("1250000","750000","0"))+
  facet_wrap(~ Organ, scales = "free_x", ncol = 3)+
  labs(title = expression(paste("ASV1 is abundant in almost all sequences of ", italic("Culex quinquefasciatus"), ", especially in ovary and pool samples")),
       caption = expression(paste("Heatmap of ", italic('Culex quinquefasciatus'), " that contains Wolbachia (log10 transformation of abundance")), x="Field", y = "ASV")
dev.off()


# Aedes aegypti
ps_aedes <- subset_samples(ps_decontam2, Species=="Aedes aegypti" & Organ!="Pool")
ps_aedes_wolbachia <- subset_taxa(ps_aedes, Genus=="Wolbachia")
ps_aedes_wolbachia <- prune_taxa(taxa_sums(ps_aedes_wolbachia) >= 1, ps_aedes_wolbachia)
ps_aedes_wolbachia <- prune_samples(sample_sums(ps_aedes_wolbachia) >= 1, ps_aedes_wolbachia)
ps_aedes_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_aedes_wolbachia),TRUE)[1:30]), ps_aedes_wolbachia)

jpeg("17-heatmap_aedes_wolbachia.jpg", width = 1080, height = 720)
plot_heatmap(ps_aedes_wolbachia, sample.label="Location", sample.order="Location", low="#000033", high="#00FF00", trans=log_trans(10))+
  #scale_fill_gradient(low="#000033", high="#FF3300",breaks=c(1250000,1000000,0), labels=c("1250000","750000","0"))+
  facet_wrap(~ Organ, scales = "free_x", ncol = 2)+
  labs(title = expression(paste("ASV1 is abundant in almost all sequences of ", italic("Aedes aegypti"),", especially in ovary samples")),
       caption = expression(paste("Heatmap of  ", italic('Aedes aegypti'), " that contains Wolbachia (log10 transformation of abundance")), x="Field", y = "ASV")
dev.off()

jpeg("17bis-heatmap_aedes_wolbachia.jpg", width = 1080, height = 720)
plot_heatmap(ps_aedes_wolbachia, sample.label="Location", sample.order="Location", low="#000033", high="#00FF00", trans=identity_trans())+
  #scale_fill_gradient(low="#000033", high="#FF3300",breaks=c(1250000,1000000,0), labels=c("1250000","750000","0"))+
  facet_wrap(~ Organ, scales = "free_x", ncol = 2)+
  labs(title = expression(paste("ASV1 is abundant in almost all sequences of ", italic("Aedes aegypti"),", especially in ovary samples")),
       caption = expression(paste("Heatmap of  ", italic('Aedes aegypti'), " that contains Wolbachia (log10 transformation of abundance")), x="Field", y = "ASV")
dev.off()



# # Whole 
# ps_whole <- subset_samples(ps_decontam2, Organ=="Whole")
# ps_whole_wolbachia <- subset_taxa(ps_whole, Genus=="Wolbachia")
# ps_whole_wolbachia <- prune_taxa(taxa_sums(ps_whole_wolbachia) >= 1, ps_whole_wolbachia)
# ps_whole_wolbachia <- prune_samples(sample_sums(ps_whole_wolbachia) >= 1, ps_whole_wolbachia)
# ps_whole_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_whole_wolbachia),TRUE)[1:30]), ps_whole_wolbachia)
# 
# jpeg("18-heatmap_whole_wolbachia.jpg", width = 1080, height = 720)
# plot_heatmap(ps_whole_wolbachia, sample.label="Location", sample.order="Field", low="#000033", high="#FF3300")+
#   labs(title = expression(paste("ASV1 is abundant in almost all whole samples")),
#        caption = expression(paste("Heatmap of whole mosquitoes that contain Wolbachia")), x="Field", y = "ASV")
# dev.off()
# 
# 
# # Intestine
# ps_intestine <- subset_samples(ps_decontam2, Organ=="Intestine")
# ps_intestine_wolbachia <- subset_taxa(ps_intestine, Genus=="Wolbachia")
# ps_intestine_wolbachia <- prune_taxa(taxa_sums(ps_intestine_wolbachia) >= 1, ps_intestine_wolbachia)
# ps_intestine_wolbachia <- prune_samples(sample_sums(ps_intestine_wolbachia) >= 1, ps_intestine_wolbachia)
# ps_intestine_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_intestine_wolbachia),TRUE)[1:30]), ps_intestine_wolbachia)
# 
# jpeg("19-heatmap_intestine_wolbachia.jpg", width = 1080, height = 720)
# plot_heatmap(ps_intestine_wolbachia, sample.label="Location", sample.order="Field", low="#000033", high="#FF3300")+
#   labs(title = expression(paste("ASV1 is abundant in almost all intestine samples")),
#        caption = expression(paste("Heatmap of itestine samples that contain Wolbachia")), x="Field", y = "ASV")
# dev.off()
# 
# # Ovary 
# ps_ovary <- subset_samples(ps_decontam2, Organ=="Ovary")
# ps_ovary_wolbachia <- subset_taxa(ps_ovary, Genus=="Wolbachia")
# ps_ovary_wolbachia <- prune_taxa(taxa_sums(ps_ovary_wolbachia) >= 1, ps_ovary_wolbachia)
# ps_ovary_wolbachia <- prune_samples(sample_sums(ps_ovary_wolbachia) >= 1, ps_ovary_wolbachia)
# ps_ovary_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_ovary_wolbachia),TRUE)[1:30]), ps_ovary_wolbachia)
# 
# jpeg("20-heatmap_ovary_wolbachia.jpg", width = 1080, height = 720)
# plot_heatmap(ps_ovary_wolbachia, sample.label="Location", sample.order="Field", low="#000033", high="#FF3300")+
#   labs(title = expression(paste("ASV1 is abundant in almost all ovary samples")),
#        caption = expression(paste("Heatmap of ovary samples that contain Wolbachia")), x="Field", y = "ASV")
# dev.off()
# 
# # Salivary gland 
# ps_gland <- subset_samples(ps_decontam2, Organ=="Salivary gland")
# ps_gland_wolbachia <- subset_taxa(ps_gland, Genus=="Wolbachia")
# ps_gland_wolbachia <- prune_taxa(taxa_sums(ps_gland_wolbachia) >= 1, ps_gland_wolbachia)
# ps_gland_wolbachia <- prune_samples(sample_sums(ps_gland_wolbachia) >= 1, ps_gland_wolbachia)
# ps_gland_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_gland_wolbachia),TRUE)[1:30]), ps_gland_wolbachia)
# 
# jpeg("21-heatmap_gland_wolbachia.jpg", width = 1080, height = 720)
# plot_heatmap(ps_gland_wolbachia, sample.label="Location", sample.order="Field", low="#000033", high="#FF3300")+
#   labs(title = expression(paste("ASV1 is abundant in almost all salivary gland samples")),
#        caption = expression(paste("Heatmap of salivary gland samples that contain Wolbachia")), x="Field", y = "ASV")
# dev.off()
# 
# # Pool 
# ps_pool <- subset_samples(ps_decontam2, Organ=="Pool")
# ps_pool_wolbachia <- subset_taxa(ps_pool, Genus=="Wolbachia")
# ps_pool_wolbachia <- prune_taxa(taxa_sums(ps_pool_wolbachia) >= 1, ps_pool_wolbachia)
# ps_pool_wolbachia <- prune_samples(sample_sums(ps_pool_wolbachia) >= 1, ps_pool_wolbachia)
# ps_pool_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_pool_wolbachia),TRUE)[1:30]), ps_pool_wolbachia)
# 
# ps_pool 
# ps_pool_wolbachia 
# 
# jpeg("22-heatmap_pool_wolbachia.jpg", width = 1080, height = 720)
# plot_heatmap(ps_pool_wolbachia, sample.label="Location", sample.order="Field", low="#000033", high="#00FF00")+
#   labs(title = expression(paste("ASV1 is abundant in almost all pooled samples")),
#        caption = expression(paste("Heatmap of pooled samples that contain Wolbachia")), x="Field", y = "ASV")
# dev.off()





#--------------------------------------------------------------------------------------------#
#--------------------------------------WOLBACHIA STATS-------------------------------------#
#--------------------------------------------------------------------------------------------#


ps_pipiens_wolbachia_field <- subset_samples(ps_pipiens_wolbachia, Field=="Field")
ps_pipiens_wolbachia_field <- prune_taxa(taxa_sums(ps_pipiens_wolbachia_field) >= 1, ps_pipiens_wolbachia_field)
ps_pipiens_wolbachia_field <- prune_samples(sample_sums(ps_pipiens_wolbachia_field) >= 1, ps_pipiens_wolbachia_field)
tax_pipiens_wolbachia_field <- as(tax_table(ps_pipiens_wolbachia_field),"matrix")
ps_pipiens_wolbachia_field
# otu_table()   OTU Table:         [ 28 taxa and 89 samples ]
# sample_data() Sample Data:       [ 89 samples by 16 sample variables ]
# tax_table()   Taxonomy Table:    [ 28 taxa by 7 taxonomic ranks ]

ps_pipiens_wolbachia_labo <- subset_samples(ps_pipiens_wolbachia, Field=="Lab ")
ps_pipiens_wolbachia_labo <- prune_taxa(taxa_sums(ps_pipiens_wolbachia_labo) >= 1, ps_pipiens_wolbachia_labo)
ps_pipiens_wolbachia_labo <- prune_samples(sample_sums(ps_pipiens_wolbachia_labo) >= 1, ps_pipiens_wolbachia_labo)
tax_labo <- as(tax_table(ps_pipiens_wolbachia_labo),"matrix")
ps_pipiens_wolbachia_labo
# otu_table()   OTU Table:         [ 8 taxa and 44 samples ]
# sample_data() Sample Data:       [ 44 samples by 16 sample variables ]
# tax_table()   Taxonomy Table:    [ 8 taxa by 7 taxonomic ranks ]


tax_pipiens <- as(tax_table(ps_pipiens),"matrix")
otu_pipiens_wolbachia <- as(otu_table(ps_pipiens_wolbachia),"matrix")
meta_pipiens_wolbachia <- as(sample_data(ps_pipiens_wolbachia),"matrix")





#--------------------------------------------------------------------------------------------#
#-------------------------------------------ADONIS-------------------------------------------#
#--------------------------------------------------------------------------------------------#

# Difference between field and labo pipiens samples ? - Effect of Field ?

### Permanova - adonis

ps_pipiens_whole_bosc_lavar <- subset_samples(ps_pipiens_whole, Location!="Camping Europe")
ps_pipiens_whole_camping_lavar <- subset_samples(ps_pipiens_whole, Location!="Bosc")

adonis(vegdist(t(otu_table(ps_pipiens_whole_bosc_lavar)), method = "bray") ~Field,
       data=as(sample_data(ps_pipiens_whole_bosc_lavar), "data.frame"), permutation = 9999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Location   1    1.1851  1.1850  8.5438 0.19621  4e-04 ***
#   Residuals 35    4.8546  0.1387         0.80379           
# Total     36    6.0397                 1.00000           
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


adonis(vegdist(t(otu_table(ps_pipiens_whole_camping_lavar)), method = "bray") ~Location,
       data=as(sample_data(ps_pipiens_whole_camping_lavar), "data.frame"), permutation = 9999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Location   1    0.5469 0.54693  2.6655 0.08985  0.036 *
#   Residuals 27    5.5400 0.20519         0.91015         
# Total     28    6.0870                 1.00000         
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


### ANOSIM - vegan
meta_pipiens_bosc_lavar <- sample_data(ps_pipiens_bosc_lavar)
meta_pipiens_camping_lavar <- sample_data(ps_pipiens_ce_lavar)

anosim(vegdist(t(otu_table(ps_pipiens_bosc_lavar))), meta_pipiens_bosc_lavar$Location, permutations=1000)
# ANOSIM statistic R: 0.1413 
# Significance: 0.026973 

anosim(vegdist(t(otu_table(ps_pipiens_ce_lavar))), meta_pipiens_camping_lavar$Location, permutations=1000)
# ANOSIM statistic R: 0.1571 
# Significance: 0.12188 

# R > 0 --> there is a difference between samples depends on they are lab or field 
# Siginifiance < 0,005 --> difference is significative 




#Difference between whole Aedes aegypti and Culex quinquefasciatus from Guadeloupe - Effect of Species ?

### Permanova - adonis

ps__gwada <- subset_samples(ps_percent, Location=="Guadeloupe" | Field=="Field")
ps_whole_gwada <- subset_samples(ps_analysis, Organ=="Whole")
meta_whole_gwada <- as(sample_data(ps_analysis),"matrix")

adonis(vegdist(t(otu_table(ps_whole_gwada)), method = "bray") ~Species,
       data=as(sample_data(ps_whole_gwada), "data.frame"), permutation = 9999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Species    2    5.6424 2.82121  16.607 0.46639  1e-04 ***
#   Residuals 38    6.4557 0.16989         0.53361           
# Total     40   12.0981                 1.00000           
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


### ANOSIM - vegan
meta_whole_gwada <- sample_data(ps_whole_gwada)

anosim(vegdist(t(otu_table(ps_whole_gwada))), metadata_analysis$Species, permutations=1000)
# ANOSIM statistic R: 0.7323 
# Significance: 0.000999 




#Difference between whole field and wolbachia- Culex quinquefasciatus - Effect of Tetracycline ?

### Permanova - adonis
ps_quinque_whole <- subset_samples(ps_quinque, Organ=="Whole")

adonis(vegdist(t(otu_table(ps_quinque_whole)), method = "bray") ~Field,
       data=as(sample_data(ps_quinque_whole), "data.frame"), permutation = 9999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Field      1    2.1141 2.11409  6.4502 0.23498  1e-04 ***
#   Residuals 21    6.8828 0.32775         0.76502           
# Total     22    8.9969                 1.00000           
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


### ANOSIM - vegan
meta_quinque_whole <- sample_data(ps_quinque_whole)
anosim(vegdist(t(otu_table(ps_quinque_whole))), metadata_quinque_whole$Field, permutations=1000)
# ANOSIM statistic R: 0.9114 
# Significance: 0.000999 




# Does location influence the structure of whole Culex microbiote ? - Effect of Location ? (confirmed effect of Species too ?)


### Permanova - adonis

adonis(vegdist(t(otu_table(ps_culex_whole_field)), method = "bray") ~Location,
       data=as(sample_data(ps_culex_whole_field), "data.frame"), permutation = 9999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Location   2    2.3422 1.17111  7.3307 0.34367  2e-04 ***
#   Residuals 28    4.4731 0.15975         0.65633           
# Total     30    6.8153                 1.00000           
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

adonis(vegdist(t(otu_table(ps_culex_whole_field)), method = "bray") ~Species,
       data=as(sample_data(ps_culex_whole_field), "data.frame"), permutation = 9999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Species    1    2.0737  2.0737  12.683 0.30427  1e-04 ***
#   Residuals 29    4.7416  0.1635         0.69573           
# Total     30    6.8153                 1.00000           
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


adonis(vegdist(t(otu_table(ps_culex_whole_field)), method = "bray") ~Species+Location,
       data=as(sample_data(ps_culex_whole_field), "data.frame"), permutation = 9999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Species    1    2.0737 2.07372 12.9807 0.30427 0.0001 ***
#   Location   1    0.2685 0.26850  1.6807 0.03940 0.1417    
# Residuals 28    4.4731 0.15975         0.65633           
# Total     30    6.8153                 1.00000           
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



adonis(vegdist(t(otu_table(ps_culex)), method = "bray") ~ Location+Species,
       data=as(sample_data(ps_culex), "data.frame"), permutation = 9999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Location    4     8.902 2.22562  6.4612 0.13063  1e-04 ***
#   Residuals 172    59.247 0.34446         0.86937           
# Total     176    68.150                 1.00000           
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


### Anosim - vegan 
meta5 <- sample_data(ps_culex_whole_field)
anosim(vegdist(t(otu_table(ps_culex_whole_field))), meta5$Location, permutations=1000)
# ANOSIM statistic R: 0.4746 
# Significance: 0.000999

anosim(vegdist(t(otu_table(ps_culex_whole_field))), meta5$Species, permutations=1000)
# ANOSIM statistic R: 0.6285 
# Significance: 0.000999 




# Difference between culex samples from Camping Europe depends on two date - Effect of Time ? 

### Permanova - adonis
ps_pipiens_camping <- subset_samples(ps_pipiens, Location=="Camping Europe")
ps_pipiens_ovary_camping <- subset_samples(ps_pipiens_camping, Organ=="Ovary")
ps_pipiens_intestine_camping <- subset_samples(ps_pipiens_camping, Organ=="Intestine")
meta_pipiens_ovary_camping <- as(sample_data(ps_pipiens_ovary_camping),"matrix")
meta_pipiens_intestine_camping <- as(sample_data(ps_pipiens_intestine_camping),"matrix")

adonis(vegdist(t(otu_table(ps_pipiens_ovary_camping)), method = "bray") ~ Date,
       data=as(sample_data(ps_pipiens_ovary_cammping), "data.frame"), permutation = 9999)
# Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)
# Date       2  0.013597 0.0067984 0.51535 0.05146  0.611
# Residuals 19  0.250645 0.0131918         0.94854       
# Total     21  0.264241                   1.00000 


adonis(vegdist(t(otu_table(ps_pipiens_intestine_camping)), method = "bray") ~ Date,
       data=as(sample_data(ps_pipiens_intestine_camping), "data.frame"), permutation = 9999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# Date       2    0.6318 0.31589 0.98822 0.09894 0.3877
# Residuals 18    5.7539 0.31966         0.90106       
# Total     20    6.3856                 1.00000 


### ANOSIM - vegan
meta_pipiens_ovary_camping <- sample_data(ps_pipiens_ovary_camping)
anosim(vegdist(t(otu_table(ps_pipiens_ovary_camping))), meta_pipiens_ovary_camping$Date, permutations=1000)
# ANOSIM statistic R: -0.01842 
# Significance: 0.51449 


meta_pipiens_intestine_camping <- sample_data(ps_pipiens_intestine_camping)
anosim(vegdist(t(otu_table(ps_pipiens_intestine_camping))), meta_pipiens_intestine_camping$Date, permutations=1000)
# ANOSIM statistic R: 0.1452 
# Significance: 0.092907 




# Difference between culex pipiens samples depends on organs - Effect of organs ? 

### Permanova - adonis
adonis(vegdist(t(otu_table(ps_pipiens_organ)), method = "bray") ~ Organ,
       data=as(sample_data(ps_pipiens_organ), "data.frame"), permutation = 9999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Organ      2    2.1702 1.08512  6.6115 0.12218  1e-04 ***
#   Residuals 95   15.5919 0.16413         0.87782           
# Total     97   17.7622                 1.00000           
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


### ANOSIM - vegan
meta_pipiens_organ <- sample_data(ps_pipiens_organ)
anosim(vegdist(t(otu_table(ps_pipiens_organ))), meta_pipiens_organ$Organ, permutations=1000)
# ANOSIM statistic R: 0.1505 
# Significance: 0.000999


### Parwise permanova 
install.packages("devtools")
devtools::install_github("leffj/mctoolsr")
library("mctoolsr")

calc_pairwise_permanovas(as(vegdist(t(otu_table(ps_pipiens_organ))), "matrix"), as(sample_data(ps_pipiens_organ), "data.frame"), "Organ")
#           X1             X2         R2  pval pvalBon pvalFDR
# 1 Intestine          Ovary 0.12644103 0.001   0.003   0.003
# 2 Intestine Salivary gland 0.01773126 0.396   1.188   0.396
# 3     Ovary Salivary gland 0.12215389 0.007   0.021   0.010




#--------------------------------------------------------------------------------------------#
#---------------------------------------SAVE SESSION-----------------------------------------#
#--------------------------------------------------------------------------------------------#

#install.packages("session", repos = "http://cran.us.r-project.org")
library(session); packageVersion("session")
save.session("plot3.Rda")
#load("phyloseq.Rda")






#--------------------------------------------------------------------------------------------#
#---------------------------------------STATS ON WHOLE---------------------------------------#
#--------------------------------------------------------------------------------------------#

# Aedes whole
ps_whole_aedes <- subset_samples(ps_whole, Species=="Aedes aegypti")
ps_whole_aedes_wolbachia <- subset_taxa(ps_whole_aedes, Genus=="Wolbachia")
ps_whole_aedes_wolbachia <- prune_taxa(taxa_sums(ps_whole_aedes_wolbachia) >= 1, ps_whole_aedes_wolbachia)
ps_whole_aedes_wolbachia <- prune_samples(sample_sums(ps_whole_aedes_wolbachia) >= 1, ps_whole_aedes_wolbachia)
ps_whole_aedes_wolbachia

otu_new_aedes_wolbachia <- as(otu_table(ps_whole_aedes_wolbachia),"matrix")
metadata_new_aedes_wolbachia <-as(sample_data(ps_whole_aedes_wolbachia),"matrix")
sum(otu_new_aedes_wolbachia)/7

# Pipiens whole 
ps_whole_pipiens_bosc <- subset_samples(ps_whole, Species=="Culex pipiens" & Location=="Bosc")
ps_whole_pipiens_bosc_wolbachia <- subset_taxa(ps_whole_pipiens_bosc, Genus=="Wolbachia")
ps_whole_pipiens_bosc_wolbachia <- prune_taxa(taxa_sums(ps_whole_pipiens_bosc_wolbachia) >= 1, ps_whole_pipiens_bosc_wolbachia)
ps_whole_pipiens_bosc_wolbachia <- prune_samples(sample_sums(ps_whole_pipiens_bosc_wolbachia) >= 1, ps_whole_pipiens_bosc_wolbachia)
ps_whole_pipiens_bosc_wolbachia
otu_new_pipiens_bosc_wolbachia <- as(otu_table(ps_whole_pipiens_bosc_wolbachia),"matrix")
sum(otu_new_pipiens_bosc_wolbachia)/15

ps_whole_pipiens_CE <- subset_samples(ps_whole, Species=="Culex pipiens" & Location=="Camping Europe")
ps_whole_pipiens_CE_wolbachia <- subset_taxa(ps_whole_pipiens_CE, Genus=="Wolbachia")
ps_whole_pipiens_CE_wolbachia <- prune_taxa(taxa_sums(ps_whole_pipiens_CE_wolbachia) >= 1, ps_whole_pipiens_CE_wolbachia)
ps_whole_pipiens_CE_wolbachia <- prune_samples(sample_sums(ps_whole_pipiens_CE_wolbachia) >= 1, ps_whole_pipiens_CE_wolbachia)
ps_whole_pipiens_CE_wolbachia
otu_new_pipiens_CE_wolbachia <- as(otu_table(ps_whole_pipiens_CE_wolbachia),"matrix")
sum(otu_new_pipiens_CE_wolbachia)/7

ps_whole_pipiens_lavar <- subset_samples(ps_whole, Species=="Culex pipiens" & Location=="Lavar (labo)")
ps_whole_pipiens_lavar_wolbachia <- subset_taxa(ps_whole_pipiens_lavar, Genus=="Wolbachia")
ps_whole_pipiens_lavar_wolbachia <- prune_taxa(taxa_sums(ps_whole_pipiens_lavar_wolbachia) >= 1, ps_whole_pipiens_lavar_wolbachia)
ps_whole_pipiens_lavar_wolbachia <- prune_samples(sample_sums(ps_whole_pipiens_lavar_wolbachia) >= 1, ps_whole_pipiens_lavar_wolbachia)
ps_whole_pipiens_lavar_wolbachia
otu_new_pipiens_lavar_wolbachia <- as(otu_table(ps_whole_pipiens_lavar_wolbachia),"matrix")
sum(otu_new_pipiens_lavar_wolbachia)/22

# Quinquefasciatus whole 
ps_whole_quinque_guada <- subset_samples(ps_whole, Species=="Culex quinquefasciatus" & Location=="Guadeloupe")
ps_whole_quinque_guada_wolbachia <- subset_taxa(ps_whole_quinque_guada, Genus=="Wolbachia")
ps_whole_quinque_guada_wolbachia <- prune_taxa(taxa_sums(ps_whole_quinque_guada_wolbachia) >= 1, ps_whole_quinque_guada_wolbachia)
ps_whole_quinque_guada_wolbachia <- prune_samples(sample_sums(ps_whole_quinque_guada_wolbachia) >= 1, ps_whole_quinque_guada_wolbachia)
ps_whole_quinque_guada_wolbachia
otu_new_quinque_guada_wolbachia <- as(otu_table(ps_whole_quinque_guada_wolbachia),"matrix")
sum(otu_new_quinque_guada_wolbachia)/7

ps_whole_quinque_TC <- subset_samples(ps_whole, Species=="Culex quinquefasciatus" & Location=="Wolbachia -")
ps_whole_quinque_TC_wolbachia <- subset_taxa(ps_whole_quinque_TC, Genus=="Wolbachia")
ps_whole_quinque_TC_wolbachia <- prune_taxa(taxa_sums(ps_whole_quinque_TC_wolbachia) >= 1, ps_whole_quinque_TC_wolbachia)
ps_whole_quinque_TC_wolbachia <- prune_samples(sample_sums(ps_whole_quinque_TC_wolbachia) >= 1, ps_whole_quinque_TC_wolbachia)
ps_whole_quinque_TC_wolbachia
otu_new_quinque_TC_wolbachia <- as(otu_table(ps_whole_quinque_TC_wolbachia),"matrix")
sum(otu_new_quinque_TC_wolbachia)/13



#--------------------------------------------------------------------------------------------#
#--------------------------------------STATS ON ALL DATA-------------------------------------#
#--------------------------------------------------------------------------------------------#
# Stats : how many Wolbachia in Culex and how many Wolbachia - ?

ps # 2454 taxa and 208 samples
sum(as(otu_table(ps),"matrix")) # 26 659 233 counts
sum(as(otu_table(ps),"matrix"))/208 # 128 169

# Raw data
ps_decontam2 # 2031 taxa and 193 samples
ps_decontam2 <- prune_taxa(taxa_sums(ps_decontam2) >= 1, ps_decontam2)
ps_decontam2 <- prune_samples(sample_sums(ps_decontam2) >= 1, ps_decontam2)

sum(as(otu_table(ps_decontam2),"matrix")) # 26 289 866 counts
sum(as(otu_table(ps_decontam2),"matrix"))/193 # 136 216
reads_by_sample <- as(colSums(as(otu_table(ps_decontam2),"matrix")), "matrix")
sqrt(mean(as(otu_table(ps_decontam2),"matrix")^2)-mean(as(otu_table(ps_decontam2),"matrix"))^2)
min(colSums(as(otu_table(ps_decontam2),"matrix"))) #10
max(colSums(as(otu_table(ps_decontam2),"matrix"))) # 1 390 286
sd(colSums(as(otu_table(ps_decontam2),"matrix"))) # 266 929
mean(colSums(as(otu_table(ps_decontam2),"matrix"))) # 136 216

# Culex
ps_culex <- subset_samples(ps_decontam2, Species!="Aedes aegypti")
ps_culex <- prune_taxa(taxa_sums(ps_culex) >= 1, ps_culex)
ps_culex <- prune_samples(sample_sums(ps_culex) >= 1, ps_culex)
ps_culex # 1394 taxa and 177 samples --> 68,6% of ASV are in Culex
sum(as(otu_table(ps_culex),"matrix")) # 19 811 598 counts --> 75,3% of counts are in Culex

# Culex wolbachia
ps_culex_wolbachia <- subset_taxa(ps_culex, Genus=="Wolbachia")
ps_culex_wolbachia <- prune_taxa(taxa_sums(ps_culex_wolbachia) >= 1, ps_culex_wolbachia)
ps_culex_wolbachia <- prune_samples(sample_sums(ps_culex_wolbachia) >= 1, ps_culex_wolbachia)
ps_culex_wolbachia # 198 taxa and 171 samples --> 14,2% of ASV are Wolbachia 
sum(as(otu_table(ps_culex_wolbachia),"matrix")) # 12 566 880 counts --> 63,4% of counts are Wolbachia within Culex 

# Culex pipiens
ps_pipiens <- subset_samples(ps_decontam2, Species=="Culex pipiens")
ps_pipiens <- prune_taxa(taxa_sums(ps_pipiens) >= 1, ps_pipiens)
ps_pipiens <- prune_samples(sample_sums(ps_pipiens) >= 1, ps_pipiens)
ps_pipiens # 1115 taxa and 142 samples --> 54,9% of ASV are in Culex
sum(as(otu_table(ps_pipiens),"matrix")) # 8 271 948 reads --> 31,5% of reads are in Culex
sum(as(otu_table(ps_pipiens),"matrix"))/142 # 58 253 reads by sample

min(colSums(as(otu_table(ps_pipiens),"matrix"))) # 12
max(colSums(as(otu_table(ps_pipiens),"matrix"))) # 873 478
mean(colSums(as(otu_table(ps_pipiens),"matrix"))) # 58 253
sd(colSums(as(otu_table(ps_pipiens),"matrix"))) # 85 417

# Culex pipiens - Wolbachia
ps_pipiens_wolbachia <- subset_taxa(ps_pipiens, Genus=="Wolbachia")
ps_pipiens_wolbachia <- prune_taxa(taxa_sums(ps_pipiens_wolbachia) >= 1, ps_pipiens_wolbachia)
ps_pipiens_wolbachia <- prune_samples(sample_sums(ps_pipiens_wolbachia) >= 1, ps_pipiens_wolbachia)
ps_pipiens_wolbachia # 109 taxa and 139 samples --> 5,37% of ASV are Wolbachia 
sum(as(otu_table(ps_pipiens_wolbachia),"matrix")) # 6 956 089 reads --> 26,46% of reads are Wolbachia within Culex 
sum(as(otu_table(ps_pipiens_wolbachia),"matrix"))/139 # 50 0043 reads of Wolbachia by sample

min(colSums(as(otu_table(ps_pipiens_wolbachia),"matrix"))) # 7
max(colSums(as(otu_table(ps_pipiens_wolbachia),"matrix"))) # 873 369
mean(colSums(as(otu_table(ps_pipiens_wolbachia),"matrix"))) # 50 043
sd(colSums(as(otu_table(ps_pipiens_wolbachia),"matrix"))) # 87 826


# Culex quinquefasciatus 
ps_quinque <- subset_samples(ps_decontam2, Species=="Culex quinquefasciatus")
ps_quinque <- prune_taxa(taxa_sums(ps_quinque) >= 1, ps_quinque)
ps_quinque <- prune_samples(sample_sums(ps_quinque) >= 1, ps_quinque)
ps_quinque # 388 taxa and 35 samples --> 19,1% of ASV are in Culex
sum(as(otu_table(ps_quinque),"matrix")) # 11 539 650 reads --> 43,9% of reads are in Culex
sum(as(otu_table(ps_quinque),"matrix"))/35 

min(colSums(as(otu_table(ps_quinque),"matrix"))) # 10
max(colSums(as(otu_table(ps_quinque),"matrix"))) # 1 390 286
mean(colSums(as(otu_table(ps_quinque),"matrix"))) # 329 704
sd(colSums(as(otu_table(ps_quinque),"matrix"))) # 440 686

# Culex quiquefasciatus - Wolbachia
ps_quinque_wolbachia <- subset_taxa(ps_quinque, Genus=="Wolbachia")
ps_quinque_wolbachia <- prune_taxa(taxa_sums(ps_quinque_wolbachia) >= 1, ps_quinque_wolbachia)
ps_quinque_wolbachia <- prune_samples(sample_sums(ps_quinque_wolbachia) >= 1, ps_quinque_wolbachia)
ps_quinque_wolbachia # 117 taxa and 32 samples --> 5,8% of ASV are Wolbachia 
sum(as(otu_table(ps_quinque_wolbachia),"matrix")) # 5 610 791 reads --> 21,3% of reads are Wolbachia within Culex

min(colSums(as(otu_table(ps_quinque_wolbachia),"matrix"))) # 5
max(colSums(as(otu_table(ps_quinque_wolbachia),"matrix"))) # 1 389 118
mean(colSums(as(otu_table(ps_quinque_wolbachia),"matrix"))) # 175 337
sd(colSums(as(otu_table(ps_quinque_wolbachia),"matrix"))) # 421 858


# Aedes aegypti 
ps_aedes <- subset_samples(ps_decontam2, Species=="Aedes aegypti")
ps_aedes <- prune_taxa(taxa_sums(ps_aedes) >= 1, ps_aedes)
ps_aedes <- prune_samples(sample_sums(ps_aedes) >= 1, ps_aedes)
ps_aedes # 744 taxa and 16 samples -->  36,6% of ASV are in Aedes aegypti
sum(as(otu_table(ps_aedes),"matrix")) # 6 478 268 reads --> 24,6% of reads are in Culex

min(colSums(as(otu_table(ps_aedes),"matrix"))) # 130
max(colSums(as(otu_table(ps_aedes),"matrix"))) # 976 860
mean(colSums(as(otu_table(ps_aedes),"matrix"))) # 404 891
sd(colSums(as(otu_table(ps_aedes),"matrix"))) # 425 845

# Aedes aegypti - Wolbachia
ps_aedes_wolbachia <- subset_taxa(ps_aedes, Genus=="Wolbachia")
ps_aedes_wolbachia <- prune_taxa(taxa_sums(ps_aedes_wolbachia) >= 1, ps_aedes_wolbachia)
ps_aedes_wolbachia <- prune_samples(sample_sums(ps_aedes_wolbachia) >= 1, ps_aedes_wolbachia)
ps_aedes_wolbachia # 10 taxa and 11 samples --> 0,49% of ASV are Wolbachia 
sum(as(otu_table(ps_aedes_wolbachia),"matrix")) # 48 059 counts --> 0,18% of counts are Wolbachia within Culex 

min(colSums(as(otu_table(ps_aedes_wolbachia),"matrix"))) # 107
max(colSums(as(otu_table(ps_aedes_wolbachia),"matrix"))) # 24 561
mean(colSums(as(otu_table(ps_aedes_wolbachia),"matrix"))) # 4369
sd(colSums(as(otu_table(ps_aedes_wolbachia),"matrix"))) # 7613

# Wolbachia -
ps_wolbachia_neg <- subset_samples(ps_decontam2, Location=="Wolbachia -")
ps_wolbachia_neg <- prune_taxa(taxa_sums(ps_wolbachia_neg) >= 1, ps_wolbachia_neg)
ps_wolbachia_neg <- prune_samples(sample_sums(ps_wolbachia_neg) >= 1, ps_wolbachia_neg)
ps_wolbachia_neg # 170 taxa and 14 samples --> 8,3% of ASV are in Wolbachia - samples
sum(as(otu_table(ps_wolbachia_neg),"matrix")) # 5 921 746 --> Wolbachia - represents 22,5% of total counts





# Ovary - Asaia 

ps_ovary <- subset_samples(ps_percent, Organ=="Ovary")
ps_ovary <- prune_taxa(taxa_sums(ps_ovary) >= 1, ps_ovary)
ps_ovary <- prune_samples(sample_sums(ps_ovary) >= 1, ps_ovary)
sum(as(otu_table(ps_ovary),"matrix")) # 4851 counts
sum(as(otu_table(ps_ovary),"matrix"))/49 # 99 counts by sample

ps_ovary_asaia <- subset_taxa(ps_ovary, Genus=="Asaia")
ps_ovary_asaia <- prune_taxa(taxa_sums(ps_ovary_asaia) >= 1, ps_ovary_asaia)
ps_ovary_asaia <- prune_samples(sample_sums(ps_ovary_asaia) >= 1, ps_ovary_asaia)
sum(as(otu_table(ps_ovary_asaia),"matrix")) # 3 counts
sum(as(otu_table(ps_ovary_asaia),"matrix"))/1 # 3 counts by sample in NP14 (ovary of Aedes aegypti)
ps_ovary_asaia


# Ovary - Wolbachia
ps_ovary_wolbachia <- subset_taxa(ps_ovary, Genus=="Wolbachia")
ps_ovary_wolbachia <- prune_taxa(taxa_sums(ps_ovary_wolbachia) >= 1, ps_ovary_wolbachia)
ps_ovary_wolbachia <- prune_samples(sample_sums(ps_ovary_wolbachia) >= 1, ps_ovary_wolbachia)
sum(as(otu_table(ps_ovary_wolbachia),"matrix")) # 4801 counts
sum(as(otu_table(ps_ovary_wolbachia),"matrix"))/49 # 97 counts by sample
ps_ovary_wolbachia

