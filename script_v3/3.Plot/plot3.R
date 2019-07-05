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





# Stats : how many Wolbachia in Culex and how many Wolbachia - ?
ps_decontam2 # 2031 taxa and 193 samples
sum(as(otu_table(ps_decontam2),"matrix")) # 26 289 866 counts

ps_culex <- subset_samples(ps_decontam2, Species!="Aedes aegypti")
ps_culex <- prune_taxa(taxa_sums(ps_culex) >= 1, ps_culex)
ps_culex <- prune_samples(sample_sums(ps_culex) >= 1, ps_culex)
ps_culex # 1394 taxa and 177 samples --> 68,6% of ASV are in Culex
sum(as(otu_table(ps_culex),"matrix")) # 19 811 598 counts --> 75,3% of counts are in Culex

ps_culex_wolbachia <- subset_taxa(ps_culex, Genus=="Wolbachia")
ps_culex_wolbachia <- prune_taxa(taxa_sums(ps_culex_wolbachia) >= 1, ps_culex_wolbachia)
ps_culex_wolbachia <- prune_samples(sample_sums(ps_culex_wolbachia) >= 1, ps_culex_wolbachia)
ps_culex_wolbachia # 198 taxa and 171 samples --> 14,2% of ASV are Wolbachia 
sum(as(otu_table(ps_culex_wolbachia),"matrix")) # 12 566 880 counts --> 63,4% of counts are Wolbachia within Culex 


ps_wolbachia_neg <- subset_samples(ps_decontam2, Location=="Wolbachia -")
ps_wolbachia_neg <- prune_taxa(taxa_sums(ps_wolbachia_neg) >= 1, ps_wolbachia_neg)
ps_wolbachia_neg <- prune_samples(sample_sums(ps_wolbachia_neg) >= 1, ps_wolbachia_neg)
ps_wolbachia_neg # 170 taxa and 14 samples --> 8,3% of ASV are in Wolbachia - samples
sum(as(otu_table(ps_wolbachia_neg),"matrix")) # 5 921 746 --> Wolbachia - represents 22,5% of total counts




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
ps_percent_nopool <- subset_samples(ps_percent, Organ!="Pool")
p <- plot_composition(ps_percent_nopool,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Phylum", 
                      numberOfTaxa = 20, 
                      fill= "Phylum") +
  facet_wrap(~ Organ + Species, scales = "free_x", ncol = 3) + 
  labs(title = "Proteobacteria is the dominant phylum",
       caption = "Taxonomic composition (20 most abundant phylum)", x="Sample", y = "Abundance")

pdf("4-taxo_phylum.pdf")
plot(p)
dev.off()

jpeg("4-taxo_phylum.jpg", width = 1080, height = 720)
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
  facet_wrap(~ Organ + Species, scales = "free_x", ncol=3) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Alphaproteobacteria is the dominant class",
       caption = "Taxonomic composition (class)", x="Sample", y = "Abundance")

pdf("5-taxo_proteo.pdf")
plot(p)
dev.off()



# Alphaproteobacteria 

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

pdf("6-taxo_alphaproteo.pdf")
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

prop.pipiens <- transform_sample_counts(full_pipiens, function(count_tab) count_tab/sum(count_tab))
jacc.pipiens <- ordinate(full_pipiens, method="NMDS", distance="jaccard")


pdf("7-NMDS_bray_full.pdf")
plot_ordination(prop.full, bray.full, color="Species", title="Bray NMDS with full body", label="Sample") +
  labs(title = "Does species influence bacterial community structure ? ",
       caption = "Bray NMDS on whole mosquitoes", x="NMDS1", y = "NMDS2")+
  geom_point(size = 4) +
  theme_gray()
dev.off()

pdf("8-NMDS_bray_full_culex.pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.full_culex, bray.full_culex, color="Field", shape="Location", title="Bray NMDS with full body - Location without Aedes aegypti", label="Sample") +
  labs(title = "Do antibiotics influence microbiote ? ",
       caption = "Bray NMDS on whole Culex mosquitoes", x="NMDS1", y = "NMDS2") +
  #stat_ellipse(geom = "polygon", level=0.6,alpha = 1/2, aes(fill = Field))+
  scale_fill_manual(values=c("yellow","green"))+
  scale_color_manual(values=c("red", "blue"))+
  geom_point(size = 5) +
  theme_gray()
dev.off()


pdf("9-NMDS_bray_full_culex_ellipse.pdf")
#jpeg("NMDS_bray_fullbody.jpg")
  plot_ordination(prop.full_culex, bray.full_culex, color="Field", shape="Location", title="Bray NMDS with full body - Location without Aedes aegypti", label="Sample") +
  labs(title = "Do antibiotics influence microbiote ? ",
       caption = "Bray NMDS on whole Culex mosquitoes", x="NMDS1", y = "NMDS2") +
  stat_ellipse(geom = "polygon", level=0.95,alpha = 1/2, aes(fill = Species))+
  scale_fill_manual(values=c("yellow","green"))+
  scale_color_manual(values=c("red", "blue"))+
  geom_point(size = 5) +
  theme_gray()
dev.off()
  
  


pdf("10-NMDS_bray_full_pipiens.pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.pipiens, bray.pipiens, color="Field", shape="Location", title="Bray NMDS with full body - Location without Aedes aegypti", label="Sample") +
  labs(title = expression(paste("Does laboratory influence microbiote of ", italic("Culex pipiens"), "?")),
       caption = "Bray NMDS", x="NMDS1", y = "NMDS2") +
  geom_point(size = 4) +
  theme_gray()
dev.off()



pdf("11-NMDS_jaccard_full_pipiens.pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.pipiens, jacc.pipiens, color="Field", shape="Location", title="Jaccard NMDS with full body - Location without Aedes aegypti", label="Sample") +
  labs(title = expression(paste("Does laboratory influence microbiote of ", italic("Culex pipiens"), "?")),
       caption = "Jaccard NMDS", x="NMDS1", y = "NMDS2") +
  geom_point(size = 4) +
  theme_gray()
dev.off()


# Organs of Culex pipiens from Cammping Europe 
ps_new <- subset_samples(ps_percent, Organ!="Full" & Species=="Culex pipiens")
ps_new <- subset_samples(ps_new, Sample!="S175")
ps_new <- subset_samples(ps_new, Location=="Camping Europe")
test <- as(sample_data(ps_new),"matrix")

prop.new <- transform_sample_counts(ps_new, function(count_tab) count_tab/sum(count_tab))
bray.new <- ordinate(ps_new, method="NMDS", distance="bray")

prop.new2 <- transform_sample_counts(ps_new, function(count_tab) count_tab/sum(count_tab))
bray.new2 <- ordinate(ps_new, method="NMDS", distance="jaccard")


pdf("11-NMDS_bray_organs_culex_CE.pdf")
plot_ordination(prop.new, bray.new, color="Organ", shape="Date", title="Bray NMDS", label="Sample") +
  labs(title = expression(paste("Organ influences the structure of ", italic("Culex pipiens"), " bacterial community")),
       caption = expression(paste("Bray NMDS on Ovary, Intestine and Salivary Gland of ", italic('Culex pipiens'), " from Camping Europe")), x="NMDS1", y = "NMDS2")+
  geom_point(size = 4) +
  theme_gray()
dev.off()

pdf("12-NMDS_jaccard_organs_culex_CE.pdf")
plot_ordination(prop.new2, bray.new2, color="Organ", shape="Date", title="Jaccard NMDS", label="Sample") +
  labs(title = expression(paste("Organ influences the structure of ", italic("Culex pipiens"), " bacterial community")),
       caption = expression(paste("Jaccard NMDS on Ovary, Intestine and Salivary Gland of ", italic('Culex pipiens'), " from Camping Europe")), x="NMDS1", y = "NMDS2")+
  geom_point(size = 4) +
  theme_gray()
dev.off()



# Intestine of Culex pipiens from Camping Europe at 2 dates
intestine <- subset_samples(ps_percent, Organ == "Intestine")
intestine_camping <- subset_samples(ps_percent, Location == "Camping Europe" & Organ =="Intestine")
intestine_camping_date <- subset_samples(intestine_camping, Date =="30/05/2017" | Date =="28/06/2017")

prop.intestine <- transform_sample_counts(intestine, function(count_tab) count_tab/sum(count_tab))
bray.intestine <- ordinate(intestine, method="NMDS", distance="bray")

prop.intestine_camping_date <- transform_sample_counts(intestine_camping_date, function(count_tab) count_tab/sum(count_tab))
bray.intestine_camping_date <- ordinate(intestine_camping_date, method="NMDS", distance="bray")

pdf("13-NMDS_bray_intestine_CE_date.pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(prop.intestine_camping_date, bray.intestine_camping_date, color="Date", title="Bray NMDS with intestine of Culex pipiens in Camping Europe", label="Sample")+
  labs(title = expression(paste("Does time of sampling influence the microbiote of ", italic("Culex pipiens"),"?")),
       caption = "Bray NMDS on intestine at Camping Europe", x="NMDS1", y = "NMDS2") +
  geom_point(size = 4) +
  theme_gray()
dev.off()



# Ovary of Culex pipiens from Camping Europe at 2 dates
ovary <- subset_samples(ps_percent, Organ == "Ovary")
ovary_culex_camping <- subset_samples(ovary, Species=="Culex pipiens" & Location == "Camping Europe")
ovary_culex_camping_date <- subset_samples(ovary_culex_camping, Date =="30/05/2017" | Date =="28/06/2017")
ovary_culex_camping_date <- subset_samples(ovary_culex_camping_date, Sample!="S103")

prop.ovary <- transform_sample_counts(ovary, function(count_tab) count_tab/sum(count_tab))
bray.ovary <- ordinate(ovary, method="NMDS", distance="bray")

prop.ovary_culex_camping_date <- transform_sample_counts(ovary_culex_camping_date, function(count_tab) count_tab/sum(count_tab))
bray.ovary_culex_camping_date <- ordinate(ovary_culex_camping_date, method="NMDS", distance="bray")

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
ps_pipiens_wolbachia <- subset_taxa(ps_pipiens, Genus=="Wolbachia")
ps_pipiens_wolbachia <- prune_taxa(taxa_sums(ps_pipiens_wolbachia) >= 1, ps_pipiens_wolbachia)
ps_pipiens_wolbachia <- prune_samples(sample_sums(ps_pipiens_wolbachia) >= 1, ps_pipiens_wolbachia)
ps_pipiens_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_pipiens_wolbachia),TRUE)[1:30]), ps_pipiens_wolbachia)

jpeg("15-heatmap_pipiens_wolbachia.jpg", width = 1080, height = 720)
plot_heatmap(ps_pipiens_wolbachia, sample.label="Field", sample.order="Field", low="#000033", high="#FF3300")+
  labs(title = expression(paste("ASV1 is abundant in almost all sequences of ", italic("Culex pipiens"))),
       caption = expression(paste("Heatmap of  ", italic('Culex pipiens'), " that contains Wolbachia")), x="Field", y = "ASV")
dev.off()

# Culex quinquefasciatus
ps_quinque <- subset_samples(ps_decontam2, Species=="Culex quinquefasciatus")
ps_quinque_wolbachia <- subset_taxa(ps_quinque, Genus=="Wolbachia")
ps_quinque_wolbachia <- prune_taxa(taxa_sums(ps_quinque_wolbachia) >= 1, ps_quinque_wolbachia)
ps_quinque_wolbachia <- prune_samples(sample_sums(ps_quinque_wolbachia) >= 1, ps_quinque_wolbachia)
ps_quinque_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_quinque_wolbachia),TRUE)[1:30]), ps_quinque_wolbachia)

jpeg("16-heatmap_quinque_wolbachia.jpg", width = 1080, height = 720)
plot_heatmap(ps_quinque_wolbachia, sample.label="Field", sample.order="Field", low="#000033", high="#FF3300")+
  labs(title = expression(paste("ASV1 is abundant in almost all sequences of ", italic("Culex quinquefasciatus"))),
       caption = expression(paste("Heatmap of  ", italic('Culex quinquefasciatus'), " that contains Wolbachia")), x="Field", y = "ASV")
dev.off()


# Aedes aegypti
ps_aedes <- subset_samples(ps_decontam2, Species=="Aedes aegypti")
ps_aedes_wolbachia <- subset_taxa(ps_aedes, Genus=="Wolbachia")
ps_aedes_wolbachia <- prune_taxa(taxa_sums(ps_aedes_wolbachia) >= 1, ps_aedes_wolbachia)
ps_aedes_wolbachia <- prune_samples(sample_sums(ps_aedes_wolbachia) >= 1, ps_aedes_wolbachia)
ps_aedes_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_aedes_wolbachia),TRUE)[1:30]), ps_aedes_wolbachia)

jpeg("17-heatmap_aedes_wolbachia.jpg", width = 1080, height = 720)
plot_heatmap(ps_aedes_wolbachia, sample.label="Field", sample.order="Field", low="#000033", high="#FF3300")+
  labs(title = expression(paste("ASV1 is abundant in almost all sequences of ", italic("Aedes aegypti"))),
       caption = expression(paste("Heatmap of  ", italic('Aedes aegypti'), " that contains Wolbachia")), x="Field", y = "ASV")
dev.off()



# Whole 
ps_whole <- subset_samples(ps_decontam2, Organ=="Full")
ps_whole_wolbachia <- subset_taxa(ps_whole, Genus=="Wolbachia")
ps_whole_wolbachia <- prune_taxa(taxa_sums(ps_whole_wolbachia) >= 1, ps_whole_wolbachia)
ps_whole_wolbachia <- prune_samples(sample_sums(ps_whole_wolbachia) >= 1, ps_whole_wolbachia)
ps_whole_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_whole_wolbachia),TRUE)[1:30]), ps_whole_wolbachia)

jpeg("18-heatmap_whole_wolbachia.jpg", width = 1080, height = 720)
plot_heatmap(ps_whole_wolbachia, sample.label="Location", sample.order="Field", low="#000033", high="#FF3300")+
  labs(title = expression(paste("ASV1 is abundant in almost all sequences of Whole mosquitoes")),
       caption = expression(paste("Heatmap of Whole mosquitoes that contain Wolbachia")), x="Field", y = "ASV")
dev.off()


# Intestine
ps_intestine <- subset_samples(ps_decontam2, Organ=="Intestine")
ps_intestine_wolbachia <- subset_taxa(ps_intestine, Genus=="Wolbachia")
ps_intestine_wolbachia <- prune_taxa(taxa_sums(ps_intestine_wolbachia) >= 1, ps_intestine_wolbachia)
ps_intestine_wolbachia <- prune_samples(sample_sums(ps_intestine_wolbachia) >= 1, ps_intestine_wolbachia)
ps_intestine_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_intestine_wolbachia),TRUE)[1:30]), ps_intestine_wolbachia)

jpeg("19-heatmap_intestine_wolbachia.jpg", width = 1080, height = 720)
plot_heatmap(ps_intestine_wolbachia, sample.label="Location", sample.order="Field", low="#000033", high="#FF3300")+
  labs(title = expression(paste("ASV1 is abundant in almost all sequences of Intestine samples")),
       caption = expression(paste("Heatmap of Intestine samples that contain Wolbachia")), x="Field", y = "ASV")
dev.off()

# Ovary 
ps_ovary <- subset_samples(ps_decontam2, Organ=="Ovary")
ps_ovary_wolbachia <- subset_taxa(ps_ovary, Genus=="Wolbachia")
ps_ovary_wolbachia <- prune_taxa(taxa_sums(ps_ovary_wolbachia) >= 1, ps_ovary_wolbachia)
ps_ovary_wolbachia <- prune_samples(sample_sums(ps_ovary_wolbachia) >= 1, ps_ovary_wolbachia)
ps_ovary_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_ovary_wolbachia),TRUE)[1:30]), ps_ovary_wolbachia)

jpeg("20-heatmap_ovary_wolbachia.jpg", width = 1080, height = 720)
plot_heatmap(ps_ovary_wolbachia, sample.label="Location", sample.order="Field", low="#000033", high="#FF3300")+
  labs(title = expression(paste("ASV1 is abundant in almost all sequences of Ovary samples")),
       caption = expression(paste("Heatmap of Ovary samples that contain Wolbachia")), x="Field", y = "ASV")
dev.off()

# Salivary gland 
ps_gland <- subset_samples(ps_decontam2, Organ=="Salivary gland")
ps_gland_wolbachia <- subset_taxa(ps_gland, Genus=="Wolbachia")
ps_gland_wolbachia <- prune_taxa(taxa_sums(ps_gland_wolbachia) >= 1, ps_gland_wolbachia)
ps_gland_wolbachia <- prune_samples(sample_sums(ps_gland_wolbachia) >= 1, ps_gland_wolbachia)
ps_gland_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_gland_wolbachia),TRUE)[1:30]), ps_gland_wolbachia)

jpeg("21-heatmap_gland_wolbachia.jpg", width = 1080, height = 720)
plot_heatmap(ps_gland_wolbachia, sample.label="Location", sample.order="Field", low="#000033", high="#FF3300")+
  labs(title = expression(paste("ASV1 is abundant in almost all sequences of Salivary Gland samples")),
       caption = expression(paste("Heatmap of Salivary Gland samples that contain Wolbachia")), x="Field", y = "ASV")
dev.off()

# Pool 
ps_pool <- subset_samples(ps_decontam2, Organ=="Pool")
ps_pool_wolbachia <- subset_taxa(ps_pool, Genus=="Wolbachia")
ps_pool_wolbachia <- prune_taxa(taxa_sums(ps_pool_wolbachia) >= 1, ps_pool_wolbachia)
ps_pool_wolbachia <- prune_samples(sample_sums(ps_pool_wolbachia) >= 1, ps_pool_wolbachia)
ps_pool_wolbachia <- prune_taxa(names(sort(taxa_sums(ps_pool_wolbachia),TRUE)[1:30]), ps_pool_wolbachia)

jpeg("22-heatmap_pool_wolbachia.jpg", width = 1080, height = 720)
plot_heatmap(ps_pool_wolbachia, sample.label="Location", sample.order="Field", low="#000033", high="#FF3300")+
  labs(title = expression(paste("ASV1 is abundant in almost all sequences of Pool samples")),
       caption = expression(paste("Heatmap of Pool samples that contain Wolbachia")), x="Field", y = "ASV")
dev.off()





#--------------------------------------------------------------------------------------------#
#--------------------------------------WOLBACHIA STATS-------------------------------------#
#--------------------------------------------------------------------------------------------#


ps_field <- subset_samples(ps_pipiens_wolbachia, Field=="Field")
ps_field <- prune_taxa(taxa_sums(ps_field) >= 1, ps_field)
ps_field <- prune_samples(sample_sums(ps_field) >= 1, ps_field)
tax_field <- as(tax_table(ps_field),"matrix")
# otu_table()   OTU Table:         [ 98 taxa and 89 samples ]
# sample_data() Sample Data:       [ 89 samples by 16 sample variables ]
# tax_table()   Taxonomy Table:    [ 98 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 98 tips and 97 internal nodes ]

ps_labo <- subset_samples(ps_pipiens_wolbachia, Field=="Lab ")
ps_labo <- prune_taxa(taxa_sums(ps_labo) >= 1, ps_labo)
ps_labo <- prune_samples(sample_sums(ps_labo) >= 1, ps_labo)
tax_labo <- as(tax_table(ps_labo),"matrix")
# otu_table()   OTU Table:         [ 25 taxa and 48 samples ]
# sample_data() Sample Data:       [ 48 samples by 16 sample variables ]
# tax_table()   Taxonomy Table:    [ 25 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 25 tips and 24 internal nodes ]


test3 <- as(tax_table(ps_pipiens),"matrix")
test4 <- as(otu_table(ps_pipiens_wolbachia),"matrix")
test5 <- as(sample_data(ps_pipiens_wolbachia),"matrix")





#--------------------------------------------------------------------------------------------#
#-------------------------------------------ADONIS-------------------------------------------#
#--------------------------------------------------------------------------------------------#

# Culex whole 
ps_culex_whole <- subset_samples(ps_percent, Organ=="Full" & Species!="Aedes aegypti")
adonis(vegdist(t(otu_table(ps_culex_whole)), method = "bray") ~ Location*Species,
       data=as(sample_data(ps_culex_whole), "data.frame"), permutation = 9999)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Location   4    8.4039 2.10098  11.316 0.42199  1e-04 ***
#   Residuals 62   11.5110 0.18566         0.57801           
# Total     66   19.9149                 1.00000           
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1





#--------------------------------------------------------------------------------------------#
#---------------------------------------SAVE SESSION-----------------------------------------#
#--------------------------------------------------------------------------------------------#

#install.packages("session", repos = "http://cran.us.r-project.org")
library(session); packageVersion("session")
save.session("plot3.Rda")
#load("phyloseq.Rda")

