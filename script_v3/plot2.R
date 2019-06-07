# SCHRIEKE Hans
# PLOT
# input : phyloseq objects from decontam step
# output : plot and stats 



#--------------------------------------------------------------------------------------------#
#--------------------------------------PATH SETTING------------------------------------------#
#--------------------------------------------------------------------------------------------#

path = "D:/stage/data/runs_new2"
setwd(path)

dir.create("plotv3") # folder for plot
path2 <- "D:/stage/data/runs_new2/plotv3"


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
library("hrbrthemes")

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

p2 <- p1 + 
  facet_wrap(~ Organ) + 
  geom_vline(xintercept = min(sample_sums(ps)), 
             color = "gray60") +
  ggtitle("Rarefaction curve") +
  xlim(0,100000) +
  ylim(0, 230) +
  theme_gray()

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
                    measures=c("Observed","Shannon","ACE", "Chao1"), 
                    nrow = 1) +
  ggtitle("Alpha diversity")
# pdf("richness_ps_decontam2_species.pdf")
# print(p1)
# dev.off()

p2 <- plot_richness(data1, 
                    x="Sample", 
                    color="Organ", 
                    measures=c("Observed","Shannon","ACE", "Chao1"), 
                    nrow = 1) +
  ggtitle("Alpha diversity")


pdf("richness_decontam2_species.pdf")
ggplot(p1$data,aes(Species,value,colour=Species,shape=Species)) +
  facet_grid(variable ~ Species, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.8, 
               position = position_dodge(width=0.9)) + 
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  ylab("Diversity index")  + xlab(NULL) + theme_gray() +
  ggtitle("Alpha diversity")
dev.off()

pdf("richness_decontam2_organ.pdf")
ggplot(p1$data,aes(Organ,value,colour=Organ,shape=Species)) +
  facet_grid(variable ~ Organ, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.8, 
               position = position_dodge(width=0.9)) + 
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  ylab("Diversity index")  + xlab(NULL) + theme_gray() +
  ggtitle("Alpha diversity ")
dev.off()



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
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Taxonomic composition - all Bacteria") +
  theme_gray()

pdf("composition_percent_species.pdf")
plot(p)
dev.off()

p <- plot_composition(ps_percent,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Class", 
                      numberOfTaxa = 20, 
                      fill= "Phylum") +
  facet_wrap(~ Organ, scales = "free_x", nrow = 3) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Taxonomic composition - all Bacteria") +
  theme_gray()

pdf("composition_percent_organ.pdf")
plot(p)
dev.off()


# ps_proteo 
p <- plot_composition(ps_proteo,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Phylum", 
                      numberOfTaxa = 20, 
                      fill= "Class") +
  facet_wrap(~ Species, scales = "free_x", nrow = 3) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Taxonomic composition - Proteobacteria") +
  theme_gray()

pdf("composition_proteobacteria_species.pdf")
plot(p)
dev.off()

p <- plot_composition(ps_proteo,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Phylum", 
                      numberOfTaxa = 20, 
                      fill= "Class") +
  facet_wrap(~ Organ, scales = "free_x", nrow = 3) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Taxonomic composition - Proteobacteria") +
  theme_gray()

pdf("composition_proteobacteria_organ.pdf")
plot(p)
dev.off()

p <- plot_composition(ps_proteo,
                      taxaRank1 = "Class",
                      taxaSet1 ="Gammaproteobacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 20, 
                      fill= "Genus") +
  facet_wrap(~ Species, scales = "free_x", nrow = 5) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Taxonomic composition - Gammaproteobacteria") +
  theme_gray()

pdf("composition_gammaproteobacteria_species.pdf")
plot(p)
dev.off()

p <- plot_composition(ps_proteo,
                      taxaRank1 = "Class",
                      taxaSet1 ="Gammaproteobacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 20, 
                      fill= "Genus") +
  facet_wrap(~ Organ, scales = "free_x", nrow = 5) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Taxonomic composition - Gammaproteobacteria") +
  theme_gray()

pdf("composition_gammaproteobacteria_organ.pdf")
plot(p)
dev.off()

p <- plot_composition(ps_proteo,
                      taxaRank1 = "Class",
                      taxaSet1 ="Alphaproteobacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 20, 
                      fill= "Genus") +
  facet_wrap(~ Species, scales = "free_x", nrow = 5) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Taxonomic composition - Alphaproteobacteria") +
  theme_gray()

pdf("composition_alphaproteobacteria_species.pdf")
plot(p)
dev.off()

p <- plot_composition(ps_proteo,
                      taxaRank1 = "Class",
                      taxaSet1 ="Alphaproteobacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 20, 
                      fill= "Genus") +
  facet_wrap(~ Organ, scales = "free_x", nrow = 5) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Taxonomic composition - Alphaproteobacteria") +
  theme_gray()

pdf("composition_alphaproteobacteria_organ.pdf")
plot(p)
dev.off()

p <- plot_composition(ps_proteo,
                      taxaRank1 = "Class",
                      taxaSet1 ="Deltaproteobacteria",
                      taxaRank2 = "Genus", 
                      numberOfTaxa = 20, 
                      fill= "Genus") +
  facet_wrap(~ Species, scales = "free_x", nrow = 5) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Taxonomic composition - Deltaproteobacteria") +
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

#--------------------------------------------------------------------------------------------#
#------------------------------------------PCOA----------------------------------------------#
#--------------------------------------------------------------------------------------------#
pdf("PCoA_percent_field1.pdf")
plot_ordination(ps_percent, ordinate(ps_percent, method ="MDS", distance = "bray"), color = "Species", shape="Field") +
  geom_point(size = 4) +
  ggtitle("Bray PCoA - Labo vs Field") +
  theme_gray()
dev.off()

pdf("PCoA_percent_field2.pdf")
plot_ordination(ps_percent, ordinate(ps_percent, method ="MDS", distance = "bray"), color = "Field", shape="Species") +
  geom_point(size = 4) +
  ggtitle("Bray PCoA - Labo vs Field") +
  theme_gray()
dev.off()


#--------------------------------------------------------------------------------------------#
#-----------------------------------NMMDS (bray distance)------------------------------------#
#--------------------------------------------------------------------------------------------#

# phyloseq objects with the different conditions
ps_percent <- subset_samples(ps_percent, Sample != "S175")


full <- subset_samples(ps_percent, Organ == "Full" | Organ == "Pool")
full_no_aedes <- subset_samples(full, Species!="Aedes aegypti")
full_no_labo <- subset_samples(full_no_aedes, Location!="Labo Tetracycline" | Location!="Lavar")

intestine <- subset_samples(ps_percent, Organ == "Intestine")
intestine_camping_date <- subset_samples(ps_percent, Location == "Camping Europe" & Organ =="Intestine" & Date =="30/05/2017" | Date =="28/06/2017")
intestine_no_lavar <- subset_samples(intestine, Location != "Lavar")
intestine_filter <- subset_samples(intestine_camping_date, Sample !="NP17" & Sample !="NP20" & Sample !="S81" & Sample != "S82")

ovary <- subset_samples(ps_percent, Organ == "Ovary")
ovary_culex_camping <- subset_samples(ovary, Species=="Culex pipiens" & Location == "Camping Europe")
ovary_culex_camping_date <- subset_samples(ovary_culex_camping, Date =="30/05/2017" | Date =="28/06/2017")
ovary_filter <- subset_samples(ovary_culex_camping_date, Sample !="S103" & Sample !="S81" & Sample != "S82")



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
  geom_point(size = 4) +
  theme_gray()
dev.off()

pdf("NMDS_bray_full(with full vs pool).pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.full, bray.full, color="Species", shape="Organ", title="Bray NMDS with full body - Full vs Pool", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()

pdf("NMDS_bray_full(with aedes).pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.full, bray.full, color="Species", shape="Location", title="Bray NMDS with full body - Location", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()

pdf("NMDS_bray_full(without aedes).pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.full_no_aedes, bray.full_no_aedes, color="Species", shape="Location", title="Bray NMDS with full body - Location without Aedes aegypti", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()

pdf("NMDS_bray_full(without aedes and with field).pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.full_no_aedes, bray.full_no_aedes, color="Field", shape="Location", title="Bray NMDS with full body - Labo vs Field", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()

pdf("NMDS_bray_full(france vs gwada with aedes and labo).pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.full, bray.full, color="Country", shape="Location", title="Bray NMDS with full body - France vs Guadeloupe", label="Sample") +
  geom_point(size = 4) +
  scale_shape_manual(values=seq(0,15)) +
  theme_gray()
dev.off()

pdf("NMDS_bray_full(france vs gwada without aedes and labo).pdf")
#jpeg("NMDS_bray_fullbody.jpg")
plot_ordination(prop.full_no_labo, bray.full_no_labo, color="Country", shape="Location", title="Bray NMDS with full body - France vs Guadeloupe without Labo and Aedes aegypti", label="Sample") +
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
plot_ordination(prop.intestine_camping_date, bray.intestine_camping_date, color="Date", title="Bray NMDS with intestine of Culex pipiens in Camping Europe", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()

#voir quels samples sont écartés et les enlever pour faire ce plot
pdf("NMDS_bray_intestine.pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(prop.intestine, bray.intestine, color="Date", shape="Location", title="Bray NMDS with intestine - Date", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off

pdf("NMDS_bray_intestine.pdf")
#jpeg("NMDS_bray_intestine.jpg")
plot_ordination(prop.intestine_no_lavar, bray.intestine_no_lavar, color="Date", shape="Location", title="Bray NMDS with intestine - Date without Lavar", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()



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

# S103, S81, S82 bordeline
pdf("NMDS_bray_ovary_camping_culex_date.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(prop.ovary_culex_camping_date, bray.ovary_culex_camping_date, color="Date", title="Bray NMDS with ovary - Culex, Camping Europe, Dates", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()

pdf("NMDS_bray_ovary_filter.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(prop.ovary_filter, bray.ovary_filter, color="Date", title="Bray NMDS with ovary - Culex, Camping Europe, Dates", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()


# NMDS of Wolbachia 
ps_wolbachia_ovary_intestine <- subset_samples(ps_wolbachia, Organ=="Intestine" | Organ=="Ovary")
ps_wolbachia_filter <- subset_samples(ps_wolbachia_ovary_intestine, Sample != "S175" & Sample != "NP38" & Sample!="S99" & Sample!="NP22" & Sample!="NP10" & Sample!="S68" &
                                        Sample!="NP2" & Sample!="NP11" & Sample!="NP8" & Sample!="NP5")

prop.wolbachia <- transform_sample_counts(ps_wolbachia, function(count_tab) count_tab/sum(count_tab))
bray.wolbachia <- ordinate(ps_wolbachia, method="NMDS", distance="bray")

prop.wolbachia_filter <- transform_sample_counts(ps_wolbachia_filter, function(count_tab) count_tab/sum(count_tab))
bray.wolbachia_filter <- ordinate(ps_wolbachia_filter, method="NMDS", distance="bray")


pdf("NMDS_bray_wolbachia.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(prop.wolbachia, bray.wolbachia, color="Individual", shape="Organ", title="Bray NMDS - Wolbachia", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off() # many samples to remove

pdf("NMDS_bray_wolbachia_filter1.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(prop.wolbachia_filter, bray.wolbachia_filter, color="Individual", shape="Organ", title="Bray NMDS - Wolbachia filter", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()

pdf("NMDS_bray_wolbachia_filter2.pdf")
#jpeg("NMDS_bray_ovary.jpg")
plot_ordination(prop.wolbachia_filter, bray.wolbachia_filter, color="Individual", shape="Location", title="Bray NMDS - Wolbachia filter", label="Sample") +
  geom_point(size = 4) +
  theme_gray()
dev.off()



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

adonis(vegdist(t(otu_table(ps_percent)), method = "bray") ~ Organ*Location*Date,
       data=as(sample_data(ps_percent), "data.frame"), permutation = 9999)

#                 Df  SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# Organ            4     6.244 1.56097  8.1461 0.10782 0.0001 ***
# Location         4    10.197 2.54920 13.3033 0.17608 0.0001 ***
# Date             4     2.744 0.68590  3.5795 0.04738 0.0001 ***
# Organ:Location   7     3.926 0.56082  2.9267 0.06779 0.0001 ***
# Organ:Date       5     0.882 0.17636  0.9204 0.01523 0.5769    
# Residuals      177    33.917 0.19162         0.58570           
# Total          201    57.909                 1.00000           
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


adonis(vegdist(t(otu_table(ps_wolbachia)), method = "bray") ~ Organ*Location*Individual,
       data=as(sample_data(ps_wolbachia), "data.frame"), permutation = 9999)

#                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Organ              4    2.7440 0.68600  3.3541 0.12551 0.0125 *
# Location           4    1.5301 0.38254  1.8704 0.06999 0.1298  
# Individual        96   10.1818 0.10606  0.5186 0.46570 0.9920  
# Organ:Location     3    0.2010 0.06700  0.3276 0.00919 0.9412  
# Organ:Individual  32    2.7069 0.08459  0.4136 0.12381 0.9956  
# Residuals         22    4.4995 0.20452         0.20580         
# Total            161   21.8634                 1.00000         
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



adonis(vegdist(t(otu_table(ps_wolbachia_ovary_intestine)), method = "bray") ~ Organ*Location*Individual,
       data=as(sample_data(ps_wolbachia_ovary_intestine), "data.frame"), permutation = 9999)


#                  Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Organ             1    0.7872 0.78720  8.1467 0.09224 0.0880 .
# Location          3    0.1664 0.05547  0.5741 0.01950 0.7193  
# Individual       43    3.9002 0.09070  0.9387 0.45698 0.5562  
# Organ:Location    3    0.2010 0.06700  0.6933 0.02355 0.6775  
# Organ:Individual 32    2.7069 0.08459  0.8754 0.31716 0.5926  
# Residuals         8    0.7730 0.09663         0.09057         
# Total            90    8.5347                 1.00000         
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1






#--------------------------------------------------------------------------------------------#
#---------------------------------------Save session-----------------------------------------#
#--------------------------------------------------------------------------------------------#

#install.packages("session", repos = "http://cran.us.r-project.org")
library(session); packageVersion("session")
save.session("plot2.Rda")
#load("phyloseq.Rda")