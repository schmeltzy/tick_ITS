##Overview: Beta diversity analyses
#load necessary packages
#beta diversity

library(phyloseq)
library(microViz)
library(dplyr)
library(ggplot2)
library(vegan)
library(pairwiseAdonis)

#set seed
set.seed(123)

#load in unrarefied phyloseq object
ps <- readRDS("phyloseq_obj.rds")

###let's subset to just look at ticks (not soil)
ticks <- subset_samples(ps, Type=="BLT")
#and ticks by habitat type
forest <- subset_samples(ticks, Habitat=="forested")
open <- subset_samples(ticks, Habitat=="open")

# does habitat, treatment, lifestage, season, or sex affect beta diversity?
#use aitchison distance for beta (clr transform followed by euclidean distance matrix)
aitch <- dist_calc(ticks, dist = "aitchison")
aitch.dist <- dist_get(aitch)

#get sample data from the matrix just to make it easier
aitch.sampledf <- data.frame(sample_data(aitch))

#permanova
adonis_habitat <- broom::tidy(adonis2(aitch.dist ~ Habitat, data = aitch.sampledf, by = "terms", perm = 999))  
#habitat p = 0.001 sig
adonis_habitat$Comparison <- "all"


#forested
aitch.forest <- dist_calc(forest, dist = "aitchison")
aitch.dist.forest <- dist_get(aitch.forest)

#get sample data from the matrix just to make it easier
forest.sampledf <- data.frame(sample_data(forest))


adonis_treatment <- broom::tidy(adonis2(aitch.dist.forest ~ Treatment, data = forest.sampledf, by = "terms", perm = 999))
#p = 0.001 sig
adonis_treatment$Comparison <- "forested"

adonis_season <- broom::tidy(adonis2(aitch.dist.forest ~ Treatment / Season, data = forest.sampledf, by = "terms", perm = 999))
#p = 0.001 sig
adonis_season$Comparison <- "forested"

adonis_lifestage <- broom::tidy(adonis2(aitch.dist.forest ~ Lifestage, data = forest.sampledf, by = "terms", perm = 999))
#p = 0.001 sig
adonis_lifestage$Comparison <- "forested"

#need to remove NAs for sex differences
adonis_sex <- broom::tidy(adonis2(aitch.dist.forest ~ Treatment / Sex, data = forest.sampledf, na.action = na.omit, by = "terms", perm = 999))  
adonis_sex$Comparison <- "Sex" #NS

rbind(adonis_habitat, adonis_treatment, adonis_season, adonis_lifestage, adonis_sex) -> adonis.all
write.csv(adonis_forest, here::here("output/adonis_forest.csv"))

## Pairwise adonis 
#pair.adonis.treat <- pairwise.adonis(aitch.dist.forest, forest.sampledf$Treatment, p.adjust.m = "fdr", perm = 999) 
#pair.adonis.treat$Comparison <- "forested"

pair.adonis.life <- pairwise.adonis(aitch.dist.forest, forest.sampledf$Lifestage, p.adjust.m = "fdr", perm = 999) 
# adult vs nymph sig diff
pair.adonis.life$Comparison <- "forest"
