##Overview: Beta diversity analyses
#load necessary packages
#beta diversity

library(phyloseq)
library(dada2)
library(ggplot2)
library(microbiome)
library(vegan)

#set seed
set.seed(123)

#load in unrarefied phyloseq object
ps <- readRDS("phyloseq_obj.rds")

#will be using Aitchinson distance to asses beta, so need to transform
#using Euclidean distances

#clr transform
clr <- microbiome::transform(ps, 'clr')

#relative abundance transform
ps_rel <- filter_taxa(ps, function(x) mean(x) > 0.1, TRUE)
ps_rel <- transform_sample_counts(ps_rel, function(x) x / sum(x) )
ps_rel
#1523 taxa and 140 samples

#filter taxa
ps_rel_filt <- filter_taxa(ps_rel, function(x) mean(x) > 1e-5, TRUE)
#1503 taxa

rel_abun <- ps_rel@otu_table
rel_abun <- as.matrix(rel_abun)
write.table(rel_abun, file = "rel_abun_new.txt", sep = "\t")


# save ps objects

save(ps_rel, file = "ps_rel.RData") #rel abundance
save(clr, file = "ps_clr.RData") #clr transformed


