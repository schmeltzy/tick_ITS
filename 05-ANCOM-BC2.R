##Overview: ancom-bc2
#load necessary packages
#calculate differential abundance

library(phyloseq)
library(dplyr)
library(ggplot2)
library(fantaxtic)
library(ANCOMBC)
library(mia)

#set seed
set.seed(123)

#load in unrarefied phyloseq object
ps <- readRDS("phyloseq_obj.rds")

#rename the NAs to higher taxa level
# Fill in names for NA taxa, not including their rank
ps <- name_na_taxa(ps, include_rank = FALSE, na_label = "Unclassified <tax>")

###let's subset to just look at ticks (not soil)
ticks <- subset_samples(ps, Type=="BLT")

# convert phylo object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ticks)

#Are there differentially abundance taxa by habitat?

tse_all_ASV_output <- ancombc2(data = tse, assay_name = "counts", tax_level = "Genus",
                               fix_formula = "Habitat + Treatment + Season", 
                           p_adj_method = "fdr", group = "Habitat", pseudo_sens = TRUE, 
                           alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = TRUE,
                           mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100)) 

res <- tse_all_ASV_output$res




