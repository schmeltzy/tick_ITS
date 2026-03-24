##Overview: Alpha diversity analyses
#load necessary packages
library(ggplot2)
library(phyloseq)
library(here)
library(ggpubr)
library(dplyr)
library(rstatix)
library(broom)
library(tibble)

#set seed for reproducibility
set.seed(123)

## load phyloseq object
ps.rare <- readRDS(here::here("output/ps.rare.rds"))

#visualize alpha diversity just to make sure everything looks normal
plot_richness(ps.rare, x="Treatment", measures=c("Shannon", "Chao1"), color="Lifestage")


####let's subset to just look at ticks (not soil)
ticks <- subset_samples(ps.rare, Type=="BLT")

#calculate shannon diversity (richness and evenness)
alphadiv <- estimate_richness(ticks, measures = c("Shannon")) %>%
  rownames_to_column(var = "Unique_Specimen_ID") %>%
  left_join(as.data.frame(sample_data(ticks)), by = "Unique_Specimen_ID")

#Look at distribution and test for normality
hist(alphadiv$Shannon) #normal-ish
shapiro.test(alphadiv$Shannon) #is normal W = 0.98223, p = 0.07


