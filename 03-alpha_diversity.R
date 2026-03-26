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
library(lme4)
library(emmeans)

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
shapiro.test(alphadiv$Shannon) #is normal W = 0.98223, p = 0.07, will move forward with linear model

# linear model shannon diversity
mod.shannon <- lm(Shannon ~ Habitat*Treatment + Season + Lifestage, data=alphadiv)
shapiro.test(residuals(mod.shannon)) # residuals normal

lm_shan <- broom::tidy(anova(mod.shannon)) # Habitat and Treatment significant 
lm_shan$Metric <- "Shannon"

# pairwise comparisons
mod_shannon_pw <- emmeans(mod.shannon, pairwise ~ Treatment, adjust = "tukey") 
mod_shannon_pw <- broom::tidy(mod_shannon_pw$contrasts)
mod_shannon_pw <- rstatix::add_significance(data = as.data.frame(mod_shannon_pw), p.col = "adj.p.value")
mod_shannon_pw$Metric <- "Shannon"

# pairwise comparisons
mod_shannon_pw2 <- emmeans(mod.shannon, pairwise ~ Treatment|Habitat, adjust = "tukey") 
mod_shannon_pw2 <- broom::tidy(mod_shannon_pw2$contrasts)
mod_shannon_pw2 <- rstatix::add_significance(data = as.data.frame(mod_shannon_pw2), p.col = "adj.p.value")
mod_shannon_pw2$Metric <- "Shannon"


