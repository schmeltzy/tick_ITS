##Overview: Beta diversity analyses
#load necessary packages
#beta diversity

library(phyloseq)
library(microViz)
library(dplyr)
library(ggplot2)
library(ggpubfigs)
library(vegan)
library(pairwiseAdonis)
library(rstatix)

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
aitch <- dist_calc(ticks, dist = "robust.aitchison")
aitch.dist <- dist_get(aitch)


#get sample data from the matrix just to make it easier
aitch.sampledf <- data.frame(sample_data(aitch))

##permanova
adonis_habitat <- broom::tidy(adonis2(aitch.dist ~ Habitat, data = aitch.sampledf, by = "terms", perm = 999))  
#habitat p = 0.001 sig
adonis_habitat$Comparison <- "all"

adonis_treatment <- broom::tidy(adonis2(aitch.dist ~ Treatment, data = aitch.sampledf, by = "terms", perm = 999))  
#p = 0.001 sig
adonis_treatment$Comparison <- "all"

adonis_lifestage <- broom::tidy(adonis2(aitch.dist ~ Lifestage, data = aitch.sampledf, by = "terms", perm = 999))  
#p = 0.001 sig
adonis_lifestage$Comparison <- "all"

adonis_season <- broom::tidy(adonis2(aitch.dist ~ Season, data = aitch.sampledf, by = "terms", perm = 999))  
#p = 0.001 sig
adonis_season$Comparison <- "all"

adonis_sex <- broom::tidy(adonis2(aitch.dist ~ Sex, data = aitch.sampledf, na.action = na.omit, by = "terms", perm = 999))  
# p = 0.001 sig
adonis_sex$Comparison <- "all"

rbind(adonis_habitat, adonis_treatment, adonis_season, adonis_lifestage, adonis_sex) -> adonis.all
write.csv(adonis.all, here::here("output/adonis_all.csv"))


#pairwise overall
pair.adonis.habitat <- pairwise.adonis(aitch.dist, aitch.sampledf$Habitat, p.adjust.m = "fdr", perm = 999) 
pair.adonis.habitat$Comparison <- "Habitat"

pair.adonis.treat <- pairwise.adonis(aitch.dist, aitch.sampledf$Treatment, p.adjust.m = "fdr", perm = 999) 
pair.adonis.treat$Comparison <- "Treatment"

pair.adonis.lifestage <- pairwise.adonis(aitch.dist, aitch.sampledf$Lifestage, p.adjust.m = "fdr", perm = 999) 
pair.adonis.lifestage$Comparison <- "Lifestage"

#save overall pairwise results
rbind(pair.adonis.habitat, pair.adonis.treat, pair.adonis.lifestage) -> adonis.pairwise
write.csv(adonis.pairwise, here::here("output/adonis_pairwise.csv"))


##forested subset
aitch.forest <- dist_calc(forest, dist = "robust.aitchison")
aitch.dist.forest <- dist_get(aitch.forest)

#get sample data from the matrix just to make it easier
forest.sampledf <- data.frame(sample_data(forest))


forest_treatment <- broom::tidy(adonis2(aitch.dist.forest ~ Treatment, data = forest.sampledf, by = "terms", perm = 999))
#p = 0.001 sig
forest_treatment$Comparison <- "forested"

forest_season <- broom::tidy(adonis2(aitch.dist.forest ~ Treatment / Season, data = forest.sampledf, by = "terms", perm = 999))
#p = 0.001 sig
forest_season$Comparison <- "forested"

forest_lifestage <- broom::tidy(adonis2(aitch.dist.forest ~ Treatment / Lifestage, data = forest.sampledf, by = "terms", perm = 999))
#p = 0.001 sig
forest_lifestage$Comparison <- "forested"

#need to remove NAs for sex differences
forest_sex <- broom::tidy(adonis2(aitch.dist.forest ~ Treatment / Sex, data = forest.sampledf, na.action = na.omit, by = "terms", perm = 999))  
forest_sex$Comparison <- "forested" #NS

rbind(forest_treatment, forest_season, forest_lifestage, forest_sex) -> adonis.forest
write.csv(adonis.forest, here::here("output/adonis_forest.csv"))

## Pairwise adonis for multi-level comparisons
pair.forest.season <- pairwise.adonis(aitch.dist.forest, forest.sampledf$Season, p.adjust.m = "fdr", perm = 999) 
pair.forest.season$Comparison <- "forest"

pair.forest.life <- pairwise.adonis(aitch.dist.forest, forest.sampledf$Lifestage, p.adjust.m = "fdr", perm = 999) 
# adult vs nymph sig diff
pair.forest.life$Comparison <- "forest"

pair.forest.all <- rbind(pair.forest.season,pair.forest.life)
write.csv(pair.forest.all, here::here("output/pairwise_forest.csv"))


##adonis on open understory subset
aitch.open <- dist_calc(open, dist = "robust.aitchison")
aitch.dist.open <- dist_get(aitch.open)

#get sample data from the matrix just to make it easier
open.sampledf <- data.frame(sample_data(open))

open_treatment <- broom::tidy(adonis2(aitch.dist.open ~ Treatment, data = open.sampledf, by = "terms", perm = 999))
#p = 0.001 sig
open_treatment$Comparison <- "open"

open_season <- broom::tidy(adonis2(aitch.dist.open ~ Treatment / Season, data = open.sampledf, by = "terms", perm = 999))
#p = 0.001 sig
open_season$Comparison <- "open"

#need to remove NAs for sex differences
open_sex <- broom::tidy(adonis2(aitch.dist.open ~ Treatment / Sex, data = open.sampledf, na.action = na.omit, by = "terms", perm = 999))  
#p = 0.03 sig
open_sex$Comparison <- "open" 

rbind(open_treatment, open_season, open_sex) -> adonis.open
write.csv(adonis.open, here::here("output/adonis_open.csv"))

#pairwise open for treatment adjust p-values
pair.open.treatment <- pairwise.adonis(aitch.dist.open, open.sampledf$Treatment, p.adjust.m = "fdr", perm = 999) 
pair.open.treatment$Comparison <- "open"

write.csv(pair.open.treatment, here:: here("output/pair_open.treatment.csv"))

##make some plots
#robust clr (same as first step of robust.aitchison) for PCA ordination

#all ticks
rclr <- tax_transform(ticks, trans = "rclr") %>%
  ord_calc(method = "PCA")

#make sure things are in order
rclr@sam_data$Lifestage <- factor(rclr@sam_data$Lifestage, c("larva", "nymph", "adult"), ordered = TRUE)
rclr@sam_data$Treatment <- factor(rclr@sam_data$Treatment, c("unmanaged", "burned", "mowed","unmowed"), ordered = TRUE)


pca.all <- rclr %>%
  ord_plot( 
    colour = "Treatment",
    shape = "Habitat",
    plot_taxa = FALSE, 
    auto_caption = NA 
  ) + stat_ellipse(aes(group = Treatment, color = Treatment), linewidth = 1.5) + theme_bw() + scale_color_manual(values = friendly_pal("nickel_five")) +
  theme_bw(base_line_size = 1.5, base_rect_size = 1) + 
  theme(axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 14)) +
  guides(color = guide_legend(title = "Treatment")) + theme(legend.position = "right")
# facet_wrap(~Lifestage)
# + theme(legend.position = c(0.8, 0.08)) + theme(legend.key.size = unit(1, "mm")) +
pca.all


##separate by habitat
#forested
rclr.forest <- tax_transform(forest, trans = "rclr") %>%
  ord_calc(method = "PCA")

#make sure things are in order
rclr.forest@sam_data$Lifestage <- factor(rclr.forest@sam_data$Lifestage, c("larva", "nymph", "adult"), ordered = TRUE)
rclr.forest@sam_data$Treatment <- factor(rclr.forest@sam_data$Treatment, c("unmanaged", "burned", "mowed","unmowed"), ordered = TRUE)

#plot to show differences by lifestage since it was significant
pca.forest <- rclr.forest %>%
  ord_plot(
    shape = "Treatment",
    colour = "Lifestage",
    plot_taxa = FALSE, 
    auto_caption = NA 
  ) + stat_ellipse(aes(group = Lifestage, color = Lifestage), linewidth = 1.5) + theme_bw() + scale_color_manual(values = friendly_pal("retro_four")) +
  theme_bw(base_line_size = 1.5, base_rect_size = 1) + 
  theme(axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 14)) +
  guides(color = guide_legend(title = "Lifestage")) + theme(legend.position = "right")
 #facet_wrap(~Season)
pca.forest


#open (differences in treatment since it was sig)
rclr.open <- tax_transform(open, trans = "rclr") %>%
  ord_calc(method = "PCA")

#make sure things are in order
rclr.open@sam_data$Lifestage <- factor(rclr.open@sam_data$Lifestage, c("larva", "nymph", "adult"), ordered = TRUE)
rclr.open@sam_data$Treatment <- factor(rclr.open@sam_data$Treatment, c("unmanaged", "burned", "mowed","unmowed"), ordered = TRUE)

pca.open <- rclr.open %>%
  ord_plot(
    shape = "Season",
    colour = "Treatment",
    plot_taxa = FALSE, 
    auto_caption = NA 
  ) + stat_ellipse(aes(group = Treatment, color = Treatment), linewidth = 1.5) + theme_bw() + scale_color_manual(values = friendly_pal("retro_four")) +
  theme_bw(base_line_size = 1.5, base_rect_size = 1) + 
  theme(axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 14)) +
  guides(color = guide_legend(title = "Treatment")) + theme(legend.position = "right")
# facet_wrap(~Lifestage)
pca.open


#=================================

  
##calculate beta dispersion
#habitat
disp_habitat <- betadisper(aitch.dist, aitch.sampledf$Habitat, bias.adjust = TRUE, type = "centroid")
broom::tidy(anova(disp_habitat)) #disp not sig diff across habitat (p=0.854)
stats_habitat <- broom::tidy(p.adjust(permutest(betadisper(aitch.dist, aitch.sampledf$Habitat, type = "centroid", 
                                                         bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_habitat <- add_significance(stats_habitat, p.col = "x")
stats_habitat$comparisons <- "Habitat"

#Life
disp_life <- betadisper(aitch.dist, aitch.sampledf$Lifestage, bias.adjust = TRUE, type = "centroid")
broom::tidy(anova(disp_life)) #disp not sig diff across life (p=0.133)
stats_life <- broom::tidy(p.adjust(permutest(betadisper(aitch.dist, aitch.sampledf$Lifestage, type = "centroid", 
                                                           bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_life <- add_significance(stats_life, p.col = "x")
stats_life$comparisons <- "Lifestage"

#treatment
disp_treatment <- betadisper(aitch.dist, aitch.sampledf$Treatment, bias.adjust = TRUE, type = "centroid")
broom::tidy(anova(disp_treatment)) #sig across treatment (p=0.0304)
stats_treatment <- broom::tidy(p.adjust(permutest(betadisper(aitch.dist, aitch.sampledf$Treatment, type = "centroid", 
                                                             bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_treatment <- add_significance(stats_treatment, p.col = "x")
stats_treatment$comparisons <- "Treatment"

#season
disp_season <- betadisper(aitch.dist, aitch.sampledf$Season, bias.adjust = TRUE, type = "centroid")
broom::tidy(anova(disp_season)) #sig across season (p=0.000227)
stats_season <- broom::tidy(p.adjust(permutest(betadisper(aitch.dist, aitch.sampledf$Season, type = "centroid", 
                                                             bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_season <- add_significance(stats_season, p.col = "x")
stats_season$comparisons <- "Season"

#sex
disp_sex <- betadisper(aitch.dist, aitch.sampledf$Sex, bias.adjust = TRUE, type = "centroid")
broom::tidy(anova(disp_sex)) #sig across sex (p=0.761)
stats_sex <- broom::tidy(p.adjust(permutest(betadisper(aitch.dist, aitch.sampledf$Sex, type = "centroid", 
                                                          bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_sex <- add_significance(stats_sex, p.col = "x")
stats_sex$comparisons <- "Sex"

#combine dispersion results

stats_all <- rbind(stats_habitat, stats_life, stats_season,stats_sex, stats_treatment)
colnames(stats_all) <- c("comparison","p.adj","p.adj.signif","variable")


##within habitat dispersion for significant main effects (just season)
##forested season
disp_forest_season <- betadisper(aitch.dist.forest, forest.sampledf$Season, bias.adjust = TRUE, type = "centroid")
stats_forest_season <- broom::tidy(anova(disp_forest_season)) #disp still sig diff across season (p=0.002)
stats_forest_season <- broom::tidy(p.adjust(permutest(betadisper(aitch.dist.forest, forest.sampledf$Season, type = "centroid", 
                                                       bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_forest_season <- add_significance(stats_forest_season, p.col = "x")
stats_forest_season$comparisons <- "forested"

##open season
disp_open_season <- betadisper(aitch.dist.open, open.sampledf$Season, bias.adjust = TRUE, type = "centroid")
stats_open_season <- broom::tidy(anova(disp_open_season)) #disp still sig diff across season (p=8.51E-06)
stats_open_season <- broom::tidy(p.adjust(permutest(betadisper(aitch.dist.open, open.sampledf$Season, type = "centroid", 
                                                                 bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_open_season <- add_significance(stats_open_season, p.col = "x")
stats_open_season$comparisons <- "open"


#combine dispersion results
stats_habitat <- rbind(stats_forest_season, stats_open_season)
colnames(stats_habitat) <- c("comparison","p.adj","p.adj.signif","variable")

write.csv(stats_all, here::here("output/dispersion_all_stats.csv"))
write.csv(stats_habitat, here::here("output/dispersion_habitat_stats.csv"))


#plot dispersion
#get dispersion distance dfs
disp_df <- as.data.frame(disp_all$distances)
colnames(disp_df) <- c("distances")
disp_df$Habitat <- paste0(aitch.sampledf$Habitat)
disp_df$Season <- paste0(aitch.sampledf$Season)
disp_df$Lifestage <- paste0(aitch.sampledf$Lifestage)
disp_df$Sex <- paste0(aitch.sampledf$Sex)
disp_df$Treatment <- paste0(aitch.sampledf$Treatment)

#get everything in order again
disp_df$Lifestage <- factor(disp_df$Lifestage, c("larva", "nymph", "adult"), ordered = TRUE)
disp_df$Treatment <- factor(disp_df$Treatment, c("unmanaged", "burned", "mowed","unmowed"), ordered = TRUE)


plot_disp <-  ggplot(disp_df, aes(x = Treatment, y = distances, color = Treatment))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + theme_bw(base_line_size = 1.5, base_rect_size = 1.75)

plot_disp <- plot_disp + geom_point(aes(color = Treatment), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) +
  ylab("Distance to Centroid") + scale_color_manual(values = friendly_pal("nickel_five"))
plot_disp <- plot_disp + theme(axis.text = element_text(face = "bold", size = 14), 
                               axis.title = element_text(face = "bold", size = 14), 
                               title = element_text(face = "bold"), axis.title.x = element_blank(), 
                               axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(title = "Treatment")) + theme(strip.text = element_text(face = "bold", size = 14))
plot_disp

#+
 # stat_pvalue_manual(disp_stats_df, label = "p.adj.signif", hide.ns = TRUE, size = 6) + 
  #theme(legend.position = c(0.85, 0.9)) + theme(legend.key.size = unit(1, "mm"))
#+ facet_grid(~factor(Location, levels=c('North', 'East', 'West'))) 
