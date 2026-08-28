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
library(ggpubr)

#set seed
set.seed(123)

#load in unrarefied phyloseq object
ps <- readRDS("phyloseq_obj.rds")

#change the variable names
ps@sam_data$Treatment <- gsub("unmanaged", "unburned", ps@sam_data$Treatment)
ps@sam_data$Treatment <- gsub("mowed", "managed", ps@sam_data$Treatment)
ps@sam_data$Treatment <- gsub("unmowed", "unmanaged", ps@sam_data$Treatment)


###let's subset to just look at ticks (not soil)
ticks <- subset_samples(ps, Type=="BLT")

##make another column with "groups"
# Extract sample data as a standard data frame
metadata <- data.frame(sample_data(ticks))

# Combine columns into a new column called 'Group'
metadata$Groups <- paste(metadata$Lifestage, metadata$Season, sep = "-")

# Assign the updated data frame back into the phyloseq object
sample_data(ticks) <- sample_data(metadata)

#and subset ticks by habitat type
forest <- subset_samples(ticks, Habitat=="forested")
open <- subset_samples(ticks, Habitat=="open")

#and also make a subset for just adults
adults <- subset_samples(ticks, Lifestage=="adult")

# and just spring forest
#and also make a subset for just adults
spring.forest <- subset_samples(forest, Season=="spring")

# does habitat, treatment, lifestage, season, or sex affect beta diversity?
#use aitchison distance for beta (clr transform followed by euclidean distance matrix)
aitch <- dist_calc(ticks, dist = "robust.aitchison")
aitch.dist <- dist_get(aitch)


#get sample data from the matrix just to make it easier
aitch.sampledf <- data.frame(sample_data(aitch))

##permanova
adonis_habitat <- broom::tidy(adonis2(aitch.dist ~ Habitat, data = aitch.sampledf, by = "terms", perm = 999))  
adonis_habitat$Comparison <- "all"

adonis_treatment <- broom::tidy(adonis2(aitch.dist ~ Treatment, data = aitch.sampledf, by = "terms", perm = 999))  
adonis_treatment$Comparison <- "all"

adonis_lifestage <- broom::tidy(adonis2(aitch.dist ~ Lifestage, data = aitch.sampledf, by = "terms", perm = 999))  
adonis_lifestage$Comparison <- "all"

adonis_season <- broom::tidy(adonis2(aitch.dist ~ Season, data = aitch.sampledf, by = "terms", perm = 999))  
adonis_season$Comparison <- "all"

adonis_sex <- broom::tidy(adonis2(aitch.dist ~ Sex, data = aitch.sampledf, na.action = na.omit, by = "terms", perm = 999))  
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
forest_treatment$Comparison <- "forested"

forest_season <- broom::tidy(adonis2(aitch.dist.forest ~ Treatment / Season, data = forest.sampledf, by = "terms", perm = 999))
forest_season$Comparison <- "forested"

forest_lifestage <- broom::tidy(adonis2(aitch.dist.forest ~ Treatment / Lifestage, data = forest.sampledf, by = "terms", perm = 999))
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
open_treatment$Comparison <- "open"

open_season <- broom::tidy(adonis2(aitch.dist.open ~ Treatment / Season, data = open.sampledf, by = "terms", perm = 999))
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

rclr.adults <- tax_transform(adults, trans = "rclr") %>%
  ord_calc(method = "PCA")

rclr.springforest <- tax_transform(spring.forest, trans = "rclr") %>%
  ord_calc(method = "PCA")

#make sure things are in order
rclr@sam_data$Lifestage <- factor(rclr@sam_data$Lifestage, c("larva", "nymph", "adult"), ordered = TRUE)
rclr@sam_data$Treatment <- factor(rclr@sam_data$Treatment, c("unburned", "burned", "managed","unmanaged"), ordered = TRUE)

rclr.adults@sam_data$Treatment <- factor(rclr.adults@sam_data$Treatment, c("unburned", "burned", "managed","unmanaged"), ordered = TRUE)

rclr.springforest@sam_data$Lifestage <- factor(rclr.springforest@sam_data$Lifestage, c("larva", "nymph", "adult"), ordered = TRUE)
rclr.springforest@sam_data$Treatment <- factor(rclr.springforest@sam_data$Treatment, c("unburned", "burned", "managed","unmanaged"), ordered = TRUE)


pca.all <- rclr %>%
  ord_plot( 
    colour = "Treatment",
    shape = "Lifestage",
    plot_taxa = FALSE, 
    auto_caption = NA 
  ) + stat_ellipse(aes(group = Treatment, color = Treatment), linewidth = 1.5) + theme_bw() +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  theme_bw(base_line_size = 1.5, base_rect_size = 1) + 
  theme(axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 14)) +
  #facet_wrap(~Season)+
  guides(color = guide_legend(title = "Treatment")) + theme(legend.position = "right")
# + theme(legend.position = c(0.8, 0.08)) + theme(legend.key.size = unit(1, "mm")) +
pca.all

### just to try

# pca.spring.forest <- rclr.springforest %>%
#   ord_plot( 
#     colour = "Lifestage",
#     shape = "Treatment",
#     plot_taxa = FALSE, 
#     auto_caption = NA 
#   ) + stat_ellipse(aes(group = Lifestage, color = Lifestage), linewidth = 1.5) + theme_bw() +
#   scale_color_manual(values = friendly_pal("retro_four")) +
#   theme_bw(base_line_size = 1.5, base_rect_size = 1) + 
#   theme(axis.text = element_text(face = "bold", size = 14),
#         axis.title = element_text(face = "bold", size = 14), 
#         title = element_text(face = "bold")) +
#   theme(strip.text = element_text(face = "bold", size = 14)) +
#   facet_wrap(~Treatment)+
#   guides(color = guide_legend(title = "Lifestage")) + theme(legend.position = "right")
# # + theme(legend.position = c(0.8, 0.08)) + theme(legend.key.size = unit(1, "mm")) +
# pca.spring.forest


###

##overlay dispersion vectors
#overlay forest dispersion vectors
# 1) Extract plotted coordinates (site scores) with the SAME scaling and axes
scores_mat.all <- vegan::scores(
  ord_get(rclr),      # ordination object from psExtra
  display = "sites",
  choices = 1:2,             # match ord_plot axes
  scaling = 2                # match ord_plot default scaling
)                            # vegan::scores is the recommended accessor for coordinates

scores_df.all <- tibble(
  .sample_name = rownames(scores_mat.all),
  Axis1 = scores_mat.all[, 1],
  Axis2 = scores_mat.all[, 2]
) %>%
  # 2) Join the sample metadata in a tidy way
  left_join(samdat_tbl(rclr), by = ".sample_name")

# 3) Compute centroids BY Treament (since colour/ellipses use Treatment)
centroids.all <- scores_df.all %>%
  group_by(Treatment) %>%
  summarise(cx = mean(Axis1, na.rm = TRUE),
            cy = mean(Axis2, na.rm = TRUE),
            .groups = "drop")

segments.all <- scores_df.all %>%
  left_join(centroids.all, by = "Treatment") %>%
  mutate(dispersion = sqrt((Axis1 - cx)^2 + (Axis2 - cy)^2))

# 4) Overlay on your existing plot (no inheritance of aesthetics)
pca.all <- pca.all +
  geom_segment(
    data = segments.all,
    aes(x = Axis1, y = Axis2, xend = cx, yend = cy, colour = Treatment),
    inherit.aes = FALSE,
    alpha = 0.6,
    linewidth = 0.5
  ) +
  geom_point(
    data = centroids.all,
    aes(x = cx, y = cy, colour = Treatment),
    inherit.aes = FALSE,
    shape = 4, size = 4, stroke = 1.5
  )
pca.all


##separate by habitat
#forested
rclr.forest <- tax_transform(forest, trans = "rclr") %>%
  ord_calc(method = "PCA")

#make sure things are in order
rclr.forest@sam_data$Lifestage <- factor(rclr.forest@sam_data$Lifestage, c("larva", "nymph", "adult"), ordered = TRUE)
rclr.forest@sam_data$Treatment <- factor(rclr.forest@sam_data$Treatment, c("unburned", "burned", "managed","unmanaged"), ordered = TRUE)


#plot to show differences by lifestage since it was significant
pca.forest.life <- rclr.forest %>%
  ord_plot(
    shape = "Season",
    colour = "Lifestage",
    plot_taxa = FALSE, 
    auto_caption = NA 
  ) + stat_ellipse(aes(group = Lifestage, color = Lifestage), linewidth = 1.5) + theme_bw() + scale_color_manual(values = friendly_pal("retro_four")) +
  theme_bw(base_line_size = 1.5, base_rect_size = 1) + 
  theme(axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 14)) +
  #facet_grid(~Season)+
  guides(color = guide_legend(title = "Lifestage", keyheight = unit(3, "mm"))) + theme(legend.position = "top")+
  theme(legend.box="vertical", legend.margin=margin())
pca.forest.life

#overlay forest dispersion vectors
# 1) Extract plotted coordinates (site scores) with the SAME scaling and axes
scores_mat.forest <- vegan::scores(
  ord_get(rclr.forest),      # ordination object from psExtra
  display = "sites",
  choices = 1:2,             # match ord_plot axes
  scaling = 2                # match ord_plot default scaling
)                            # vegan::scores is the recommended accessor for coordinates

scores_df.forest <- tibble(
  .sample_name = rownames(scores_mat.forest),
  Axis1 = scores_mat.forest[, 1],
  Axis2 = scores_mat.forest[, 2]
) %>%
  # 2) Join the sample metadata in a tidy way
  left_join(samdat_tbl(rclr.forest), by = ".sample_name")

# 3) Compute centroids BY Lifestage (since colour/ellipses use Lifestage)
centroids.forest.life <- scores_df.forest %>%
  group_by(Lifestage) %>%
  summarise(cx = mean(Axis1, na.rm = TRUE),
            cy = mean(Axis2, na.rm = TRUE),
            .groups = "drop")

segments.forest.life <- scores_df.forest %>%
  left_join(centroids.forest.life, by = "Lifestage") %>%
  mutate(dispersion = sqrt((Axis1 - cx)^2 + (Axis2 - cy)^2))

# 4) Overlay on your existing plot (no inheritance of aesthetics)
pca.forest.life <- pca.forest.life +
  geom_segment(
    data = segments.forest.life,
    aes(x = Axis1, y = Axis2, xend = cx, yend = cy, colour = Lifestage),
    inherit.aes = FALSE,
    alpha = 0.6,
    linewidth = 0.5
  ) +
  geom_point(
    data = centroids.forest.life,
    aes(x = cx, y = cy, colour = Lifestage),
    inherit.aes = FALSE,
    shape = 4, size = 4, stroke = 1.5
  )
pca.forest.life


##and plot for season since it was sig
pca.forest.season <- rclr.forest %>%
  ord_plot(
    shape = "Treatment",
    colour = "Season",
    plot_taxa = FALSE, 
    auto_caption = NA 
  ) + stat_ellipse(aes(group = Season, color = Season), linewidth = 1.5) + theme_bw() + scale_color_manual(values = friendly_pal("zesty_four")) +
  theme_bw(base_line_size = 1.5, base_rect_size = 1) + 
  theme(axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  guides(color = guide_legend(title = "Season"), fill = guide_legend(label.theme = element_text(size = 10, lineheight = 0.8)))+
  theme(strip.text = element_text(face = "bold", size = 14)) +
  guides(color = guide_legend(title = "Season", keyheight = unit(3, "mm"))) + theme(legend.position = "top")+
  theme(legend.box="vertical", legend.margin=margin())
pca.forest.season

#overlay forest dispersion vectors
# 1) Extract plotted coordinates (site scores) with the SAME scaling and axes
scores_mat.forest <- vegan::scores(
  ord_get(rclr.forest),      # ordination object from psExtra
  display = "sites",
  choices = 1:2,             # match ord_plot axes
  scaling = 2                # match ord_plot default scaling
)                            # vegan::scores is the recommended accessor for coordinates

scores_df.forest <- tibble(
  .sample_name = rownames(scores_mat.forest),
  Axis1 = scores_mat.forest[, 1],
  Axis2 = scores_mat.forest[, 2]
) %>%
  # 2) Join the sample metadata in a tidy way
  left_join(samdat_tbl(rclr.forest), by = ".sample_name")

# 3) Compute centroids BY Season (since colour/ellipses use Season)
centroids.forest.season <- scores_df.forest %>%
  group_by(Season) %>%
  summarise(cx = mean(Axis1, na.rm = TRUE),
            cy = mean(Axis2, na.rm = TRUE),
            .groups = "drop")

segments.forest.season <- scores_df.forest %>%
  left_join(centroids.forest.season, by = "Season") %>%
  mutate(dispersion = sqrt((Axis1 - cx)^2 + (Axis2 - cy)^2))

# 4) Overlay on your existing plot (no inheritance of aesthetics)
pca.forest.season <- pca.forest.season +
  geom_segment(
    data = segments.forest.season,
    aes(x = Axis1, y = Axis2, xend = cx, yend = cy, colour = Season),
    inherit.aes = FALSE,
    alpha = 0.6,
    linewidth = 0.5
  ) +
  geom_point(
    data = centroids.forest.season,
    aes(x = cx, y = cy, colour = Season),
    inherit.aes = FALSE,
    shape = 4, size = 4, stroke = 1.5
  )
pca.forest.season



###open (differences in treatment since it was sig)
rclr.open <- tax_transform(open, trans = "rclr") %>%
  ord_calc(method = "PCA")

#make sure things are in order
rclr.open@sam_data$Lifestage <- factor(rclr.open@sam_data$Lifestage, c("larva", "nymph", "adult"), ordered = TRUE)
rclr.open@sam_data$Treatment <- factor(rclr.open@sam_data$Treatment, c("unburned", "burned", "managed","unmanaged"), ordered = TRUE)


pca.open <- rclr.open %>%
  ord_plot(
    shape = "Season",
    colour = "Treatment",
    plot_taxa = FALSE, 
    auto_caption = NA 
  ) + stat_ellipse(aes(group = Treatment, color = Treatment), linewidth = 1.5) + theme_bw() + scale_color_manual(values = friendly_pal("ito_seven")) +
  theme_bw(base_line_size = 1.5, base_rect_size = 1) + 
  theme(axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 14)) +
  guides(color = guide_legend(title = "Treatment", keyheight = unit(3, "mm")))+ theme(legend.position = "top")+
  theme(legend.box="vertical", legend.margin=margin())
pca.open


##overlay dispersion vectors
#overlay forest dispersion vectors
# 1) Extract plotted coordinates (site scores) with the SAME scaling and axes
scores_mat.open <- vegan::scores(
  ord_get(rclr.open),      # ordination object from psExtra
  display = "sites",
  choices = 1:2,             # match ord_plot axes
  scaling = 2                # match ord_plot default scaling
)                            # vegan::scores is the recommended accessor for coordinates

scores_df.open <- tibble(
  .sample_name = rownames(scores_mat.open),
  Axis1 = scores_mat.open[, 1],
  Axis2 = scores_mat.open[, 2]
) %>%
  # 2) Join the sample metadata in a tidy way
  left_join(samdat_tbl(rclr.open), by = ".sample_name")

# 3) Compute centroids BY Treament (since colour/ellipses use Treatment)
centroids.open <- scores_df.open %>%
  group_by(Treatment) %>%
  summarise(cx = mean(Axis1, na.rm = TRUE),
            cy = mean(Axis2, na.rm = TRUE),
            .groups = "drop")

segments.open <- scores_df.open %>%
  left_join(centroids.open, by = "Treatment") %>%
  mutate(dispersion = sqrt((Axis1 - cx)^2 + (Axis2 - cy)^2))

# 4) Overlay on your existing plot (no inheritance of aesthetics)
pca.open <- pca.open +
  geom_segment(
    data = segments.open,
    aes(x = Axis1, y = Axis2, xend = cx, yend = cy, colour = Treatment),
    inherit.aes = FALSE,
    alpha = 0.6,
    linewidth = 0.5
  ) +
  geom_point(
    data = centroids.open,
    aes(x = cx, y = cy, colour = Treatment),
    inherit.aes = FALSE,
    shape = 4, size = 4, stroke = 1.5
  )
pca.open

##save our individual plots
ggplot2::ggsave(here::here("output/betas_all.png"), pca.all,
                height = 400, width = 600, units = "mm",
                scale = 0.5, dpi = 1000)

ggplot2::ggsave(here::here("output/betas_forest.life.png"), pca.forest.life,
                height = 400, width = 600, units = "mm",
                scale = 0.5, dpi = 1000)

ggplot2::ggsave(here::here("output/betas_forest.season.png"), pca.forest.season,
                height = 400, width = 600, units = "mm",
                scale = 0.5, dpi = 1000)

ggplot2::ggsave(here::here("output/betas_open.png"), pca.open,
                height = 400, width = 600, units = "mm",
                scale = 0.5, dpi = 1000)


#=================================

  
##calculate beta dispersion
#all
disp_all <- betadisper(aitch.dist, aitch.sampledf$Habitat, bias.adjust = TRUE, type = "centroid")
broom::tidy(anova(disp_all)) #disp not sig diff across all (p=0.854)
stats_all <- tibble::enframe(p.adjust(permutest(betadisper(aitch.dist, aitch.sampledf$Habitat, type = "centroid", 
                                                           bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_all <- add_significance(stats_all, p.col = "value")
stats_all$comparisons <- "all"

#habitat
disp_habitat <- betadisper(aitch.dist, aitch.sampledf$Habitat, bias.adjust = TRUE, type = "centroid")
broom::tidy(anova(disp_habitat)) #disp not sig diff across habitat (p=0.854)
stats_habitat <- tibble::enframe(p.adjust(permutest(betadisper(aitch.dist, aitch.sampledf$Habitat, type = "centroid", 
                                                         bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_habitat <- add_significance(stats_habitat, p.col = "value")
stats_habitat$comparisons <- "Habitat"

stats_habitat$group1 <- sapply(strsplit(stats_habitat$name, "-"), `[`,1)
stats_habitat$group2 <- sapply(strsplit(stats_habitat$name, "-"), `[`,2)

#Life
disp_life <- betadisper(aitch.dist, aitch.sampledf$Lifestage, bias.adjust = TRUE, type = "centroid")
broom::tidy(anova(disp_life)) #disp not sig diff across life (p=0.133)
stats_life <- tibble::enframe(p.adjust(permutest(betadisper(aitch.dist, aitch.sampledf$Lifestage, type = "centroid", 
                                                           bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_life <- add_significance(stats_life, p.col = "value")
stats_life$comparisons <- "Lifestage"

stats_life$group1 <- sapply(strsplit(stats_life$name, "-"), `[`,1)
stats_life$group2 <- sapply(strsplit(stats_life$name, "-"), `[`,2)

#treatment
disp_treatment <- betadisper(aitch.dist, aitch.sampledf$Treatment, bias.adjust = TRUE, type = "centroid")
broom::tidy(anova(disp_treatment)) #sig across treatment (p=0.0304)
stats_treatment <- tibble::enframe(p.adjust(permutest(betadisper(aitch.dist, aitch.sampledf$Treatment, type = "centroid", 
                                                             bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_treatment <- add_significance(stats_treatment, p.col = "value")
stats_treatment$comparisons <- "Treatment"

stats_treatment$group1 <- sapply(strsplit(stats_treatment$name, "-"), `[`,1)
stats_treatment$group2 <- sapply(strsplit(stats_treatment$name, "-"), `[`,2)


#season
disp_season <- betadisper(aitch.dist, aitch.sampledf$Season, bias.adjust = TRUE, type = "centroid")
broom::tidy(anova(disp_season)) #sig across season (p=0.000227)
stats_season <- tibble::enframe(p.adjust(permutest(betadisper(aitch.dist, aitch.sampledf$Season, type = "centroid", 
                                                             bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_season <- add_significance(stats_season, p.col = "value")
stats_season$comparisons <- "Season"

stats_season$group1 <- sapply(strsplit(stats_season$name, "-"), `[`,1)
stats_season$group2 <- sapply(strsplit(stats_season$name, "-"), `[`,2)

#sex
disp_sex <- betadisper(aitch.dist, aitch.sampledf$Sex, bias.adjust = TRUE, type = "centroid")
broom::tidy(anova(disp_sex)) #sig across sex (p=0.761)
stats_sex <- tibble::enframe(p.adjust(permutest(betadisper(aitch.dist, aitch.sampledf$Sex, type = "centroid", 
                                                          bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_sex <- add_significance(stats_sex, p.col = "value")
stats_sex$comparisons <- "Sex"

stats_sex$group1 <- sapply(strsplit(stats_sex$name, "-"), `[`,1)
stats_sex$group2 <- sapply(strsplit(stats_sex$name, "-"), `[`,2)


#combine dispersion results
stats_all <- rbind(stats_habitat, stats_life, stats_season, stats_sex, stats_treatment)
colnames(stats_all) <- c("comparison","p.adj","p.adj.signif","variable","group1","group2")
write.csv(stats_all, here::here("output/dispersion_stats_all.csv"))

stats_all$y.position <- 22

##within habitat dispersion for significant main effects (just season since t)
##forested season
disp_forest_season <- betadisper(aitch.dist.forest, forest.sampledf$Season, bias.adjust = TRUE, type = "centroid")
stats_forest_season <- broom::tidy(anova(disp_forest_season)) #disp still sig diff across season (p=0.002)
stats_forest_season <- tibble::enframe(p.adjust(permutest(betadisper(aitch.dist.forest, forest.sampledf$Season, type = "centroid", 
                                                       bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_forest_season <- add_significance(stats_forest_season, p.col = "value")
stats_forest_season$comparisons <- "forested"

stats_forest_season$group1 <- sapply(strsplit(stats_forest_season$name, "-"), `[`,1)
stats_forest_season$group2 <- sapply(strsplit(stats_forest_season$name, "-"), `[`,2)


##open season
disp_open_season <- betadisper(aitch.dist.open, open.sampledf$Season, bias.adjust = TRUE, type = "centroid")
stats_open_season <- broom::tidy(anova(disp_open_season)) #disp still sig diff across season (p=8.51E-06)
stats_open_season <- tibble::enframe(p.adjust(permutest(betadisper(aitch.dist.open, open.sampledf$Season, type = "centroid", 
                                                                 bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_open_season <- add_significance(stats_open_season, p.col = "value")
stats_open_season$comparisons <- "open"

stats_open_season$group1 <- sapply(strsplit(stats_open_season$name, "-"), `[`,1)
stats_open_season$group2 <- sapply(strsplit(stats_open_season$name, "-"), `[`,2)


#combine dispersion results
stats_habitat <- rbind(stats_forest_season, stats_open_season)
colnames(stats_habitat) <- c("comparison","p.adj","p.adj.signif","habitat","group1","group2")

write.csv(stats_habitat, here::here("output/dispersion_habitat_stats.csv"))

stats_habitat$y.position <- 23

#and for open_treatment
disp_open_treat <- betadisper(aitch.dist.open, open.sampledf$Treatment, bias.adjust = TRUE, type = "centroid")
stats_open_treat <- broom::tidy(anova(disp_open_treat)) 
stats_open_treat <- tibble::enframe(p.adjust(permutest(betadisper(aitch.dist.open, open.sampledf$Treatment, type = "centroid", 
                                                               bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_open_treat <- add_significance(stats_open_treat, p.col = "value")
stats_open_treat$comparisons <- "open"

stats_open_treat$group1 <- sapply(strsplit(stats_open_treat$name, "-"), `[`,1)
stats_open_treat$group2 <- sapply(strsplit(stats_open_treat$name, "-"), `[`,2)

#write dispersion open stats
colnames(stats_open_treat) <- c("comparison","p.adj","p.adj.signif","treatment","group1","group2")
write.csv(stats_open_treat, here::here("output/dispersion_open_treat_stats.csv"))


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
disp_df$Treatment <- factor(disp_df$Treatment, c("unburned", "burned", "managed", "unmanaged"), ordered = TRUE)

#plot dispersion
#do disp treatment first to match betas
colnames(stats_treatment) <- c("comparison","p.adj","p.adj.signif","variable","group1","group2")
stats_treatment$y.position <- 22

#then plot
plot_disp_treat <-  ggplot(disp_df, aes(x = Treatment, y = distances, color = Treatment))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + theme_bw(base_line_size = 1.5, base_rect_size = 1.75)

plot_disp_treat <- plot_disp_treat + geom_point(aes(color = Treatment), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) +
  ylab("Distance to Centroid") + scale_color_manual(values = friendly_pal("ito_seven"))
plot_disp_treat <- plot_disp_treat + theme(axis.text = element_text(face = "bold", size = 14), 
                               axis.title = element_text(face = "bold", size = 14), 
                               title = element_text(face = "bold"), axis.title.x = element_text(face = "bold", color = "white"), 
                               axis.text.x = element_text(face = "bold", color = "white"), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(title = "Treatment")) + theme(strip.text = element_text(face = "bold", size = 14))+
stat_pvalue_manual(stats_treatment, label = "p.adj.signif", hide.ns = TRUE, size = 6)+
theme(legend.position = c(0.15, 0.9)) + theme(legend.key.size = unit(1, "mm"))
plot_disp_treat



disp_df_forested <- subset(disp_df, Habitat == "forested")

#then for lifestage
colnames(stats_life) <- c("comparison","p.adj","p.adj.signif","variable","group1","group2")
stats_life$y.position <- 22

#plot lifestage disp
plot_disp_life <-  ggplot(disp_df_forested, aes(x = Lifestage, y = distances, color = Lifestage))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + theme_bw(base_line_size = 1.5, base_rect_size = 1.75)

plot_disp_life <- plot_disp_life + geom_point(aes(color = Lifestage), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) +
  ylab("Distance to Centroid") + scale_color_manual(values = friendly_pal("retro_four"))
plot_disp_life <- plot_disp_life + theme(axis.text = element_text(face = "bold", size = 14), 
                                         axis.title = element_text(face = "bold", size = 14), 
                                         title = element_text(face = "bold"), axis.title.x = element_blank(), 
                                         axis.text.x = element_text(face = "bold", size = 12), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(title = "Lifestage")) + theme(strip.text = element_text(face = "bold", size = 14))+
  #stat_pvalue_manual(stats_life, label = "p.adj.signif", hide.ns = TRUE, size = 6)+
  guides(color = guide_legend(title = "Lifestage", keyheight = unit(3, "mm"))) + theme(legend.position = "top")+
  theme(legend.box="vertical", legend.margin=margin())
  #facet_grid(~Season)
plot_disp_life


#then for season (both open and forested)
colnames(stats_season) <- c("comparison","p.adj","p.adj.signif","variable","group1","group2")
stats_season$y.position <- 23

#forested first
colnames(stats_forest_season) <- c("comparison","p.adj","p.adj.signif","variable","group1","group2")
stats_forest_season$y.position <- 22


#then plot
plot_disp_season.forest <-  ggplot(disp_df_forested, aes(x = Season, y = distances, color = Season))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + theme_bw(base_line_size = 1.5, base_rect_size = 1.75)+
  stat_pvalue_manual(stats_forest_season, label = "p.adj.signif", hide.ns = TRUE, size = 6, inherit.aes = FALSE)

plot_disp_season.forest <- plot_disp_season.forest + geom_point(aes(color = Season), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) +
  ylab("Distance to Centroid") + scale_color_manual(values = friendly_pal("zesty_four"))
plot_disp_season.forest <- plot_disp_season.forest + theme(axis.text = element_text(face = "bold", size = 14), 
                                         axis.title = element_text(face = "bold", size = 14), 
                                         title = element_text(face = "bold"), axis.title.x = element_blank(), 
                                         axis.text.x = element_text(face = "bold", size = 12), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(title = "Season"), fill = guide_legend(label.theme = element_text(size = 10, lineheight = 0.8)))+
  theme(strip.text = element_text(face = "bold", size = 14))+
  stat_pvalue_manual(stats_forest_season, label = "p.adj.signif", hide.ns = TRUE, size = 6, inherit.aes = FALSE)+
  guides(color = guide_legend(title = "Season", keyheight = unit(3, "mm"))) +
  theme(legend.box="vertical", legend.margin=margin())
plot_disp_season.forest


#subset to "open" habitat
disp_df_open <- subset(disp_df, Habitat == "open")


#then disp open treatment
colnames(stats_open_treat) <- c("comparison","p.adj","p.adj.signif","variable","group1","group2")
stats_open_treat$y.position <- 21

#then plot
disp_df_open <- subset(disp_df, Habitat == "open")

plot_disp_treat.open <-  ggplot(disp_df_open, aes(x = Treatment, y = distances, color = Treatment))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + theme_bw(base_line_size = 1.5, base_rect_size = 1.75)

plot_disp_treat.open <- plot_disp_treat.open + geom_point(aes(color = Treatment), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) +
  ylab("Distance to Centroid") + scale_color_manual(values = friendly_pal("ito_seven"))
plot_disp_open_treat <- plot_disp_treat.open + theme(axis.text = element_text(face = "bold", size = 14), 
                                                  axis.title = element_text(face = "bold", size = 14), 
                                                  title = element_text(face = "bold"), axis.title.x = element_blank(), 
                                                  axis.text.x = element_text(face = "bold", size = 12), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(title = "Treatment")) + theme(strip.text = element_text(face = "bold", size = 14))+
 stat_pvalue_manual(stats_open_treat, label = "p.adj.signif", hide.ns = TRUE, size = 6)+
  guides(color = guide_legend(title = "Treatment", keyheight = unit(3, "mm"))) +
           theme(legend.position = "top")+
  theme(legend.box="vertical", legend.margin=margin())
plot_disp_open_treat


##figure out how we want to arrange our plots
#Dispersion 1
plot_pca_disp_all <- ggarrange(pca.all, plot_disp_treat, nrow = 1, ncol = 2, common.legend = TRUE, legend = "left") -> pca_disp_all
plot_pca_disp_all
ggplot2::ggsave(here::here("output/pca_disp_all.png"), pca_disp_all,
                height = 450, width = 600, units = "mm",
                scale = 0.5, dpi = 1000) #this is done


#Dispersion 2 a 
pca_disp_forest_life <- ggarrange(pca.forest.life, plot_disp_life, nrow = 2, ncol = 1, common.legend = TRUE, legend = "top")
pca_disp_forest_life <-  annotate_figure(
  pca_disp_forest_life,
  top = ggpubr::text_grob("Forested\n", face = "bold", size = 14, lineheight = 0.9))

plot(pca_disp_forest_life)
ggplot2::ggsave(here::here("output/pca_disp_forest_life.png"), pca_disp_forest_life,
                height = 450, width = 600, units = "mm",
                scale = 0.5, dpi = 1000) #a and b


#Dispersion 2 b
ggarrange(pca.forest.season, plot_disp_season.forest, nrow = 2, ncol = 1, common.legend = TRUE, legend = "top") -> pca_disp_forest_season
pca_disp_forest_season <-  annotate_figure(
  pca_disp_forest_season,
  top = ggpubr::text_grob("Forested\n", face = "bold", size = 14, lineheight = 0.9))
plot(pca_disp_forest_season)
ggplot2::ggsave(here::here("output/pca_disp_forest_season.png"), pca_disp_forest_season,
                height = 450, width = 600, units = "mm",
                scale = 0.5, dpi = 1000) #c and d


#Dispersion 2 c
ggarrange(pca.open, plot_disp_open_treat, nrow = 2, ncol = 1, common.legend = TRUE, legend = "top") -> pca_disp_open_treat
pca_disp_open_treat <-  annotate_figure(
  pca_disp_open_treat,
  top = ggpubr::text_grob("Open\n", face = "bold", size = 14, lineheight = 0.9))

plot(pca_disp_open_treat)
ggplot2::ggsave(here::here("output/pca_disp_open_treat.png"), pca_disp_open_treat,
                height = 450, width = 600, units = "mm",
                scale = 0.5, dpi = 1000) #e and f

#put all together
spacer <- ggplot() + theme_void()
pca_dispersion_combo <- ggarrange(pca_disp_forest_life, spacer, pca_disp_forest_season, spacer, 
                                  pca_disp_open_treat, nrow = 1, ncol = 5,
                                  widths = c(1, 0.12, 1, 0.12, 1),  # <- increase 0.06 for larger gaps
                                  common.legend = FALSE)
plot(pca_dispersion_combo)
ggplot2::ggsave(here::here("output/pca_dispersion_combo.png"), pca_dispersion_combo,
                height = 450, width = 600, units = "mm",
                scale = 0.5, dpi = 1000) #combo plot of PCAs and dispersion results for significant betas

