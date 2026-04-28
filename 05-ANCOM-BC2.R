##Overview: ancom-bc2
#load necessary packages
#calculate differential abundance

library(phyloseq)
library(dplyr)
library(ggplot2)
library(fantaxtic)
library(ANCOMBC)
library(mia)
library(ggpubfigs)

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

#Are there differentially abundant taxa by habitat?

tse_all_ASV_output <- ancombc2(data = tse, assay_name = "counts", tax_level = "Genus",
                               fix_formula = "Habitat", 
                           p_adj_method = "fdr", group = "Habitat", pseudo_sens = TRUE, 
                           alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = FALSE,
                           mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100)) 

res <- tse_all_ASV_output$res

#make another df of the output we want (non-intercept columns)
df_habitat <- data.frame(c(res[1], res[3], res[5], res[7]), res[9], res[11], res[13], res[15], res[17])
df_habitat$comparison <- "Habitat"
colnames(df_habitat)[2] <- "LFC"
colnames(df_habitat)[3] <- "SE"
colnames(df_habitat)[4] <- "Wstat"
colnames(df_habitat)[5] <- "p"
colnames(df_habitat)[6] <- "p.adj"
colnames(df_habitat)[7] <- "diff"
colnames(df_habitat)[8] <- "pass_sens"
colnames(df_habitat)[9] <- "diff_robust"

#add another column for plot coloring for our non-sig taxa
df_habitat$Genus <- paste0(df_habitat$taxon)
df_habitat$Genus[df_habitat$p.adj >= 0.05] <- "Other"

#subset to keep diff abundant
#let's keep the ones that don't pass sensitivity first
df_habitat <- subset(df_habitat, !(diff == "TRUE" & pass_sens == "FALSE"))

#set some colors (colorblind friendly of course)
colors <- c(friendly_pal("glasbey_twelve"), "lightgrey")
#put genera in order (need to put "other" last so it's greyed out)
df_habitat$Genus <- factor(df_habitat$Genus, levels = c("Alternaria", "Pilidium", "Piskurozyma", "Plectosphaerella", "Ramularia", "Thyridium", "Unclassified Amphisphaeriales", "Unclassified Basidiomycota", "Unclassified Leotiomycetes", "Unclassified Mycosphaerellaceae", "Unclassified Tremellales", "Vishniacozyma", "Other"))

#and make our plot
vol_plot_all <- df_habitat %>%
  ggplot(aes(x = LFC,
             y = -log10(p.adj),
             color = Genus)) + 
  geom_point(size = 3.5, alpha = 0.8)+
  scale_color_manual(values = colors)

vol_plot_all <- vol_plot_all + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")

vol_plot_all + theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p-value") + labs(color = 'Genera') + theme(strip.text = element_text(face = "bold", size = 12)) 

##############

  
#now do separated by habitat
#are there differentially abundant taxa within habitat type based on lifestage, season, and sex?

#forested first
forest <- subset_samples(ticks, Habitat=="forested")
#convert phylo object to tree summarized experiment
tse_forest <- mia::convertFromPhyloseq(forest)

tse_forest_output <- ancombc2(data = tse_forest, assay_name = "counts", tax_level = "Genus",
                               fix_formula = "Treatment", 
                               p_adj_method = "fdr", group = "Treatment", pseudo_sens = TRUE, 
                               alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = FALSE,
                               mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100)) 

res_forest <- tse_forest_output$res

#make another df of the output we want (non-intercept columns)
df_forest_treat <- data.frame(c(res_forest[1], res_forest[3], res_forest[5], res_forest[7]), res_forest[9], res_forest[11], res_forest[13], res_forest[15], res_forest[17])
df_forest_treat$comparison <- "Treatment"
colnames(df_forest_treat)[2] <- "LFC"
colnames(df_forest_treat)[3] <- "SE"
colnames(df_forest_treat)[4] <- "Wstat"
colnames(df_forest_treat)[5] <- "p"
colnames(df_forest_treat)[6] <- "p.adj"
colnames(df_forest_treat)[7] <- "diff"
colnames(df_forest_treat)[8] <- "pass_sens"
colnames(df_forest_treat)[9] <- "diff_robust"

#add another column for plot coloring for our non-sig taxa
df_forest_treat$Genus <- paste0(df_forest_treat$taxon)
df_forest_treat$Genus[df_forest_treat$p.adj >= 0.05] <- "Other"

#subset to keep diff abundant
#let's keep the ones that don't pass sensitivity first
df_forest_treat <- subset(df_forest_treat, !(diff == "TRUE" & pass_sens == "FALSE"))

#set some colors (colorblind friendly of course)
#colors <- c(friendly_pal("glasbey_twelve"), "lightgrey")
#put genera in order (need to put "other" last so it's greyed out)
df_forest_treat$Genus <- factor(df_forest_treat$Genus, levels = c("Alternaria", "Pilidium", "Piskurozyma", "Plectosphaerella", "Ramularia", "Thyridium", "Unclassified Amphisphaeriales", "Unclassified Basidiomycota", "Unclassified Leotiomycetes", "Unclassified Mycosphaerellaceae", "Unclassified Tremellales", "Vishniacozyma", "Other"))

#and make our plot
vol_plot_forest_treat <- df_forest_treat %>%
  ggplot(aes(x = LFC,
             y = -log10(p.adj),
             color = Genus)) + 
  geom_point(size = 3.5, alpha = 0.8)+
  scale_color_manual(values = "lightgrey")

vol_plot_forest_treat <- vol_plot_forest_treat + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")

vol_plot_forest_treat + theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p-value") + labs(color = 'Genera') + theme(strip.text = element_text(face = "bold", size = 12)) 

### none pass the p-value adjustment for forested between treatments

###############
#move on to lifestage forested
tse_forest_life_output <- ancombc2(data = tse_forest, assay_name = "counts", tax_level = "Genus",
                              fix_formula = "Lifestage", 
                              p_adj_method = "fdr", group = "Lifestage", pseudo_sens = TRUE, 
                              alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = TRUE,
                              mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100)) 

res_forest_life <- tse_forest_life_output$res_pair

#make another df of the output we want (non-intercept columns with intercept being adult in this case)
df_forest_AvL <- data.frame(c(res_forest_life[1], res_forest_life[2], res_forest_life[5], res_forest_life[8]), res_forest_life[11], res_forest_life[14], res_forest_life[17], res_forest_life[20], res_forest_life[23])
df_forest_AvL$Habitat <- "forested"
df_forest_AvL$comparison <- "adult v. larva"
colnames(df_forest_AvL)[2] <- "LFC"
colnames(df_forest_AvL)[3] <- "SE"
colnames(df_forest_AvL)[4] <- "Wstat"
colnames(df_forest_AvL)[5] <- "p"
colnames(df_forest_AvL)[6] <- "p.adj"
colnames(df_forest_AvL)[7] <- "diff"
colnames(df_forest_AvL)[8] <- "pass_sens"
colnames(df_forest_AvL)[9] <- "diff_robust"

df_forest_AvN <- data.frame(c(res_forest_life[1], res_forest_life[3], res_forest_life[6], res_forest_life[9]), res_forest_life[12], res_forest_life[15], res_forest_life[18], res_forest_life[21], res_forest_life[24])
df_forest_AvN$Habitat <- "forested"
df_forest_AvN$comparison <- "adult v. nymph"
colnames(df_forest_AvN)[2] <- "LFC"
colnames(df_forest_AvN)[3] <- "SE"
colnames(df_forest_AvN)[4] <- "Wstat"
colnames(df_forest_AvN)[5] <- "p"
colnames(df_forest_AvN)[6] <- "p.adj"
colnames(df_forest_AvN)[7] <- "diff"
colnames(df_forest_AvN)[8] <- "pass_sens"
colnames(df_forest_AvN)[9] <- "diff_robust"

df_forest_NvL <- data.frame(c(res_forest_life[1], res_forest_life[4], res_forest_life[7], res_forest_life[10]), res_forest_life[13], res_forest_life[16], res_forest_life[19], res_forest_life[22], res_forest_life[25])
df_forest_NvL$Habitat <- "forested"
df_forest_NvL$comparison <- "nymph v. larva"
colnames(df_forest_NvL)[2] <- "LFC"
colnames(df_forest_NvL)[3] <- "SE"
colnames(df_forest_NvL)[4] <- "Wstat"
colnames(df_forest_NvL)[5] <- "p"
colnames(df_forest_NvL)[6] <- "p.adj"
colnames(df_forest_NvL)[7] <- "diff"
colnames(df_forest_NvL)[8] <- "pass_sens"
colnames(df_forest_NvL)[9] <- "diff_robust"

## Merge dataframes
df_forest_life_pair <- dplyr::bind_rows(df_forest_AvL, df_forest_AvN, df_forest_NvL)

#add another column for plot coloring for our non-sig taxa
df_forest_life_pair$Genus <- paste0(df_forest_life_pair$taxon)
df_forest_life_pair$Genus[df_forest_life_pair$p.adj >= 0.05] <- "Other"

#subset to keep diff abundant
#let's keep the ones that don't pass sensitivity first
df_forest_life_pair <- subset(df_forest_life_pair, !(diff == "TRUE" & pass_sens == "FALSE"))

#put genera in order (need to put "other" last so it's greyed out)
df_forest_life_pair$Genus <- factor(df_forest_life_pair$Genus, levels = c("Papiliotrema", "Unclassified Amphisphaeriales", "Vishniacozyma", "Other"))
colors <- c("#0000FF", "#FF0000", "#00FF00", "lightgrey")


#and make our plot
vol_plot_forest_life_pair <- df_forest_life_pair %>%
  ggplot(aes(x = LFC,
             y = -log10(p.adj),
             color = Genus)) +
  geom_point(size = 3.5, alpha = 0.8)+
  scale_color_manual(values = colors)+
  facet_grid(~comparison)

vol_plot_forest_life_pair <- vol_plot_forest_life_pair +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")

vol_plot_forest_life_pair + theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14),
        title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p-value") + labs(color = 'Genera') + theme(strip.text = element_text(face = "bold", size = 12))

##only sig changes in adults v. nymphs

##okay let's do season forested


