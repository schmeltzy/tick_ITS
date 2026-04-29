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
tse_habitat_output <- ancombc2(data = tse, assay_name = "counts", tax_level = "Genus",
                               fix_formula = "Habitat", 
                           p_adj_method = "fdr", group = "Habitat", pseudo_sens = TRUE, 
                           alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 

res <- tse_habitat_output$res

#make another df of the output we want (non-intercept columns)
df_habitat <- data.frame(c(res[1], res[3], res[5], res[7]), res[9], res[11], res[13], res[15], res[17])
df_habitat$comparison <- "forested v. open"
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
  scale_color_manual(values = colors)+
  facet_grid(~comparison)

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

#should probably save the figure

#### then do for treatments
tse_treatment_output <- ancombc2(data = tse, assay_name = "counts", tax_level = "Genus",
                               fix_formula = "Treatment", 
                               p_adj_method = "fdr", group = "Treatment", pseudo_sens = TRUE, 
                               alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = TRUE,
                               mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100)) 

res.pair <- tse_treatment_output$res_pair

#make another df of the output we want (non-intercept columns)
df_BvMowed <- data.frame(c(res.pair[1], res.pair[2], res.pair[8], res.pair[14]), res.pair[20], res.pair[26], res.pair[32], res.pair[38], res.pair[44])
df_BvMowed$comparison <- "burned v. mowed"
colnames(df_BvMowed)[2] <- "LFC"
colnames(df_BvMowed)[3] <- "SE"
colnames(df_BvMowed)[4] <- "Wstat"
colnames(df_BvMowed)[5] <- "p"
colnames(df_BvMowed)[6] <- "p.adj"
colnames(df_BvMowed)[7] <- "diff"
colnames(df_BvMowed)[8] <- "pass_sens"
colnames(df_BvMowed)[9] <- "diff_robust"

df_BvUnmanaged <- data.frame(c(res.pair[1], res.pair[3], res.pair[9], res.pair[15]), res.pair[21], res.pair[27], res.pair[33], res.pair[39], res.pair[45])
df_BvUnmanaged$comparison <- "burned v. unmanaged"
colnames(df_BvUnmanaged)[2] <- "LFC"
colnames(df_BvUnmanaged)[3] <- "SE"
colnames(df_BvUnmanaged)[4] <- "Wstat"
colnames(df_BvUnmanaged)[5] <- "p"
colnames(df_BvUnmanaged)[6] <- "p.adj"
colnames(df_BvUnmanaged)[7] <- "diff"
colnames(df_BvUnmanaged)[8] <- "pass_sens"
colnames(df_BvUnmanaged)[9] <- "diff_robust"

df_BvUnmowed <- data.frame(c(res.pair[1], res.pair[4], res.pair[10], res.pair[16]), res.pair[22], res.pair[28], res.pair[34], res.pair[40], res.pair[46])
df_BvUnmowed$comparison <- "burned v. unmowed"
colnames(df_BvUnmowed)[2] <- "LFC"
colnames(df_BvUnmowed)[3] <- "SE"
colnames(df_BvUnmowed)[4] <- "Wstat"
colnames(df_BvUnmowed)[5] <- "p"
colnames(df_BvUnmowed)[6] <- "p.adj"
colnames(df_BvUnmowed)[7] <- "diff"
colnames(df_BvUnmowed)[8] <- "pass_sens"
colnames(df_BvUnmowed)[9] <- "diff_robust"

df_UnmanagedvMowed <- data.frame(c(res.pair[1], res.pair[5], res.pair[11], res.pair[17]), res.pair[23], res.pair[29], res.pair[35], res.pair[41], res.pair[47])
df_UnmanagedvMowed$comparison <- "unmanaged v. mowed"
colnames(df_UnmanagedvMowed)[2] <- "LFC"
colnames(df_UnmanagedvMowed)[3] <- "SE"
colnames(df_UnmanagedvMowed)[4] <- "Wstat"
colnames(df_UnmanagedvMowed)[5] <- "p"
colnames(df_UnmanagedvMowed)[6] <- "p.adj"
colnames(df_UnmanagedvMowed)[7] <- "diff"
colnames(df_UnmanagedvMowed)[8] <- "pass_sens"
colnames(df_UnmanagedvMowed)[9] <- "diff_robust"

df_UnmowedvMowed <- data.frame(c(res.pair[1], res.pair[6], res.pair[12], res.pair[18]), res.pair[24], res.pair[30], res.pair[36], res.pair[42], res.pair[48])
df_UnmowedvMowed$comparison <- "unmowed v. mowed"
colnames(df_UnmowedvMowed)[2] <- "LFC"
colnames(df_UnmowedvMowed)[3] <- "SE"
colnames(df_UnmowedvMowed)[4] <- "Wstat"
colnames(df_UnmowedvMowed)[5] <- "p"
colnames(df_UnmowedvMowed)[6] <- "p.adj"
colnames(df_UnmowedvMowed)[7] <- "diff"
colnames(df_UnmowedvMowed)[8] <- "pass_sens"
colnames(df_UnmowedvMowed)[9] <- "diff_robust"

df_UnmowedvUnmanaged <- data.frame(c(res.pair[1], res.pair[7], res.pair[13], res.pair[19]), res.pair[25], res.pair[31], res.pair[37], res.pair[43], res.pair[49])
df_UnmowedvUnmanaged$comparison <- "unmowed v. unmanaged"
colnames(df_UnmowedvUnmanaged)[2] <- "LFC"
colnames(df_UnmowedvUnmanaged)[3] <- "SE"
colnames(df_UnmowedvUnmanaged)[4] <- "Wstat"
colnames(df_UnmowedvUnmanaged)[5] <- "p"
colnames(df_UnmowedvUnmanaged)[6] <- "p.adj"
colnames(df_UnmowedvUnmanaged)[7] <- "diff"
colnames(df_UnmowedvUnmanaged)[8] <- "pass_sens"
colnames(df_UnmowedvUnmanaged)[9] <- "diff_robust"

#merge dataframes
df_treatment_pair <- dplyr::bind_rows(df_BvMowed, df_BvUnmanaged, df_BvUnmowed, df_UnmanagedvMowed, df_UnmowedvMowed, df_UnmowedvUnmanaged)

#add another column for plot coloring for our non-sig taxa
df_treatment_pair$Genus <- paste0(df_treatment_pair$taxon)
df_treatment_pair$Genus[df_treatment_pair$p.adj >= 0.05] <- "Other"

#subset to keep diff abundant
#let's keep the ones that don't pass sensitivity first
df_treatment_pair <- subset(df_treatment_pair, !(diff == "TRUE" & pass_sens == "FALSE"))

#set some colors (colorblind friendly of course)
colors <- c("#00FFBE", "lightgrey")
#put genera in order (need to put "other" last so it's greyed out)
df_treatment_pair$Genus <- factor(df_treatment_pair$Genus, levels = c("Unclassified Mycosphaerellaceae", "Other"))
#put treatments in the right order so they don't look so disorganized
df_treatment_pair$comparison <- factor(df_treatment_pair$comparison, levels = c("unmowed v. mowed", "unmowed v. unmanaged", "unmanaged v. mowed", "burned v. unmowed", "burned v. mowed", "burned v. unmanaged"))


#and make our plot
vol_plot_treatment <- df_treatment_pair %>%
  ggplot(aes(x = LFC,
             y = -log10(p.adj),
             color = Genus)) + 
  geom_point(size = 3.5, alpha = 0.8)+
  facet_grid(~comparison, labeller = label_wrap_gen(width = 10))+
  scale_color_manual(values = colors)

vol_plot_treatment <- vol_plot_treatment + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")

vol_plot_treatment + theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p-value") + labs(color = 'Taxa') +
  theme(strip.text = element_text(face = "bold", size = 12)) +
  theme(legend.position = "top", legend.text = element_text(face = "bold", size = 12))

##only diff in Mycosphaerellaceae in sugarbush treatments vs. burned
#save fig



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

# merge dataframes
df_forest_life_pair <- dplyr::bind_rows(df_forest_AvL, df_forest_AvN, df_forest_NvL)

#add another column for plot coloring for our non-sig taxa
df_forest_life_pair$Genus <- paste0(df_forest_life_pair$taxon)
df_forest_life_pair$Genus[df_forest_life_pair$p.adj >= 0.05] <- "Other"

#subset to keep diff abundant
#let's keep the ones that don't pass sensitivity first
df_forest_life_pair <- subset(df_forest_life_pair, !(diff == "TRUE" & pass_sens == "FALSE"))

#put genera in order (need to put "other" last so it's greyed out)
df_forest_life_pair$Genus <- factor(df_forest_life_pair$Genus, levels = c("Papiliotrema", "Unclassified Amphisphaeriales", "Vishniacozyma", "Other"))
colors <- c("#DC267F", "#FFD300", "#1F9698", "lightgrey")


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
tse_forest_season_output <- ancombc2(data = tse_forest, assay_name = "counts", tax_level = "Genus",
                                   fix_formula = "Season", 
                                   p_adj_method = "fdr", group = "Season", pseudo_sens = TRUE, 
                                   alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = FALSE,
                                   ) 

res_forest_season <- tse_forest_season_output$res

#make another df of the output we want (non-intercept columns with intercept being adult in this case)
df_forest_season <- data.frame(c(res_forest_season[1], res_forest_season[3], res_forest_season[5], res_forest_season[7]), res_forest_season[9], res_forest_season[11], res_forest_season[13], res_forest_season[15], res_forest_season[17])
df_forest_season$Habitat <- "forested"
df_forest_season$comparison <- "fall v. spring"
colnames(df_forest_season)[2] <- "LFC"
colnames(df_forest_season)[3] <- "SE"
colnames(df_forest_season)[4] <- "Wstat"
colnames(df_forest_season)[5] <- "p"
colnames(df_forest_season)[6] <- "p.adj"
colnames(df_forest_season)[7] <- "diff"
colnames(df_forest_season)[8] <- "pass_sens"
colnames(df_forest_season)[9] <- "diff_robust"

#add another column for plot coloring for our non-sig taxa
df_forest_season$Genus <- paste0(df_forest_season$taxon)
df_forest_season$Genus[df_forest_season$p.adj >= 0.05] <- "Other"

#subset to keep diff abundant
#let's keep the ones that don't pass sensitivity first
df_forest_season <- subset(df_forest_season, !(diff == "TRUE" & pass_sens == "FALSE"))

#put genera in order (need to put "other" last so it's greyed out)
df_forest_season$Genus <- factor(df_forest_season$Genus, levels = c("Ramularia", "Other"))
colors <- c("#FF00B6", "lightgrey")


#and make our plot
vol_plot_forest_season <- df_forest_season %>%
  ggplot(aes(x = LFC,
             y = -log10(p.adj),
             color = Genus)) +
  geom_point(size = 3.5, alpha = 0.8)+
  scale_color_manual(values = colors)+
  facet_grid(~comparison)

vol_plot_forest_season <- vol_plot_forest_season +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")

vol_plot_forest_season + theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14),
        title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p-value") + labs(color = 'Genera') + theme(strip.text = element_text(face = "bold", size = 12))

##only Ramularia




####now do open habitats
open <- subset_samples(ticks, Habitat=="open")
#convert phylo object to tree summarized experiment
tse_open <- mia::convertFromPhyloseq(open)

tse_open_output <- ancombc2(data = tse_open, assay_name = "counts", tax_level = "Genus",
                              fix_formula = "Treatment", 
                              p_adj_method = "fdr", group = "Treatment", pseudo_sens = TRUE, 
                              alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = FALSE,
                              mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100)) 

res_open <- tse_open_output$res

#make another df of the output we want (non-intercept columns)
df_open_treat <- data.frame(c(res_open[1], res_open[3], res_open[5], res_open[7]), res_open[9], res_open[11], res_open[13], res_open[15], res_open[17])
df_open_treat$comparison <- "Treatment"
colnames(df_open_treat)[2] <- "LFC"
colnames(df_open_treat)[3] <- "SE"
colnames(df_open_treat)[4] <- "Wstat"
colnames(df_open_treat)[5] <- "p"
colnames(df_open_treat)[6] <- "p.adj"
colnames(df_open_treat)[7] <- "diff"
colnames(df_open_treat)[8] <- "pass_sens"
colnames(df_open_treat)[9] <- "diff_robust"

#add another column for plot coloring for our non-sig taxa
df_open_treat$Genus <- paste0(df_open_treat$taxon)
df_open_treat$Genus[df_open_treat$p.adj >= 0.05] <- "Other"

#subset to keep diff abundant
#let's keep the ones that don't pass sensitivity first
df_open_treat <- subset(df_open_treat, !(diff == "TRUE" & pass_sens == "FALSE"))

#set some colors (colorblind friendly of course)
#colors <- c(friendly_pal("glasbey_twelve"), "lightgrey")
#put genera in order (need to put "other" last so it's greyed out)
df_open_treat$Genus <- factor(df_open_treat$Genus, levels = c("Alternaria", "Pilidium", "Piskurozyma", "Plectosphaerella", "Ramularia", "Thyridium", "Unclassified Amphisphaeriales", "Unclassified Basidiomycota", "Unclassified Leotiomycetes", "Unclassified Mycosphaerellaceae", "Unclassified Tremellales", "Vishniacozyma", "Other"))

#and make our plot
vol_plot_open_treat <- df_open_treat %>%
  ggplot(aes(x = LFC,
             y = -log10(p.adj),
             color = Genus)) + 
  geom_point(size = 3.5, alpha = 0.8)+
  scale_color_manual(values = "lightgrey")

vol_plot_open_treat <- vol_plot_open_treat + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")

vol_plot_open_treat + theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p-value") + labs(color = 'Genera') + theme(strip.text = element_text(face = "bold", size = 12)) 

