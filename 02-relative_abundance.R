##Overview: using rarefied phyloseq objects to look at relative abundance measures
## relative abundance and visualization
#load necessary packages
library(ggplot2)
library(phyloseq)
library(here)
library(fantaxtic)
library(ggpubr)
library(ggpubfigs)
library(dplyr)
library(rstatix)
library(broom)

#set seed for reproducibility
set.seed(123)


## load phyloseq object
ps.rare <- readRDS(here::here("output/ps.rare.rds"))


####let's subset to just look at ticks (not soil)

ticks <- subset_samples(ps.rare, Type=="BLT")


## Transform to relative abundance ====
ps.trans <- transform_sample_counts(ticks, function(OTU) (OTU/sum(OTU)*100))
#check sample sums to 100
sample_sums(ps.trans)

#rename the NAs to higher taxa level
# Fill in names for NA taxa, not including their rank
ps.trans <- name_na_taxa(ps.trans, include_rank = FALSE, na_label = "Unclassified <tax>")
tax_table(ps.trans)

#make a quick lil plot of things glommed by Phylum

#sample colors
sample_colors <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","darkblue","#BAB0AC")
tol12 <- c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", 
           "#DDCC77", "#661100", "#CC6677", "#AA4499", "#882255", "#CCBB44","black") 

ps.phy <- tax_glom(ps.trans, taxrank = "Phylum", NArm = FALSE)
plot_bar(ps.phy, fill = "Phylum") +
  scale_fill_manual(values = sample_colors) 

##prune to get most abundant taxa 
top30 <- names(sort(taxa_sums(ps.trans), decreasing = TRUE))[1:30] 
ps.top30 <- prune_taxa(top30, ps.trans)

## Plot relative abundance ====
ps.top30 <- tax_glom(ps.top30, taxrank = "Family", NArm = FALSE)

plot.relAbund <- plot_bar(ps.top30, x="Sample", fill="Family") +
  #facet_grid(vars(Habitat), scales = "free_x") + 
  scale_fill_manual(values = tol12) 
plot.relAbund <- plot.relAbund + theme_bw(base_line_size = 1, base_rect_size = 1.5) + 
  theme(axis.text.x = element_text(size = 0, angle = 90), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + xlab("Sample") + ylab("Relative Abundance (%)") +
  labs(title = "Top 30 most abundant tick-associated fungi") + theme(strip.text = element_text(face = "bold", size = 11)) +
  labs(fill = "Family") 
plot.relAbund

#we see some Hypocreales in there

#melt at genus level
ps_trans_genus_melt <- ps.trans %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt()

#get genera with mean relative abundance >1% across all samples
genus_sum <- ps_trans_genus_melt %>% group_by(Genus) %>% dplyr::summarise(Average = mean(Abundance), sd = sd(Abundance), min = min(Abundance), max = max(Abundance))
genus_sub <- genus_sum[which(genus_sum$Average > 1),]
gen_names <- genus_sub$Genus
gen_names

write.csv(genus_sub, here::here("output/top_relabund_genera.csv"))

#replace genus with <1% abundance with NAs
ps_trans_genus_melt$genus <- ps_trans_genus_melt$Genus

ps_trans_genus_melt$genus[ps_trans_genus_melt$genus != "Alternaria" & 
                          ps_trans_genus_melt$genus != "Mycosphaerella" &
                          ps_trans_genus_melt$genus != "Paraphaeosphaeria" &
                          ps_trans_genus_melt$genus != "Ramularia" &
                          ps_trans_genus_melt$genus != "Unclassified Basidiomycota" &
                          ps_trans_genus_melt$genus != "Unclassified Fungi" &
                          ps_trans_genus_melt$genus != "Unclassified Mycosphaerellaceae" &
                          ps_trans_genus_melt$genus != "Cladosporium" &
                          ps_trans_genus_melt$genus != "Papiliotrema" & 
                          ps_trans_genus_melt$genus != "Plectosphaerella" &
                          ps_trans_genus_melt$genus != "Unclassified Amphisphaeriales" &
                          ps_trans_genus_melt$genus != "Unclassified Didymellaceae" &
                          ps_trans_genus_melt$genus != "Unclassified Leotiomycetes" & 
                          ps_trans_genus_melt$genus != "Vishniacozyma"  ] <- NA

#make our plot
colors <- c(friendly_pal("glasbey_twelve"),"lightgrey","maroon")
#reorder colors to match genera of other plots for consistency
colors <- c("#9A4D42","#FF0000","#000033","#00FF00","maroon","#009FFF","#FFD300","#005300","#0000FF","orange","#00FFBE","#FF00B6","#1F9698","lightgrey")

plot_genus = ggplot(ps_trans_genus_melt, aes(x = Treatment, y=Abundance)) + 
  geom_bar(stat="identity", position="fill", aes(fill = reorder(genus, Abundance))) +
  scale_fill_manual(values= colors, 
                   na.value = "transparent")  +
  facet_grid(~Habitat, scales = "free_x", space = "free_x") +
  theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  ylab("Relative Abundance (%)") +
  xlab("Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  theme(strip.text = element_text(face = "bold", size = 11)) +
  theme(legend.position = "right", legend.text = element_text(face = "bold", size = 10))+
  labs(fill = "Genus")
plot_genus

#save fig
ggplot2::ggsave(here::here("output/plot_rel_abund.png"), plot_genus,
                height = 500, width = 500, units = "mm",
                scale = 0.5, dpi = 1000)


#================

#what if we do the same thing but separated by habitat and look at top 30 most abundant
forest <- subset_samples(ps.trans, Habitat=="forested")
open <- subset_samples(ps.trans, Habitat=="open")


###maybe 1% abundance is better
#melt at genus level
ps_trans_forest_melt <- forest %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt()

#get genera with mean relative abundance >1% across all samples
forest_sum <- ps_trans_forest_melt %>% group_by(Genus) %>% dplyr::summarise(Average = mean(Abundance), sd = sd(Abundance), min = min(Abundance), max = max(Abundance))
forest_sub <- forest_sum[which(forest_sum$Average > 1),]
forest_names <- forest_sub$Genus
forest_names

write.csv(forest_sub, here::here("output/top_forest_relabund_genera.csv"))

#replace genus with <1% abundance with NAs
ps_trans_forest_melt$genus <- ps_trans_forest_melt$Genus

ps_trans_forest_melt$genus[ps_trans_forest_melt$genus != "Cladosporium" & 
                            ps_trans_forest_melt$genus != "Fusarium" &
                            ps_trans_forest_melt$genus != "Papiliotrema" &
                            ps_trans_forest_melt$genus != "Plectosphaerella" &
                            ps_trans_forest_melt$genus != "Unclassified Basidiomycota" &
                            ps_trans_forest_melt$genus != "Unclassified Fungi" &
                            ps_trans_forest_melt$genus != "Unclassified Mycosphaerellaceae" &
                            ps_trans_forest_melt$genus != "Curvibasidium" &
                            ps_trans_forest_melt$genus != "Mycosphaerella" &
                            ps_trans_forest_melt$genus != "Ramularia" &
                            ps_trans_forest_melt$genus != "Unclassified Helotiales" &
                            ps_trans_forest_melt$genus != "Pilidium" &
                            ps_trans_forest_melt$genus != "Unclassified Amphisphaeriales" &
                            ps_trans_forest_melt$genus != "Unclassified Didymellaceae" &
                            ps_trans_forest_melt$genus != "Unclassified Leotiomycetes" & 
                            ps_trans_forest_melt$genus != "Vishniacozyma"  ] <- NA

#make our plot
#reorder colors to match genera of other plots for consistency
colors <- c("#FF00B6","lightgrey","#00FFBE","orange","#1F9698","#FFD300","#009FFF","#00FF00","#005300","#000033","#9A4D42","maroon","#FF0000","#FF9DA7","#440154","#0000FF")
#ugh these are reversed!
colors <- rev(c("#FF00B6","lightgrey","#00FFBE","orange","#1F9698","#FFD300","#009FFF","#00FF00","#005300","#000033","#9A4D42","maroon","#FF0000","#FF9DA7","#440154","#0000FF"))
#also reorder lifestage levels
ps_trans_forest_melt$Lifestage <- factor(ps_trans_forest_melt$Lifestage, levels = c("larva", "nymph","adult"))
#and season levels so they group together nicely
ps_trans_forest_melt$Season <- factor(ps_trans_forest_melt$Season, levels = c("spring", "fall"))

plot_forest = ggplot(ps_trans_forest_melt, aes(x = Lifestage, y=Abundance)) + 
  geom_bar(stat="identity", position="fill", aes(fill = reorder(genus, Abundance))) +
  scale_fill_manual(values= colors, 
                    na.value = "transparent")  +
  facet_grid(Treatment~Season, scales = "free_x", space = "free_x") +
  theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  ylab("Relative Abundance (%)") +
  #xlab("Life Stage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.title.x = element_blank()) +
  theme(strip.text = element_text(face = "bold", size = 11)) +
  theme(legend.position = "right", legend.text = element_text(face = "bold", size = 10))+
  labs(fill = "Genus")
plot_forest

#save fig
ggplot2::ggsave(here::here("output/plot_forest_relabund.png"), plot_forest,
                height = 500, width = 500, units = "mm",
                scale = 0.5, dpi = 1000)

#==============
#now with open data
#melt at genus level
ps_trans_open_melt <- open %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt()

#get genera with mean relative abundance >1% across all samples
open_sum <- ps_trans_open_melt %>% group_by(Genus) %>% dplyr::summarise(Average = mean(Abundance), sd = sd(Abundance), min = min(Abundance), max = max(Abundance))
open_sub <- open_sum[which(open_sum$Average > 1),]
open_names <- open_sub$Genus
open_names

write.csv(open_sub, here::here("output/top_open_relabund_genera.csv"))

#replace genus with <1% abundance with NAs
ps_trans_open_melt$genus <- ps_trans_open_melt$Genus

ps_trans_open_melt$genus[ps_trans_open_melt$genus != "Cladosporium" & 
                             ps_trans_open_melt$genus != "Alternaria" &
                             ps_trans_open_melt$genus != "Papiliotrema" &
                             ps_trans_open_melt$genus != "Dioszegia" &
                             ps_trans_open_melt$genus != "Hannaella" &
                             ps_trans_open_melt$genus != "Unclassified Fungi" &
                             ps_trans_open_melt$genus != "Mrakia" &
                             ps_trans_open_melt$genus != "Paraphaeosphaeria" &
                             ps_trans_open_melt$genus != "Parastagonospora" &
                             ps_trans_open_melt$genus != "Septoria" &
                             ps_trans_open_melt$genus != "Stagonosporopsis" &
                             ps_trans_open_melt$genus != "Tricellula" &
                             ps_trans_open_melt$genus != "Unclassified Didymellaceae" &
                             ps_trans_open_melt$genus != "Vishniacozyma"  ] <- NA

#make our plot
#reorder colors to match genera of other plots for consistency
colors <- c("lightgrey","#1F9698","#0000FF","orange","#FF0000","#005300","#FF9DA7","maroon","#440154","#9A4D42","#00FF00","#000033","#009FFF","#FFD300")
#reverse the order
colors <- rev(c("lightgrey","#1F9698","#0000FF","orange","#FF0000","#005300","#FF9DA7","maroon","#440154","#9A4D42","#00FF00","#000033","#009FFF","#FFD300"))

ps_trans_open_melt$Season <- factor(ps_trans_open_melt$Season, levels = c("spring", "fall"))

plot_open = ggplot(ps_trans_open_melt, aes(x = Lifestage, y=Abundance)) + 
  geom_bar(stat="identity", position="fill", aes(fill = reorder(genus, Abundance))) +
  scale_fill_manual(values= colors, 
                    na.value = "transparent")  +
  facet_grid(~Treatment+Season, scales = "free_x", space = "free_x") +
  theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  ylab("Relative Abundance (%)") +
  xlab("Season") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  theme(strip.text = element_text(face = "bold", size = 11)) +
  theme(legend.position = "right", legend.text = element_text(face = "bold", size = 10))+
  labs(fill = "Genus")
plot_open

#save fig
ggplot2::ggsave(here::here("output/plot_open_relabund.png"), plot_open,
                height = 500, width = 500, units = "mm",
                scale = 0.5, dpi = 1000)


#===============================

## Let's try to just look at some of the entomopathogens
## Subset into the 4 we know about (mostly)
ps.pathogens <- subset_taxa(ps.trans, Genus == "Metarhizium" | 
                              Genus == "Beauveria" | 
                              Genus == "Lecanicillium" |
                              Genus == "Cordyceps" |
                              Genus == "Paecilomyces" )

#trim samples to just get the ones with entos
path.pruned <- prune_samples(sample_sums(ps.pathogens@otu_table) > 0, ps.pathogens)

#put treatments in order for the plot
path.pruned@sam_data$Treatment <- factor(path.pruned@sam_data$Treatment, levels = c("unmowed", "mowed","unmanaged"))
#add column for taxon ID with genus+species info
tax_df <- as.data.frame(tax_table(path.pruned))

tax_df$Species <- paste(tax_df$Genus, tax_df$Species, sep = " ")

tax_table(path.pruned) <- tax_table(as.matrix(tax_df))

## entomopathogen relative abundance figure
plot.entos <- plot_bar(path.pruned, x="Lifestage", fill = "Species") +
  facet_wrap(~Treatment) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = sample_colors)
plot.entos <- plot.entos + theme_bw(base_line_size = 1, base_rect_size = 1) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(face = "bold", size = 11.5, angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + xlab("Life Stage") + ylab("Relative Abundance (%)") +
  theme(strip.text = element_text(face = "bold", size = 11)) +
  theme(legend.position = "right", legend.text = element_text(face = "bold", size = 10))
plot.entos
ggplot2::ggsave(here::here("output/plot.entos.png"), plot.entos,
                height = 400, width = 600, units = "mm",
                scale = 0.5, dpi = 1000)

#we don't have a lot of replication but we can see if we can do stats
melt <- psmelt(ps.pathogens)
melt %>%
  select(Sample, Abundance, Treatment, Genus, Season, Sex, Habitat, Lifestage) -> melt

#look at distrbution
hist(melt$Abundance) #definitely not normal

#take average relative abundance between lifestages in forested habitat
melt %>%
  filter(Habitat == "forested") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = Lifestage) %>%
  arrange()

## let's compare mean relative abundance of entomopathogens by Lifestage, Habitat, Treatment, Sex, Season

kw_Season <- broom::tidy(kruskal.test(Abundance ~ Season, data = melt)) #sig
kw_Season$var <- "season"

kw_Treatment <- broom::tidy(kruskal.test(Abundance ~ Treatment, data = melt)) #NS
kw_Treatment$var <- "treatment"

kw_Habitat <- broom::tidy(kruskal.test(Abundance ~ Habitat, data = melt)) #NS
kw_Habitat$var <- "Habitat"

kw_Sex <- broom::tidy(kruskal.test(Abundance ~ Sex, data = melt)) #NS
kw_Sex$var <- "sex"

kw_Lifestage <- broom::tidy(kruskal.test(Abundance ~ Lifestage, data = melt)) #sig
kw_Lifestage$var <- "lifestage"
dunn_Lifestage <- as.data.frame(dunn_test(Abundance ~ Lifestage, data = melt, p.adjust.method = "fdr")) # all vs. all sig
dunn_Lifestage$var <- "lifestage"


means.life <- melt %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = Lifestage) %>%
  arrange()
means.season <- melt %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = Season) %>%
  arrange()
means.hab <- melt %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = Habitat) %>%
  arrange()
means.treatment<- melt %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = Treatment) %>%
  arrange()
means.sex<- melt %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = Sex) %>%
  arrange()

#total relative abundance by lifestage by treatment
total.abund.life.mowed <- melt %>%
  filter(Treatment == "mowed") %>%
  summarise(total = sum(Abundance), .by = Lifestage, Treatment = "mowed") %>%
  arrange()

total.abund.life.unmowed <- melt %>%
  filter(Treatment == "unmowed") %>%
  summarise(total = sum(Abundance), .by = Lifestage, Treatment = "unmowed") %>%
  arrange()

total.abund.life.unmanaged <- melt %>%
  filter(Treatment == "unmanaged") %>%
  summarise(total = sum(Abundance), .by = Lifestage, Treatment = "unmanaged") %>%
  arrange()

#save output
rbind(kw_Season, kw_Habitat, kw_Lifestage, kw_Sex, kw_Treatment) -> kw_output
write.csv(kw_output, here::here("output/kw_output.csv"))

rbind(dunn_Lifestage) -> dunn_output
write.csv(dunn_output, here::here("output/dunn_output.csv"))

rbind(total.abund.life.mowed, total.abund.life.unmowed, total.abund.life.unmanaged) -> total.abund.life.output
write.csv(total.abund.life.output, here::here("output/total_abundance_lifestage_pathogens.csv"))

bind_rows(means.hab, means.life, means.treatment, means.season, means.sex) -> means.all.output
write.csv(means.all.output, here::here("output/means_all_pathogens.csv"))





