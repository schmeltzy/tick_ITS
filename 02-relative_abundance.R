##Overview: using rarefied phyloseq objects to look at relative abundance measures
## relative abundance and visualization
#load necessary packages
library(ggplot2)
library(phyloseq)
library(here)
library(fantaxtic)
#library(speedyseq)
#library(microViz)
#library(patchwork)
library(ggpubr)
library(dplyr)
library(rstatix)
library(broom)
#library(effsize)

#set seed for reproducibility
set.seed(123)


## load phyloseq object
ps.rare <- readRDS(here::here("output/ps.rare.rds"))


####let's subset to just look at ticks (not soil)

ticks <- subset_samples(ps.rare, Type=="BLT")


## Transform to relative abundance ====
ps.trans <- transform_sample_counts(ticks, function(OTU) (OTU/sum(OTU)*100))

#rename the NAs to higher taxa level
# Fill in names for NA taxa, not including their rank
ps.trans <- name_na_taxa(ps.trans, include_rank = FALSE, na_label = "Unclassified <tax>")
tax_table(ps.trans)

#make a quick lil plot of things glommed by Phylum

#sample colors
sample_colors <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","darkblue","#BAB0AC")
tol12 <- c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", 
           "#DDCC77", "#661100", "#CC6677", "#AA4499", "#882255", "#CCBB44","black") 
tol.extras <- c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", 
                "#DDCC77", "#661100", "#CC6677", "#AA4499", "#882255", "#CCBB44","black","#F28E2B","#E15759")

ps.phy <- tax_glom(ps.trans, taxrank = "Phylum", NArm = FALSE)
plot_bar(ps.phy, fill = "Phylum") +
  scale_fill_manual(values = sample_colors) 

## Prune to get most abundant taxa ====
top100 <- names(sort(taxa_sums(ps.trans), decreasing = TRUE))[1:100] 
ps.top100 <- prune_taxa(top100, ps.trans)

ps.top100@sam_data$Site <- factor(ps.top100@sam_data$Site)
ps.top100@sam_data$Habitat <- factor(ps.top100@sam_data$Habitat)

#trim samples to just get the ones with taxa remaining
ps.top100 <- prune_samples(sample_sums(ps.top100@otu_table) > 0, ps.top100)

## Plot top100
plot.relAbund <- plot_bar(ps.top100, x="Sample", fill="Class") + 
  #facet_grid(vars(Habitat), scales = "free_x") + 
  scale_fill_manual(values = sample_colors) 
plot.relAbund <- plot.relAbund + theme_bw(base_line_size = 1, base_rect_size = 1) + 
  theme(axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + xlab("Sample") + ylab("Relative Abundance (%)") +
  labs(title = "Tick-associated mycobiome") + theme(strip.text = element_text(face = "bold", size = 11)) + labs(fill = "Class") 
plot.relAbund

## okay what if we just do top 30 taxa, do we see Hypocreales?
## Prune most abundant taxa ====
top30 <- names(sort(taxa_sums(ps.trans), decreasing = TRUE))[1:30] 
ps.top30 <- prune_taxa(top30, ps.trans)

## Plot relative abundance ====
ps.top30 <- tax_glom(ps.top30, taxrank = "Family", NArm = FALSE)

plot.relAbund <- plot_bar(ps.top30, x="Sample", fill="Family") +
  #+facet_grid(vars(Lifestage), scales = "free_x") + 
  scale_fill_manual(values = tol12) 
plot.relAbund <- plot.relAbund + theme_bw(base_line_size = 1, base_rect_size = 1.5) + 
  theme(axis.text.x = element_text(size = 0, angle = 90), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + xlab("Sample") + ylab("Relative Abundance (%)") +
  labs(title = "Top 30 most abundant tick-associated fungi") + theme(strip.text = element_text(face = "bold", size = 11)) + labs(fill = "Family") 
plot.relAbund

## Yes they're in there 

## Let's try to just look at some of the entomopathogens
## Subset into the 4 we know about (mostly)
ps.gen <- tax_glom(ps.trans, taxrank = "Genus", NArm = FALSE)
ps.pathogens <- subset_taxa(ps.gen, Genus == "Metarhizium" | 
              Genus == "Beauveria" | 
              Genus == "Lecanicillium" |
              Genus == "Cordyceps" )

#trim samples to just get the ones with entos
path.pruned <- prune_samples(sample_sums(ps.pathogens@otu_table) > 0, ps.pathogens)

plot.entos <- plot_bar(path.pruned, x="Lifestage", fill = "Genus") +
  facet_wrap(~Treatment) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = sample_colors)
plot.entos <- plot.entos + theme_bw(base_line_size = 1, base_rect_size = 1) + 
  theme(axis.text.x = element_text(face = "bold", size = 11.5, angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + xlab("") + ylab("Relative Abundance (%)") +
  labs(title = "Relative abundance of tick-associated entomopathogens in forested habitat") + theme(strip.text = element_text(face = "bold", size = 11)) + labs(fill = "Genus") 
plot.entos
ggplot2::ggsave(here::here("output/plot.entos.png"), plot.entos,
                height = 447, width = 520, units = "mm",
                scale = 0.5, dpi = 1000)

#we don't have a lot of replication but we can see if we can do stats
melt <- psmelt(ps.pathogens)
melt %>%
  select(Sample, Abundance, Treatment, Genus, Season, Sex, Habitat, Lifestage) -> melt

## let's compare relative abundance of entomopathogens by Lifestage, Habitat, Treatment, Sex, Season

wil_Season <- broom::tidy(wilcox.test(Abundance ~ Season, data = melt)) #sig
wil_Season$var <- "season"

kw_Treatment <- broom::tidy(kruskal.test(Abundance ~ Treatment, data = melt)) #overall NS but p = 0.06
kw_Treatment$var <- "treatment"

wil_Habitat <- broom::tidy(wilcox.test(Abundance ~ Habitat, data = melt)) #sig
wil_Habitat$var <- "Habitat"

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

#save output
rbind(wil_Season, wil_Habitat) -> wil_output
write.csv(wil_output, here::here("output/wil_output.csv"))

rbind(kw_Lifestage, kw_Sex, kw_Treatment) -> kw_output
write.csv(kw_output, here::here("output/kw_output.csv"))

rbind(dunn_Lifestage) -> dunn_output
write.csv(dunn_output, here::here("output/dunn_output.csv"))
