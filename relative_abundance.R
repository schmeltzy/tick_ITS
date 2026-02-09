## relative abundance and figures

## Set seed ====
set.seed(123)

## Load libraries ====
library(ggplot2)
library(phyloseq)
#library(speedyseq)
library(microViz)
library(ggpubr)
library(dplyr)
library(rstatix)
library(broom)
#library(effsize)

## load phyloseq object
ps <- readRDS("phyloseq_obj.rds")

## validate phyloseq object and fix if you need to
ps <- phyloseq_validate(ps)

ps <- tax_fix(
  ps,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
)

####let's subset to just look at ticks (not soil)

ticks <- subset_samples(ps, Species=="BLT")


## Transform to relative abundance ====
ps.trans <- transform_sample_counts(ticks, function(OTU) OTU/sum(OTU))

#make a quick lil plot of things glommed by Phylum
ps.phy <- tax_glom(ps.trans, taxrank = "Phylum", NArm = FALSE)
plot_bar(ps.phy, fill = "Phylum")

## Prune most abundant taxa ====
top100 <- names(sort(taxa_sums(ps.trans), decreasing = TRUE))[1:100] 
ps.top100 <- prune_taxa(top100, ps.trans)

## Plot relative abundance ====
ps.top100 <- tax_glom(ps.top100, taxrank = "Class", NArm = FALSE)


## make tickID a factor
ps.top100@sam_data$Specimen_ID <- as.factor(ps.top100@sam_data$Specimen_ID)

#sample colors
sample_colors <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
                    "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                    "#920000","#924900","#db6d00","#24ff24","#ffff6d","#E69F00")


ps.top100@sam_data$Site <- factor(ps.top100@sam_data$Site)
ps.top100@sam_data$Habitat <- factor(ps.top100@sam_data$Habitat)

## Plot by treatment and season
plot.relAbund <- plot_bar(ps.top100, x="Specimen_ID", fill="Class") + 
  facet_grid(vars(Habitat), scales = "free_x") + scale_fill_manual(values = sample_colors) 
plot.relAbund <- plot.relAbund + theme_bw(base_line_size = 1, base_rect_size = 1.5) + 
  theme(axis.text.x = element_text(face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + xlab("SampleID") + ylab("Relative Abundance (%)") +
  labs(title = "Tick-associated mycobiome") + theme(strip.text = element_text(face = "bold", size = 11)) + labs(fill = "Class") 
plot.relAbund

## okay what if we just do top 50 taxa, do we see our Hypocreales friends?
## Prune most abundant taxa ====
top50 <- names(sort(taxa_sums(ps.trans), decreasing = TRUE))[1:50] 
ps.top50 <- prune_taxa(top50, ps.trans)

## Plot relative abundance ====
ps.top50 <- tax_glom(ps.top50, taxrank = "Order", NArm = FALSE)

plot.relAbund <- plot_bar(ps.top50, x="Specimen_ID", fill="Order") 
  #+facet_grid(vars(Lifestage), scales = "free_x") + scale_fill_manual(values = sample_colors) 
plot.relAbund <- plot.relAbund + theme_bw(base_line_size = 1, base_rect_size = 1.5) + 
  theme(axis.text.x = element_text(face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + xlab("") + ylab("Relative Abundance (%)") +
  labs(title = "Tick-associated mycobiome") + theme(strip.text = element_text(face = "bold", size = 11)) + labs(fill = "Order") 
plot.relAbund

## Yay they're in there 

## Let's try to just look at some of the entomopathogens
## Subset into the 4 we care about (mostly)
ps.pathogens <- subset_taxa(ps.trans, Genus == "Metarhizium" | 
              Genus == "Beauveria" | 
              Genus == "Lecanicillium" |
              Genus == "Cordyceps" ) 

plot_bar(ps.pathogens, x="Genus", fill = "Species") + facet_grid(~Treatment)

#GP.chl <- subset_taxa(ps, Genus=="Metarhizium")
#plot_tree(GP.chl, color="Treatment", shape="Genus", label.tips="Genus", size="abundance", plot.margin=0.6)
