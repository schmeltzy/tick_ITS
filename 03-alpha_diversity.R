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
library(emmeans)
library(ggpubfigs)
library(svglite)

#set seed for reproducibility
set.seed(123)

## load phyloseq object
ps.rare <- readRDS(here::here("output/ps.rare.rds"))

#visualize alpha diversity just to make sure everything looks normal
plot_richness(ps.rare, x="Treatment", measures=c("Shannon","Observed"), color="Lifestage")
#+facet_grid(~Season)

####let's subset to just look at ticks (not soil)
ticks <- subset_samples(ps.rare, Type=="BLT")
plot_richness(ticks, x="Lifestage", measures=c("Shannon"), color="Treatment") +
  facet_grid(~Habitat+Season)

#calculate shannon diversity (richness and evenness)
alphadiv <- phyloseq::estimate_richness(ticks, measures = c("Shannon")) %>%
  rownames_to_column(var = "Unique_Specimen_ID") %>%
  left_join(as.data.frame(sample_data(ticks)), by = "Unique_Specimen_ID")

#Look at distribution and test for normality
hist(alphadiv$Shannon) #normal-ish
shapiro.test(alphadiv$Shannon) #is normal W = 0.98223, p = 0.07, will move forward with linear model


#=========================

# 1. Broad Landscape Scale 
global_model <- lm(Shannon ~ Habitat + Habitat:Treatment + Habitat:Season + Habitat:Treatment*Lifestage, data = alphadiv)
shapiro.test(residuals(global_model)) #residuals normal

global_lm_shan <- broom::tidy(anova(global_model)) # Habitat and Treatment and Lifestage interactions significant 
global_lm_shan$Metric <- "Shannon"

# 2. Micro Scale (To isolate effects within that single life stage)
sex_model <- lm(Shannon ~ Habitat + Habitat:Treatment + Season + Sex, 
                data = subset(alphadiv, Lifestage == "adult")) # Replace 'adult' if it's a different stage
shapiro.test(residuals(sex_model)) #residuals normal

sex_lm_shan <- broom::tidy(anova(sex_model)) # Sex NS
sex_lm_shan$Metric <- "Shannon"

#check out means across treatments and seasons
emmip(global_model, Habitat ~ Treatment | Lifestage*Season)

# pairwise comparisons
global_mod_shannon_pw <- emmeans(global_model, pairwise ~ Treatment, adjust = "tukey") 
global_mod_shannon_pw <- broom::tidy(global_mod_shannon_pw$contrasts)
global_mod_shannon_pw <- rstatix::add_significance(data = as.data.frame(global_mod_shannon_pw), p.col = "adj.p.value")
global_mod_shannon_pw$Metric <- "Shannon"

global_mod_shannon_pw2 <- emmeans(global_model, pairwise ~ Treatment | Lifestage, adjust = "tukey") 
global_mod_shannon_pw2 <- broom::tidy(global_mod_shannon_pw2$contrasts)
global_mod_shannon_pw2 <- rstatix::add_significance(data = as.data.frame(global_mod_shannon_pw2), p.col = "adj.p.value")
global_mod_shannon_pw2$Metric <- "Shannon"

global_mod_shannon_pw3 <- emmeans(global_model,pairwise ~ Treatment + Season | Lifestage, adjust = "tukey") 
global_mod_shannon_pw3 <- broom::tidy(global_mod_shannon_pw3$contrasts)
global_mod_shannon_pw3 <- rstatix::add_significance(data = as.data.frame(global_mod_shannon_pw3), p.col = "adj.p.value")
global_mod_shannon_pw3$Metric <- "Shannon"

#========================


# linear model shannon diversity
mod.shannon <- glm(Shannon ~ Habitat + Treatment*Lifestage + Season, data=alphadiv)
shapiro.test(residuals(mod.shannon)) # residuals normal

glm_shan <- broom::tidy(anova(mod.shannon)) # Habitat and Treatment significant 
glm_shan$Metric <- "Shannon"

#check out means across treatments
emmip(mod.shannon, Habitat ~ Treatment | Lifestage)

# pairwise comparisons
mod_shannon_pw <- emmeans(mod.shannon, pairwise ~ Treatment, adjust = "tukey") 
mod_shannon_pw <- broom::tidy(mod_shannon_pw$contrasts)
mod_shannon_pw <- rstatix::add_significance(data = as.data.frame(mod_shannon_pw), p.col = "adj.p.value")
mod_shannon_pw$Metric <- "Shannon"

mod_shannon_pw2 <- emmeans(mod.shannon, pairwise ~ Treatment | Lifestage, adjust = "tukey") 
mod_shannon_pw2 <- broom::tidy(mod_shannon_pw2$contrasts)
mod_shannon_pw2 <- rstatix::add_significance(data = as.data.frame(mod_shannon_pw2), p.col = "adj.p.value")
mod_shannon_pw2$Metric <- "Shannon"

# save
write.csv(glm_shan, here::here("output/glm_alpha_shannon.csv"))
write.csv(mod_shannon_pw, here::here("output/pairwise_alpha_treatment.csv"))
write.csv(mod_shannon_pw2, here::here("output/pairwise_alpha_lifestage.csv"))

##plot alphas
#make sure things are in order
alphadiv$Lifestage <- factor(alphadiv$Lifestage, c("larva", "nymph", "adult"), ordered = TRUE)
alphadiv$Treatment <- factor(alphadiv$Treatment, c("unmanaged", "burned", "mowed","unmowed"), ordered = TRUE)

#change the variable names
alphadiv$Treatment <- gsub("unmanaged", "unburned", alphadiv$Treatment)
alphadiv$Treatment <- gsub("mowed", "managed", alphadiv$Treatment)
alphadiv$Treatment <- gsub("unmowed", "unmanaged", alphadiv$Treatment)


#make sig df
sig.df <- mod_shannon_pw2

sig.df$group1 <- paste0(sig.df$contrast)
sig.df$group2 <- paste0(sig.df$contrast)

sig.df$group1 <- sapply(strsplit(basename(sig.df$group1), "-"), `[`,1)
sig.df$group1  <- gsub("forested", "", sig.df$group1) 
sig.df$group1  <- gsub("open", "", sig.df$group1) 

sig.df$group2 <- sapply(strsplit(basename(sig.df$group2), "-"), `[`,2)
sig.df$group2  <- gsub("forested", "", sig.df$group2) 
sig.df$group2  <- gsub("open", "", sig.df$group2) 


sig.df <- subset(sig.df, adj.p.value.signif != "ns" & adj.p.value != "NA")

#set y position of bars on plot
sig.df$y.position <- 3.9

#rename column
colnames(sig.df)[10] <- "p.adj.signif"

sig.df <- sig.df %>%
  add_xy_position(
    x = "Lifestage", 
    group = "Treatment",    # Crucial for dodging
    data = alphadiv, 
    formula = Shannon ~ Lifestage
  )

sig.df$y.position[1] <- 0
sig.df$y.position[2] <- 0
sig.df$y.position[3] <- 3.5
sig.df$y.position[4] <- 3.5


#add Habitat facet
sig.df$Habitat[1] <- "forested"
sig.df$Habitat[2] <- "forested"
sig.df$Habitat[3] <- "open"
sig.df$Habitat[4] <- "forested"

#the x.positions are wrong? so need to manually input
#sig.df$xmin[1] <- 0
#sig.df$xmax[1] <- 0 
#sig.df$xmin[2] <- 0
#sig.df$xmax[2] <- 0
sig.df$xmin[3] <- 0.8
sig.df$xmax[3] <- 1.2
sig.df$xmin[4] <- 1.8
sig.df$xmax[4] <- 2.2


as.factor(alphadiv$Lifestage) -> alphadiv$Lifestage
as.factor(alphadiv$Treatment) -> alphadiv$Treatment
as.factor(sig.df$group1) -> sig.df$group1
as.factor(sig.df$group2) -> sig.df$group2

#skip first two rows since we can't span facets
plot.sig.df <- tail(sig.df, -2)

#relevel factors
alphadiv$Treatment <- factor(alphadiv$Treatment, levels = c("unburned", "burned", "managed", "unmanaged"))


plot.alpha <-  ggplot(alphadiv, aes(x = Lifestage, y = Shannon, color=Treatment)) + 
  geom_boxplot(lwd = 1.1, outlier.colour = "NA") + 
  scale_color_manual(values = friendly_pal("ito_seven"))+
  stat_pvalue_manual(plot.sig.df, label = "p.adj.signif", hide.ns = TRUE, inherit.aes = FALSE, size = 6)+
  geom_point(aes(color=Treatment), size = 1.5, alpha = 0.5, 
           position = position_jitterdodge(jitter.width = 0.1)) + 
  theme_bw(base_line_size = 1.5, base_rect_size = 1.75) +
  theme(axis.text = element_text(face = "bold", size = 11.5),
        axis.title = element_text(face = "bold", size = 12),
        title = element_text(face = "bold"))
plot.alpha

plot.alpha <- plot.alpha + 
  facet_grid(~Habitat, scales = "free", space = "free")+
  theme(strip.text = element_text(size = 12, face = "bold"))+
  labs(color = "Treatment") + xlab("Lifestage")
plot.alpha

#save plot
ggplot2::ggsave(here::here("output/alphas.plot.png"), plot.alpha,
                height = 300, width = 500, units = "mm",
                scale = 0.5, dpi = 1000)

