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

#set seed for reproducibility
set.seed(123)

## load phyloseq object
ps.rare <- readRDS(here::here("output/ps.rare.rds"))

#visualize alpha diversity just to make sure everything looks normal
plot_richness(ps.rare, x="Treatment", measures=c("Shannon", "Chao1"), color="Lifestage")
#+facet_grid(~Season)

####let's subset to just look at ticks (not soil)
ticks <- subset_samples(ps.rare, Type=="BLT")
plot_richness(ticks, x="Lifestage", measures=c("Shannon"), color="Treatment") +
  facet_grid(~Habitat+Season)

#calculate shannon diversity (richness and evenness)
alphadiv <- estimate_richness(ticks, measures = c("Shannon")) %>%
  rownames_to_column(var = "Unique_Specimen_ID") %>%
  left_join(as.data.frame(sample_data(ticks)), by = "Unique_Specimen_ID")

#Look at distribution and test for normality
hist(alphadiv$Shannon) #normal-ish
shapiro.test(alphadiv$Shannon) #is normal W = 0.98223, p = 0.07, will move forward with linear model


# linear model shannon diversity
mod.shannon <- glm(Shannon ~ Habitat + Treatment*Lifestage + Season, data=alphadiv)
shapiro.test(residuals(mod.shannon)) # residuals normal

glm_shan <- broom::tidy(anova(mod.shannon)) # Habitat and Treatment significant 
glm_shan$Metric <- "Shannon"

#check out means across lifestages
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

sig.df$y.position[1] <- 3.7
sig.df$y.position[2] <- 3.9
sig.df$y.position[3] <- 3.5
sig.df$y.position[4] <- 3.4


#the x.positions are wrong? so need to manually input
sig.df$xmin[1] <- 2.7
sig.df$xmax[1] <- 3.1 
sig.df$xmin[2] <- 2.7
sig.df$xmax[2] <- 3.3
sig.df$xmin[3] <- 2.7
sig.df$xmax[3] <- 2.9
sig.df$xmin[4] <- 1.8
sig.df$xmax[4] <- 2.2



as.factor(alphadiv$Lifestage) -> alphadiv$Lifestage
as.factor(alphadiv$Treatment) -> alphadiv$Treatment
as.factor(sig.df$group1) -> sig.df$group1
as.factor(sig.df$group2) -> sig.df$group2


plot.alpha <-  ggplot(alphadiv, aes(x = Lifestage, y = Shannon, color=Treatment)) + 
  geom_boxplot(lwd = 1.1, outlier.colour = "NA") + 
  scale_color_manual(values = friendly_pal("nickel_five"))+
  stat_pvalue_manual(sig.df, label = "p.adj.signif", hide.ns = TRUE, inherit.aes = TRUE, size = 6)

plot.alpha <- plot.alpha + 
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
#+ stat_pvalue_manual(sig.df, label = "p.adj.signif", hide.ns = TRUE, size = 6)
plot.alpha

