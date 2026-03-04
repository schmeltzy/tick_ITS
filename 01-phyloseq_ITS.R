#Overview: bringing ITS ASV tables and taxonomy (and phylogenetic tree if you have it) into phyloseq

library(phyloseq)
library(ggplot2)
library(readxl)
library(tidyverse)
library(vegan)
library(FSA)


#call up taxa and asv table from dada2 output
taxa <- readRDS("taxa_both.rds")
seqtab.nochim <- readRDS("ASV_both.rds")
#fitGTR <- readRDS("fitGTR.rds")

#import sample data
mapfile <- read_excel("Tick_amps_metadata_for_R.xlsx")
samdf <- as.data.frame(mapfile)

#create phyloseq object
theme_set(theme_bw())
set.seed(123)

#get rid of the extra stuff at the end of rownames
rownames(seqtab.nochim) <- sub("\\_.*", "", rownames(seqtab.nochim))
samples.out <- rownames(seqtab.nochim)
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxa))

ps <- prune_samples(sample_names(ps) != "NonLibraryControl", ps) # Remove control sample
ps <- prune_samples(sample_names(ps) != "NonTemplateControl", ps) # Remove control sample

##cleaning taxonomy table
tax_table(ps)[, colnames(tax_table(ps))] <- gsub(tax_table(ps)[, colnames(tax_table(ps))],     pattern = "[a-z]__", replacement = "")

## Remove non-fungal ASVs
ps <- subset_taxa(ps, Kingdom=="Fungi") # 
ntaxa(ps) #1582 ASVs left


#shorten ASV names but retain DNA sequence
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

## Plot sample library size 
df <- as.data.frame(sample_data(ps)) 
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
plot_df <- ggplot(data=df, aes(x = Index, y = LibrarySize)) + geom_point()
plot_df

#rarefaction curve to see if we need to drop samples or rarefy
asvdf <- as.matrix(as.data.frame(otu_table(ps)))
sort(sample_sums(ps))
asv.rarecurve = rarecurve(asvdf, step = 100, label = FALSE, abline(v= 1121, col = "blue"), col = "black") #1121 lowest # of reads after removing samples less than 1000
asv.rarecurve = rarecurve(asvdf, step = 100, label = FALSE, abline(v= 2008, col = "red"), col = "black") #would lose 20 samples and capture majority of diversity but not all
asv.rarecurve = rarecurve(asvdf, step = 100, label = FALSE, abline(v= 3004, col = "purple"), col = "black") #would lose 34 samples and capture almost all diversity

#let's zoom in a bit
asv.rarecurve = rarecurve(asvdf, step = 20, label = FALSE, xlim=c(0, 4000), abline(v= 1121, col = "blue"), col = "black")
asv.rarecurve = rarecurve(asvdf, step = 100, label = FALSE, xlim=c(0, 4000), abline(v= 2008, col = "red"), col = "black")
asv.rarecurve = rarecurve(asvdf, step = 100, label = FALSE, xlim=c(0, 4000), abline(v= 3004, col = "purple"), col = "black")

# so if we choose to rarefy to 2008, we will lose 20 samples and will capture most of the diversity but not all.

#remove any samples with less than 1000 reads
ps <- prune_samples(sample_sums(ps) >= 1000, ps) #154 samples remaining, lost 4
ps

#check taxa table is still correct
ps <- prune_taxa(taxa_sums(ps@otu_table) > 0, ps) # 1575 total taxa
sum(sample_sums(ps)) # 1165175 total reads

#save phyloseq object
saveRDS(ps, "phyloseq_obj.rds")

saveRDS(ps, here::here("output/og.ps.rds")) 

#visualize alpha diversity just to make sure everything looks normal
plot_richness(ps, x="Treatment", measures=c("Shannon", "Simpson"), color="Lifestage")


# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Lifestage", shape ="Treatment", title="Bray NMDS", label = "Specimen_ID")
#SBadultMale5fall looks really weird, might remove later after exploration
#ps <- prune_samples(sample_names(ps) != "SBTadult5Mfall", ps) # Remove weird sample



###let's subset to just look at ticks (not soil)

ticks <- subset_samples(ps, Species=="BLT")

asvdf_ticks <- as.data.frame(otu_table(ticks))
asv.rarecurve.ticks = rarecurve(asvdf_ticks, step = 10, label = F)


### all ticks
#alphas
plot_richness(ticks, x="Season", measures=c("Shannon"), color="Treatment")+
  facet_grid(~Lifestage)

#ordinate OTUs for all ticks
ticks.prop <- transform_sample_counts(ticks, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ticks.prop, method="NMDS", distance="bray")
p1 <- plot_ordination(ticks.prop, ord.nmds.bray, label = "Specimen_ID", color="Treatment", shape ="Lifestage", title="Bray NMDS")
p1
p1 + facet_wrap(~Season)
#sample 128 looks pretty weird but it's my only nymph from that site for the fall season so I kinda wanna leave it for now

ordu_ticks = ordinate(ticks.prop, "PCoA", "jaccard")
plot_ordination(ticks.prop, ordu_rr, color="Understory", shape="Treatment")+
  facet_grid(~Lifestage)

#all ticks
ticks %>%                                           #phyloseq object
  plot_richness(
    x = "Treatment",                                 #change to the variable you want to compare
    measures = c("Shannon")) +           #choose diversity measure
  facet_wrap(~Understory)+
  geom_boxplot(aes(fill = Lifestage), show.legend = TRUE)

#sample colors
sample_colors <- c("lightblue", "orchid","lightgreen","orange")

#violin plot
ticks %>%                                                              #phyloseq object
  plot_richness(
    x = "Season",                                                    #compare diversity of datatype
    measures = c("Shannon")) +                              #choose diversity measures
  geom_boxplot(aes(fill = Treatment), show.legend = TRUE )+             #make violin plot, set fill aes to sampletype
  facet_wrap(~Lifestage+Sex)+
  theme_linedraw()+                                                     #change theme to classic
  xlab(NULL)+                                                           #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),                     #adjust y-axis text
        axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 1, angle = 90),           #adjust x-axis label position
        axis.title.y = element_text(size = 10))+                        #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 10))+           #adjust headings
  scale_fill_manual(values = sample_colors)+                            #set fill colors
  theme(plot.title=element_text(size = 10, face = "bold", hjust = 0.5)) #change title size, face and position


alphadiv <- estimate_richness(ticks, measures = c("Shannon")) %>%
  rownames_to_column(var = "Unique_Specimen_ID") %>%
  left_join(as.data.frame(sample_data(ticks)), by = "Unique_Specimen_ID") 

#test for differences in alpha
kruskal.test(Shannon ~ Season, data = alphadiv) #NS
kruskal.test(Shannon ~ Understory, data = alphadiv) #S
#data:  Shannon by Understory
#Kruskal-Wallis chi-squared = 12.201, df = 1, p-value = 0.0004778
kruskal.test(Shannon ~ Treatment, data = alphadiv) #S
#data:  Shannon by Treatment
#Kruskal-Wallis chi-squared = 20.015, df = 3, p-value = 0.0001686
kruskal.test(Shannon ~ Lifestage, data = alphadiv) #NS
kruskal.test(Shannon ~ Sex, data = alphadiv) #NS
kruskal.test(Shannon ~ Site, data = alphadiv) #S 
#Kruskal-Wallis chi-squared = 11.271, df = 1, p-value = 0.0007872


#do post-hoc tests
pairwise.wilcox.test(alphadiv$Shannon, alphadiv$Treatment, p.adjust.method = 'fdr') 
# data:  alphadiv$Shannon and alphadiv$Treatment 
# 
#             burned mowed   unmanaged
#   mowed     0.647  -       -        
#   unmanaged 0.060  5.8e-05 -        
#   unmowed   0.671  0.070   0.012   
# 
# P value adjustment method: fdr 





##and now to separate sites
richville <- subset_samples(ticks, Site=="RR")
sugarbush <- subset_samples(ticks, Site=="SB")

#alphas
plot_richness(richville, x="Season", measures=c("Shannon", "Simpson"), color="Treatment")
plot_richness(sugarbush, x="Season", measures=c("Shannon", "Simpson"), color="Treatment")

#ordinate OTUs for burn/unburned
rr.prop <- transform_sample_counts(richville, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(rr.prop, method="NMDS", distance="bray")
p1 <- plot_ordination(rr.prop, ord.nmds.bray, color="Treatment", shape ="Season", title="Bray NMDS")
p1
p1 + facet_wrap(~Sex)

ordu_rr = ordinate(rr.prop, "PCoA", "jaccard")
plot_ordination(rr.prop, ordu_rr, color="Treatment", shape="Season")


#ordinate OTUs for sugarbush
sb.prop <- transform_sample_counts(sugarbush, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(sb.prop, method="NMDS", distance="bray")
p2 <- plot_ordination(sb.prop, ord.nmds.bray, color="Treatment", shape ="Season", label = "specimen_id", title="Bray NMDS")
p2
p2 + facet_wrap(~Lifestage)

ordu_sb = ordinate(sb.prop, "PCoA", "jaccard")
plot_ordination(sb.prop, ordu_sb, color="Season", shape="Treatment")

## okay there are no super obvious patterns? Annoying but okay

#sugarbush
sugarbush %>%                                           #phyloseq object
  plot_richness(
    x = "Treatment",                                 #change to the variable you want to compare
    measures = c("Shannon")) +           #choose diversity measure
  facet_wrap(~Lifestage)+
  geom_boxplot(aes(fill = Treatment), show.legend = FALSE)

#sample colors
sample_colors <- c("lightblue", "orchid","lightgreen")

#violin plot
sugarbush %>%                                                              #phyloseq object
  plot_richness(
    x = "Treatment",                                                    #compare diversity of datatype
    measures = c("Shannon")) +                              #choose diversity measures
  geom_boxplot(aes(fill = Treatment), show.legend = FALSE)+             #make violin plot, set fill aes to sampletype
  facet_wrap(~Season)+
  theme_linedraw()+                                                     #change theme to classic
  xlab(NULL)+                                                           #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),                     #adjust y-axis text
        axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 10))+                        #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 10))+           #adjust headings
  scale_fill_manual(values = sample_colors)+                            #set fill colors
  theme(plot.title=element_text(size = 10, face = "bold", hjust = 0.5)) #change title size, face and position


alphadiv <- estimate_richness(sugarbush, measures = c("Shannon")) %>%
  rownames_to_column(var = "Unique_Specimen_ID") %>%
  left_join(as.data.frame(sample_data(sugarbush)), by = "Unique_Specimen_ID") 

#test for differences in alpha
kruskal.test(Shannon ~ Season, data = alphadiv) #NS
kruskal.test(Shannon ~ Treatment, data = alphadiv) #S
# data:  Shannon by Treatment
# Kruskal-Wallis chi-squared = 3.9658, df = 1, p-value = 0.04643
kruskal.test(Shannon ~ Lifestage, data = alphadiv) #NS
kruskal.test(Shannon ~ Sex, data = alphadiv) #NS ooo so close! 0.06




##burned plots
#richville
richville %>%                                           #phyloseq object
  plot_richness(
    x = "Treatment",                                 #change to the variable you want to compare
    measures = c("Shannon")) +           #choose diversity measure
 # facet_wrap(~Lifestage)+
  geom_boxplot(aes(fill = Treatment), show.legend = FALSE)

#sample colors
#sample_colors <- c("lightblue", "orchid","lightgreen")

#violin plot
richville %>%                                                              #phyloseq object
  plot_richness(
    x = "Treatment",                                                    #compare diversity of datatype
    measures = c("Shannon")) +                              #choose diversity measures
  geom_boxplot(aes(fill = Treatment), show.legend = FALSE)+             #make violin plot, set fill aes to sampletype
  #facet_wrap(~Season)+
  theme_linedraw()+                                                     #change theme to classic
  xlab(NULL)+                                                           #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),                     #adjust y-axis text
        axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 10))+                        #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 10))+           #adjust headings
  scale_fill_manual(values = sample_colors)+                            #set fill colors
  theme(plot.title=element_text(size = 10, face = "bold", hjust = 0.5)) #change title size, face and position


alphadiv <- estimate_richness(richville, measures = c("Shannon")) %>%
  rownames_to_column(var = "Unique_Specimen_ID") %>%
  left_join(as.data.frame(sample_data(richville)), by = "Unique_Specimen_ID") 

#test for differences in alpha
kruskal.test(Shannon ~ Season, data = alphadiv) #NS
kruskal.test(Shannon ~ Treatment, data = alphadiv) #S
# data:  Shannon by Treatment
# Kruskal-Wallis chi-squared = 4.6922, df = 1, p-value = 0.0303
kruskal.test(Shannon ~ Sex, data = alphadiv) #NS



