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
asv.rarecurve = rarecurve(asvdf, step = 20, label = FALSE, xlim=c(0, 4000), abline(v= 2008, col = "red"), col = "black")
asv.rarecurve = rarecurve(asvdf, step = 20, label = FALSE, xlim=c(0, 4000), abline(v= 3004, col = "purple"), col = "black")

# so if we choose to rarefy to 2008, we will lose 20 samples total and will capture most of the diversity but not all.

#remove any samples with less than 1000 reads
ps <- prune_samples(sample_sums(ps) >= 1000, ps) #154 samples remaining, lost 4
ps

#check taxa table is still correct
ps <- prune_taxa(taxa_sums(ps@otu_table) > 0, ps) # 1575 total taxa
sum(sample_sums(ps)) # 1165175 total reads

#save phyloseq object
saveRDS(ps, "phyloseq_obj.rds")

saveRDS(ps, here::here("output/og.ps.rds")) 


## Let's see how rarefying goes and we can compare
#can really only use this for observed richness estimates and not for beta diversity analyses
#rarefy_even_depth automatically removes ASVs that become 0 in abundance
ps.rare <- rarefy_even_depth(ps, sample.size = 2008, rngseed = 999) #97 ASVs removed, 16 samples removed

ps.rare #1478 ASVs, 138 samples
sum(sample_sums(ps.rare)) #277104
sample_sums(ps.rare)
rarefied_names <- sample_names(ps.rare)
names <- sample_names(ps)
write.table(rarefied_names, file = "rarefied_names.txt", sep = "\t")
write.table(names, file = "full_names.txt", sep = "\t")

summary(taxa_sums(ps.rare)) # minimum times a taxa appears: 1.0, 1st quartile: 9.0, median: 29.0, mean: 187.5     3rd quantile: 100.0 , max: 32418.0 

saveRDS(ps.rare, here::here("output/ps.rare.rds"))
