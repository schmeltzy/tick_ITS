## dealing with ITS amplicon data where reads may have read into the primer on the other side of too-short amplicons
## 2 x 301 bp MiSeq
## Overview: preprocessing reads using cutadapt and dada2 to remove primers, adapters, quality trimming and filtering, and chimera removal after merging. This is based on the tutorial on the dada2 website.

#load necessary packages
library(dada2)
library(ShortRead)
library(Biostrings)
library(DECIPHER)
#library(phangorn)
library(here)

#set seed for reproducibility
set.seed(123)

#change path to directory where the fastq files live
path <- "~/Documents/Computing/Tick_ITS_redo_7.31.25/Fastq/amplicons"
list.files(path)

#generate matched lists of forward and reverse files and parse out the sample names
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#identify primers
FWD <- "GGCTTGGTCATTTAGAGGAAGTAA" 
REV <- "CGGCTGCGTTCTTCATCGATGC" 

#verify presence and orientation of primers
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

#check for presence of ambiguous bases
fnFs.filtN <- file.path(path, "filtN_both", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN_both", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) #this will take awhile

#count primers in forward and reverse reads
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), FWD.ReverseReads = sapply(FWD.orients,
    primerHits, fn = fnRs.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
    fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#remove primers with cutadapt
cutadapt <- "/Users/eschmeltzer/miniforge3/envs/cutadapt/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

#cutadapt has a lot of output just fyi
path.cut <- file.path(path, "cutadapt_both")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-m", 20, # -m 20 to remove zero-length reads and specify min length =20 bp
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
} ## this will take awhile too

##sanity check! count the presence of primers in the first cutadapt-ed sample
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients,
                                     primerHits, fn = fnRs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
                                    fn = fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

list.files(path.cut)

#generate new matched lists of forward and reverse files and parse out the sample names
primersremFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
primersremRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(primersremFs, get.sample.name))
head(sample.names)

#inspect read quality profiles
#forward
plotQualityProfile(primersremFs[9:11]) #drop in quality around 230 where we trimmed primers
#reverse
plotQualityProfile(primersremRs[9:11])


#try without truncating
filt_bf_bowFs <- file.path(path.cut, "filtered_bf_bowtie", basename(primersremFs))
filt_bf_bowRs <- file.path(path.cut, "filtered_bf_bowtie", basename(primersremRs))

out <- filterAndTrim(primersremFs, filt_bf_bowFs, primersremRs, filt_bf_bowRs, maxN = 0, maxEE = c(8, 8), truncQ = 8,
                           minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out) #lose ~70% of reads


#let's learn error rates
errF <- learnErrors(filt_bf_bowFs, multithread = TRUE)
errR <- learnErrors(filt_bf_bowRs, multithread = TRUE)

#visualize error rates
plotErrors(errF, nominalQ = TRUE)


##core sample inference
dadaFs <- dada(filt_bf_bowFs, err = errF, multithread = TRUE)
dadaRs <- dada(filt_bf_bowRs, err = errR, multithread = TRUE)

#merge pairs
mergers <- mergePairs(dadaFs, filt_bf_bowFs, dadaRs, filt_bf_bowRs, verbose=TRUE)

#construct ASV table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#calculate proportion of non-chimeric seqs
sum(seqtab.nochim)/sum(seqtab)
sum(seqtab.nochim) # 1960531 reads


#inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.nochim)))

#track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN),
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "dada_merged.csv")


#assign taxonomy using UNITE database with non-fungal Eukaryotic outgroups
unite.ref <- "~/Documents/Computing/Tick_ITS_and_Meta_6.18.25/unite_db/sh_general_release_all_19.02.2025/sh_general_release_dynamic_all_19.02.2025.fasta"
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)
#this takes a long time

#inspect taxonomic assignments
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
dim(taxa) # 3285 ASVs


# #make phylogenetic tree
# sequences<-getSequences(seqtab.nochim)
# names(sequences)<-sequences
# alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA) #this will take awhile
# phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
# dm <- dist.ml(phang.align) #this will take awhile
# treeNJ <- NJ(dm) # Note, tip order != sequence order
# fit = pml(treeNJ, data=phang.align)
# fitGTR <- update(fit, k=4, inv=0.2)
# fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                     rearrangement = "NNI", control = pml.control(trace = 0)) #this takes awhile
# #this last command will take many many many many many hours
# write.tree(fitGTR$tree, file = "GTR.phy")

#save files for later use in phyloseq
saveRDS(taxa, "taxa_both.rds")
saveRDS(seqtab.nochim, "ASV_both.rds")
#saveRDS(fitGTR, "fitGTR_both.rds")







