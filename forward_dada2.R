## REDOING DADA2 only using forward reads to see if we have better results
## 2 x 301 bp MiSeq
#load necessary packages
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(DECIPHER)
library(phangorn)

#change path to directory where the fastq files live
path <- "~/Documents/Computing/Tick_ITS_redo_7.31.25/Fastq/amplicons"
list.files(path)

#generate matched lists of forward and reverse files and parse out the sample names
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
#fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

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
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
#fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE)

#count primers in forrward and reverse reads
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
                                              fn = fnFs.filtN[[1]]))

#remove primers with cutadapt
cutadapt <- "/Users/eschmeltzer/miniforge3/envs/cutadapt/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

#cutadapt has a lot of output just fyi
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
#fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
#R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-m", 20, # -m 20 to remove zero-length reads and specify min length =20 bp
                             "-o", fnFs.cut[i], # output files
                             fnFs.filtN[i])) # input files
}

##sanity check! count the presence of primers in the first cutadapt-ed sample
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
                                                                                                      fn = fnFs.cut[[1]]))

#uhhhh okay there are still primers in there, bummer. Just a few though (single digits) so we'll move on

####alright let's put in dada2!
#generate new matched lists of forward and reverse files and parse out the sample names
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
#cutRs <- sort(list.files(path.cut, pattern = "_2.fastq.gz", full.names = TRUE))


# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)


#inspect read quality profiles
#forward
plotQualityProfile(cutFs[5:7]) #drop in quality around 200 
#reverse
#plotQualityProfile(cutRs[1:2])

#filter and trim
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
#filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = c(2), truncQ = 2,
                     minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out)

#Still losing over half of reads...

#let's learn error rates
errF <- learnErrors(filtFs, multithread = TRUE)

#errR <- learnErrors(filtRs, multithread = TRUE)

#visualize error rates
plotErrors(errF, nominalQ = TRUE)


##core sample inference
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
#dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

#merge pairs
#mergers <- mergePairs(dadaFs, filtFs, verbose=TRUE)

#construct ASV table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#calculate proportion of non-chimeric seqs
sum(seqtab.nochim)/sum(seqtab)

sum(seqtab.nochim) # 1273010 reads

#inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.nochim)))

#track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN),
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "dada_forwards.csv")

#assign taxonomy using UNITE database with non-fungal Eukaryotic outgroups
unite.ref <- "~/Documents/Computing/Tick_ITS_and_Meta_6.18.25/unite_db/sh_general_release_all_19.02.2025/sh_general_release_dynamic_all_19.02.2025.fasta"
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

saveRDS(taxa, "taxa.rds")
saveRDS(seqtab.nochim, "ASV.rds")

#inspect taxonomic assignments
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
dim(taxa) # 4325 ASVs


#make phylogenetic tree
sequences<-getSequences(seqtab.nochim)
names(sequences)<-sequences
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA) #this will take awhile
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align) #this will take awhile
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "NNI", control = pml.control(trace = 0)) #this takes awhile
#this last command will take many many many many many hours
write.tree(fitGTR$tree, file = "GTR.phy")

#save files for later use in phyloseq
saveRDS(taxa, "taxa.rds")
saveRDS(seqtab.nochim, "ASV.rds")
saveRDS(fitGTR, "fitGTR.rds")





