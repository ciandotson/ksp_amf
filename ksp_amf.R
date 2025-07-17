# This is an R script on how to use the dada2 pipeline to denoise Illumina amplicon sequencing data # 
# into analyzed sequence variants (ASVs). This tutorial is based on the published #
# tutorial (https://benjjneb.github.io/dada2/tutorial_1_8.html) and is adapted to fit # 
# the Kankakee Sand Prairie fungal LSU sequencing data #

# After cloning this repository, the script will automatically change the working directory to the #
# cloned repository as long as it is cloned directly into your home directory #
setwd("~/ksp_amf")

#### Argument Parsing ####
# This portion makes it possible for the entire R script to be run from the command line. This just basically #
# tells the R script where to find the forward reads to take as input, making it very easy to run the whole #
# script in one go #

# If required, the R package `optparse` is downloaded #
if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse")
library(optparse); packageVersion("optparse")

# A special list of all the parsable objects is made #
option_list <- list(
  make_option("--forward", type = "character", help = "filepath containing the raw, untrimmed reads"))

# The list is changed to a format that the parasable objects are able to be read #
opt <- parse_args(OptionParser(option_list=option_list))

# The forward reads are saved as the `for_reads` object #
for_reads <- opt$forward

# This string of commands allows the Rscript to find the forward reads, which, for those seeking to reproduce #
# our results, will have to download the reads from Sequence Read Archive (SRA). All of the other relevant information #
# for this pipeline, such as the metadata and the non-R scripts, are found in the github and are automatically downloaded #
# when the repository is cloned. As such, the remaining relevant files can be called from the cloned repository #
# from within the script without having to do any additional downloads #

#### Installing and Initializing the pipeline ####
# this takes about a half an hour to both install and update all of your packages. Basically, you 
# you use BiocManager (which, if you haven't installed, will be installed with this command)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.18", force = TRUE)
library(dada2); packageVersion('dada2')

BiocManager::install("ShortRead")
library(ShortRead); packageVersion('ShortRead')

BiocManager::install("Biostrings")
library(Biostrings); packageVersion('Biostrings')

# Read in the metadata from the cloned repository #
all.met <- read.csv2(file = './metadata/Fungal_Community_Soil_Samples_KSP.csv', sep = ',')

# The metadata is reorganized such that the samples are ordered based on the number of sample they are #
rownames(all.met) <- all.met$Sample
all.met$Order <- as.numeric(gsub("JW_", "", all.met$Number))
all.met <- all.met[order(all.met$Order),]

# Make sure input value from the command line outputs the file paths for the raw forward #
for.fp <- for_reads
list.files(for.fp)

# As a sanity check, add the file paths to the metadata data.frame and make sure the filepaths align to the sample names #
all.met$Forward <- for.fp
sample.names <- rownames(all.met)

#### Trimming Primers Using cutadapt ####
## Though cutadapt is a bash (command line) tool, we can actually use bash tools directly from the R command line as both languages are Unix based ## 
## As long as cutadapt is installed prior to running the R.script, then it will work! This can be performed by following the dircetions of the README.md ##

# Create a string for the forward and reverse primers that will be manipulated to tell cutadapt what sequence motifs to cut from the reads #
for.primer <- "CAGCCGCGGTAATTCCAGCT"

# We need to know if our primers are in the correct orientation within our reads, so we need to be able to see if the original primer, its reverse, its complement, or its reverse complement is found in the reads #
# We can do this by first creating a function that takes a primer and outputs these four orientations as below #
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer) # converts primer string to DNAString object #
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna), RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString)) # outputs the string of the four orientations #
}

# Now we can actually find the orientations for the forward primer #
for.ornt <- allOrients(for.primer)
for.ornt

## Ambiguous base calls (base call N) make this process difficult, so we just remove any reads that have them! ##
# Create a new directory with associated filepaths to output the "pre-filtered" forward reads #
if(!dir.exists('./reads'))
  dir.create('./reads')
if(!dir.exists('./reads/filtN'))
  dir.create(('./reads/filtN'))
forfilt.fp <- file.path('./reads/filtN', paste0(sample.names, '_filtN_R1.fastq.gz'))

# Perform the filtering with the only unqiue parameter being maxN = 0 #
prefilt.track <- filterAndTrim(for.fp, forfilt.fp, rev.fp, revfilt.fp, maxN = 0, multithread = FALSE)
# The results are written directly to the new files we denoted in the last two commands, so nothing new needs to be saved #

# To compare all of the orientations of our primers to the newly pre-filtered reads, we can make a function that lines up the orientations to these reads and tells us the number of hits #
primer.hits <- function(primer, fp){
  # function that outputs the number of hits a primer orientation matches to a read #
  nhits <- vcountPattern(primer, sread(readFastq(fp)), fixed = FALSE)
  return(sum(nhits>0))
}

# now we can make a for loop to tell us how many times the different orientations match the reads #
n.hits <- matrix(data = c(0,0,0,0),nrow=1, ncol =4)
n.hits <- as.data.frame(n.hits)

for(i in 1:length(forfilt.fp)){
  n.hits[1,1] <- n.hits[1,1] + primer.hits(for.ornt[1], forfilt.fp[[i]])
  n.hits[1,2] <- n.hits[1,2] + primer.hits(for.ornt[2], forfilt.fp[[i]])
  n.hits[1,3] <- n.hits[1,3] + primer.hits(for.ornt[3], forfilt.fp[[i]])
  n.hits[1,4] <- n.hits[1,4] + primer.hits(for.ornt[4], forfilt.fp[[i]])
}

## Now that we know that the primers are (mostly) where they are supposed to be, we can actually use cutadapt to trim them off ##
# Make a directory and filepaths for the reads after the primers have been trimmed # 
path.cut <- "./reads/ptrim"
if(!dir.exists(path.cut)) dir.create(path.cut)
forcut.fp <- file.path(path.cut, paste0(sample.names, '_ptrim_R1.fastq.gz'))

# Take the reverse complement of the primers #
for.rc <- dada2::rc(for.primer)

# Estbalish the paramters for cutadapt #
R1.flags <- paste("-g", for.primer) # tells cutadapt the sequence of the forward primer to cut #

# Actual cutadapt command that loops through all files #
for(i in seq_along(for.fp)){
  system2(cutadapt, args = c(R1.flags, "-n", 2, "-o", forcut.fp[i], forfilt.fp[i]))
}

## Checking to make sure all of the primers were trimmed ##
t.hits <- matrix(data = c(0,0,0,0),nrow=1, ncol =4)
t.hits <- as.data.frame(t.hits)

for(i in 1:nrow(for.fp)){
  t.hits[1,1] <- t.hits[1,1] + primer.hits(for.ornt[1], forcut.fp[[i]])
  t.hits[1,2] <- t.hits[1,2] + primer.hits(for.ornt[2], forcut.fp[[i]])
  t.hits[1,3] <- t.hits[1,3] + primer.hits(for.ornt[3], forcut.fp[[i]])
  t.hits[1,4] <- t.hits[1,4] + primer.hits(for.ornt[4], forcut.fp[[i]])
}

#### dada2 Pipeline Implementation ####
# Here is where we start using dada2 for what it was meant to do: denoise the reads into ASVs. We start by doing some additional filtering to make the pipeline go faster #
# Make the directory and file paths for the final filtering of the reads #
path.filt <- "./reads/filtered"
if(!dir.exists(path.filt)) dir.create(path.filt)
forpost.fp <- file.path("./reads/filtered", paste0(sample.names, '_F_filt.fastq.gz'))

# Filter out reads that have any ambiguous base calls (maxN = 0) and reads with expected errors above 2 (maxEE = 2), while truncating reads to 200 bp (truncLen = 200, as described in Lekberg et al., 2023) or at first instance in which expected error is above 2 (truncQ = 2) # 
postfilt.track <- filterAndTrim(forcut.fp, forpost.fp, maxN = 0, maxEE = 2, truncLen = 200,
                     truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = FALSE, verbose = TRUE)

# With filtering finished, dada2 can learn the error rates that are specific to the given dataset #
for.er <- learnErrors(forpost.fp, multithread=FALSE, verbose = TRUE)

# Another way to help speed up the process is by dereplication, or only saving one copy of each each unique sequence, which is what is done below #
for.derep <- derepFastq(forpost.fp, verbose = TRUE)

# Finally, dada2 can take the error model to denoise the derepelicated sequneces into Analyzed Sequence Variants (ASVs) #
for.dada <- dada(for.derep, err = for.er, multithread = FALSE, verbose = TRUE)

save.image("./ksp_amf.RData")

# The denoised reads can be more simply represented with an ASV table, which is produced below #
ksp.st <- makeSequenceTable(for.dada)

# Chimeras, or artifacts of DNA amplification and/or sequencing, can be simply removed using the below command
nochim_ksp.st <- removeBimeraDenovo(ksp.st, method = 'consensus', multithread = FALSE, verbose = TRUE)

# Finally, we can track how many reads passed each step of the pipeline with the code below #
getN <- function(x) sum(getUniques(x))
final.track <- cbind(prefilt.track[,1], prefilt.track[,2], postfilt.track[,2], sapply(for.dada, getN), colSums(nochim_ksp.st))
colnames(final.track) <- c("pre-cutadapt", "post-cutadapt", "filtered", "denoised", "nonchim")
final.track <- as.data.frame(final.track)

#### Assigning Taxonomy ####
# To assign taxonomy, we use the MaarJAM database as a reference that has been adapted to be read by dada2 #
ksp.taxa <- assignTaxonomy(nochim_ksp.st, "./reference/maarjam_dada2.fasta", multithread = FALSE, verbose = TRUE)
ksp.taxa <- as.matrix(ksp.taxa)

save.image("./ksp_amf.RData")