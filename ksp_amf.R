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

save.image("./ksp_amf.RData")

