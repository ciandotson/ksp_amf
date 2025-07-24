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
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse")
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

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

if (!requireNamespace("dada2", quietly = TRUE)) BiocManager::install("dada2")
library(dada2); packageVersion('dada2')

if (!requireNamespace("ShortRead", quietly = TRUE)) BiocManager::install("ShortRead")
library(ShortRead); packageVersion('ShortRead')

if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
library(Biostrings); packageVersion('Biostrings')

if(!requireNamespace("cgwtools", quietly = TRUE)) install.packages("cgwtools")
library(cgwtools); packageVersion("cgwtools")

# Read in the metadata from the cloned repository #
ksp.met <- read.csv2(file = './metadata/Fungal_Community_Soil_Samples_KSP.csv', sep = ',')

# The metadata is reorganized such that the samples are ordered based on the number of sample they are #
rownames(ksp.met) <- ksp.met$Sample
ksp.met$Order <- as.numeric(gsub("JW_", "", ksp.met$Number))
ksp.met <- ksp.met[order(ksp.met$Order),]

# Add data entries for site location and amount of innoculant #
for(i in 1:nrow(ksp.met)){ 
  ksp.met$Site[i] <- substr(ksp.met$Sample[i], 1,1)
  if(substr(ksp.met$Sample[i], 2,2) == 'C'){
    ksp.met$Treatment[i] <- 'Control'
  }else{
    ksp.met$Treatment[i] <- substr(ksp.met$Sample[i], 2,3) 
  }
}
ksp.met$Treatment <- gsub('HI', 'High', ksp.met$Treatment)
ksp.met$Treatment <- gsub('LO', 'Low', ksp.met$Treatment)
ksp.met$Treatment <- gsub('yc', 'MycoBloom', ksp.met$Treatment)

# Make sure input value from the command line outputs the file paths for the raw forward #
for.fp <- for_reads
list.files(for.fp)

# As a sanity check, add the file paths to the metadata data.frame and make sure the filepaths align to the sample names #
ksp.met$Forward <- for.fp
sample.names <- rownames(ksp.met)

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
prefilt.track <- filterAndTrim(for.fp, forfilt.fp, maxN = 0, multithread = TRUE, verbose = TRUE)
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

save.image("./ksp_amf.RData")

# Actual cutadapt command that loops through all files #
for(i in seq_along(list.files(for.fp))){
  system(paste0("cutadapt -g ", for.primer, " -o ", forcut.fp[i], " ", forfilt.fp[i]))
}

## Checking to make sure all of the primers were trimmed ##
t.hits <- matrix(data = c(0,0,0,0),nrow=1, ncol =4)
t.hits <- as.data.frame(t.hits)

for(i in 1:length(forcut.fp)){
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
                     truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE, verbose = TRUE)

# With filtering finished, dada2 can learn the error rates that are specific to the given dataset #
for.er <- learnErrors(forpost.fp, multithread=TRUE, verbose = TRUE)

# Another way to help speed up the process is by dereplication, or only saving one copy of each each unique sequence, which is what is done below #
for.derep <- derepFastq(forpost.fp,verbose = TRUE)

# Finally, dada2 can take the error model to denoise the derepelicated sequneces into Analyzed Sequence Variants (ASVs) #
for.dada <- dada(for.derep, err = for.er, multithread = TRUE, verbose = TRUE)

save.image("./ksp_amf.RData")

# The denoised reads can be more simply represented with an ASV table, which is produced below #
ksp.st <- makeSequenceTable(for.dada)

# Chimeras, or artifacts of DNA amplification and/or sequencing, can be simply removed using the below command
nochim_ksp.st <- removeBimeraDenovo(ksp.st, method = 'consensus', multithread = TRUE, verbose = TRUE)

save(nochim_ksp.st, file = "./abridged.RData")

if(!requireNamespace('cgwtools')) install.packages('cgwtools')
library(cgwtools); packageVersion('cgwtools')
resave(ksp.met, file = "./abridged.RData")

# Finally, we can track how many reads passed each step of the pipeline with the code below #
getN <- function(x) sum(getUniques(x))
final.track <- cbind(prefilt.track[,1], prefilt.track[,2], postfilt.track[,2], sapply(for.dada, getN), colSums(nochim_ksp.st))
colnames(final.track) <- c("pre-cutadapt", "post-cutadapt", "filtered", "denoised", "nonchim")
final.track <- as.data.frame(final.track)

#### Assigning Taxonomy ####
# To assign taxonomy, we use the MaarJAM database as a reference that has been adapted to be read by dada2 #
ksp.taxa <- assignTaxonomy(nochim_ksp.st, "./reference/maarjam_dada2.fasta", multithread = TRUE, verbose = TRUE)
ksp.taxa <- as.data.frame(ksp.taxa)
colnames(ksp.taxa) <- c('Family', "Genus", 'Species')

# Finally, we can change the format of our ASV table ("nochim_ksp.st") and taxonomy tables ("ksp.taxa") into more convenient formats #
ksp.taxa <- as.matrix(ksp.taxa)

# make sure the rownames of the metadata tables matches the column names of the ASV table #
rownames(nochim_ksp.st) <- rownames(ksp.met)

resave(ksp.taxa, file = "./abridged.RData")
save.image("./ksp_amf.RData")

#### Constructing phyloseq object ####
# The phyloseq object is a means of combining and summarizing all relevant features of a microbiome dataset through an ASV table, #
# a metadata table, a taxonomy table, a phylogenetic tree, and a DNAStringSet Object. We will be focusing on the first three objects, #
# which we have already made. #

if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
library(phyloseq); packageVersion("phyloseq")

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr); packageVersion('dplyr')

# Here we actually make the phyloseq object (raw_ksp.ps) that contains the ASV table (nochim_ksp.st), taxonomy table (ksp.taxa), #
# and metadata table (ksp.met) #
nochim_ksp.st <- t(nochim_ksp.st)
raw_ksp.ps <- phyloseq(otu_table(nochim_ksp.st, taxa_are_rows = TRUE),
                       tax_table(ksp.taxa),
                       sample_data(ksp.met))

resave(raw_ksp.ps, file = "./abridged.RData")

# Filter out samples that did not have any matches to the MaarJAM database #
raw_ksp.ps <- subset_taxa(raw_ksp.ps, !is.na(Family))

# We can save the DNA sequences of each ASV as they are saved as the row names of the ASV and taxonomy tables and then change the taxa names to "ASVn", #
# where n is the number ASV the sequence corresponds to, which is ranked by total abundance across all samples #
raw_ksp.dna <- Biostrings::DNAStringSet(taxa_names(raw_ksp.ps))
names(raw_ksp.dna) <- taxa_names(raw_ksp.ps)
raw_ksp.ps <- merge_phyloseq(raw_ksp.ps, raw_ksp.dna)
taxa_names(raw_ksp.ps) <- paste0("ASV", seq(ntaxa(raw_ksp.ps)))

# We can take this a step further and add the lowest level classification to the ASV in parentheses so we don't have to flip between taxonomy tables for exploratory # 
raw_ksp.tax <- as.data.frame(tax_table(raw_ksp.ps))
for(i in 1:nrow(raw_ksp.tax)){
  if(!is.na(raw_ksp.tax$Species[i])){
    taxa_names(raw_ksp.ps)[i] <- paste0(taxa_names(raw_ksp.ps)[i], '(', raw_ksp.tax$Species[i], ')')
  } else if(!is.na(raw_ksp.tax$Genus[i])){
    taxa_names(raw_ksp.ps)[i] <- paste0(taxa_names(raw_ksp.ps)[i], '(', raw_ksp.tax$Genus[i], ')')
  } else if(!is.na(raw_ksp.tax$Family[i])){
    taxa_names(raw_ksp.ps)[i] <- paste0(taxa_names(raw_ksp.ps)[i], '(', raw_ksp.tax$Family[i], ')')
  }
}

# Since we are interested in seeing if the inoculant communities makes it into the incumbent community, we can make a separate phyloseq object that just has the MycoBloom Community #
myc.ps <- subset_samples(raw_ksp.ps, Treatment == "MycoBloom")
myc.ps <- subset_taxa(myc.ps, taxa_sums(myc.ps) > 0)

resave(myc.ps, file = "./abridged.RData")

# This function below is a handy tool that decomposes phyloseq objects to a list of data.frames that make the data more easily manipulated and available, especially for seeing #
# how taxonomy may be related to abundance when we look at the "fra" data.frames (short for Frankenstein as it is just the ASV table appended to the taxonomy table) #
decompose_ps <- function(ps, label){
  # function that decomposes a phyloseq object into separate data.frame and refseq objects (does not include tree) #
  tax.tab <- as.data.frame(tax_table(ps))
  otu.tab <- as.data.frame(otu_table(ps))
  met.tab <- as(sample_data(ps), 'data.frame')
  dna.tab <- refseq(ps)
  fra.tab <- cbind(tax.tab, otu.tab)
  decomposed = list(
    tax = tax.tab,
    otu = otu.tab,
    met = met.tab,
    dna = dna.tab,
    fra = fra.tab
  )
  assign(label, decomposed, envir = .GlobalEnv)
  invisible(decomposed)
}

# To make this decomposed phyloseq object, we just supply the phyloseq object and what we want to call it #
decompose_ps(raw_ksp.ps, "raw_ksp")

# Each component can be accessed by typing the name of the decomposed phyloseq object ("raw_ksp") followed by an accessor ("$") and the three letter abbreviation of the object #
# ("tax" for taxonomy table, "otu" for asv table, "met" for metadata table, "dna" for DNAStringSet, and "fra" for the Frankenstein table) #
raw_ksp$fra

resave(raw_ksp, file = "./abridged.RData")

#### Taxonomy Validation ####
# Here we blast all of the sequences we have returned and make sure that the asvs we get back actually correspond to something that is at least mostly fungal #
if(!requireNamespace("rBLAST")) BiocManager::install("rBLAST")
library(rBLAST); packageVersion('rBLAST')

# Create local blast database from the 16S rRNA database using rBLAST #
blast.tar <- blast_db_get("SSU_eukaryote_rRNA.tar.gz", baseURL = 'https://ftp.ncbi.nlm.nih.gov/blast/db/', check_update = TRUE)
untar(blast.tar, exdir = './reference/SSU_database')
list.files('./reference/SSU_database')
blast.db <- blast(db = './reference/SSU_database/SSU_eukaryote_rRNA')

# Performs the blast for each read and returns the best hit #
ksp.hits <- matrix(nrow = nrow(raw_ksp$tax), ncol = 12)
ksp.hits <- as.data.frame(ksp.hits) 
hold <- c()
for(i in 1:length(raw_ksp$dna)){
  hold <- predict(blast.db, raw_ksp$dna[i])
  ksp.hits[i,] <- hold[1,]
  raw_ksp$tax$Best_Hit[i] <- hold[1, 2]
}

# Filter out reads that do not correspond to a NCBI entry #
filt_ksp.tax <- filter(raw_ksp$tax, !is.na(raw_ksp$tax$Best_Hit))

# Output the resulting NCBI entry names to a list #
if(!dir.exists("./blast_hits")){
  dir.create('./blast_hits')
}
write.table(filt_ksp.tax$Best_Hit, './blast_hits/ksp_blast_hits.txt')

# Call the python script to retrieve the taxonomies of the matched entries #
system('python3 ~/ksp_amf/SSU_BLAST.py -i ./blast_hits/ksp_blast_hits.txt -o ./blast_hits/ksp_ncbi_hits.csv')

# Read in the output from the python script and make new taxonomy table "s_ncbi_fin.tax" #
ksp_ncbi.taxa <- read.csv2('./blast_hits/ksp_ncbi_hits.csv', header = FALSE, fill = TRUE)
ksp_ncbi.int <- strsplit(as.character(ksp_ncbi.taxa$V1), ",")
ksp_ncbi_fin.tax <- do.call(rbind, lapply(ksp_ncbi.int, function(x) { length(x) <- max(sapply(ksp_ncbi.int, length)); x }))

# Join the two taxa tables into `double_ksp.tax` (double because it is validated by both BLAST and the MaarJAM databse) and use it to filter the old phyloseq object; #
# make the concatenated taxa table the taxa table for the phyloseq object #
double_ksp.tax <- cbind(filt_ksp.tax, ksp_ncbi_fin.tax)
ksp.ps <- subset_taxa(raw_ksp.ps, taxa_names(raw_ksp.ps) %in% rownames(double_ksp.tax))
tax_table(ksp.ps) <- as.matrix(double_ksp.tax)

# Filter out ASVs that correspond to non-fungal targets #
final_ksp.ps <- subset_taxa(ksp.ps, X2 == "Fungi")

#### Phylogenetic Tree Construction for Soils ####
# Add Outgroups to the data to try and catch any non-AMF reads that have made it through #
decompose_ps(final_ksp.ps, 'final_ksp')
out.rna <- readRNAStringSet('./reference/outgroup.fasta')
out.dna <- DNAStringSet(out.rna)
names(out.dna) <- c("Outgroup1", "Outgroup2")
tree.dna <- c(final_ksp$dna, out.dna)

# Output the reads into a fasta file #
writeXStringSet(tree.dna, "./reads/ksp_input.fasta")

# Perform a multiple sequence alignment using MAFFT #
system('mafft --auto --thread -1 ./reads/ksp_input.fasta > ./reads/ksp_aligned.fasta')

# Construct a tree using IQTree with a general time reversible model with a gamma distribution and invariant site copies #
system('iqtree -s ./reads/ksp_aligned.fasta -m GTR+G+I -nt AUTO')

if(!requireNamespace('ape')) install.packages('ape')
library(ape); packageVersion('ape')

if(!requireNamespace('phytools')) install.packages('phytools')
library(phytools); packageVersion('phytools')

# Read the tree into R and change the tip labels such that they match what was originally input # 
ksp.tre <- read.tree("./reads/ksp_aligned.fasta.treefile")
ksp.tre$tip.label <- sub("^(ASV[0-9]+)_([^_]+)_$", "\\1(\\2)", ksp.tre$tip.label)

# Save the tree into the phyloseq object #
phy_tree(final_ksp.ps) <- ksp.tre

# Now we can actually do the filtering using the tree by first labelling our outgroups #
out.nam <- c("Outgroup1", "Outgroup2")

# Find the most recent common ancestor node for the two outgroups #
out.mrca <- getMRCA(ksp.tre, out.nam)

# Find all of the descendants of this MRCA #
out.des <- getDescendants(ksp.tre, out.mrca)
out.tip <- ksp.tre$tip.label[out.des[out.des <= length(ksp.tre$tip.label)]]

# Finally, remove the taxa denoted in the tips denoted in out.tip from the phyloseq object #
final_ksp.ps <- subset_taxa(final_ksp.ps, !taxa_names(final_ksp.ps) %in% out.tip)

# Remove the control samples #
final_ksp.ps <- subset_samples(final_ksp.ps, Treatment != "TC")
final_ksp.ps <- subset_taxa(final_ksp.ps, taxa_sums(final_ksp.ps) > 0)

# Reorganize the ASV table in decresaing order of taxa abundance #
decompose_ps(final_ksp.ps, 'final_ksp')
final_ksp$otu <- arrange(final_ksp$otu, desc(rowSums(final_ksp$otu)))
final_ksp$tax <- final_ksp$tax[rownames(final_ksp$otu),]

# Resave the new final phyloseq object #
final_ksp.ps <- phyloseq(otu_table(final_ksp$otu, taxa_are_rows = TRUE),
                         tax_table(as.matrix(final_ksp$tax)),
                         sample_data(final_ksp$met),
                         refseq(final_ksp$dna),
                         phy_tree(phy_tree(final_ksp.ps)))

# Rename taxa so they correspond to their new ordering based on abundance #
taxa_names(final_ksp.ps) <- paste0("ASV", seq(ntaxa(final_ksp.ps)))
final_ksp.tax <- as.data.frame(tax_table(final_ksp.ps))
for(i in 1:nrow(final_ksp.tax)){
  if(!is.na(final_ksp.tax$Species[i])){
    if(final_ksp.tax$Species[i] != "sp.") {
      taxa_names(final_ksp.ps)[i] <- paste0(taxa_names(final_ksp.ps)[i], '(', final_ksp.tax$Species[i], ')')}
      else(taxa_names(final_ksp.ps)[i] <- paste0(taxa_names(final_ksp.ps)[i], '(', final_ksp.tax$Genus[i], ')'))
  } else if(!is.na(final_ksp.tax$Genus[i])){
    taxa_names(final_ksp.ps)[i] <- paste0(taxa_names(final_ksp.ps)[i], '(', final_ksp.tax$Genus[i], ')')
  } else if(!is.na(final_ksp.tax$Family[i])){
    taxa_names(final_ksp.ps)[i] <- paste0(taxa_names(final_ksp.ps)[i], '(', final_ksp.tax$Family[i], ')')
  }
}

decompose_ps(final_ksp.ps, 'final_ksp')
resave(final_ksp.ps, file = "./hold.RData") 

# Since we are interested in seeing if the inoculant communities makes it into the incumbent community, we can make a separate phyloseq object that just has the MycoBloom Community #
myc.ps <- subset_samples(final_ksp.ps, Treatment == "MycoBloom")
myc.ps <- subset_taxa(myc.ps, taxa_sums(myc.ps) > 0)
decompose_ps(myc.ps, 'myc')
myc$fra <- arrange(myc$fra, desc(myc$fra$MycoBloom1))

save.image("ksp_amf.RData")

#### Alpha Diversity Measurement and Visualization ####
ksp.rich <- estimate_richness(final_ksp.ps) # automatically performs alpha diversity calculations #
ksp.rich <- cbind(final_ksp$met, ksp.rich)
myc.rich <- filter(ksp.rich, Treatment == "MycoBloom")
ksp.rich <- filter(ksp.rich, Treatment != "MycoBloom")

# ANOVA for each kind of alpha diversity 
ksp.sha <- aov(Shannon ~ Site*Treatment, ksp.rich) 
summary(ksp.sha)
ksp_sha.hsd <- TukeyHSD(ksp.sha)

ksp.sim <- aov(Simpson ~ Site*Treatment, ksp.rich) 
summary(ksp.sim)
ksp_sim.hsd <- TukeyHSD(ksp.sim)

ksp.cha <- aov(Chao1 ~ Site*Treatment, ksp.rich) 
summary(ksp.cha)
ksp_cha.hsd <- TukeyHSD(ksp.cha)

if(!requireNamespace('ggprism')) install.packages('ggprism')
library(ggprism); packageVersion('ggprism')

if(!requireNamespace('ggplot2')) install.packages('ggplot2')
library(ggplot2); packageVersion('ggplot2')

ksp_rich.mnsd <- ksp.rich %>%
  group_by(Site, Treatment) %>%
  summarize(
    sha.mean = mean(Shannon),
    evn.mean = mean(Simpson),
    cha.mean = mean(Chao1),
    sha.sd = sd(Shannon),
    evn.sd = sd(Simpson),
    cha.sd = sd(Chao1),
    .groups = "drop" # Prevent grouping in the result
  )

ksp_rich.mnsd <- as.data.frame(ksp_rich.mnsd)

# Individual Site  Diversity # 
a.rich <- filter(ksp.rich, Site == 'A')
b.rich <- filter(ksp.rich, Site == 'B')
c.rich <- filter(ksp.rich, Site == 'C')
d.rich <- filter(ksp.rich, Site == 'D')
e.rich <- filter(ksp.rich, Site == 'E')
f.rich <- filter(ksp.rich, Site == 'F')
g.rich <- filter(ksp.rich, Site == 'G')
h.rich <- filter(ksp.rich, Site == 'H')

a.rich <- arrange(a.rich, Treatment)
b.rich <- arrange(b.rich, Treatment)
c.rich <- arrange(c.rich, Treatment)
d.rich <- arrange(d.rich, Treatment)
e.rich <- arrange(e.rich, Treatment)
f.rich <- arrange(f.rich, Treatment)
g.rich <- arrange(g.rich, Treatment)
h.rich <- arrange(h.rich, Treatment)

# Shannon Diversity #
if(!requireNamespace("multcompView")) install.packages("multcompView")
library(multcompView); packageVersion('multcompView')

a.sha <- aov(Shannon~Treatment, a.rich)
summary(a.sha)
a_sha.hsd <- TukeyHSD(a.sha)
a_sha.hsd
a_sha.let <- multcompLetters4(a.sha, a_sha.hsd)

b.sha <- aov(Shannon~Treatment, b.rich)
summary(b.sha)
b_sha.hsd <- TukeyHSD(b.sha)
b_sha.hsd
b_sha.let <- multcompLetters4(b.sha, b_sha.hsd)

c.sha <- aov(Shannon~Treatment, c.rich)
summary(c.sha)
c_sha.hsd <- TukeyHSD(c.sha)
c_sha.hsd
c_sha.let <- multcompLetters4(c.sha, c_sha.hsd)

d.sha <- aov(Shannon~Treatment, d.rich)
summary(d.sha)
d_sha.hsd <- TukeyHSD(d.sha)
d_sha.hsd
d_sha.let <- multcompLetters4(d.sha, d_sha.hsd)

e.sha <- aov(Shannon~Treatment, e.rich)
summary(e.sha)
e_sha.hsd <- TukeyHSD(e.sha)
e_sha.hsd
e_sha.let <- multcompLetters4(e.sha, e_sha.hsd)

f.sha <- aov(Shannon~Treatment, f.rich)
summary(f.sha)
f_sha.hsd <- TukeyHSD(f.sha)
f_sha.hsd
f_sha.let <- multcompLetters4(f.sha, f_sha.hsd)

g.sha <- aov(Shannon~Treatment, g.rich)
summary(g.sha)
g_sha.hsd <- TukeyHSD(g.sha)
g_sha.hsd
g_sha.let <- multcompLetters4(g.sha, g_sha.hsd)

h.sha <- aov(Shannon~Treatment, h.rich)
summary(h.sha)
h_sha.hsd <- TukeyHSD(h.sha)
h_sha.hsd
h_sha.let <- multcompLetters4(h.sha, h_sha.hsd)

sha.let <- c(a_sha.let$Treatment$Letters,
             b_sha.let$Treatment$Letters,
             c_sha.let$Treatment$Letters,
             d_sha.let$Treatment$Letters,
             e_sha.let$Treatment$Letters,
             f_sha.let$Treatment$Letters,
             g_sha.let$Treatment$Letters,
             h_sha.let$Treatment$Letters)
plot.rich <- cbind(ksp_rich.mnsd, sha.let)

plot.rich$Treats <- factor(plot.rich$Treatment, levels = c("Control", "Low", "High"))
plot.rich$sha.let[7] <- 'ab'
plot.rich$sha.let[8] <- 'b'
plot.rich$sha.let[9] <- 'a'
plot.rich$sha.let[19] <- 'b'
plot.rich$sha.let[20] <- 'a'
plot.rich$sha.let[22] <- 'b'
plot.rich$sha.let[23] <- 'a'
plot.rich$sha.let[24] <- 'b'

sha.plot <- ggplot(plot.rich, aes(x = Site, y = sha.mean, fill = Treats, color = Treats)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  theme_prism() +
  ylab('Shannon Diversity') +
  scale_fill_manual(values = c("white", "gray", "#4D4D4D", 'black')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,6)) +
  ggtitle('Shannon Diversity') +
  geom_text(aes(label = sha.let, y = sha.mean + sha.sd + 0.1), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 12) +
  geom_errorbar(aes(ymin = sha.mean - sha.sd, ymax = sha.mean + sha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2)

# Simpson Diversity #
a.sim <- aov(Simpson~Treatment, a.rich)
summary(a.sim)
a_sim.hsd <- TukeyHSD(a.sim)
a_sim.hsd
a_sim.let <- multcompLetters4(a.sim, a_sim.hsd)

b.sim <- aov(Simpson~Treatment, b.rich)
summary(b.sim)
b_sim.hsd <- TukeyHSD(b.sim)
b_sim.hsd
b_sim.let <- multcompLetters4(b.sim, b_sim.hsd)

c.sim <- aov(Simpson~Treatment, c.rich)
summary(c.sim)
c_sim.hsd <- TukeyHSD(c.sim)
c_sim.hsd
c_sim.let <- multcompLetters4(c.sim, c_sim.hsd)

d.sim <- aov(Simpson~Treatment, d.rich)
summary(d.sim)
d_sim.hsd <- TukeyHSD(d.sim)
d_sim.hsd 
d_sim.let <- multcompLetters4(d.sim, d_sim.hsd)

e.sim <- aov(Simpson~Treatment, e.rich)
summary(e.sim)
e_sim.hsd <- TukeyHSD(e.sim)
e_sim.hsd
e_sim.let <- multcompLetters4(e.sim, e_sim.hsd)

f.sim <- aov(Simpson~Treatment, f.rich)
summary(f.sim)
f_sim.hsd <- TukeyHSD(f.sim)
f_sim.hsd
f_sim.let <- multcompLetters4(f.sim, f_sim.hsd)

g.sim <- aov(Simpson~Treatment, g.rich)
summary(g.sim)
g_sim.hsd <- TukeyHSD(g.sim)
g_sim.hsd
g_sim.let <- multcompLetters4(g.sim, g_sim.hsd)

h.sim <- aov(Simpson~Treatment, h.rich)
summary(h.sim)
h_sim.hsd <- TukeyHSD(h.sim)
h_sim.hsd
h_sim.let <- multcompLetters4(h.sim, h_sim.hsd)

sim.let <- c(a_sim.let$Treatment$Letters,
             b_sim.let$Treatment$Letters,
             c_sim.let$Treatment$Letters,
             d_sim.let$Treatment$Letters,
             e_sim.let$Treatment$Letters,
             f_sim.let$Treatment$Letters,
             g_sim.let$Treatment$Letters,
             h_sim.let$Treatment$Letters)
plot.rich <- cbind(plot.rich, sim.let)

plot.rich$sim.let[20] <- 'a'
plot.rich$sim.let[21] <- 'ab'
plot.rich$sim.let[23] <- 'a'
plot.rich$sim.let[24] <- 'ab'


evn.plot <- ggplot(plot.rich, aes(x = Site, y = evn.mean, fill = Treats, color = Treats)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  theme_prism() +
  ylab("Simpson's Diversity") +
  scale_fill_manual(values = c("white", "gray", "#4D4D4D", 'black')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,1.1)) +
  ggtitle("Simpson's Diversity") +
  geom_text(aes(label = sim.let, y = evn.mean + evn.sd + 0.01), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 12) +
  geom_errorbar(aes(ymin = evn.mean - evn.sd, ymax = evn.mean + evn.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2)

# Otu Richness #
a.cha <- aov(Chao1~Treatment, a.rich)
summary(a.cha)
a_cha.hsd <- TukeyHSD(a.cha)
a_cha.hsd
a_cha.let <- multcompLetters4(a.cha, a_cha.hsd)

b.cha <- aov(Chao1~Treatment, b.rich)
summary(b.cha)
b_cha.hsd <- TukeyHSD(b.cha)
b_cha.hsd
b_cha.let <- multcompLetters4(b.cha, b_cha.hsd)

c.cha <- aov(Chao1~Treatment, c.rich)
summary(c.cha)
c_cha.hsd <- TukeyHSD(c.cha)
c_cha.hsd
c_cha.let <- multcompLetters4(c.cha, c_cha.hsd)

d.cha <- aov(Chao1~Treatment, d.rich)
summary(d.cha)
d_cha.hsd <- TukeyHSD(d.cha)
d_cha.hsd # Low vs. control #
d_cha.let <- multcompLetters4(d.cha, d_cha.hsd)

e.cha <- aov(Chao1~Treatment, e.rich)
summary(e.cha)
e_cha.hsd <- TukeyHSD(e.cha)
e_cha.hsd
e_cha.let <- multcompLetters4(e.cha, e_cha.hsd)

f.cha <- aov(Chao1~Treatment, f.rich)
summary(f.cha)
f_cha.hsd <- TukeyHSD(f.cha)
f_cha.hsd
f_cha.let <- multcompLetters4(f.cha, f_cha.hsd)

g.cha <- aov(Chao1~Treatment, g.rich)
summary(g.cha)
g_cha.hsd <- TukeyHSD(g.cha)
g_cha.hsd
g_cha.let <- multcompLetters4(g.cha, g_cha.hsd)

h.cha <- aov(Chao1~Treatment, h.rich)
summary(h.cha)
h_cha.hsd <- TukeyHSD(h.cha)
h_cha.hsd
h_cha.let <- multcompLetters4(h.cha, h_cha.hsd)

cha.let <- c(a_cha.let$Treatment$Letters,
             b_cha.let$Treatment$Letters,
             c_cha.let$Treatment$Letters,
             d_cha.let$Treatment$Letters,
             e_cha.let$Treatment$Letters,
             f_cha.let$Treatment$Letters,
             g_cha.let$Treatment$Letters,
             h_cha.let$Treatment$Letters)
plot.rich <- cbind(plot.rich, cha.let)

cha.plot <- ggplot(plot.rich, aes(x = Site, y = cha.mean, fill = Treats, color = Treats)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  theme_prism() +
  ylab("ASV Richness") +
  scale_fill_manual(values = c("white", "gray", "#4D4D4D", 'black')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,250)) +
  ggtitle("Chao1 Diversity") +
  geom_text(aes(label = cha.let, y = cha.mean + cha.sd + 0.1), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 12) +
  geom_errorbar(aes(ymin = cha.mean - cha.sd, ymax = cha.mean + cha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2)

if(!requireNamespace('patchwork')) install.packages('patchwork')
library(patchwork); packageVersion('patchwork')

alpha.plot <- (cha.plot) /
(evn.plot) /
(sha.plot) +
  plot_layout(guides = 'collect')
  

#### Beta Diversity ####
set.seed(248)
final_ksp_prop.ps <- transform_sample_counts(final_ksp.ps, function(otu) otu/sum(otu))
ord.nmds.wuni <- ordinate(final_ksp_prop.ps, method="NMDS", distance="wunifrac")
ksp.bdist <- phyloseq::distance(final_ksp.ps, method = "wunifrac")

if(!requireNamespace("vegan")) install.package("vegan")
library(vegan); packageVersion('vegan')
ksp.perm <- adonis2(ksp.bdist~Site*Treatment, data = final_ksp$met)
ksp.perm

plot_ordination(final_ksp_prop.ps, ord.nmds.wuni, color="Site", shape = 'Treatment', title="NMDS") +
  theme_prism() +
  geom_point(size = 6) 

#### Per Sample Analysis ####
if(!requireNamespace("devtools")) install.packages('devtools')
library(devtools); packageVersion('devtools')

if(!requireNamespace("pairwiseAdonis")) devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis); packageVersion("pairwiseAdonis")

# A #
a.ps <- subset_samples(final_ksp.ps, Site == 'A')
a.met <- as(sample_data(a.ps), 'data.frame')
a_prop.ps <- transform_sample_counts(a.ps, function(otu) otu/sum(otu))
a_ord.nmds.wuni <- ordinate(a_prop.ps, method="NMDS", distance="wunifrac")
a.bdist <- phyloseq::distance(a.ps, method = "wunifrac")

a.perm <- adonis2(a.bdist~Treatment, data = a.met)
a.perm

a.pair <- pairwise.adonis2(a.bdist~Treatment, data = a.met)
a.pair

plot_ordination(a_prop.ps, a_ord.nmds.wuni, color="Treatment", title="A Samples NMDS") +
  theme_prism() +
  geom_point(size = 6)

# B #
b.ps <- subset_samples(final_ksp.ps, Site == 'B')
b.met <- as(sample_data(b.ps), 'data.frame')
b_prop.ps <- transform_sample_counts(b.ps, function(otu) otu/sum(otu))
b_ord.nmds.wuni <- ordinate(b_prop.ps, method="NMDS", distance="wunifrac")
b.bdist <- phyloseq::distance(b.ps, method = "wunifrac")

b.perm <- adonis2(b.bdist~Treatment, data = b.met)
b.perm

b.pair <- pairwise.adonis2(b.bdist~Treatment, data = b.met)
b.pair

plot_ordination(b_prop.ps, b_ord.nmds.wuni, color="Treatment", title="B Samples NMDS") +
  theme_prism() +
  geom_point(size = 6)

# C #
c.ps <- subset_samples(final_ksp.ps, Site == 'C')
c.met <- as(sample_data(c.ps), 'data.frame')
c_prop.ps <- transform_sample_counts(c.ps, function(otu) otu/sum(otu))
c_ord.nmds.wuni <- ordinate(c_prop.ps, method="NMDS", distance="wunifrac")
c.bdist <- phyloseq::distance(c.ps, method = "wunifrac")

c.perm <- adonis2(c.bdist~Treatment, data = c.met)
c.perm

c.pair <- pairwise.adonis2(c.bdist~Treatment, data = c.met)
c.pair

plot_ordination(c_prop.ps, c_ord.nmds.wuni, color="Treatment", title="C Samples NMDS") +
  theme_prism() +
  geom_point(size = 6) 

# D #
d.ps <- subset_samples(final_ksp.ps, Site == 'D')
d.met <- as(sample_data(d.ps), 'data.frame')
d_prop.ps <- transform_sample_counts(d.ps, function(otu) otu/sum(otu))
d_ord.nmds.wuni <- ordinate(d_prop.ps, method="NMDS", distance="wunifrac")
d.bdist <- phyloseq::distance(d.ps, method = "wunifrac")

d.perm <- adonis2(d.bdist~Treatment, data = d.met)
d.perm

d.pair <- pairwise.adonis2(d.bdist~Treatment, data = d.met)
d.pair

plot_ordination(d_prop.ps, d_ord.nmds.wuni, color="Treatment", title="D Samples NMDS") +
  theme_prism() +
  geom_point(size = 6) 

# E #
e.ps <- subset_samples(final_ksp.ps, Site == 'E')
e.met <- as(sample_data(e.ps), 'data.frame')
e_prop.ps <- transform_sample_counts(e.ps, function(otu) otu/sum(otu))
e_ord.nmds.wuni <- ordinate(e_prop.ps, method="NMDS", distance="wunifrac")
e.bdist <- phyloseq::distance(e.ps, method = "wunifrac")

e.perm <- adonis2(e.bdist~Treatment, data = e.met)
e.perm

e.pair <- pairwise.adonis2(e.bdist~Treatment, data = e.met)
e.pair

plot_ordination(e_prop.ps, e_ord.nmds.wuni, color="Treatment", title="E Samples NMDS") +
  theme_prism() +
  geom_point(size = 6) 

# F #
f.ps <- subset_samples(final_ksp.ps, Site == 'F')
f.met <- as(sample_data(f.ps), 'data.frame')
f_prop.ps <- transform_sample_counts(f.ps, function(otu) otu/sum(otu))
f_ord.nmds.wuni <- ordinate(f_prop.ps, method="NMDS", distance="wunifrac")
f.bdist <- phyloseq::distance(f.ps, method = "wunifrac")

f.perm <- adonis2(f.bdist~Treatment, data = f.met)
f.perm

f.pair <- pairwise.adonis2(f.bdist~Treatment, data = f.met)
f.pair

plot_ordination(f_prop.ps, f_ord.nmds.wuni, color="Treatment", title="F Samples NMDS") +
  theme_prism() +
  geom_point(size = 6) 

# G #
g.ps <- subset_samples(final_ksp.ps, Site == 'G')
g.met <- as(sample_data(g.ps), 'data.frame')
g_prop.ps <- transform_sample_counts(g.ps, function(otu) otu/sum(otu))
g_ord.nmds.wuni <- ordinate(g_prop.ps, method="NMDS", distance="wunifrac")
g.bdist <- phyloseq::distance(g.ps, method = "wunifrac")

g.perm <- adonis2(g.bdist~Treatment, data = g.met)
g.perm

g.pair <- pairwise.adonis2(g.bdist~Treatment, data = g.met)
g.pair

plot_ordination(g_prop.ps, g_ord.nmds.wuni, color="Treatment", title="G Samples NMDS") +
  theme_prism() +
  geom_point(size = 6) 

# H #
h.ps <- subset_samples(final_ksp.ps, Site == 'G')
h.met <- as(sample_data(h.ps), 'data.frame')
h_prop.ps <- transform_sample_counts(h.ps, function(otu) otu/sum(otu))
h_ord.nmds.wuni <- ordinate(h_prop.ps, method="NMDS", distance="wunifrac")
h.bdist <- phyloseq::distance(h.ps, method = "wunifrac")

h.perm <- adonis2(h.bdist~Treatment, data = h.met)
h.perm

h.pair <- pairwise.adonis2(h.bdist~Treatment, data = b.met)
h.pair

plot_ordination(h_prop.ps, h_ord.nmds.wuni, color="Treatment", title="H Samples NMDS") +
  theme_prism() +
  geom_point(size = 6)

#### Stacked Histograms ####
# First we start by making a color pallette for each unique ASV #
if(!requireNamespace("Polychrome")) install.packages("Polychrome")
library(Polychrome); packageVersion("Polychrome")
ksp.colr <- createPalette(ntaxa(final_ksp.ps),  c("#ff0000", "#00ff00", "#0000ff"))
ksp.colr <- as.data.frame(ksp.colr)
rownames(ksp.colr) <- taxa_names(final_ksp.ps)

# Add a gray color for "Other"
ksp.colr[2003,] <- "#D4D4D4" 
rownames(ksp.colr)[2003] <- "Other" 

# Save 'ASV' as it's own unique level of taxonomy #
final_ksp$tax$ASV <- taxa_names(final_ksp.ps)
tax_table(final_ksp.ps) <- as.matrix(final_ksp$tax)

# Save a phyloseq Object that conatins only the samples of interest #
myc.ps <- subset_samples(final_ksp.ps, Treatment == "MycoBloom")
myc.ps <- subset_taxa(myc.ps, taxa_sums(myc.ps) > 0)

# Collect only the top 19 ASVs, and group the remaining ASVs into "Other" and save their names in their orders of abundance to make a subsetted palette for plotting #
if(!requireNamespace("microbiomeutilities")) devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities); packageVersion('microbiomeutilities')

if(!requireNamespace("microbiome")) BiocManager::install("microbiome")
library(microbiome); packageVersion("microbiome")

hg_myc.ps <- aggregate_top_taxa2(myc.ps, top = 9, level = "ASV")
hg_myc.name <- names(sort(taxa_sums(hg_myc.ps), decreasing = TRUE))
hg_myc.colr <- hg_myc.colr[hg_myc.name,]

# Save the phyloseq object as a data.frame and make factors that guide the plot what to plot #
hg_myc.df <- psmelt(hg_myc.ps)
hg_myc.df$ASVs <- factor(hg_myc.df$ASV, levels = hg_myc.name)
hg_myc.df$Soil <- factor(hg_myc.df$Treatment, levels = c("Control", "Low", "High"))

# Plot the histogram #
hg_myc.plot <- ggplot(hg_myc.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(name = "Fungal ASV", values = hg_myc.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "MycoBloom")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(color = "black", size = 18, family = "Liberation Sans", angle = -45, vjust = 0.6, hjust = 0.1),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right')




save.image("./ksp_amf.RData")