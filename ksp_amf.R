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

# A finalial list of all the parsable objects is made #
option_list <- list(
  make_option("--forward", type = "character", help = "filepath containing the raw, untrimmed reads"))

# The list is changed to a format that the parasable objects are able to be read #
opt <- parse_args(OptionParser(option_list=option_list))

# The forward reads are saved as the `for_reads` object #
for_reads <- opt$forward

# This command moves all of the output to a log file 'ksp_amf.log 
sink(file = './ksp_amf.log', append = TRUE, type = c("output", "message")) # Redirects stdout (e.g., print, cat) and stderr (e.g., warnings, message) #
cat("## Script started at", Sys.time(), "\n\n")
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
ksp.met$Forward <- list.files(for.fp)
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
t.hits <- matrix(data = c(0,0,0,0),nrow=1, ncol = 4)
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

# Save the ASV table to a new .RData object that will serve as the abridged data produced by this pipeline #
save(nochim_ksp.st, file = "./abridged.RData")

# Do the same thing for the metadata #
if(!requireNamespace('cgwtools')) install.packages('cgwtools')
library(cgwtools); packageVersion('cgwtools')

resave(ksp.met, file = "./abridged.RData")

# make sure the rownames of the metadata tables matches the column names of the ASV table #
rownames(nochim_ksp.st) <- rownames(ksp.met)

# Finally, we can track how many reads passed each step of the pipeline with the code below #
getN <- function(x) sum(getUniques(x))
final.track <- cbind(prefilt.track[,1], prefilt.track[,2], postfilt.track[,2], sapply(for.dada, getN), rowSums(nochim_ksp.st))
rownames(final.track) <- rownames(ksp.met)
colnames(final.track) <- c("pre-cutadapt", "post-cutadapt", "filtered", "denoised", "nonchim")
final.track <- as.data.frame(final.track)

resave(final.track, file = 'abridged.RData')

#### Assigning Taxonomy ####
# To assign taxonomy, we use the MaarJAM database as a reference that has been adapted to be read by dada2 #
ksp.taxa <- assignTaxonomy(nochim_ksp.st, "./reference/maarjam_dada2.fasta", minBoot = 0, outputBootstraps = TRUE, multithread = TRUE, verbose = TRUE)
ksp.taxa <- as.data.frame(ksp.taxa)
colnames(ksp.taxa) <- c('Family', "Genus", 'Species', 'Family_Boot', 'Genus_Boot', 'Species_Boot')

# Finally, we can change the format of our ASV table ("nochim_ksp.st") and taxonomy tables ("ksp.taxa") into more convenient formats #
ksp.taxa <- as.matrix(ksp.taxa)

# Save the taxonomy table to the abridged .RData file #
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

# This function below is a handy tool that decomposes phyloseq objects to a list of data.frames that make the data more easily manipulated and available, especially for seeing #
# how taxonomy may be related to abundance when we look at the "fra" data.frames (short for Frankenstein as it is just the ASV table appended to the taxonomy table) #
decompose_ps <- function(ps, label){
  # function that decomposes a phyloseq object into separate data.frame and refseq objects (does not include tree) #
  tax.tab <- as.data.frame(tax_table(ps))
  otu.tab <- as.data.frame(otu_table(ps))
  met.tab <- as(sample_data(ps), 'data.frame')
  if(!is.null(access(ps, "refseq"))){
    dna.tab <- refseq(ps)
  }
  fra.tab <- cbind(tax.tab, otu.tab)
  if(!is.null(access(ps, "refseq"))){
    decomposed = list(
    tax = tax.tab,
    otu = otu.tab,
    met = met.tab,
    dna = dna.tab,
    fra = fra.tab)
  }else {
    decomposed = list(
      tax = tax.tab,
      otu = otu.tab,
      met = met.tab,
      fra = fra.tab) 
  }
  assign(label, decomposed, envir = .GlobalEnv)
  invisible(decomposed)
}

# To make this decomposed phyloseq object, we just supply the phyloseq object and what we want to call it #
decompose_ps(raw_ksp.ps, "raw_ksp")

# Each component can be accessed by typing the name of the decomposed phyloseq object ("raw_ksp") followed by an accessor ("$") and the three letter abbreviation of the object #
# ("tax" for taxonomy table, "otu" for asv table, "met" for metadata table, "dna" for DNAStringSet, and "fra" for the Frankenstein table) #
raw_ksp$fra

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

# Read in the output from the python script and make new taxonomy table "ksp_ncbi_fin.tax" #
ksp_ncbi.taxa <- read.csv2('./blast_hits/ksp_ncbi_hits.csv', header = FALSE, fill = TRUE)
ksp_ncbi.int <- strsplit(as.character(ksp_ncbi.taxa$V1), ",")
ksp_ncbi_fin.tax <- do.call(rbind, lapply(ksp_ncbi.int, function(x) { length(x) <- max(sapply(ksp_ncbi.int, length)); x }))

# Join the two taxa tables into `double_ksp.tax` (double because it is validated by both BLAST and the MaarJAM databse) and use it to filter the old phyloseq object; #
# make the concatenated taxa table the taxa table for the phyloseq object #
double_ksp.tax <- cbind(filt_ksp.tax, ksp_ncbi_fin.tax)
ksp.ps <- subset_taxa(raw_ksp.ps, taxa_names(raw_ksp.ps) %in% rownames(double_ksp.tax))
tax_table(ksp.ps) <- as.matrix(double_ksp.tax)

# Filter out ASVs that correspond to non-fungal targets #
fun_ksp.ps <- subset_taxa(ksp.ps, X2 == "Fungi")

#### Phylogenetic Tree Construction for Soils ####
# Add Outgroups to the data to try and catch any non-AMF reads that have made it through #
decompose_ps(fun_ksp.ps, 'fun_ksp')
out.rna <- readRNAStringSet('./reference/outgroup.fasta')
out.dna <- DNAStringSet(out.rna)
names(out.dna) <- c("Outgroup1", "Outgroup2")
tree.dna <- c(fun_ksp$dna, out.dna)

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
phy_tree(fun_ksp.ps) <- ksp.tre

# Now we can actually do the filtering using the tree by first labelling our outgroups #
out.nam <- c("Outgroup1", "Outgroup2")

# Find the most recent common ancestor node for the two outgroups #
out.mrca <- getMRCA(ksp.tre, out.nam)

# Find all of the descendants of this MRCA #
out.des <- getDescendants(ksp.tre, out.mrca)
out.tip <- ksp.tre$tip.label[out.des[out.des <= length(ksp.tre$tip.label)]]

# Finally, remove the taxa denoted in the tips denoted in out.tip from the phyloseq object #
fun_ksp.ps <- subset_taxa(fun_ksp.ps, !taxa_names(fun_ksp.ps) %in% out.tip)

# Remove the control samples #
fun_ksp.ps <- subset_samples(fun_ksp.ps, Treatment != "TC")
fun_ksp.ps <- subset_taxa(fun_ksp.ps, taxa_sums(fun_ksp.ps) > 0)

# To give us the option to remove half of the control samples ofr the sake of having balanced statitsical tests, I have made this function to remove one of each replicate. #
# The decsion for a replicate to be selected is based purely on the number of reads the replicate has, with the sample with higher reads being selected. #
con.ps <- subset_samples(fun_ksp.ps, Treatment == "Control")
decompose_ps(con.ps, 'con')
con$met$Reads <- colSums(con$otu)
remn.met <- c()
for(i in LETTERS[1:8]){
  group <- filter(con$met, Site == i)
  for(j in 1:3){
    grouped <- filter(group, Replicate == j)
    if(grouped$Reads[1] <= grouped$Reads[2]) {
      remn.met <- rbind(remn.met, grouped[1,])
    }else(remn.met <- rbind(remn.met, grouped[2,]))
  }
}

# Save only the control samples that are NOT found in the list #
filt_ksp.ps <- subset_samples(fun_ksp.ps, !sample_names(fun_ksp.ps) %in% rownames(remn.met))
filt_ksp.ps <- subset_taxa(filt_ksp.ps, taxa_sums(final_ksp.ps) > 0)

## One last filtering step: only save ASVs with at least 0.1% prevalence within each group of samples #
# Subset the data by Site (8 phyloseq) objects #
a.ps <- subset_samples(filt_ksp.ps, Site == "A")
b.ps <- subset_samples(filt_ksp.ps, Site == "B")
c.ps <- subset_samples(filt_ksp.ps, Site == "C")
d.ps <- subset_samples(filt_ksp.ps, Site == "D")
e.ps <- subset_samples(filt_ksp.ps, Site == "E")
f.ps <- subset_samples(filt_ksp.ps, Site == "F")
g.ps <- subset_samples(filt_ksp.ps, Site == "G")
h.ps <- subset_samples(filt_ksp.ps, Site == "H")

# Transform each subset of the data using Total Sum Scaling $
a.ps <- transform_sample_counts(a.ps, function(x) x/sum(x))
b.ps <- transform_sample_counts(b.ps, function(x) x/sum(x))
c.ps <- transform_sample_counts(c.ps, function(x) x/sum(x))
d.ps <- transform_sample_counts(d.ps, function(x) x/sum(x))
e.ps <- transform_sample_counts(e.ps, function(x) x/sum(x))
f.ps <- transform_sample_counts(f.ps, function(x) x/sum(x))
g.ps <- transform_sample_counts(g.ps, function(x) x/sum(x))
h.ps <- transform_sample_counts(h.ps, function(x) x/sum(x))

# Only save ASVs that correspond to at least 0.5% of all reads (0.005) #
a.ps <- subset_taxa(a.ps, taxa_sums(a.ps) > 0.005)
b.ps <- subset_taxa(b.ps, taxa_sums(b.ps) > 0.005)
c.ps <- subset_taxa(c.ps, taxa_sums(c.ps) > 0.005)
d.ps <- subset_taxa(d.ps, taxa_sums(d.ps) > 0.005)
e.ps <- subset_taxa(e.ps, taxa_sums(e.ps) > 0.005)
f.ps <- subset_taxa(f.ps, taxa_sums(f.ps) > 0.005)
g.ps <- subset_taxa(g.ps, taxa_sums(g.ps) > 0.005)
h.ps <- subset_taxa(h.ps, taxa_sums(h.ps) > 0.005)

# Save the names of each individual data set as a character #
final_asv.names <- c()
for(i in letters[1:8]){
  final_asv.names <- c(final_asv.names, taxa_names(get(paste0(i, '.ps'))))
}

# Create a phyloseq object with the saved names based on prevalence #
final_ksp.ps <- subset_taxa(filt_ksp.ps, taxa_names(filt_ksp.ps) %in% final_asv.names)

# Make a data frame to be used to filter taxa based on bootstrap values #
final_ksp.tax <- as.data.frame(tax_table(final_ksp.ps))
for(i in 4:6){
  final_ksp.tax[,i] <- as.numeric(final_ksp.tax[,i])
}
# Remove taxa that have bootstrap values below 50 #
ksp.rem <- filter(final_ksp.tax, Family_Boot < 50)
final_ksp.ps <- subset_taxa(final_ksp.ps, !taxa_names(final_ksp.ps) %in% rownames(ksp.rem))
final_ksp.tax <- final_ksp.tax[taxa_names(final_ksp.ps),]

# Assign final taxonomy based on Bootstrap Values #
taxa_names(final_ksp.ps) <- paste0('ASV', seq(ntaxa(final_ksp.ps)))
for(i in 1:nrow(final_ksp.tax)){
  if(final_ksp.tax$Genus_Boot[i] >= 50){
    taxa_names(final_ksp.ps)[i] <- paste0(taxa_names(final_ksp.ps)[i], '(', final_ksp.tax$Genus[i], ')')
  } else if(final_ksp.tax$Family_Boot[i] >= 50){
    taxa_names(final_ksp.ps)[i] <- paste0(taxa_names(final_ksp.ps)[i], '(', final_ksp.tax$Family[i], ')')
    } else{
      taxa_names(final_ksp.ps)[i] <- paste0(taxa_names(final_ksp.ps)[i], '(NA)')
  }
}

# Remake the taxa table to only include MaarJAM databse hits #
rownames(final_ksp.tax) <- taxa_names(final_ksp.ps)
final_ksp.tax$ASV <- rownames(final_ksp.tax)
tax_table(final_ksp.ps) <- as.matrix(final_ksp.tax)[,c("Family", "Genus", "Species", "ASV", "Family_Boot", "Genus_Boot", "Species_Boot")]

# Save Treatment and Site as factors in the metadata table, as well as their pairwise combinations #
sample_data(final_ksp.ps)$Treats <- factor(sample_data(final_ksp.ps)$Treatment, levels = c("Control", "Low", "High", "MycoBloom"))
sample_data(final_ksp.ps)$Sites <- factor(sample_data(final_ksp.ps)$Site, levels = c("A", "B", "C", "D", "E", "F", "G", "H", "M"))
sample_data(final_ksp.ps)$Both <- substring(sample_names(final_ksp.ps), 1, 3)
for(i in 1:nsamples(final_ksp.ps)){
  if(substring(sample_data(final_ksp.ps)$Both[i],3,3) %in% c("1", "2", "3", "4", "5", "6")){
    sample_data(final_ksp.ps)$Both[i] <- substring(sample_data(final_ksp.ps)$Both[i],1,2)
  }
}
sample_data(final_ksp.ps)$Group <- factor(sample_data(final_ksp.ps)$Both, levels = unique(sort(sample_data(final_ksp.ps)$Both)))

# Since we are interested in seeing if the inoculant communities makes it into the incumbent community, we can make a separate phyloseq object that just has the MycoBloom Community #
myc.ps <- subset_samples(final_ksp.ps, Treatment == "MycoBloom")
myc.ps <- subset_taxa(myc.ps, taxa_sums(myc.ps) > 0)
decompose_ps(myc.ps, 'myc')
myc$fra <- arrange(myc$fra, desc(myc$fra$MycoBloom1))
decompose_ps(final_ksp.ps, 'final_ksp')

# Save the final and mycbloom phyloseq and decomposed phyloseq objects to abridged.RData #
resave(myc.ps, file = './abridged.RData')
resave(myc, file = './abridged.RData')
resave(final_ksp.ps, file = './abridged.RData')
resave(final_ksp, file = './abridged.RData')

#### Alpha Diversity Measurement and Visualization ####
# Create a data.frame of the alpha diversity values and their associated metadata values #
ksp.rich <- estimate_richness(final_ksp.ps) # automatically performs alpha diversity calculations #
ksp.rich <- cbind(final_ksp$met, ksp.rich)
ksp.rich <- arrange(ksp.rich, Sites, Treats)

# Calculate the Shannon Evenness of each Sample (Shannon diversity divided by the natural log of the observed species) #
for(i in 1:nrow(ksp.rich)){
  ksp.rich$ShaEvn[i] <- ksp.rich$Shannon[i] / log(ksp.rich$Chao1[i])
}

# Save and remove the values for the MycoBloom community #
myc.rich <- filter(ksp.rich, Treatment == "MycoBloom")
ksp.rich <- filter(ksp.rich, Treatment != "MycoBloom")
ksp.rich$Treats <- factor(ksp.rich$Treatment, levels = c("Control", "Low", "High"))

resave(ksp.rich, file = './abridged.RData')

# ANOVA for each kind of alpha diversity #
if(!requireNamespace("car")) install.packages("car")
library(car); packageVersion("car")

if(!requireNamespace("lmerTest")) install.packages("lmerTest")
library(lmerTest); packageVersion("lmerTest")

if(!requireNamespace("emmeans")) install.packages("emmeans")
library(emmeans); packageVersion("emmeans")

if(!requireNamespace("performance")) install.packages("performance")
library(performance); packageVersion("performance")

## Shannon Diversity ##
# Perform the ANOVA and normality and variance tests for a multiple linear regression #
ksp_sha.fit <- lm(Shannon ~ Treatment * Site, data = ksp.rich)
ksp_sha.sum <- Anova(ksp_sha.fit, type = 3) 
ksp_sha.sum
shapiro.test(residuals(ksp_sha.fit))
leveneTest(ksp_sha.fit)

# Perform ANOVA and diagnostics for the mixed linear model #
ksp_mix.sha <- lmer(Shannon ~ Treats + (1 | Sites), data = ksp.rich)
ksp_mix_sha.dia <- check_model(ksp_mix.sha, panel = TRUE)
ksp_mix_sha.dia
shapiro.test(resid(ksp_mix.sha))
leveneTest(resid(ksp_mix.sha) ~ ksp.rich$Treats)
ksp_mix_sha.sum <- summary(ksp_mix.sha)
ksp_mix_sha.sum
ksp_mix_sha.aov <- Anova(ksp_mix.sha, type = 3)
ksp_mix_sha.aov
ksp_sha.emn <- emmeans(ksp_mix.sha,pairwise~ Treats, adjust = "none")
ksp_sha.emn
resave(ksp_sha.sum, file = './abridged.RData')
resave(ksp_mix_sha.sum, file = './abridged.RData')
resave(ksp_mix_sha.aov, file = 'abridged.RData')
resave(ksp_sha.emn, file = './abridged.RData')

## Shannon Evenness ##
# Perform the ANOVA and normality and variance tests #
ksp_sev.fit <- lm(ShaEvn ~ Treats * Sites, data = ksp.rich)
ksp_sev.sum <- Anova(ksp_sev.fit, type = 3)
ksp_sev.sum
shapiro.test(residuals(ksp_sev.fit))
leveneTest(ksp_sev.fit)

# Perform ANOVA and diagnostics for the mixed linear model #
ksp_mix.sev <- lmer(ShaEvn ~ Treats + (1 | Sites), data = ksp.rich)
ksp_mix_sev.dia <- check_model(ksp_mix.sev, panel = TRUE)
ksp_mix_sev.dia
shapiro.test(resid(ksp_mix.sev))
leveneTest(resid(ksp_mix.sev) ~ ksp.rich$Treatment)
ksp_mix_sev.sum <- summary(ksp_mix.sev)
ksp_mix_sev.sum
ksp_mix_sev.aov <- Anova(ksp_mix.sev, type = 3)
ksp_mix_sev.aov
ksp_sev.emn <- emmeans(ksp_mix.sev,pairwise~ Treats, adjust = "none")
ksp_sev.emn
resave(ksp_sev.sum, file = './abridged.RData')
resave(ksp_mix_sev.sum, file = './abridged.RData')
resave(ksp_mix_sev.aov, file = 'abridged.RData')
resave(ksp_sev.emn, file = './abridged.RData')


# Perform the ANOVA and normality and variance tests #
ksp_cha.fit <- lm(Chao1 ~ Treats * Sites, ksp.rich) 
ksp_cha.sum <- Anova(ksp_cha.fit, type = 3)
ksp_cha.sum
shapiro.test(residuals(ksp_cha.fit))
leveneTest(ksp_cha.fit)

# Perform ANOVA and diagnostics for the mixed linear model #
ksp_mix.cha <- lmer(Chao1 ~ Treats + (1 | Sites), data = ksp.rich)
ksp_mix_cha.dia <- check_model(ksp_mix.cha, panel = TRUE)
shapiro.test(resid(ksp_mix.cha))
leveneTest(resid(ksp_mix.cha) ~ ksp.rich$Treatment)
ksp_mix_cha.sum <- summary(ksp_mix.cha)
ksp_mix_cha.sum
ksp_mix_cha.aov <- Anova(ksp_mix.cha, type = 3)
ksp_mix_cha.aov
ksp_cha.emn <- emmeans(ksp_mix.cha,pairwise~ Treats, adjust = "none")
ksp_cha.emn
resave(ksp_cha.sum, file = './abridged.RData')
resave(ksp_mix_cha.sum, file = './abridged.RData')
resave(ksp_mix_cha.aov, file = 'abridged.RData')
resave(ksp_cha.emn, file = './abridged.RData')

# Find the mean and standard deviation of each grouping #
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
ksp_rich.mnsd$Treats <- factor(ksp_rich.mnsd$Treatment, levels = c("Control", "Low", "High"))
ksp_rich.mnsd$Sites <- factor(ksp_rich.mnsd$Site, levels = c("A", "B", "C", "D", "E", "F", "G", "H"))
ksp_rich.mnsd <- arrange(ksp_rich.mnsd, Sites, Treats)

## Individual Site  Diversity ## 
# Subset the dataset by site # 
a.rich <- filter(ksp.rich, Site == 'A')
b.rich <- filter(ksp.rich, Site == 'B')
c.rich <- filter(ksp.rich, Site == 'C')
d.rich <- filter(ksp.rich, Site == 'D')
e.rich <- filter(ksp.rich, Site == 'E')
f.rich <- filter(ksp.rich, Site == 'F')
g.rich <- filter(ksp.rich, Site == 'G')
h.rich <- filter(ksp.rich, Site == 'H')

## Shannon Diversity ##
if(!requireNamespace("multcompView")) install.packages("multcompView")
library(multcompView); packageVersion('multcompView')

# Perform the ANOVA within each site across Treatments #
a.sha <- aov(Shannon~Treats, a.rich)
a_sha.sum <- summary(a.sha)
a_sha.sum
a_sha.hsd <- TukeyHSD(a.sha)
a_sha.hsd
a_sha.let <- multcompLetters4(a.sha, a_sha.hsd)
a_sha.let <- a_sha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(a_sha.sum, file = './abridged.RData')
resave(a_sha.hsd, file = './abridged.RData')

b.sha <- aov(Shannon~Treats, b.rich)
b.sha_sum <- summary(b.sha)
b.sha_sum
b_sha.hsd <- TukeyHSD(b.sha)
b_sha.hsd
b_sha.let <- multcompLetters4(b.sha, b_sha.hsd)
b_sha.let <- b_sha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(b_sha.sum, file = './abridged.RData')
resave(b_sha.hsd, file = './abridged.RData')

c.sha <- aov(Shannon~Treats, c.rich)
c_sha.sum <- summary(c.sha)
c_sha.sum
c_sha.hsd <- TukeyHSD(c.sha)
c_sha.hsd
c_sha.let <- multcompLetters4(c.sha, c_sha.hsd)
c_sha.let <- c_sha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(c_sha.sum, file = './abridged.RData')
resave(c_sha.hsd, file = './abridged.RData')

d.sha <- aov(Shannon~Treats, d.rich)
d_sha.sum <- summary(d.sha)
d_sha.sum
d_sha.hsd <- TukeyHSD(d.sha)
d_sha.hsd
d_sha.let <- multcompLetters4(d.sha, d_sha.hsd)
d_sha.let <- d_sha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(d_sha.sum, file = './abridged.RData')
resave(d_sha.hsd, file = './abridged.RData')

e.sha <- aov(Shannon~Treats, e.rich)
e_sha.sum <- summary(e.sha)
e_sha.sum
e_sha.hsd <- TukeyHSD(e.sha)
e_sha.hsd
e_sha.let <- multcompLetters4(e.sha, e_sha.hsd)
e_sha.let <- e_sha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(e_sha.sum, file = './abridged.RData')
resave(e_sha.hsd, file = './abridged.RData')

f.sha <- aov(Shannon~Treats, f.rich)
f_sha.sum <- summary(f.sha)
f_sha.sum
f_sha.hsd <- TukeyHSD(f.sha)
f_sha.hsd
f_sha.let <- multcompLetters4(f.sha, f_sha.hsd)
f_sha.let <- f_sha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(f_sha.sum, file = './abridged.RData')
resave(f_sha.hsd, file = './abridged.RData')

g.sha <- aov(Shannon~Treats, g.rich)
g_sha.sum <- summary(g.sha)
g_sha.sum
g_sha.hsd <- TukeyHSD(g.sha)
g_sha.hsd
g_sha.let <- multcompLetters4(g.sha, g_sha.hsd)
g_sha.let <- g_sha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(g_sha.sum, file = './abridged.RData')
resave(g_sha.hsd, file = './abridged.RData')

h.sha <- aov(Shannon~Treats, h.rich)
h_sha.sum <- summary(h.sha)
h_sha.sum
h_sha.hsd <- TukeyHSD(h.sha)
h_sha.hsd
h_sha.let <- multcompLetters4(h.sha, h_sha.hsd)
h_sha.let <- h_sha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(h_sha.sum, file = './abridged.RData')
resave(h_sha.hsd, file = './abridged.RData')

# Save the letters and create a data frame used for the plot #
sha.let <- c(a_sha.let,
             b_sha.let,
             c_sha.let,
             d_sha.let,
             e_sha.let,
             f_sha.let,
             g_sha.let,
             h_sha.let)
plot.rich <- cbind(ksp_rich.mnsd, sha.let)

# Plot the Results for Shannon Diversity #
sha.plot <- ggplot(plot.rich, aes(x = Sites, y = sha.mean, fill = Treats, color = Treats)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  theme_prism() +
  ylab('Shannon Diversity') +
  scale_fill_manual(values = c("white", "gray", "#4D4D4D", 'black'), name = "Inoculant Level") +
  scale_color_manual(values = c('black', 'black', 'black'), name = "Inoculant Level") +
  scale_y_continuous(limits = c(0,5), breaks = seq(0,5,1), expand = expansion(mult = c(0, 0.05))) +
  geom_text(aes(label = sha.let, y = sha.mean + sha.sd + 0.1), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 8) +
  geom_errorbar(aes(ymin = sha.mean - sha.sd, ymax = sha.mean + sha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  theme(legend.title = element_text(size = 20, family = "Liberation Sans", face = 'bold'),
        legend.text = element_text(size = 16, family = "Liberation Sans", face = 'bold'),
        legend.key.spacing.y = unit(3, 'mm'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, family = 'Liberation Sans', face = 'bold')) +
  labs(tag = 'C.')

resave(sha.plot, file = './abridged.RData')

# Shannon Diversity #
a.sim <- aov(ShaEvn~Treats, a.rich)
a_sim.sum <- summary(a.sim)
a_sim.sum
a_sim.hsd <- TukeyHSD(a.sim)
a_sim.hsd
a_sim.let <- multcompLetters4(a.sim, a_sim.hsd)
a_sim.let <- a_sim.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(a_sim.sum, file = './abridged.RData')
resave(a_sim.hsd, file = './abridged.RData')

b.sim <- aov(ShaEvn~Treats, b.rich)
b_sim.sum <- summary(b.sim)
b_sim.sum
b_sim.hsd <- TukeyHSD(b.sim)
b_sim.hsd
b_sim.let <- multcompLetters4(b.sim, b_sim.hsd)
b_sim.let <- b_sim.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(b_sim.sum, file = './abridged.RData')
resave(b_sim.hsd, file = './abridged.RData')

c.sim <- aov(ShaEvn~Treats, c.rich)
c_sim.sum <- summary(c.sim)
c_sim.sum
c_sim.hsd <- TukeyHSD(c.sim)
c_sim.hsd
c_sim.let <- multcompLetters4(c.sim, c_sim.hsd)
c_sim.let <- c_sim.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(c_sim.sum, file = './abridged.RData')
resave(c_sim.hsd, file = './abridged.RData')

d.sim <- aov(ShaEvn~Treats, d.rich)
d_sim.sum <- summary(d.sim)
d_sim.sum
d_sim.hsd <- TukeyHSD(d.sim)
d_sim.hsd 
d_sim.let <- multcompLetters4(d.sim, d_sim.hsd)
d_sim.let <- d_sim.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(d_sim.sum, file = './abridged.RData')
resave(d_sim.hsd, file = './abridged.RData')

e.sim <- aov(ShaEvn~Treats, e.rich)
e_sim.sum <- summary(e.sim)
e_sim.sum
e_sim.hsd <- TukeyHSD(e.sim)
e_sim.hsd
e_sim.let <- multcompLetters4(e.sim, e_sim.hsd)
e_sim.let <- e_sim.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(e_sim.sum, file = './abridged.RData')
resave(e_sim.hsd, file = './abridged.RData')

f.sim <- aov(ShaEvn~Treats, f.rich)
f_sim.sum <- summary(f.sim)
f_sim.sum
f_sim.hsd <- TukeyHSD(f.sim)
f_sim.hsd
f_sim.let <- multcompLetters4(f.sim, f_sim.hsd)
f_sim.let <- f_sim.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(f_sim.sum, file = './abridged.RData')
resave(f_sim.hsd, file = './abridged.RData')

g.sim <- aov(ShaEvn~Treats, g.rich)
g_sim.sum <- summary(g.sim)
g_sim.sum
g_sim.hsd <- TukeyHSD(g.sim)
g_sim.hsd
g_sim.let <- multcompLetters4(g.sim, g_sim.hsd)
g_sim.let <- g_sim.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(g_sim.sum, file = './abridged.RData')
resave(g_sim.hsd, file = './abridged.RData')

h.sim <- aov(ShaEvn~Treats, h.rich)
h_sim.sum <- summary(h.sim)
h_sim.sum
h_sim.hsd <- TukeyHSD(h.sim)
h_sim.hsd
h_sim.let <- multcompLetters4(h.sim, h_sim.hsd)
h_sim.let <- h_sim.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(h_sim.sum, file = './abridged.RData')
resave(h_sim.hsd, file = './abridged.RData')

sim.let <- c(a_sim.let,
             b_sim.let,
             c_sim.let,
             d_sim.let,
             e_sim.let,
             f_sim.let,
             g_sim.let,
             h_sim.let)
plot.rich <- cbind(plot.rich, sim.let)

evn.plot <- ggplot(plot.rich, aes(x = Sites, y = evn.mean, fill = Treats, color = Treats)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  theme_prism() +
  ylab("Shannon Evenness") +
  scale_fill_manual(values = c("white", "gray", "#4D4D4D", 'black'), name = 'Inoculant Level') +
  scale_color_manual(values = c('black', 'black', 'black'), name = 'Inoculant Level') +
  scale_y_continuous(limits = c(0,1.1), breaks = seq(0,1,0.2),expand = expansion(mult = c(0, 0.05))) +
  geom_text(aes(label = sim.let, y = evn.mean + evn.sd + 0.01), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 8) +
  geom_errorbar(aes(ymin = evn.mean - evn.sd, ymax = evn.mean + evn.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  theme(legend.title = element_text(size = 20, family = "Liberation Sans", face = 'bold'),
        legend.text = element_text(size = 16, family = "Liberation Sans", face = 'bold'),
        legend.key.spacing.y = unit(3, 'mm'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  labs(tag = 'B.')

resave(evn.plot, file = './abridged.RData')

# Otu Richness #
a.cha <- aov(Chao1~Treats, a.rich)
a_cha.sum <- summary(a.cha)
a_cha.sum
a_cha.hsd <- TukeyHSD(a.cha)
a_cha.hsd
a_cha.let <- multcompLetters4(a.cha, a_cha.hsd)
a_cha.let <- a_cha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(a_cha.sum, file = './abridged.RData')
resave(a_cha.hsd, file = './abridged.RData')

b.cha <- aov(Chao1~Treats, b.rich)
b_cha.sum <- summary(b.cha)
b_cha.sum
b_cha.hsd <- TukeyHSD(b.cha)
b_cha.hsd
b_cha.let <- multcompLetters4(b.cha, b_cha.hsd)
b_cha.let <- b_cha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(b_cha.sum, file = './abridged.RData')
resave(b_cha.hsd, file = './abridged.RData')

c.cha <- aov(Chao1~Treats, c.rich)
c_cha.sum <- summary(c.cha)
c_cha.sum
c_cha.hsd <- TukeyHSD(c.cha)
c_cha.hsd
c_cha.let <- multcompLetters4(c.cha, c_cha.hsd)
c_cha.let <- c_cha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(c_cha.sum, file = './abridged.RData')
resave(c_cha.hsd, file = './abridged.RData')

d.cha <- aov(Chao1~Treats, d.rich)
d_cha.sum <- summary(d.cha)
d_cha.sum
d_cha.hsd <- TukeyHSD(d.cha)
d_cha.hsd
d_cha.let <- multcompLetters4(d.cha, d_cha.hsd)
d_cha.let <- d_cha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(d_cha.sum, file = './abridged.RData')
resave(d_cha.hsd, file = './abridged.RData')

e.cha <- aov(Chao1~Treats, e.rich)
e_cha.sum <- summary(e.cha)
e_cha.sum
e_cha.hsd <- TukeyHSD(e.cha)
e_cha.hsd
e_cha.let <- multcompLetters4(e.cha, e_cha.hsd)
e_cha.let <- e_cha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(e_cha.sum, file = './abridged.RData')
resave(e_cha.hsd, file = './abridged.RData')

f.cha <- aov(Chao1~Treats, f.rich)
f_cha.sum <- summary(f.cha)
f_cha.sum
f_cha.hsd <- TukeyHSD(f.cha)
f_cha.hsd
f_cha.let <- multcompLetters4(f.cha, f_cha.hsd)
f_cha.let <- f_cha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(f_cha.sum, file = './abridged.RData')
resave(f_cha.hsd, file = './abridged.RData')

g.cha <- aov(Chao1~Treats, g.rich)
g_cha.sum <- summary(g.cha)
g_cha.sum
g_cha.hsd <- TukeyHSD(g.cha)
g_cha.hsd
g_cha.let <- multcompLetters4(g.cha, g_cha.hsd)
g_cha.let <- g_cha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(g_cha.sum, file = './abridged.RData')
resave(g_cha.hsd, file = './abridged.RData')

h.cha <- aov(Chao1~Treats, h.rich)
h_cha.sum <- summary(h.cha)
h_cha.sum
h_cha.hsd <- TukeyHSD(h.cha)
h_cha.hsd
h_cha.let <- multcompLetters4(h.cha, h_cha.hsd)
h_cha.let <- h_cha.let$Treats$Letters[levels(ksp.rich$Treats)]
resave(h_cha.sum, file = './abridged.RData')
resave(h_cha.hsd, file = './abridged.RData')

cha.let <- c(a_cha.let,
             b_cha.let,
             c_cha.let,
             d_cha.let,
             e_cha.let,
             f_cha.let,
             g_cha.let,
             h_cha.let)
plot.rich <- cbind(plot.rich, cha.let)

resave(plot.rich, file = './abridged.RData')

cha.plot <- ggplot(plot.rich, aes(x = Sites, y = cha.mean, fill = Treats, color = Treats)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  theme_prism() +
  ylab("ASV Richness") +
  scale_fill_manual(values = c("white", "gray", "#4D4D4D", 'black'), name = "Inoculant Level") +
  scale_color_manual(values = c('black', 'black', 'black'), name = "Inoculant Level") +
  scale_y_continuous(limits = c(0,220), breaks = seq(0, 200, 40), expand = expansion(mult = c(0,0.05))) +
  geom_text(aes(label = cha.let, y = cha.mean + cha.sd + 1), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 8) +
  geom_errorbar(aes(ymin = cha.mean - cha.sd, ymax = cha.mean + cha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  theme(legend.title = element_text(size = 20, family = "Liberation Sans", face = 'bold'),
        legend.text = element_text(size = 16, family = "Liberation Sans", face = 'bold'),
        legend.key.spacing.y = unit(3, 'mm'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(tag = 'A.')

resave(cha.plot, file = "./abridged.RData")

if(!requireNamespace('patchwork')) install.packages('patchwork')
library(patchwork); packageVersion('patchwork')

alpha.plot <- (cha.plot) /
(evn.plot) /
(sha.plot) +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom',
        legend.key.spacing.x = unit(3, 'cm'),
        legend.text = element_text(size = 18, color = 'black', face = 'bold', family = "Liberation Sans"),
        legend.title = element_blank())

resave(alpha.plot, file = "./abridged.RData")

#### Beta Diversity ####
# For the entire dataset, calculate the weighted unifrac distances from the total sum scaled (TSS) data #
if(!requireNamespace("vegan")) install.package("vegan")
library(vegan); packageVersion('vegan')
ksp_prop.ps <- transform_sample_counts(final_ksp.ps, function(otu) otu/sum(otu))
set.seed(248)
ksp.bdist <- phyloseq::distance(ksp_prop.ps, method = "unifrac")
ksp_unif.pcoa <- phyloseq::ordinate(ksp_prop.ps, method = "PCoA", distance = ksp.bdist)
ord.nmds.wuni <- metaMDS(ksp.bdist,
                         k = 2,
                         try = 100,
                         trymax = 1000,
                         previous.best = ksp_unif.pcoa$vectors)

# Perform PermANOVA for the entire dataset
decompose_ps(final_ksp.ps, 'final_ksp')

ksp_by.perm <- adonis2(ksp.bdist~Treats*Sites, data = final_ksp$met, by = "terms")
ksp_by.perm
ksp.perm <- adonis2(ksp.bdist~Treats*Sites, data = final_ksp$met)
ksp.perm
ksp.drda <- dbrda(ksp.bdist ~ Treats * Sites, data = final_ksp$met)
ksp_drda.sum <- summary(ksp.drda)
ksp_drda.sum

resave(ksp_by.perm, file = './abridged.RData')
resave(ksp.perm, file = './abridged.RData')

beta.plot<- plot_ordination(ksp_prop.ps, ord.nmds.wuni, shape = "Treatment", color="Site", title="NMDS") +
  theme_prism() +
  geom_point(size = 6) 

resave(beta.plot, file = './abridged.RData')

#### Per Sample Analysis ####
if(!requireNamespace("devtools")) install.packages('devtools')
library(devtools); packageVersion('devtools')

if(!requireNamespace("pairwiseAdonis")) devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis); packageVersion("pairwiseAdonis")

## A ##
# Create a phyloseq object with just the sites of interest and calculate weighted unifrac distance #
a.ps <- subset_samples(final_ksp.ps, Site == 'A')
a.ps <- subset_taxa(a.ps, taxa_sums(a.ps) > 0)
a.met <- as(sample_data(a.ps), 'data.frame')
a_prop.ps <- transform_sample_counts(a.ps, function(otu) otu/sum(otu))
a.dist <- as.dist(as.matrix(ksp.bdist)[sample_names(a.ps), sample_names(a.ps)])
a_ord.nmds.wuni <- ordinate(a_prop.ps, method="NMDS", distance= a.dist)

# Perform PermANOVA #
a.perm <- adonis2(a.bdist~Treatment, data = a.met)
a.perm
resave(a.perm, file = 'abridged.RData')

# Plot the results #
a_beta.plot <- plot_ordination(a_prop.ps, a_ord.nmds.wuni, color="Treatment", title="A Samples NMDS") +
  theme_prism() +
  geom_point(size = 6)

resave(a_beta.plot, file = 'abridged.RData')

## B ##
# Create a phyloseq object with just the sites of interest and calculate weighted unifrac distance #
b.ps <- subset_samples(final_ksp.ps, Site == 'B')
b.ps <- subset_taxa(b.ps, taxa_sums(b.ps) > 0)
b.met <- as(sample_data(b.ps), 'data.frame')
b_prop.ps <- transform_sample_counts(b.ps, function(otu) otu/sum(otu))
b.dist <- as.dist(as.matrix(ksp.bdist)[sample_names(b.ps), sample_names(b.ps)])
b_ord.nmds.wuni <- ordinate(b_prop.ps, method="NMDS", distance=b.dist)

# Perform PermANOVA #
b.perm <- adonis2(b.bdist~Treatment, data = b.met)
b.perm
resave(b.perm, file = './abridged.RData')

# Plot the results #
b_beta.plot <- plot_ordination(b_prop.ps, b_ord.nmds.wuni, color="Treatment", title="B Samples NMDS") +
  theme_prism() +
  geom_point(size = 6)
resave(b_beta.plot, file = './abridged.RData')

## C ##
# Create a phyloseq object with just the sites of interest and calculate weighted unifrac distance #
c.ps <- subset_samples(final_ksp.ps, Site == 'C')
c.ps <- subset_taxa(c.ps, taxa_sums(c.ps) > 0)
c.met <- as(sample_data(c.ps), 'data.frame')
c_prop.ps <- transform_sample_counts(c.ps, function(otu) otu/sum(otu))
c.dist <- as.dist(as.matrix(ksp.bdist)[sample_names(c.ps), sample_names(c.ps)])
c_ord.nmds.wuni <- ordinate(c_prop.ps, method="NMDS", distance=c.dist)

# Perform PermANOVA #
c.perm <- adonis2(c.bdist~Treatment, data = c.met)
c.perm
resave(c.perm, file = './abridged.RData')

# Plot the results #
c_beta.plot <- plot_ordination(c_prop.ps, c_ord.nmds.wuni, color="Treatment", title="C Samples NMDS") +
  theme_prism() +
  geom_point(size = 6)
resave(c_beta.plot, './abridged.RData')

## D ##
# Create a phyloseq object with just the sites of interest and calculate weighted unifrac distance #
d.ps <- subset_samples(final_ksp.ps, Site == 'D')
d.ps <- subset_taxa(d.ps, taxa_sums(d.ps) > 0)
d.met <- as(sample_data(d.ps), 'data.frame')
d_prop.ps <- transform_sample_counts(d.ps, function(otu) otu/sum(otu))
d.dist <- as.dist(as.matrix(ksp.bdist)[sample_names(d.ps), sample_names(d.ps)])
d_ord.nmds.wuni <- ordinate(d_prop.ps, method="NMDS", distance=d.dist)

# Perform PermANOVA #
d.perm <- adonis2(d.bdist~Treatment, data = d.met)
d.perm
resave(d.perm, file = './abridged.RData')

# Plot the results #
d_beta.plot <- plot_ordination(d_prop.ps, d_ord.nmds.wuni, color="Treatment", title="D Samples NMDS") +
  theme_prism() +
  geom_point(size = 6)
resave(d_beta.plot, file = './abridged.RData')

## E ##
# Create a phyloseq object with just the sites of interest and calculate weighted unifrac distance #
e.ps <- subset_samples(final_ksp.ps, Site == 'E')
e.ps <- subset_taxa(e.ps, taxa_sums(e.ps) > 0)
e.met <- as(sample_data(e.ps), 'data.frame')
e_prop.ps <- transform_sample_counts(e.ps, function(otu) otu/sum(otu))
e.dist <- as.dist(as.matrix(ksp.bdist)[sample_names(e.ps), sample_names(e.ps)])
e_ord.nmds.wuni <- ordinate(e_prop.ps, method="NMDS", distance=e.dist)

# Perform PermANOVA #
e.perm <- adonis2(e.bdist~Treatment, data = e.met)
e.perm
resave(e.perm, file = 'abridged.RData')

# Plot the results #
e_beta.plot <- plot_ordination(e_prop.ps, e_ord.nmds.wuni, color="Treatment", title="E Samples NMDS") +
  theme_prism() +
  geom_point(size = 6)
resave(e_beta.plot, file = './abridged.RData')

## F ##
# Create a phyloseq object with just the sites of interest and calculate weighted unifrac distance #
f.ps <- subset_samples(final_ksp.ps, Site == 'F')
f.ps <- subset_taxa(f.ps, taxa_sums(f.ps) > 0)
f.met <- as(sample_data(f.ps), 'data.frame')
f_prop.ps <- transform_sample_counts(f.ps, function(otu) otu/sum(otu))
f.dist <- as.dist(as.matrix(ksp.bdist)[sample_names(f.ps), sample_names(f.ps)])
f_ord.nmds.wuni <- ordinate(f_prop.ps, method="NMDS", distance=f.dist)

# Perform PermANOVA #
f.perm <- adonis2(f.bdist~Treatment, data = f.met)
f.perm
resave(f.perm, file = './abridged.RData')

# Plot the results #
f_beta.plot <- plot_ordination(f_prop.ps, f_ord.nmds.wuni, color="Treatment", title="F Samples NMDS") +
  theme_prism() +
  geom_point(size = 6)
resave(f_beta.plot, file = 'abridged.RData', file = './abrdiged.RData')

## G ##
# Create a phyloseq object with just the sites of interest and calculate weighted unifrac distance #
g.ps <- subset_samples(final_ksp.ps, Site == 'G')
g.ps <- subset_taxa(g.ps, taxa_sums(g.ps) > 0)
g.met <- as(sample_data(g.ps), 'data.frame')
g_prop.ps <- transform_sample_counts(g.ps, function(otu) otu/sum(otu))
g.dist <- as.dist(as.matrix(ksp.bdist)[sample_names(g.ps), sample_names(g.ps)])
g_ord.nmds.wuni <- ordinate(g_prop.ps, method="NMDS", distance=g.dist)

# Perform PermANOVA #
g.perm <- adonis2(g.bdist~Treatment, data = g.met)
g.perm
resave(g.perm, file = './abridged.RData')

# Plot the results #
g_beta.plot <- plot_ordination(g_prop.ps, g_ord.nmds.wuni, color="Treatment", title="G Samples NMDS") +
  theme_prism() +
  geom_point(size = 6)
resave(g_beta.plot, file = './abridged.RData')

## H ##
# Create a phyloseq object with just the sites of interest and calculate weighted unifrac distance #
h.ps <- subset_samples(final_ksp.ps, Site == 'H')
h.ps <- subset_taxa(h.ps, taxa_sums(h.ps) > 0)
h.met <- as(sample_data(h.ps), 'data.frame')
h_prop.ps <- transform_sample_counts(h.ps, function(otu) otu/sum(otu))
h.dist <- as.dist(as.matrix(ksp.bdist)[sample_names(h.ps), sample_names(h.ps)])
h_ord.nmds.wuni <- ordinate(h_prop.ps, method="NMDS", distance=h.dist)

# Perform PermANOVA #
h.perm <- adonis2(h.bdist~Treatment, data = h.met)
h.perm
resave(h.perm, file = './abridged.RData')

# Plot the results #
h_beta.plot <- plot_ordination(h_prop.ps, h_ord.nmds.wuni, color="Treatment", title="H Samples NMDS") +
  theme_prism() +
  geom_point(size = 6)
resave(h_beta.plot, file = './abridged.RData')

#### Stacked Histograms ####
# First we start by making a color pallette for each unique ASV #
if(!requireNamespace("Polychrome")) install.packages("Polychrome")
library(Polychrome); packageVersion("Polychrome")
ksp.colr <- createPalette(ntaxa(final_ksp.ps),  c("#ff0000", "#00ff00", "#0000ff"))
ksp.colr <- as.data.frame(ksp.colr)
rownames(ksp.colr) <- taxa_names(final_ksp.ps)

# Add a gray color for "Other"
ksp.colr[ntaxa(final_ksp.ps) + 1,] <- "#D4D4D4" 
rownames(ksp.colr)[ntaxa(final_ksp.ps) + 1] <- "Other" 

## MycoBloom Sample ##
# Save a phyloseq Object that contains only the samples of interest #
myc.ps <- subset_samples(final_ksp.ps, Treatment == "MycoBloom")
myc.ps <- subset_taxa(myc.ps, taxa_sums(myc.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_myc.name <- names(sort(taxa_sums(myc.ps), decreasing = TRUE))[1:9]
hg_myc.name <- c(hg_myc.name, "Other")
hg_myc.colr <- ksp.colr[hg_myc.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_myc.ps <- merge_taxa(myc.ps, taxa_names(myc.ps)[!taxa_names(myc.ps) %in% hg_myc.name])
taxa_names(hg_myc.ps)[5] <- "Other"
tax_table(hg_myc.ps)[5, "ASV"] <- "Other"

# Save the phyloseq object as a data.frame and make factors that guide the plot what to plot #
hg_myc.df <- psmelt(hg_myc.ps)
hg_myc.df$ASVs <- factor(hg_myc.df$ASV, levels = hg_myc.name)
hg_myc.df$Soil <- factor(hg_myc.df$Treatment, levels = c("Control", "Low", "High", "MycoBloom"))

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
resave(hg_myc.plot, file = './abridged.RData')

## Site A Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_a.ps <- subset_samples(final_ksp.ps, Site == "A" | Treatment == "MycoBloom")
hg_a.ps <- subset_taxa(hg_a.ps, taxa_sums(hg_a.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_a.name <- names(sort(taxa_sums(hg_a.ps), decreasing = TRUE))[1:9]
hg_a.name <- c(hg_a.name, "Other")
hg_a.colr <- ksp.colr[hg_a.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_a.ps <- merge_taxa(hg_a.ps, taxa_names(hg_a.ps)[!taxa_names(hg_a.ps) %in% hg_a.name])
taxa_names(hg_a.ps)[6] <- "Other"
tax_table(hg_a.ps)[6, "ASV"] <- "Other"

# Save the phyloseq object as a data.frame and make factors that guide the plot what to plot #
hg_a.df <- psmelt(hg_a.ps)
hg_a.df$ASVs <- factor(hg_a.df$ASV, levels = hg_a.name)
hg_a.df$Soil <- factor(hg_a.df$Treatment, levels = c("Control", "Low", "High", "MycoBloom"))

# Plot the histogram #
hg_a.plot <- ggplot(hg_a.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(name = "Fungal ASV", values = hg_a.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "Site A")) +
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

resave(hg_a.plot, file = './abridged.RData')

## Site B Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_b.ps <- subset_samples(final_ksp.ps, Site == "B" | Treatment == "MycoBloom")
hg_b.ps <- subset_taxa(hg_b.ps, taxa_sums(hg_b.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_b.name <- names(sort(taxa_sums(hg_b.ps), decreasing = TRUE))[1:9]
hg_b.name <- c(hg_b.name, "Other")
hg_b.colr <- ksp.colr[hg_b.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_b.ps <- merge_taxa(hg_b.ps, taxa_names(hg_b.ps)[!taxa_names(hg_b.ps) %in% hg_b.name])
taxa_names(hg_b.ps)[8] <- "Other"
tax_table(hg_b.ps)[8, "ASV"] <- "Other"

# Save the phyloseq object as a data.frame and make factors that guide the plot what to plot #
hg_b.df <- psmelt(hg_b.ps)
hg_b.df$ASVs <- factor(hg_b.df$ASV, levels = hg_b.name)
hg_b.df$Soil <- factor(hg_b.df$Treatment, levels = c("Control", "Low", "High", "MycoBloom"))

# Plot the histogram #
hg_b.plot <- ggplot(hg_b.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(name = "Fungal ASV", values = hg_b.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "Site B")) +
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
resave(hg_b.plot, file = './abridged.RData')

## Site C Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_c.ps <- subset_samples(final_ksp.ps, Site == "C" | Treatment == "MycoBloom")
hg_c.ps <- subset_taxa(hg_c.ps, taxa_sums(hg_c.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_c.name <- names(sort(taxa_sums(hg_c.ps), decreasing = TRUE))[1:9]
hg_c.name <- c(hg_c.name, "Other")
hg_c.colr <- ksp.colr[hg_c.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_c.ps <- merge_taxa(hg_c.ps, taxa_names(hg_c.ps)[!taxa_names(hg_c.ps) %in% hg_c.name])
taxa_names(hg_c.ps)[9] <- "Other"
tax_table(hg_c.ps)[9, "ASV"] <- "Other"

# Save the phyloseq object as a data.frame and make factors that guide the plot what to plot #
hg_c.df <- psmelt(hg_c.ps)
hg_c.df$ASVs <- factor(hg_c.df$ASV, levels = hg_c.name)
hg_c.df$Soil <- factor(hg_c.df$Treatment, levels = c("Control", "Low", "High", "MycoBloom"))

# Plot the histogram #
hg_c.plot <- ggplot(hg_c.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(name = "Fungal ASV", values = hg_c.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "Site C")) +
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
resave(hg_c.plot, file = './abridged.RData')

## Site D Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_d.ps <- subset_samples(final_ksp.ps, Site == "D" | Treatment == "MycoBloom")
hg_d.ps <- subset_taxa(hg_d.ps, taxa_sums(hg_d.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_d.name <- names(sort(taxa_sums(hg_d.ps), decreasing = TRUE))[1:9]
hg_d.name <- c(hg_d.name, "Other")
hg_d.colr <- ksp.colr[hg_d.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_d.ps <- merge_taxa(hg_d.ps, taxa_names(hg_d.ps)[!taxa_names(hg_d.ps) %in% hg_d.name])
taxa_names(hg_d.ps)[7] <- "Other"
tax_table(hg_d.ps)[7, "ASV"] <- "Other"

# Save the phyloseq object as a data.frame and make factors that guide the plot what to plot #
hg_d.df <- psmelt(hg_d.ps)
hg_d.df$ASVs <- factor(hg_d.df$ASV, levels = hg_d.name)
hg_d.df$Soil <- factor(hg_d.df$Treatment, levels = c("Control", "Low", "High", "MycoBloom"))

# Plot the histogram #
hg_d.plot <- ggplot(hg_d.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(name = "Fungal ASV", values = hg_d.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "Site D")) +
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
resave(hg_d.plot, file = './abridged.RData')

## Site E Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_e.ps <- subset_samples(final_ksp.ps, Site == "E" | Treatment == "MycoBloom")
hg_e.ps <- subset_taxa(hg_e.ps, taxa_sums(hg_e.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_e.name <- names(sort(taxa_sums(hg_e.ps), decreasing = TRUE))[1:9]
hg_e.name <- c(hg_e.name, "Other")
hg_e.colr <- ksp.colr[hg_e.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_e.ps <- merge_taxa(hg_e.ps, taxa_names(hg_e.ps)[!taxa_names(hg_e.ps) %in% hg_e.name])
taxa_names(hg_e.ps)[9] <- "Other"
tax_table(hg_e.ps)[9, "ASV"] <- "Other"

# Save the phyloseq object as a data.frame and make factors that guide the plot what to plot #
hg_e.df <- psmelt(hg_e.ps)
hg_e.df$ASVs <- factor(hg_e.df$ASV, levels = hg_e.name)
hg_e.df$Soil <- factor(hg_e.df$Treatment, levels = c("Control", "Low", "High", "MycoBloom"))

# Plot the histogram #
hg_e.plot <- ggplot(hg_e.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(name = "Fungal ASV", values = hg_e.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "Site E")) +
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
resave(hg_e.plot, file = './abridged.RData')

## Site F Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_f.ps <- subset_samples(final_ksp.ps, Site == "F" | Treatment == "MycoBloom")
hg_f.ps <- subset_taxa(hg_f.ps, taxa_sums(hg_f.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_f.name <- names(sort(taxa_sums(hg_f.ps), decreasing = TRUE))[1:9]
hg_f.name <- c(hg_f.name, "Other")
hg_f.colr <- ksp.colr[hg_f.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_f.ps <- merge_taxa(hg_f.ps, taxa_names(hg_f.ps)[!taxa_names(hg_f.ps) %in% hg_f.name])
taxa_names(hg_f.ps)[4] <- "Other"
tax_table(hg_f.ps)[4, "ASV"] <- "Other"

# Save the phyloseq object as a data.frame and make factors that guide the plot what to plot #
hg_f.df <- psmelt(hg_f.ps)
hg_f.df$ASVs <- factor(hg_f.df$ASV, levels = hg_f.name)
hg_f.df$Soil <- factor(hg_f.df$Treatment, levels = c("Control", "Low", "High", "MycoBloom"))

# Plot the histogram #
hg_f.plot <- ggplot(hg_f.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(name = "Fungal ASV", values = hg_f.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "Site F")) +
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
resave(hg_f.plot, file = './abridged.RData')

## Site G Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_g.ps <- subset_samples(final_ksp.ps, Site == "G" | Treatment == "MycoBloom")
hg_g.ps <- subset_taxa(hg_g.ps, taxa_sums(hg_g.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_g.name <- names(sort(taxa_sums(hg_g.ps), decreasing = TRUE))[1:9]
hg_g.name <- c(hg_g.name, "Other")
hg_g.colr <- ksp.colr[hg_g.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_g.ps <- merge_taxa(hg_g.ps, taxa_names(hg_g.ps)[!taxa_names(hg_g.ps) %in% hg_g.name])
taxa_names(hg_g.ps)[7] <- "Other"
tax_table(hg_g.ps)[7, "ASV"] <- "Other"

# Save the phyloseq object as a data.frame and make factors that guide the plot what to plot #
hg_g.df <- psmelt(hg_g.ps)
hg_g.df$ASVs <- factor(hg_g.df$ASV, levels = hg_g.name)
hg_g.df$Soil <- factor(hg_g.df$Treatment, levels = c("Control", "Low", "High", "MycoBloom"))

# Plot the histogram #
hg_g.plot <- ggplot(hg_g.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(name = "Fungal ASV", values = hg_g.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "Site G")) +
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
resave(hg_g.plot, file = './abridged.RData')

## Site H Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_h.ps <- subset_samples(final_ksp.ps, Site == "H" | Treatment == "MycoBloom")
hg_h.ps <- subset_taxa(hg_h.ps, taxa_sums(hg_h.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_h.name <- names(sort(taxa_sums(hg_h.ps), decreasing = TRUE))[1:9]
hg_h.name <- c(hg_h.name, "Other")
hg_h.colr <- ksp.colr[hg_h.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_h.ps <- merge_taxa(hg_h.ps, taxa_names(hg_h.ps)[!taxa_names(hg_h.ps) %in% hg_h.name])
taxa_names(hg_h.ps)[10] <- "Other"
tax_table(hg_h.ps)[10, "ASV"] <- "Other"

# Save the phyloseq object as a data.frame and make factors that guide the plot what to plot #
hg_h.df <- psmelt(hg_h.ps)
hg_h.df$ASVs <- factor(hg_h.df$ASV, levels = hg_h.name)
hg_h.df$Soil <- factor(hg_h.df$Treatment, levels = c("Control", "Low", "High", "MycoBloom"))

# Plot the histogram #
hg_h.plot <- ggplot(hg_h.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(name = "Fungal ASV", values = hg_h.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "Site H")) +
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
resave(hg_h.plot, file = './abridged.RData')

## All Samples Across Sites Pooled ##
hg_ksp.colr <- ksp.colr[hg_myc.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_ksp.ps <- merge_taxa(final_ksp.ps, taxa_names(final_ksp.ps)[!taxa_names(final_ksp.ps) %in% hg_myc.name])
taxa_names(hg_ksp.ps)[1] <- "Other"
tax_table(hg_ksp.ps)[1, "ASV"] <- "Other"

# Save the phyloseq object as a data.frame and make factors that guide the plot what to plot #
hg_ksp.df <- psmelt(hg_ksp.ps)
hg_ksp.df$ASVs <- factor(hg_ksp.df$ASV, levels = hg_myc.name)
hg_ksp.df$Soil <- factor(hg_ksp.df$Treatment, levels = c("Control", "Low", "High", "MycoBloom"))

# Plot the histogram #
hg_ksp.plot <- ggplot(hg_ksp.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(name = "Fungal ASV", values = hg_ksp.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "Across All Sites")) +
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
resave(hg_ksp.plot, file = './abridged.RData')

#### Set Plot for overlapping taxa with MycoBloom and Treatments ####
# Save the names of the taxa found in the mycobloom sample #
myc.name <- taxa_names(myc.ps)

# Subset the whole phyloseq object by treatments and filter to only contain taxa found in the MycoBloom #
con.ps <- subset_samples(final_ksp.ps, Treatment == "Control")
con.ps <- subset_taxa(con.ps, taxa_names(con.ps) %in% myc.name)
con.ps <- subset_taxa(con.ps, taxa_sums(con.ps) > 0)
con.name <- taxa_names(con.ps)

hig.ps <- subset_samples(final_ksp.ps, Treatment == "High")
hig.ps <- subset_taxa(hig.ps, taxa_names(hig.ps) %in% myc.name)
hig.ps <- subset_taxa(hig.ps, taxa_sums(hig.ps) > 0)
hig.name <- taxa_names(hig.ps)

low.ps <- subset_samples(final_ksp.ps, Treatment == "Low")
low.ps <- subset_taxa(low.ps, taxa_names(low.ps) %in% myc.name)
low.ps <- subset_taxa(low.ps, taxa_sums(low.ps) > 0)
low.name <- taxa_names(low.ps)

# Construct a data.frame that contains TRUE and FALSE values for taxa presence #
if(!requireNamespace("ComplexUpset")) install.packages("ComplexUpset")
library(ComplexUpset); packageVersion("ComplexUpset")
all.name <- list(MycoBloom = myc.name, Control = con.name, High = hig.name, Low = low.name)
all.asv <- data.frame(ASV = myc.name)
for (group in names(all.name)) {
  all.asv[[group]] <- all.asv$ASV %in% all.name[[group]]
}

myc_nc.name <- filter(all.asv, High == "TRUE" | Low == "TRUE")
myc_nc.name <- filter(myc_nc.name, Control == "FALSE")
resave(myc_nc.name, file = 'abridged.RData')

myc_myc.name <- filter(all.asv, Control == "FALSE" & High == "FALSE" & Low == "FALSE")
resave(myc_myc.name, file = './abridged.RData')

# Create the set plot #
myc.set <- upset(
  all.asv,
  intersect = c("MycoBloom", "Control", "High", "Low"),
  base_annotations = list(
    'Intersection Size' = intersection_size(width = 0.9)),
  name = "Shared Inoculant ASVs",
  queries = list(
    upset_query(intersect = c("MycoBloom", "High", "Low"), color = "black", fill = "blue"),
    upset_query(intersect = c("MycoBloom", 'High'), color = "black", fill = "blue")),
  )
resave(myc.set, file = './abridged.RData')

## Site A Set Plot
a_con.ps <- subset_samples(a.ps, Treatment == "Control")
a_con.ps <- subset_taxa(a_con.ps, taxa_names(a_con.ps) %in% myc.name)
a_con.ps <- subset_taxa(a_con.ps, taxa_sums(a_con.ps) > 0)
a_con.name <- taxa_names(a_con.ps)

a_hig.ps <- subset_samples(a.ps, Treatment == "High")
a_hig.ps <- subset_taxa(a_hig.ps, taxa_names(a_hig.ps) %in% myc.name)
a_hig.ps <- subset_taxa(a_hig.ps, taxa_sums(a_hig.ps) > 0)
a_hig.name <- taxa_names(a_hig.ps)

a_low.ps <- subset_samples(a.ps, Treatment == "Low")
a_low.ps <- subset_taxa(a_low.ps, taxa_names(a_low.ps) %in% myc.name)
a_low.ps <- subset_taxa(a_low.ps, taxa_sums(a_low.ps) > 0)
a_low.name <- taxa_names(a_low.ps)

# Construct a data.frame that contains TRUE and FALSE values for taxa presence #
a.name <- list(MycoBloom = myc.name, Control = a_con.name, High = a_hig.name, Low = a_low.name)
a.asv <- data.frame(ASV = myc.name)
for (group in names(a.name)) {
  a.asv[[group]] <- a.asv$ASV %in% a.name[[group]]}

# Save the names of the taxa found in just the Mycobloom and ones found in everything but the control #
a_nc.name <- filter(a.asv, High == "TRUE" | Low == "TRUE")
a_nc.name <- filter(a_nc.name, Control == "FALSE")
resave(a_nc.name, file = 'abridged.RData')

a_myc.name <- filter(a.asv, Control == "FALSE" & High == "FALSE" & Low == "FALSE")
resave(a_myc.name, file = './abridged.RData')

# Create the set plot #
a.set <- upset(
  a.asv,
  intersect = c("MycoBloom", "Control", "High", "Low"),
  base_annotations = list(
    'Intersection size' = intersection_size(width = 0.9)
  ),
  queries = list(
    upset_query(intersect = "MycoBloom", color = "black", fill = "orange"),
    upset_query(intersect = c("MycoBloom", "High", "Low"), color = "black", fill = "blue"),
    upset_query(intersect = c("MycoBloom", "Low"), color = "black", fill = "blue"))
  )
resave(a.set, file = './abridged.RData')

## Site B Set Plot
b_con.ps <- subset_samples(b.ps, Treatment == "Control")
b_con.ps <- subset_taxa(b_con.ps, taxa_names(b_con.ps) %in% myc.name)
b_con.ps <- subset_taxa(b_con.ps, taxa_sums(b_con.ps) > 0)
b_con.name <- taxa_names(b_con.ps)

b_hig.ps <- subset_samples(b.ps, Treatment == "High")
b_hig.ps <- subset_taxa(b_hig.ps, taxa_names(b_hig.ps) %in% myc.name)
b_hig.ps <- subset_taxa(b_hig.ps, taxa_sums(b_hig.ps) > 0)
b_hig.name <- taxa_names(b_hig.ps)

b_low.ps <- subset_samples(b.ps, Treatment == "Low")
b_low.ps <- subset_taxa(b_low.ps, taxa_names(b_low.ps) %in% myc.name)
b_low.ps <- subset_taxa(b_low.ps, taxa_sums(b_low.ps) > 0)
b_low.name <- taxa_names(b_low.ps)

# Construct a data.frame that contains TRUE and FALSE values for taxa presence #
b.name <- list(MycoBloom = myc.name, Control = b_con.name, High = b_hig.name, Low = b_low.name)
b.asv <- data.frame(ASV = myc.name)
for (group in names(b.name)) {
  b.asv[[group]] <- b.asv$ASV %in% b.name[[group]]}

# Save the names of the taxa found in just the Mycobloom and ones found in everything but the control #
b_nc.name <- filter(b.asv, High == "TRUE" | Low == "TRUE")
b_nc.name <- filter(b_nc.name, Control == "FALSE")
resave(b_nc.name, file = 'abridged.RData')

b_myc.name <- filter(b.asv, Control == "FALSE" & High == "FALSE" & Low == "FALSE")
resave(b_myc.name, file = './abridged.RData')

# Create the set plot #
b.set <- upset(
  b.asv,
  intersect = c("MycoBloom", "Control", "High", "Low"),
  base_annotations = list(
    'Intersection size' = intersection_size(width = 0.9)
  ),
  queries = list(
    upset_query(intersect = "MycoBloom", color = "black", fill = "orange"),
    upset_query(intersect = c("MycoBloom", "Low"), color = "black", fill = "blue"),
    upset_query(intersect = c("MycoBloom", "High", "Low"), color = "black", fill = "blue"))
  )
resave(b.set, file = './abridged.RData')

## Site C Set Plot
c_con.ps <- subset_samples(c.ps, Treatment == "Control")
c_con.ps <- subset_taxa(c_con.ps, taxa_names(c_con.ps) %in% myc.name)
c_con.ps <- subset_taxa(c_con.ps, taxa_sums(c_con.ps) > 0)
c_con.name <- taxa_names(c_con.ps)

c_hig.ps <- subset_samples(c.ps, Treatment == "High")
c_hig.ps <- subset_taxa(c_hig.ps, taxa_names(c_hig.ps) %in% myc.name)
c_hig.ps <- subset_taxa(c_hig.ps, taxa_sums(c_hig.ps) > 0)
c_hig.name <- taxa_names(c_hig.ps)

c_low.ps <- subset_samples(c.ps, Treatment == "Low")
c_low.ps <- subset_taxa(c_low.ps, taxa_names(c_low.ps) %in% myc.name)
c_low.ps <- subset_taxa(c_low.ps, taxa_sums(c_low.ps) > 0)
c_low.name <- taxa_names(c_low.ps)

# Construct a data.frame that contains TRUE and FALSE values for taxa presence #
c.name <- list(MycoBloom = myc.name, Control = c_con.name, High = c_hig.name, Low = c_low.name)
c.asv <- data.frame(ASV = myc.name)
for (group in names(c.name)) {
  c.asv[[group]] <- c.asv$ASV %in% c.name[[group]]}

# Save the names of the taxa found in just the Mycobloom and ones found in everything but the control #
c_nc.name <- filter(c.asv, High == "TRUE" | Low == "TRUE")
c_nc.name <- filter(c_nc.name, Control == "FALSE")
resave(c_nc.name, file = 'abridged.RData')

c_myc.name <- filter(c.asv, Control == "FALSE" & High == "FALSE" & Low == "FALSE")
resave(c_myc.name, file = './abridged.RData')

# Create the set plot #
c.set <- upset(
  c.asv,
  intersect = c("MycoBloom", "Control", "High", "Low"),
  base_annotations = list(
    'Intersection size' = intersection_size(width = 0.9)
  ),
  queries = list(
    upset_query(intersect = "MycoBloom", color = "black", fill = "orange"),
    upset_query(intersect = c("MycoBloom", "High", "Low"), color = "black", fill = "blue"),
    upset_query(intersect = c("MycoBloom", "High"), color = "black", fill = "blue"),
    upset_query(intersect = c("MycoBloom", "Low"), color = "black", fill = "blue"))
)
resave(c.set, file = './abridged.RData')

## Site D Set Plot
d_con.ps <- subset_samples(d.ps, Treatment == "Control")
d_con.ps <- subset_taxa(d_con.ps, taxa_names(d_con.ps) %in% myc.name)
d_con.ps <- subset_taxa(d_con.ps, taxa_sums(d_con.ps) > 0)
d_con.name <- taxa_names(d_con.ps)

d_hig.ps <- subset_samples(d.ps, Treatment == "High")
d_hig.ps <- subset_taxa(d_hig.ps, taxa_names(d_hig.ps) %in% myc.name)
d_hig.ps <- subset_taxa(d_hig.ps, taxa_sums(d_hig.ps) > 0)
d_hig.name <- taxa_names(d_hig.ps)

d_low.ps <- subset_samples(d.ps, Treatment == "Low")
d_low.ps <- subset_taxa(d_low.ps, taxa_names(d_low.ps) %in% myc.name)
d_low.ps <- subset_taxa(d_low.ps, taxa_sums(d_low.ps) > 0)
d_low.name <- taxa_names(d_low.ps)

# Construct a data.frame that contains TRUE and FALSE values for taxa presence #
d.name <- list(MycoBloom = myc.name, Control = d_con.name, High = d_hig.name, Low = d_low.name)
d.asv <- data.frame(ASV = myc.name)
for (group in names(d.name)) {
  d.asv[[group]] <- d.asv$ASV %in% d.name[[group]]}

# Save the names of the taxa found in just the Mycobloom and ones found in everything but the control #
d_nc.name <- filter(d.asv, High == "TRUE" | Low == "TRUE")
d_nc.name <- filter(d_nc.name, Control == "FALSE")
resave(d_nc.name, file = 'abridged.RData')

d_myc.name <- filter(d.asv, Control == "FALSE" & High == "FALSE" & Low == "FALSE")
resave(d_myc.name, file = './abridged.RData')

# Create the set plot #
d.set <- upset(
  d.asv,
  intersect = c("MycoBloom", "Control", "Low", "High"),
  base_annotations = list(
    'Intersection size' = intersection_size(width = 0.9)
  ),
  queries = list(
    upset_query(intersect = "MycoBloom", color = "black", fill = "orange"),
    upset_query(intersect = c("MycoBloom", "High", "Low"), color = "black", fill = "blue"),
    upset_query(intersect = c("MycoBloom", "Low"), color = "black", fill = "blue"),
    upset_query(intersect = c("MycoBloom", "High"), color = "black", fill = "blue"))
)
resave(d.set, file = './abridged.RData')

## Site E Set Plot
e_con.ps <- subset_samples(e.ps, Treatment == "Control")
e_con.ps <- subset_taxa(e_con.ps, taxa_names(e_con.ps) %in% myc.name)
e_con.ps <- subset_taxa(e_con.ps, taxa_sums(e_con.ps) > 0)
e_con.name <- taxa_names(e_con.ps)

e_hig.ps <- subset_samples(e.ps, Treatment == "High")
e_hig.ps <- subset_taxa(e_hig.ps, taxa_names(e_hig.ps) %in% myc.name)
e_hig.ps <- subset_taxa(e_hig.ps, taxa_sums(e_hig.ps) > 0)
e_hig.name <- taxa_names(e_hig.ps)

e_low.ps <- subset_samples(e.ps, Treatment == "Low")
e_low.ps <- subset_taxa(e_low.ps, taxa_names(e_low.ps) %in% myc.name)
e_low.ps <- subset_taxa(e_low.ps, taxa_sums(e_low.ps) > 0)
e_low.name <- taxa_names(e_low.ps)

# Construct a data.frame that contains TRUE and FALSE values for taxa presence #
e.name <- list(MycoBloom = myc.name, Control = e_con.name, High = e_hig.name, Low = e_low.name)
e.asv <- data.frame(ASV = myc.name)
for (group in names(e.name)) {
  e.asv[[group]] <- e.asv$ASV %in% e.name[[group]]}

# Save the names of the taxa found in just the Mycobloom and ones found in everything but the control #
e_nc.name <- filter(e.asv, High == "TRUE" | Low == "TRUE")
e_nc.name <- filter(e_nc.name, Control == "FALSE")
resave(e_nc.name, file = 'abridged.RData')

e_myc.name <- filter(e.asv, Control == "FALSE" & High == "FALSE" & Low == "FALSE")
resave(e_myc.name, file = './abridged.RData')

# Create the set plot #
e.set <- upset(
  e.asv,
  intersect = c("MycoBloom", "Control", "High", "Low"),
  base_annotations = list(
    'Intersection size' = intersection_size(width = 0.9)
  ),
  queries = list(
    upset_query(intersect = "MycoBloom", color = "black", fill = "orange"),
    upset_query(intersect = c("MycoBloom", "High", "Low"), color = "black", fill = "blue"))
  )
resave(e.set, file = './abridged.RData')

## Site F Set Plot
f_con.ps <- subset_samples(f.ps, Treatment == "Control")
f_con.ps <- subset_taxa(f_con.ps, taxa_names(f_con.ps) %in% myc.name)
f_con.ps <- subset_taxa(f_con.ps, taxa_sums(f_con.ps) > 0)
f_con.name <- taxa_names(f_con.ps)

f_hig.ps <- subset_samples(f.ps, Treatment == "High")
f_hig.ps <- subset_taxa(f_hig.ps, taxa_names(f_hig.ps) %in% myc.name)
f_hig.ps <- subset_taxa(f_hig.ps, taxa_sums(f_hig.ps) > 0)
f_hig.name <- taxa_names(f_hig.ps)

f_low.ps <- subset_samples(f.ps, Treatment == "Low")
f_low.ps <- subset_taxa(f_low.ps, taxa_names(f_low.ps) %in% myc.name)
f_low.ps <- subset_taxa(f_low.ps, taxa_sums(f_low.ps) > 0)
f_low.name <- taxa_names(f_low.ps)

# Construct a data.frame that contains TRUE and FALSE values for taxa presence #
f.name <- list(MycoBloom = myc.name, Control = f_con.name, High = f_hig.name, Low = f_low.name)
f.asv <- data.frame(ASV = myc.name)
for (group in names(e.name)) {
  f.asv[[group]] <- f.asv$ASV %in% f.name[[group]]}

# Save the names of the taxa found in just the Mycobloom and ones found in everything but the control #
f_nc.name <- filter(f.asv, High == "TRUE" | Low == "TRUE")
f_nc.name <- filter(f_nc.name, Control == "FALSE")
resave(f_nc.name, file = 'abridged.RData')

f_myc.name <- filter(f.asv, Control == "FALSE" & High == "FALSE" & Low == "FALSE")
resave(f_myc.name, file = './abridged.RData')

# Create the set plot #
f.set <- upset(
  f.asv,
  intersect = c("MycoBloom", "Control", "High", "Low"),
  base_annotations = list(
    'Intersection size' = intersection_size(width = 0.9)
  ),
  queries = list(
    upset_query(intersect = "MycoBloom", color = "black", fill = "orange"),
    upset_query(intersect = c("MycoBloom", "High", "Low"), color = "black", fill = "blue"),
    upset_query(intersect = c("MycoBloom", "High"), color = "black", fill = "blue"))
)
resave(f.set, file = './abridged.RData')

## Site G Set Plot
g_con.ps <- subset_samples(g.ps, Treatment == "Control")
g_con.ps <- subset_taxa(g_con.ps, taxa_names(g_con.ps) %in% myc.name)
g_con.ps <- subset_taxa(g_con.ps, taxa_sums(g_con.ps) > 0)
g_con.name <- taxa_names(g_con.ps)

g_hig.ps <- subset_samples(g.ps, Treatment == "High")
g_hig.ps <- subset_taxa(g_hig.ps, taxa_names(g_hig.ps) %in% myc.name)
g_hig.ps <- subset_taxa(g_hig.ps, taxa_sums(g_hig.ps) > 0)
g_hig.name <- taxa_names(g_hig.ps)

g_low.ps <- subset_samples(g.ps, Treatment == "Low")
g_low.ps <- subset_taxa(g_low.ps, taxa_names(g_low.ps) %in% myc.name)
g_low.ps <- subset_taxa(g_low.ps, taxa_sums(g_low.ps) > 0)
g_low.name <- taxa_names(g_low.ps)

# Construct a data.frame that contains TRUE and FALSE values for taxa presence #
g.name <- list(MycoBloom = myc.name, Control = g_con.name, High = g_hig.name, Low = g_low.name)
g.asv <- data.frame(ASV = myc.name)
for (group in names(g.name)) {
  g.asv[[group]] <- g.asv$ASV %in% g.name[[group]]}

# Save the names of the taxa found in just the Mycobloom and ones found in everything but the control #
g_nc.name <- filter(g.asv, High == "TRUE" | Low == "TRUE")
g_nc.name <- filter(g_nc.name, Control == "FALSE")
resave(a_nc.name, file = 'abridged.RData')

g_myc.name <- filter(g.asv, Control == "FALSE" & High == "FALSE" & Low == "FALSE")
resave(g_myc.name, file = './abridged.RData')

# Create the set plot #
g.set <- upset(
  g.asv,
  intersect = c("MycoBloom", "Control", "High", "Low"),
  base_annotations = list(
    'Intersection size' = intersection_size(width = 0.9)
  ),
  queries = list(
    upset_query(intersect = "MycoBloom", color = "black", fill = "orange"),
    upset_query(intersect = c("MycoBloom", "High", "Low"), color = "black", fill = "blue"),
    upset_query(intersect = c("MycoBloom", "Low"), color = "black", fill = "blue"),
    upset_query(intersect = c("MycoBloom", "High"), color = "black", fill = "blue"))
)
resave(g.set, file = './abridged.RData')

## Site H Set Plot
h_con.ps <- subset_samples(h.ps, Treatment == "Control")
h_con.ps <- subset_taxa(h_con.ps, taxa_names(h_con.ps) %in% myc.name)
h_con.ps <- subset_taxa(h_con.ps, taxa_sums(h_con.ps) > 0)
h_con.name <- taxa_names(h_con.ps)

h_hig.ps <- subset_samples(h.ps, Treatment == "High")
h_hig.ps <- subset_taxa(h_hig.ps, taxa_names(h_hig.ps) %in% myc.name)
h_hig.ps <- subset_taxa(h_hig.ps, taxa_sums(h_hig.ps) > 0)
h_hig.name <- taxa_names(h_hig.ps)

h_low.ps <- subset_samples(h.ps, Treatment == "Low")
h_low.ps <- subset_taxa(h_low.ps, taxa_names(h_low.ps) %in% myc.name)
h_low.ps <- subset_taxa(h_low.ps, taxa_sums(h_low.ps) > 0)
h_low.name <- taxa_names(h_low.ps)

# Construct a data.frame that contains TRUE and FALSE values for taxa presence #
h.name <- list(MycoBloom = myc.name, Control = h_con.name, High = h_hig.name, Low = h_low.name)
h.asv <- data.frame(ASV = myc.name)
for (group in names(h.name)) {
  h.asv[[group]] <- h.asv$ASV %in% h.name[[group]]}

# Save the names of the taxa found in just the Mycobloom and ones found in everything but the control #
h_nc.name <- filter(h.asv, High == "TRUE" | Low == "TRUE")
h_nc.name <- filter(h_nc.name, Control == "FALSE")
resave(h_nc.name, file = 'abridged.RData')

h_myc.name <- filter(h.asv, Control == "FALSE" & High == "FALSE" & Low == "FALSE")
resave(h_myc.name, file = './abridged.RData')

# Create the set plot #
h.set <- upset(
  h.asv,
  intersect = c("MycoBloom", "Control", "High", "Low"),
  base_annotations = list(
    'Intersection size' = intersection_size(width = 0.9)
  ),
  queries = list(
    upset_query(intersect = "MycoBloom", color = "black", fill = "orange"),
    upset_query(intersect = c("MycoBloom", "High", "Low"), color = "black", fill = "blue"),
    upset_query(intersect = c("MycoBloom", "Low"), color = "black", fill = "blue"),
    upset_query(intersect = c("MycoBloom", "High"), color = "black", fill = "blue"))
  )
resave(h.set, file = './abridged.RData')

#### Percentage of Community Made Up of Things found in the Inoculant ####
# Make a phyloseq object that contains just the total counts of ASVs found in MycoBloom Community #
tax_check.ps <- subset_taxa(final_ksp.ps, taxa_names(final_ksp.ps) %in% taxa_names(myc.ps))
tax_check.ps <- subset_samples(tax_check.ps, Treatment != 'MycoBloom')

# Filter such that only samples with High inoculant levels #
tax_check_hi.ps <- subset_samples(tax_check.ps, Treatment == "High")
hig.ps <- subset_samples(final_ksp.ps, Treatment == "High")

# Filter such that only samples with Low inoculant levels #
tax_check_lo.ps <- subset_samples(tax_check.ps, Treatment == "Low")
low.ps <- subset_samples(final_ksp.ps, Treatment == "Low")

# Filter such that only samples with Control inoculant levels #
tax_check_co.ps <- subset_samples(tax_check.ps, Treatment == "Control")
con.ps <- subset_samples(final_ksp.ps, Treatment == "Control")

# Calculate and save the percentages in which the community is made up of MycoBloom Members #
ksp_perc_co.df <- data.frame(Percent = sample_sums(tax_check_co.ps)/ sample_sums(con.ps),
                             Treatment = "Control")
ksp_perc_hi.df <- data.frame(Percent = sample_sums(tax_check_hi.ps)/ sample_sums(hig.ps),
                             Treatment = "High")
ksp_perc_lo.df <- data.frame(Percent = sample_sums(tax_check_lo.ps)/ sample_sums(low.ps),
                             Treatment = "Low")
ksp_perc.df <- rbind(ksp_perc_co.df, ksp_perc_hi.df, ksp_perc_lo.df)
ksp_perc.df$Site <- substring(rownames(ksp_perc.df), 1,1)

ksp_perc.df <- arrange(ksp_perc.df, Site, Treatment)

# Perform a One way Anova to see if there are differences #
ksp_perc.df$Prop <- ksp_perc.df$Percent / 100
ksp_perc.df$Trans <- asin(sqrt(ksp_perc.df$Prop))
ksp_perc.fit <- lm(Trans ~ Treatment * Site, data = ksp_perc.df)
ksp_perc.aov <- aov(ksp_perc.fit)
summary(ksp_perc.aov)
leveneTest(ksp_perc.fit)
shapiro.test(resid(ksp_perc.fit))
ksp_perc.hsd <- as.data.frame(TukeyHSD(ksp_perc.aov, 'Treatment:Site', ordered = TRUE)$`Treatment:Site`)
ksp_perc.sig <- filter(ksp_perc.hsd, `p adj` < 0.05)
# Treatment does not significantly affect the percentage of community composition made up of MycoBloom ASVs # 

#### Differential Abundance ####
if(!requireNamespace("maaslin3")) BiocManager::install("biobakery/maaslin3")
library(maaslin3); packageVersion("maaslin3")

if(!dir.exists('./maas')){
  dir.create('./maas')
}

# Make a phyloseq object that excludes the MycoBloom Community #
nomy.ps <- subset_samples(final_ksp.ps, Treatment != "MycoBloom")
decompose_ps(nomy.ps, 'nomy')

### Site A ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
a.maas <- maaslin3(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/a.maas',
                   formula = '~ Group',
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Group,AC'))
# Save the results tables for abundance and prevalence, respectively #  
a_abun.res <- a.maas$fit_data_abundance$results
a_prev.res <- a.maas$fit_data_prevalence$results

## Abundance Associations ##
# Save only the results of with values corresponding to communities found within Site A #
a_alla.res <- c()
for(i in 1:nrow(a_abun.res)){
  if(a_abun.res$value[i] == "AHI" | a_abun.res$value[i] == "ALO"){
    a_alla.res <- rbind(a_alla.res, a_abun.res[i,])
  }
}

# Save all of the significant results to a new data frame #
a_siga.res <- c()
for(i in 1:nrow(a_alla.res)){
  if(!is.na(a_alla.res$qval_individual[i]) && a_alla.res$qval_individual[i] < 0.1){
    a_siga.res <- rbind(a_siga.res, a_alla.res[i,])
  }
}

## Prevalence Associations ##
# Save only the results of with values corresponding to communities found within Site A #
a_allp.res <- c()
for(i in 1:nrow(a_prev.res)){
  if(a_prev.res$value[i] == "AHI" | a_prev.res$value[i] == "ALO"){
    a_allp.res <- rbind(a_allp.res, a_prev.res[i,])
  }
}

# Save all of the significant results to a new data frame #
a_sigp.res <- c()
for(i in 1:nrow(a_allp.res)){
  if(!is.na(a_allp.res$qval_individual[i]) && a_allp.res$qval_individual[i] < 0.1){
    a_sigp.res <- rbind(a_sigp.res, a_allp.res[i,])
  }
}

### Site B ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
b.maas <- maaslin3(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/b.maas',
                   formula = '~ Group',
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Group,BC'))
# Save the results tables for abundance and prevalence, respectively #  
b_abun.res <- b.maas$fit_data_abundance$results
b_prev.res <- b.maas$fit_data_prevalence$results

## Abundance Associations ##
# Save only the results of with values corresponding to communities found within Site B #
b_alla.res <- c()
for(i in 1:nrow(b_abun.res)){
  if(b_abun.res$value[i] == "BHI" | b_abun.res$value[i] == "BLO"){
    b_alla.res <- rbind(b_alla.res, b_abun.res[i,])
  }
}

# Save all of the significant results to a new data frame #
b_siga.res <- c()
for(i in 1:nrow(b_alla.res)){
  if(!is.na(b_alla.res$qval_individual[i]) && b_alla.res$qval_individual[i] < 0.1){
    b_siga.res <- rbind(b_siga.res, b_alla.res[i,])
  }
}

## Prevalence Associations ##
# Save only the results of with values corresponding to communities found within Site B #
b_allp.res <- c()
for(i in 1:nrow(b_prev.res)){
  if(b_prev.res$value[i] == "BHI" | b_prev.res$value[i] == "BLO"){
    b_allp.res <- rbind(b_allp.res, b_prev.res[i,])
  }
}

# Save all of the significant results to a new data frame #
b_sigp.res <- c()
for(i in 1:nrow(b_allp.res)){
  if(!is.na(b_allp.res$qval_individual[i]) && b_allp.res$qval_individual[i] < 0.1){
    b_sigp.res <- rbind(b_sigp.res, b_allp.res[i,])
  }
}

### Site C ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
c.maas <- maaslin3(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/c.maas',
                   formula = '~ Group',
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Group,CC'))
# Save the results tables for abundance and prevalence, respectively #  
c_abun.res <- c.maas$fit_data_abundance$results
c_prev.res <- c.maas$fit_data_prevalence$results

## Abundance Associations ##
# Save only the results of with values corresponding to communities found within Site B #
c_alla.res <- c()
for(i in 1:nrow(c_abun.res)){
  if(c_abun.res$value[i] == "CHI" | c_abun.res$value[i] == "CLO"){
    c_alla.res <- rbind(c_alla.res, c_abun.res[i,])
  }
}

# Save all of the significant results to a new data frame #
c_siga.res <- c()
for(i in 1:nrow(c_alla.res)){
  if(!is.na(c_alla.res$qval_individual[i]) && c_alla.res$qval_individual[i] < 0.1){
    c_siga.res <- rbind(c_siga.res, c_alla.res[i,])
  }
}

## Prevalence Associations ##
# Save only the results of with values corresponding to communities found within Site B #
c_allp.res <- c()
for(i in 1:nrow(c_prev.res)){
  if(c_prev.res$value[i] == "CHI" | c_prev.res$value[i] == "CLO"){
    c_allp.res <- rbind(c_allp.res, c_prev.res[i,])
  }
}

# Save all of the significant results to a new data frame #
c_sigp.res <- c()
for(i in 1:nrow(c_allp.res)){
  if(!is.na(c_allp.res$qval_individual[i]) && c_allp.res$qval_individual[i] < 0.1){
    c_sigp.res <- rbind(c_sigp.res, c_allp.res[i,])
  }
}

### Site D ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
d.maas <- maaslin3(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/d.maas',
                   formula = '~ Group',
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Group,DC'))
# Save the results tables for abundance and prevalence, respectively #  
d_abun.res <- d.maas$fit_data_abundance$results
d_prev.res <- d.maas$fit_data_prevalence$results

## Abundance Associations ##
# Save only the results of with values corresponding to communities found within Site B #
d_alla.res <- c()
for(i in 1:nrow(d_abun.res)){
  if(d_abun.res$value[i] == "DHI" | d_abun.res$value[i] == "DLO"){
    d_alla.res <- rbind(d_alla.res, d_abun.res[i,])
  }
}

# Save all of the significant results to a new data frame #
d_siga.res <- c()
for(i in 1:nrow(d_alla.res)){
  if(!is.na(d_alla.res$qval_individual[i]) && d_alla.res$qval_individual[i] < 0.1){
    d_siga.res <- rbind(d_siga.res, d_alla.res[i,])
  }
}

## Prevalence Associations ##
# Save only the results of with values corresponding to communities found within Site B #
d_allp.res <- c()
for(i in 1:nrow(d_prev.res)){
  if(d_prev.res$value[i] == "DHI" | d_prev.res$value[i] == "DLO"){
    d_allp.res <- rbind(d_allp.res, d_prev.res[i,])
  }
}

# Save all of the significant results to a new data frame #
d_sigp.res <- c()
for(i in 1:nrow(d_allp.res)){
  if(!is.na(d_allp.res$qval_individual[i]) && d_allp.res$qval_individual[i] < 0.1){
    d_sigp.res <- rbind(d_sigp.res, d_allp.res[i,])
  }
}

### Site E ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
e.maas <- maaslin3(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/e.maas',
                   formula = '~ Group',
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Group,EC'))
# Save the results tables for abundance and prevalence, respectively #  
e_abun.res <- e.maas$fit_data_abundance$results
e_prev.res <- e.maas$fit_data_prevalence$results

## Abundance Associations ##
# Save only the results of with values corresponding to communities found within Site B #
e_alla.res <- c()
for(i in 1:nrow(e_abun.res)){
  if(e_abun.res$value[i] == "EHI" | e_abun.res$value[i] == "ELO"){
    e_alla.res <- rbind(e_alla.res, e_abun.res[i,])
  }
}

# Save all of the significant results to a new data frame #
e_siga.res <- c()
for(i in 1:nrow(e_alla.res)){
  if(!is.na(e_alla.res$qval_individual[i]) && e_alla.res$qval_individual[i] < 0.1){
    e_siga.res <- rbind(e_siga.res, e_alla.res[i,])
  }
}

## Prevalence Associations ##
# Save only the results of with values corresponding to communities found within Site B #
e_allp.res <- c()
for(i in 1:nrow(e_prev.res)){
  if(e_prev.res$value[i] == "EHI" | e_prev.res$value[i] == "ELO"){
    e_allp.res <- rbind(e_allp.res, e_prev.res[i,])
  }
}

# Save all of the significant results to a new data frame #
e_sigp.res <- c()
for(i in 1:nrow(e_allp.res)){
  if(!is.na(e_allp.res$qval_individual[i]) && e_allp.res$qval_individual[i] < 0.1){
    e_sigp.res <- rbind(e_sigp.res, e_allp.res[i,])
  }
}

### Site F ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
f.maas <- maaslin3(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/f.maas',
                   formula = '~ Group',
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Group,FC'))
# Save the results tables for abundance and prevalence, respectively #  
f_abun.res <- f.maas$fit_data_abundance$results
f_prev.res <- f.maas$fit_data_prevalence$results

## Abundance Associations ##
# Save only the results of with values corresponding to communities found within Site B #
f_alla.res <- c()
for(i in 1:nrow(f_abun.res)){
  if(f_abun.res$value[i] == "FHI" | f_abun.res$value[i] == "FLO"){
    f_alla.res <- rbind(f_alla.res, f_abun.res[i,])
  }
}

# Save all of the significant results to a new data frame #
f_siga.res <- c()
for(i in 1:nrow(f_alla.res)){
  if(!is.na(f_alla.res$qval_individual[i]) && f_alla.res$qval_individual[i] < 0.1){
    f_siga.res <- rbind(f_siga.res, f_alla.res[i,])
  }
}

## Prevalence Associations ##
# Save only the results of with values corresponding to communities found within Site B #
f_allp.res <- c()
for(i in 1:nrow(f_prev.res)){
  if(f_prev.res$value[i] == "FHI" | f_prev.res$value[i] == "FLO"){
    f_allp.res <- rbind(f_allp.res, f_prev.res[i,])
  }
}

# Save all of the significant results to a new data frame #
f_sigp.res <- c()
for(i in 1:nrow(f_allp.res)){
  if(!is.na(f_allp.res$qval_individual[i]) && f_allp.res$qval_individual[i] < 0.1){
    f_sigp.res <- rbind(f_sigp.res, f_allp.res[i,])
  }
}

### Site G ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
g.maas <- maaslin3(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/g.maas',
                   formula = '~ Group',
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Group,GC'))
# Save the results tables for abundance and prevalence, respectively #  
g_abun.res <- g.maas$fit_data_abundance$results
g_prev.res <- g.maas$fit_data_prevalence$results

## Abundance Associations ##
# Save only the results of with values corresponding to communities found within Site B #
g_alla.res <- c()
for(i in 1:nrow(g_abun.res)){
  if(g_abun.res$value[i] == "GHI" | g_abun.res$value[i] == "GLO"){
    g_alla.res <- rbind(g_alla.res, g_abun.res[i,])
  }
}

# Save all of the significant results to a new data frame #
g_siga.res <- c()
for(i in 1:nrow(g_alla.res)){
  if(!is.na(g_alla.res$qval_individual[i]) && g_alla.res$qval_individual[i] < 0.1){
    g_siga.res <- rbind(g_siga.res, g_alla.res[i,])
  }
}

## Prevalence Associations ##
# Save only the results of with values corresponding to communities found within Site B #
g_allp.res <- c()
for(i in 1:nrow(g_prev.res)){
  if(g_prev.res$value[i] == "GHI" | g_prev.res$value[i] == "GLO"){
    g_allp.res <- rbind(g_allp.res, g_prev.res[i,])
  }
}

# Save all of the significant results to a new data frame #
g_sigp.res <- c()
for(i in 1:nrow(g_allp.res)){
  if(!is.na(g_allp.res$qval_individual[i]) && g_allp.res$qval_individual[i] < 0.1){
    g_sigp.res <- rbind(g_sigp.res, g_allp.res[i,])
  }
}

### Site H ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
h.maas <- maaslin3(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/h.maas',
                   formula = '~ Group',
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Group,HC'))
# Save the results tables for abundance and prevalence, respectively #  
h_abun.res <- h.maas$fit_data_abundance$results
h_prev.res <- h.maas$fit_data_prevalence$results

## Abundance Associations ##
# Save only the results of with values corresponding to communities found within Site B #
h_alla.res <- c()
for(i in 1:nrow(h_abun.res)){
  if(h_abun.res$value[i] == "HHI" | h_abun.res$value[i] == "HLO"){
    h_alla.res <- rbind(h_alla.res, h_abun.res[i,])
  }
}

# Save all of the significant results to a new data frame #
h_siga.res <- c()
for(i in 1:nrow(h_alla.res)){
  if(!is.na(h_alla.res$qval_individual[i]) && h_alla.res$qval_individual[i] < 0.1){
    h_siga.res <- rbind(h_siga.res, h_alla.res[i,])
  }
}

## Prevalence Associations ##
# Save only the results of with values corresponding to communities found within Site B #
h_allp.res <- c()
for(i in 1:nrow(h_prev.res)){
  if(h_prev.res$value[i] == "HHI" | h_prev.res$value[i] == "HLO"){
    h_allp.res <- rbind(h_allp.res, h_prev.res[i,])
  }
}

# Save all of the significant results to a new data frame #
h_sigp.res <- c()
for(i in 1:nrow(h_allp.res)){
  if(!is.na(h_allp.res$qval_individual[i]) && h_allp.res$qval_individual[i] < 0.1){
    h_sigp.res <- rbind(h_sigp.res, h_allp.res[i,])
  }
}

save.image("./ksp_amf.RData")

# This just outputs the final time in which the Rscript finishes running #
cat("\n## Script finished at", Sys.time(), "\n")
sink(type = "message")
sink()