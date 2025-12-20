# This is an R script on how to use the dada2 pipeline to denoise Illumina amplicon sequencing data # 
# into analyzed sequence variants (ASVs). This tutorial is based on the published #
# tutorial (https://benjjneb.github.io/dada2/tutorial_1_8.html) and is adapted to fit # 
# the Kankakee Sand Prairie fungal SSU sequencing data #

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
    fra.tab <- cbind(tax.tab, as.character(dna.tab),otu.tab)
  } else {
    fra.tab <- cbind(tax.tab, otu.tab)
  }
  for(j in 1:ncol(fra.tab)){
    if(colnames(fra.tab)[j] == "as.character(dna.tab)"){
      colnames(fra.tab)[j] = "DNA Sequence"
    }
  }
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
      fra = fra.tab,
      ) 
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
filt_ksp.ps <- subset_taxa(filt_ksp.ps, taxa_sums(filt_ksp.ps) > 0)

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

# As a final taxa limiting step, join taxa that are 95% related (According to the phylogentic distances between ASVs) # 
final_ksp.ps <- tip_glom(final_ksp.ps, h = 0.05)

# Fix the taxonomy table using the original tax table #
decompose_ps(final_ksp.ps, 'final_ksp')
fin_ksp.tax <- c()
for(i in rownames(final_ksp.tax)){
  if(i %in% taxa_names(final_ksp.ps)){
    fin_ksp.tax <- rbind(fin_ksp.tax, final_ksp.tax[i,c("Family", "Genus", "Species", "ASV")])
  }
}
tax_table(final_ksp.ps) <- as.matrix(fin_ksp.tax)

# Fix the ASVs such that they are numbered in order of abundance #
fix_tax_names <- function(ps, label){
  # function that fixes the names of taxa such that they represent their order of abundance #
  tax.list <- names(sort(taxa_sums(ps), decreasing = TRUE))
  tax.df <- as.data.frame(tax_table(ps))
  new_tax.df <- tax.df[tax.list,]
  for(i in 1:nrow(new_tax.df)){
    new_tax.df$Abundance_Rank[i] <- i
    new_tax.df$Lowest[i] <- sub(".*\\((.*)\\)", "\\1", new_tax.df$ASV[i])
    new_tax.df$ASV[i] <- paste0('ASV', as.character(new_tax.df$Abundance_Rank[i]), '(', new_tax.df$Lowest[i], ')')}
  new_tax.df <- new_tax.df[taxa_names(ps),]
  tax_table(ps) <- as.matrix(new_tax.df)
  taxa_names(ps) <- new_tax.df$ASV
  assign(label, ps, envir = .GlobalEnv)
}

fix_tax_names(final_ksp.ps, 'final_ksp.ps')

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

## Shannon Diversity ##
# Perform the ANOVA and normality and variance tests for a multiple linear regression #
ksp_sha.fit <- lm(Shannon ~ Treatment * Site, data = ksp.rich)
ksp_sha.sum <- Anova(ksp_sha.fit, type = 3) 
ksp_sha.sum
shapiro.test(residuals(ksp_sha.fit))
leveneTest(ksp_sha.fit)
AIC(ksp_sha.fit)
resave(ksp_sha.sum, file = './abridged.RData')

## Shannon Evenness ##
# Perform the ANOVA and normality and variance tests #
ksp_sev.fit <- lm(ShaEvn ~ Treats * Sites, data = ksp.rich)
ksp_sev.sum <- Anova(ksp_sev.fit, type = 3)
ksp_sev.sum
shapiro.test(residuals(ksp_sev.fit))
leveneTest(ksp_sev.fit)
AIC(ksp_sev.fit)
resave(ksp_sev.emn, file = './abridged.RData')


# Perform the ANOVA and normality and variance tests #
ksp_cha.fit <- lm(Chao1 ~ Treats * Sites, ksp.rich) 
ksp_cha.sum <- Anova(ksp_cha.fit, type = 3)
ksp_cha.sum
shapiro.test(residuals(ksp_cha.fit))
leveneTest(ksp_cha.fit)
AIC(ksp_cha.fit)
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
  geom_text(aes(label = sim.let, y = evn.mean + evn.sd + 0.025), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 8) +
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
  scale_y_continuous(limits = c(0,160), breaks = seq(0, 150, 30), expand = expansion(mult = c(0,0.05))) +
  geom_text(aes(label = cha.let, y = cha.mean + cha.sd + 5), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 8) +
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
### Make an upset plot to determine whether to use weighted or unweighted unifrac ###

## Subset the whole phyloseq object by treatments and filter to only contain taxa found in the MycoBloom ##
# Control #
con.ps <- subset_samples(final_ksp.ps, Treatment == "Control")
con.ps <- subset_taxa(con.ps, taxa_sums(con.ps) > 0)
con.name <- taxa_names(con.ps)

# High #
hig.ps <- subset_samples(final_ksp.ps, Treatment == "High")
hig.ps <- subset_taxa(hig.ps, taxa_sums(hig.ps) > 0)
hig.name <- taxa_names(hig.ps)

# Low #
low.ps <- subset_samples(final_ksp.ps, Treatment == "Low")
low.ps <- subset_taxa(low.ps, taxa_sums(low.ps) > 0)
low.name <- taxa_names(low.ps)

all.usdf <- list(Control = con.name, High = hig.name, Low = low.name)
all.asv <- data.frame(ASV = taxa_names(final_ksp.ps))
for (group in names(all.usdf)) {
  all.asv[[group]] <- all.asv$ASV %in% all.usdf[[group]]
}

for(i in 2:4){
  for(j in 1:nrow(all.asv)){
    all.asv[j,i] <- as.integer(all.asv[j,i])
  }
}

if(!requireNamespace('UpSetR')) install.packages('UpSetR')
library(UpSetR); packageVersion('UpSetR')

all.set <- upset(all.asv, order.by = 'freq')
resave(all.set, file = './abridged.RData')

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

resave(ksp_by.perm, file = './abridged.RData')
resave(ksp.perm, file = './abridged.RData')

beta.plot<- plot_ordination(ksp_prop.ps, ksp_unif.pcoa, shape = "Treats", color="Site", title="All Samples") +
  theme_prism() +
  geom_point(size = 10) +
  scale_x_continuous(name = "PCoA 1 (13.45%)") +
  scale_y_continuous(name = "PCoA 2 (11.24%)") +
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  annotate("richtext", x = 0.05, y = -0.29,
           label = paste0("R<sup>2</sup> = ", round(ksp.perm$R2[1], 3), "; pseudo-<em>F</em><sub>24,48</sub> = ", round(ksp.perm$F[1], 3),'; <em>p</em> < 0.001'),
           size = 10, fontface = 'bold', family = 'Liberation Sans')
  
TukeyHSD(betadisper(ksp.bdist, group = final_ksp$met$Treats, type = 'median'))

resave(beta.plot, file = './abridged.RData')

#### Per Sample Analysis ####
## A ##
# Create a phyloseq object with just the sites of interest and calculate weighted unifrac distance #
a.ps <- subset_samples(final_ksp.ps, Site == 'A')
a.ps <- subset_taxa(a.ps, taxa_sums(a.ps) > 0)
a.met <- as(sample_data(a.ps), 'data.frame')
a_prop.ps <- transform_sample_counts(a.ps, function(otu) otu/sum(otu))
a.dist <- as.dist(as.matrix(ksp.bdist)[sample_names(a.ps), sample_names(a.ps)])
a_ord.nmds.wuni <- ordinate(a_prop.ps, method="NMDS", distance= a.dist)

# Perform PermANOVA #
a.perm <- adonis2(a.dist~Treatment, data = a.met)
a.perm
resave(a.perm, file = 'abridged.RData')

# Plot the results #
a_beta.plot <- plot_ordination(a_prop.ps, a_ord.nmds.wuni, color="Treatment", title="A Samples NMDS") +
  theme_prism() +
  geom_point(size = 12)

TukeyHSD(betadisper(a.dist, group = sample_data(a.ps)$Treats, type = 'median'))

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
b.perm <- adonis2(b.dist~Treatment, data = b.met)
b.perm
resave(b.perm, file = './abridged.RData')

# Plot the results #
b_beta.plot <- plot_ordination(b_prop.ps, b_ord.nmds.wuni, color="Treatment", title="B Samples NMDS") +
  theme_prism() +
  geom_point(size = 12)
resave(b_beta.plot, file = './abridged.RData')

TukeyHSD(betadisper(b.dist, group = sample_data(b.ps)$Treats, type = 'median'))

## C ##
# Create a phyloseq object with just the sites of interest and calculate weighted unifrac distance #
c.ps <- subset_samples(final_ksp.ps, Site == 'C')
c.ps <- subset_taxa(c.ps, taxa_sums(c.ps) > 0)
c.met <- as(sample_data(c.ps), 'data.frame')
c_prop.ps <- transform_sample_counts(c.ps, function(otu) otu/sum(otu))
c.dist <- as.dist(as.matrix(ksp.bdist)[sample_names(c.ps), sample_names(c.ps)])
c_ord.nmds.wuni <- ordinate(c_prop.ps, method="NMDS", distance=c.dist)

# Perform PermANOVA #
c.perm <- adonis2(c.dist~Treatment, data = c.met)
c.perm
resave(c.perm, file = './abridged.RData')

# Plot the results #
c_beta.plot <- plot_ordination(c_prop.ps, c_ord.nmds.wuni, color="Treatment", title="C Samples NMDS") +
  theme_prism() +
  geom_point(size = 12)
resave(c_beta.plot, file = './abridged.RData')

TukeyHSD(betadisper(c.dist, group = sample_data(c.ps)$Treats, type = 'median'))

## D ##
# Create a phyloseq object with just the sites of interest and calculate weighted unifrac distance #
d.ps <- subset_samples(final_ksp.ps, Site == 'D')
d.ps <- subset_taxa(d.ps, taxa_sums(d.ps) > 0)
d.met <- as(sample_data(d.ps), 'data.frame')
d_prop.ps <- transform_sample_counts(d.ps, function(otu) otu/sum(otu))
d.dist <- as.dist(as.matrix(ksp.bdist)[sample_names(d.ps), sample_names(d.ps)])
d_ord.nmds.wuni <- ordinate(d_prop.ps, method="NMDS", distance=d.dist)

# Perform PermANOVA #
d.perm <- adonis2(d.dist~Treatment, data = d.met)
d.perm
resave(d.perm, file = './abridged.RData')

# Plot the results #
d_beta.plot <- plot_ordination(d_prop.ps, d_ord.nmds.wuni, color="Treatment", title="D Samples NMDS") +
  theme_prism() +
  geom_point(size = 12)
resave(d_beta.plot, file = './abridged.RData')

TukeyHSD(betadisper(b.dist, group = sample_data(b.ps)$Treats, type = 'median'))

## E ##
# Create a phyloseq object with just the sites of interest and calculate weighted unifrac distance #
e.ps <- subset_samples(final_ksp.ps, Site == 'E')
e.ps <- subset_taxa(e.ps, taxa_sums(e.ps) > 0)
e.met <- as(sample_data(e.ps), 'data.frame')
e_prop.ps <- transform_sample_counts(e.ps, function(otu) otu/sum(otu))
e.dist <- as.dist(as.matrix(ksp.bdist)[sample_names(e.ps), sample_names(e.ps)])
e_ord.nmds.wuni <- ordinate(e_prop.ps, method="NMDS", distance=e.dist)

# Perform PermANOVA #
e.perm <- adonis2(e.dist~Treatment, data = e.met)
e.perm
resave(e.perm, file = 'abridged.RData')

# Plot the results #
e_beta.plot <- plot_ordination(e_prop.ps, e_ord.nmds.wuni, color="Treatment", title="E Samples NMDS") +
  theme_prism() +
  geom_point(size = 12)
resave(e_beta.plot, file = './abridged.RData')

TukeyHSD(betadisper(e.dist, group = sample_data(e.ps)$Treats, type = 'median'))

## F ##
# Create a phyloseq object with just the sites of interest and calculate weighted unifrac distance #
f.ps <- subset_samples(final_ksp.ps, Site == 'F')
f.ps <- subset_taxa(f.ps, taxa_sums(f.ps) > 0)
f.met <- as(sample_data(f.ps), 'data.frame')
f_prop.ps <- transform_sample_counts(f.ps, function(otu) otu/sum(otu))
f.dist <- as.dist(as.matrix(ksp.bdist)[sample_names(f.ps), sample_names(f.ps)])
f_ord.nmds.wuni <- ordinate(f_prop.ps, method="NMDS", distance=f.dist)

# Perform PermANOVA #
f.perm <- adonis2(f.dist~Treatment, data = f.met)
f.perm
resave(f.perm, file = './abridged.RData')

# Plot the results #
f_beta.plot <- plot_ordination(f_prop.ps, f_ord.nmds.wuni, color="Treatment", title="F Samples NMDS") +
  theme_prism() +
  geom_point(size = 12)
resave(f_beta.plot, file = './abridged.RData')

TukeyHSD(betadisper(f.dist, group = sample_data(f.ps)$Treats, type = 'median'))

## G ##
# Create a phyloseq object with just the sites of interest and calculate weighted unifrac distance #
g.ps <- subset_samples(final_ksp.ps, Site == 'G')
g.ps <- subset_taxa(g.ps, taxa_sums(g.ps) > 0)
g.met <- as(sample_data(g.ps), 'data.frame')
g_prop.ps <- transform_sample_counts(g.ps, function(otu) otu/sum(otu))
g.dist <- as.dist(as.matrix(ksp.bdist)[sample_names(g.ps), sample_names(g.ps)])
g_ord.nmds.wuni <- ordinate(g_prop.ps, method="NMDS", distance=g.dist)

# Perform PermANOVA #
g.perm <- adonis2(g.dist~Treatment, data = g.met)
g.perm
resave(g.perm, file = './abridged.RData')

# Plot the results #
g_beta.plot <- plot_ordination(g_prop.ps, g_ord.nmds.wuni, color="Treatment", title="G Samples NMDS") +
  theme_prism() +
  geom_point(size = 12)
resave(g_beta.plot, file = './abridged.RData')

TukeyHSD(betadisper(g.dist, group = sample_data(g.ps)$Treats, type = 'median'))

## H ##
# Create a phyloseq object with just the sites of interest and calculate weighted unifrac distance #
h.ps <- subset_samples(final_ksp.ps, Site == 'H')
h.ps <- subset_taxa(h.ps, taxa_sums(h.ps) > 0)
h.met <- as(sample_data(h.ps), 'data.frame')
h_prop.ps <- transform_sample_counts(h.ps, function(otu) otu/sum(otu))
h.dist <- as.dist(as.matrix(ksp.bdist)[sample_names(h.ps), sample_names(h.ps)])
h_ord.nmds.wuni <- ordinate(h_prop.ps, method="NMDS", distance=h.dist)

# Perform PermANOVA #
h.perm <- adonis2(h.dist~Treatment, data = h.met)
h.perm
resave(h.perm, file = './abridged.RData')

# Plot the results #
h_beta.plot <- plot_ordination(h_prop.ps, h_ord.nmds.wuni, color="Treatment", title="H Samples NMDS") +
  theme_prism() +
  geom_point(size = 12)
resave(h_beta.plot, file = './abridged.RData')

TukeyHSD(betadisper(h.dist, group = sample_data(h.ps)$Treats, type = 'median'))

#### Stacked Histograms ####
# First we start by making a color pallette for each unique ASV #
if(!requireNamespace("Polychrome")) install.packages("Polychrome")
library(Polychrome); packageVersion("Polychrome")
set.seed(248)
ksp.colr <- createPalette(ntaxa(final_ksp.ps),  c("#ff0000", "#00ff00", "#0000ff"))
ksp.colr <- as.data.frame(ksp.colr)
rownames(ksp.colr) <- taxa_names(final_ksp.ps)

# Add a gray color for "Other"
ksp.colr[ntaxa(final_ksp.ps) + 1,] <- "#D4D4D4" 
rownames(ksp.colr)[ntaxa(final_ksp.ps) + 1] <- "Other" 

## All non-mycobloom samples ##
hg_all.ps <- subset_samples(final_ksp.ps, Treatment != 'MycoBloom')
hg_all.ps <- subset_samples(hg_all.ps, taxa_sums(hg_all.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_all.name <- names(sort(taxa_sums(hg_all.ps), decreasing = TRUE))[1:19]
hg_all.name <- c(hg_all.name, "Other")
hg_all.colr <- ksp.colr[hg_all.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_all.ps <- merge_taxa(hg_all.ps, taxa_names(hg_all.ps)[!taxa_names(hg_all.ps) %in% hg_all.name])
taxa_names(hg_all.ps)[5] <- "Other"
tax_table(hg_all.ps)[5, "ASV"] <- "Other"

# Save the phyloseq object as a data.frame and make factors that guide the plot what to plot #
hg_all.df <- psmelt(hg_all.ps)
hg_all.df$ASVs <- factor(hg_all.df$ASV, levels = hg_all.name)
hg_all.df$Soil <- factor(hg_all.df$Treatment, levels = c("Control", "Low", "High"))

# Plot the histogram #
hg_all.plot <- ggplot(hg_all.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(name = "Fungal ASV", values = hg_all.colr) +
  scale_y_continuous(breaks = seq(0,1, 0.20), expand = expansion(mult = c(0.002, 0.02))) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(color = "black", size = 28, family = "Liberation Sans", face = 'bold'),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        legend.background = element_rect(color = 'black'),
        axis.text.y.left  = element_text(size = 24, family = "Liberation Sans", face = 'bold', vjust = 0.45),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        axis.ticks.y.left = element_line(linewidth = 1.5),
        legend.key.spacing.y = unit(2, 'mm'),
        legend.position = 'right')

## MycoBloom Sample ##
# Save a phyloseq Object that contains only the samples of interest #
myc.ps <- subset_samples(final_ksp.ps, Treatment == "MycoBloom")
myc.ps <- subset_taxa(myc.ps, taxa_sums(myc.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_myc.name <- names(sort(taxa_sums(myc.ps), decreasing = TRUE))[1:19]
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
  scale_y_continuous(breaks = seq(0,1, 0.20), expand = expansion(mult = c(0.002, 0.02)), sec.axis = dup_axis(name = "MycoBloom")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(color = "black", size = 28, family = "Liberation Sans", face = 'bold'),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        legend.background = element_rect(color = 'black'),
        axis.text.y.left  = element_text(size = 24, family = "Liberation Sans", face = 'bold', vjust = 0.45),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        axis.ticks.y.left = element_line(linewidth = 1.5),
        legend.key.spacing.y = unit(2, 'mm'),
        legend.position = 'right')
resave(hg_myc.plot, file = './abridged.RData')

## Site A Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_a.ps <- subset_samples(final_ksp.ps, Site == "A")
hg_a.ps <- subset_taxa(hg_a.ps, taxa_sums(hg_a.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_a.name <- names(sort(taxa_sums(hg_a.ps), decreasing = TRUE))[1:19]
hg_a.name <- c(hg_a.name, "Other")
hg_a.colr <- ksp.colr[hg_a.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_a.ps <- merge_taxa(hg_a.ps, taxa_names(hg_a.ps)[!taxa_names(hg_a.ps) %in% hg_a.name])
taxa_names(hg_a.ps)[16] <- "Other"
tax_table(hg_a.ps)[16, "ASV"] <- "Other"

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
  scale_y_continuous(breaks = seq(0,1, 0.20), expand = expansion(mult = c(0.002, 0.02)), sec.axis = dup_axis(name = "Site A")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(color = "black", size = 28, family = "Liberation Sans", face = 'bold'),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        legend.background = element_rect(color = 'black'),
        axis.text.y.left  = element_text(size = 24, family = "Liberation Sans", face = 'bold', vjust = 0.45),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        axis.ticks.y.left = element_line(linewidth = 1.5),
        legend.key.spacing.y = unit(2, 'mm'),
        legend.position = 'right')

resave(hg_a.plot, file = './abridged.RData')

## Site B Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_b.ps <- subset_samples(final_ksp.ps, Site == "B")
hg_b.ps <- subset_taxa(hg_b.ps, taxa_sums(hg_b.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_b.name <- names(sort(taxa_sums(hg_b.ps), decreasing = TRUE))[1:19]
hg_b.name <- c(hg_b.name, "Other")
hg_b.colr <- ksp.colr[hg_b.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_b.ps <- merge_taxa(hg_b.ps, taxa_names(hg_b.ps)[!taxa_names(hg_b.ps) %in% hg_b.name])
taxa_names(hg_b.ps)[16] <- "Other"
tax_table(hg_b.ps)[16, "ASV"] <- "Other"

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
  scale_y_continuous(breaks = seq(0,1, 0.20), expand = expansion(mult = c(0.002, 0.02)), sec.axis = dup_axis(name = "Site B")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(color = "black", size = 28, family = "Liberation Sans", face = 'bold'),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        legend.background = element_rect(color = 'black'),
        axis.text.y.left  = element_text(size = 24, family = "Liberation Sans", face = 'bold', vjust = 0.45),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        axis.ticks.y.left = element_line(linewidth = 1.5),
        legend.key.spacing.y = unit(2, 'mm'),
        legend.position = 'right')

resave(hg_b.plot, file = './abridged.RData')

## Site C Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_c.ps <- subset_samples(final_ksp.ps, Site == "C")
hg_c.ps <- subset_taxa(hg_c.ps, taxa_sums(hg_c.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_c.name <- names(sort(taxa_sums(hg_c.ps), decreasing = TRUE))[1:19]
hg_c.name <- c(hg_c.name, "Other")
hg_c.colr <- ksp.colr[hg_c.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_c.ps <- merge_taxa(hg_c.ps, taxa_names(hg_c.ps)[!taxa_names(hg_c.ps) %in% hg_c.name])
taxa_names(hg_c.ps)[20] <- "Other"
tax_table(hg_c.ps)[20, "ASV"] <- "Other"

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
  scale_y_continuous(breaks = seq(0,1, 0.20), expand = expansion(mult = c(0.002, 0.02)), sec.axis = dup_axis(name = "Site C")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(color = "black", size = 28, family = "Liberation Sans", face = 'bold'),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        legend.background = element_rect(color = 'black'),
        axis.text.y.left  = element_text(size = 24, family = "Liberation Sans", face = 'bold', vjust = 0.45),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        axis.ticks.y.left = element_line(linewidth = 1.5),
        legend.key.spacing.y = unit(2, 'mm'),
        legend.position = 'right')
resave(hg_c.plot, file = './abridged.RData')

## Site D Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_d.ps <- subset_samples(final_ksp.ps, Site == "D")
hg_d.ps <- subset_taxa(hg_d.ps, taxa_sums(hg_d.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_d.name <- names(sort(taxa_sums(hg_d.ps), decreasing = TRUE))[1:19]
hg_d.name <- c(hg_d.name, "Other")
hg_d.colr <- ksp.colr[hg_d.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_d.ps <- merge_taxa(hg_d.ps, taxa_names(hg_d.ps)[!taxa_names(hg_d.ps) %in% hg_d.name])
taxa_names(hg_d.ps)[17] <- "Other"
tax_table(hg_d.ps)[17, "ASV"] <- "Other"

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
  scale_y_continuous(breaks = seq(0,1, 0.20), expand = expansion(mult = c(0.002, 0.02)), sec.axis = dup_axis(name = "Site D")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(color = "black", size = 28, family = "Liberation Sans", face = 'bold'),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        legend.background = element_rect(color = 'black'),
        axis.text.y.left  = element_text(size = 24, family = "Liberation Sans", face = 'bold', vjust = 0.45),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        axis.ticks.y.left = element_line(linewidth = 1.5),
        legend.key.spacing.y = unit(2, 'mm'),
        legend.position = 'right')
resave(hg_d.plot, file = './abridged.RData')

## Site E Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_e.ps <- subset_samples(final_ksp.ps, Site == "E")
hg_e.ps <- subset_taxa(hg_e.ps, taxa_sums(hg_e.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_e.name <- names(sort(taxa_sums(hg_e.ps), decreasing = TRUE))[1:19]
hg_e.name <- c(hg_e.name, "Other")
hg_e.colr <- ksp.colr[hg_e.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_e.ps <- merge_taxa(hg_e.ps, taxa_names(hg_e.ps)[!taxa_names(hg_e.ps) %in% hg_e.name])
taxa_names(hg_e.ps)[3] <- "Other"
tax_table(hg_e.ps)[3, "ASV"] <- "Other"

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
  scale_y_continuous(breaks = seq(0,1, 0.20), expand = expansion(mult = c(0.002, 0.02)), sec.axis = dup_axis(name = "Site E")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(color = "black", size = 28, family = "Liberation Sans", face = 'bold'),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        legend.background = element_rect(color = 'black'),
        axis.text.y.left  = element_text(size = 24, family = "Liberation Sans", face = 'bold', vjust = 0.45),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        axis.ticks.y.left = element_line(linewidth = 1.5),
        legend.key.spacing.y = unit(2, 'mm'),
        legend.position = 'right')
resave(hg_e.plot, file = './abridged.RData')

## Site F Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_f.ps <- subset_samples(final_ksp.ps, Site == "F")
hg_f.ps <- subset_taxa(hg_f.ps, taxa_sums(hg_f.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_f.name <- names(sort(taxa_sums(hg_f.ps), decreasing = TRUE))[1:19]
hg_f.name <- c(hg_f.name, "Other")
hg_f.colr <- ksp.colr[hg_f.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_f.ps <- merge_taxa(hg_f.ps, taxa_names(hg_f.ps)[!taxa_names(hg_f.ps) %in% hg_f.name])
taxa_names(hg_f.ps)[10] <- "Other"
tax_table(hg_f.ps)[10, "ASV"] <- "Other"

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
  scale_y_continuous(breaks = seq(0,1, 0.20), expand = expansion(mult = c(0.002, 0.02)), sec.axis = dup_axis(name = "Site F")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(color = "black", size = 28, family = "Liberation Sans", face = 'bold'),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        legend.background = element_rect(color = 'black'),
        axis.text.y.left  = element_text(size = 24, family = "Liberation Sans", face = 'bold', vjust = 0.45),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        axis.ticks.y.left = element_line(linewidth = 1.5),
        legend.key.spacing.y = unit(2, 'mm'),
        legend.position = 'right')
resave(hg_f.plot, file = './abridged.RData')

## Site G Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_g.ps <- subset_samples(final_ksp.ps, Site == "G")
hg_g.ps <- subset_taxa(hg_g.ps, taxa_sums(hg_g.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_g.name <- names(sort(taxa_sums(hg_g.ps), decreasing = TRUE))[1:19]
hg_g.name <- c(hg_g.name, "Other")
hg_g.colr <- ksp.colr[hg_g.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_g.ps <- merge_taxa(hg_g.ps, taxa_names(hg_g.ps)[!taxa_names(hg_g.ps) %in% hg_g.name])
taxa_names(hg_g.ps)[5] <- "Other"
tax_table(hg_g.ps)[5, "ASV"] <- "Other"

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
  scale_y_continuous(breaks = seq(0,1, 0.20), expand = expansion(mult = c(0.002, 0.02)), sec.axis = dup_axis(name = "Site G")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(color = "black", size = 28, family = "Liberation Sans", face = 'bold'),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        legend.background = element_rect(color = 'black'),
        axis.text.y.left  = element_text(size = 24, family = "Liberation Sans", face = 'bold', vjust = 0.45),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        axis.ticks.y.left = element_line(linewidth = 1.5),
        legend.key.spacing.y = unit(2, 'mm'),
        legend.position = 'right')
resave(hg_g.plot, file = './abridged.RData')

## Site H Samples ##
# Save a phyloseq Object that contains only the samples of interest #
hg_h.ps <- subset_samples(final_ksp.ps, Site == "H")
hg_h.ps <- subset_taxa(hg_h.ps, taxa_sums(hg_h.ps) > 0)

# Save the taxa names of the top 9 taxa and add "Other" to the end. #
hg_h.name <- names(sort(taxa_sums(hg_h.ps), decreasing = TRUE))[1:19]
hg_h.name <- c(hg_h.name, "Other")
hg_h.colr <- ksp.colr[hg_h.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_h.ps <- merge_taxa(hg_h.ps, taxa_names(hg_h.ps)[!taxa_names(hg_h.ps) %in% hg_h.name])
taxa_names(hg_h.ps)[15] <- "Other"
tax_table(hg_h.ps)[15, "ASV"] <- "Other"

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
  scale_y_continuous(breaks = seq(0,1, 0.20), expand = expansion(mult = c(0.002, 0.02)), sec.axis = dup_axis(name = "Site H")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(color = "black", size = 28, family = "Liberation Sans", face = 'bold'),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        legend.background = element_rect(color = 'black'),
        axis.text.y.left  = element_text(size = 24, family = "Liberation Sans", face = 'bold', vjust = 0.45),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        axis.ticks.y.left = element_line(linewidth = 1.5),
        legend.key.spacing.y = unit(2, 'mm'),
        legend.position = 'right')
resave(hg_h.plot, file = './abridged.RData')

## All Samples Across Sites Pooled ##
hg_ksp.colr <- ksp.colr[hg_myc.name,]

# Create a phyloseq object that saves the abundances of the top 9 ASVs and groups the remaining taxa into "Other" #
hg_ksp.ps <- merge_taxa(final_ksp.ps, taxa_names(final_ksp.ps)[!taxa_names(final_ksp.ps) %in% hg_myc.name])
hg_ksp.ps <- subset_samples(hg_ksp.ps, Treatment != "MycoBloom")
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
  scale_y_continuous(breaks = seq(0,1, 0.20), expand = expansion(mult = c(0.002, 0.02)), sec.axis = dup_axis(name = "MycoBloom ASVs")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(color = "black", size = 28, family = "Liberation Sans", face = 'bold'),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        legend.background = element_rect(color = 'black'),
        axis.text.y.left  = element_text(size = 24, family = "Liberation Sans", face = 'bold', vjust = 0.45),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        axis.ticks.y.left = element_line(linewidth = 1.5),
        legend.key.spacing.y = unit(2, 'mm'),
        legend.position = 'right')
resave(hg_ksp.plot, file = './abridged.RData')

#### Differential Abundance ####
if(!requireNamespace("Maaslin2")) BiocManager::install("biobakery/maaslin2")
library(Maaslin2); packageVersion("Maaslin2")

if(!dir.exists('./maas')){
  dir.create('./maas')
}

# Make a phyloseq object that excludes the MycoBloom Community #
nomy.ps <- subset_samples(final_ksp.ps, Treatment != "MycoBloom")
decompose_ps(nomy.ps, 'nomy')

# Save the taxa names of the MycoBloom Community #
myc.name <- taxa_names(myc.ps)

### Site A ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
a.maas <- Maaslin2(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/a.maas',
                   fixed_effects = "Both",
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Both,AC'),
                   plot_scatter = FALSE,
                   plot_heatmap = FALSE,)

# Save the results tables for abundance and prevalence, respectively #  
a_maas.res <- a.maas$results

# Save only the results of with values corresponding to communities found within Site A #
a_comp.res <- c()
for(i in 1:nrow(a_maas.res)){
  if(a_maas.res$value[i] == "AHI" | a_maas.res$value[i] == "ALO"){
    a_comp.res <- rbind(a_comp.res, a_maas.res[i,])
  }
}

# Save only the significant results #
a_sigs.res <- c()
for(i in 1:nrow(a_comp.res)){
  if(a_comp.res$qval[i] < 0.05){
    a_sigs.res <- rbind(a_sigs.res, a_comp.res[i,])
  }
}

# Save the significant results that correspond to taxa found in the mycobloom community #
a_mycs.res <- c()
for(i in 1:nrow(a_sigs.res)){
  if(a_sigs.res$feature[i] %in% myc.name){
    a_mycs.res <- rbind(a_mycs.res, a_sigs.res[i,])
  }
}

### Site B ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
b.maas <- Maaslin2(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/b.maas',
                   fixed_effects = "Both",
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Both,BC'),
                   plot_scatter = FALSE,
                   plot_heatmap = FALSE,)

# Save the results tables for abundance and prevalence, respectively #  
b_maas.res <- b.maas$results

# Save only the results of with values corresponding to communities found within Site A #
b_comp.res <- c()
for(i in 1:nrow(b_maas.res)){
  if(b_maas.res$value[i] == "BHI" | b_maas.res$value[i] == "BLO"){
    b_comp.res <- rbind(b_comp.res, b_maas.res[i,])
  }
}

# Save only the significant results #
b_sigs.res <- c()
for(i in 1:nrow(b_comp.res)){
  if(b_comp.res$qval[i] < 0.05){
    b_sigs.res <- rbind(b_sigs.res, b_comp.res[i,])
  }
}

# Save the significant results that correspond to taxa found in the mycobloom community #
b_mycs.res <- c()
for(i in 1:nrow(b_sigs.res)){
  if(b_sigs.res$feature[i] %in% myc.name){
    b_mycs.res <- rbind(b_mycs.res, b_sigs.res[i,])
  }
}

### Site C ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
c.maas <- Maaslin2(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/c.maas',
                   fixed_effects = "Both",
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Both,CC'),
                   plot_scatter = FALSE,
                   plot_heatmap = FALSE,)

# Save the results tables for abundance and prevalence, respectively #  
c_maas.res <- c.maas$results

# Save only the results of with values corresponding to communities found within Site A #
c_comp.res <- c()
for(i in 1:nrow(c_maas.res)){
  if(c_maas.res$value[i] == "CHI" | c_maas.res$value[i] == "CLO"){
    c_comp.res <- rbind(c_comp.res, c_maas.res[i,])
  }
}

# Save only the significant results #
c_sigs.res <- c()
for(i in 1:nrow(c_comp.res)){
  if(c_comp.res$qval[i] < 0.05){
    c_sigs.res <- rbind(c_sigs.res, c_comp.res[i,])
  }
}

# Save the significant results that correspond to taxa found in the mycobloom community #
c_mycs.res <- c()
for(i in 1:nrow(c_sigs.res)){
  if(c_sigs.res$feature[i] %in% myc.name){
    c_mycs.res <- rbind(c_mycs.res, c_sigs.res[i,])
  }
}

### Site D ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
d.maas <- Maaslin2(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/d.maas',
                   fixed_effects = "Both",
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Both,DC'),
                   plot_scatter = FALSE,
                   plot_heatmap = FALSE,)

# Save the results tables for abundance and prevalence, respectively #  
d_maas.res <- d.maas$results

# Save only the results of with values corresponding to communities found within Site A #
d_comp.res <- c()
for(i in 1:nrow(d_maas.res)){
  if(d_maas.res$value[i] == "DHI" | d_maas.res$value[i] == "DLO"){
    d_comp.res <- rbind(d_comp.res, d_maas.res[i,])
  }
}

# Save only the significant results #
d_sigs.res <- c()
for(i in 1:nrow(d_comp.res)){
  if(d_comp.res$qval[i] < 0.05){
    d_sigs.res <- rbind(d_sigs.res, d_comp.res[i,])
  }
}

# Save the significant results that correspond to taxa found in the mycobloom community #
d_mycs.res <- c()
for(i in 1:nrow(d_sigs.res)){
  if(d_sigs.res$feature[i] %in% myc.name){
    d_mycs.res <- rbind(d_mycs.res, d_sigs.res[i,])
  }
}

### Site E ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
e.maas <- Maaslin2(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/e.maas',
                   fixed_effects = "Both",
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Both,EC'),
                   plot_scatter = FALSE,
                   plot_heatmap = FALSE,)

# Save the results tables for abundance and prevalence, respectively #  
e_maas.res <- e.maas$results

# Save only the results of with values corresponding to communities found within Site A #
e_comp.res <- c()
for(i in 1:nrow(e_maas.res)){
  if(e_maas.res$value[i] == "EHI" | e_maas.res$value[i] == "ELO"){
    e_comp.res <- rbind(e_comp.res, e_maas.res[i,])
  }
}

# Save only the significant results #
e_sigs.res <- c()
for(i in 1:nrow(e_comp.res)){
  if(e_comp.res$qval[i] < 0.05){
    e_sigs.res <- rbind(e_sigs.res, e_comp.res[i,])
  }
}

# Save the significant results that correspond to taxa found in the mycobloom community #
e_mycs.res <- c()
for(i in 1:nrow(e_sigs.res)){
  if(e_sigs.res$feature[i] %in% myc.name){
    e_mycs.res <- rbind(e_mycs.res, e_sigs.res[i,])
  }
}

### Site F ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
f.maas <- Maaslin2(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/f.maas',
                   fixed_effects = "Both",
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Both,FC'),
                   plot_scatter = FALSE,
                   plot_heatmap = FALSE,)

# Save the results tables for abundance and prevalence, respectively #  
f_maas.res <- f.maas$results

# Save only the results of with values corresponding to communities found within Site A #
f_comp.res <- c()
for(i in 1:nrow(f_maas.res)){
  if(f_maas.res$value[i] == "FHI" | f_maas.res$value[i] == "FLO"){
    f_comp.res <- rbind(f_comp.res, f_maas.res[i,])
  }
}

# Save only the significant results #
f_sigs.res <- c()
for(i in 1:nrow(f_comp.res)){
  if(f_comp.res$qval[i] < 0.05){
    f_sigs.res <- rbind(f_sigs.res, f_comp.res[i,])
  }
}

# Save the significant results that correspond to taxa found in the mycobloom community #
f_mycs.res <- c()
for(i in 1:nrow(f_sigs.res)){
  if(f_sigs.res$feature[i] %in% myc.name){
    f_mycs.res <- rbind(f_mycs.res, f_sigs.res[i,])
  }
}


### Site G ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
g.maas <- Maaslin2(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/g.maas',
                   fixed_effects = "Both",
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Both,GC'),
                   plot_scatter = FALSE,
                   plot_heatmap = FALSE,)

# Save the results tables for abundance and prevalence, respectively #  
g_maas.res <- g.maas$results

# Save only the results of with values corresponding to communities found within Site A #
g_comp.res <- c()
for(i in 1:nrow(g_maas.res)){
  if(g_maas.res$value[i] == "GHI" | g_maas.res$value[i] == "GLO"){
    g_comp.res <- rbind(g_comp.res, g_maas.res[i,])
  }
}

# Save only the significant results #
g_sigs.res <- c()
for(i in 1:nrow(g_comp.res)){
  if(g_comp.res$qval[i] < 0.05){
    g_sigs.res <- rbind(g_sigs.res, g_comp.res[i,])
  }
}

# Save the significant results that correspond to taxa found in the mycobloom community #
g_mycs.res <- c()
for(i in 1:nrow(g_sigs.res)){
  if(g_sigs.res$feature[i] %in% myc.name){
    g_mycs.res <- rbind(g_mycs.res, g_sigs.res[i,])
  }
}

### Site H ###
# Calculate differential abundance using the mixed linear model with treatment as a fixed effect and site as a random effect #
h.maas <- Maaslin2(input_data = nomy$otu,
                   input_metadata = nomy$met,
                   output = './maas/h.maas',
                   fixed_effects = "Both",
                   normalization = "TSS",
                   transform = "LOG",
                   correction = "BH",
                   reference = c('Both,HC'),
                   plot_scatter = FALSE,
                   plot_heatmap = FALSE,)

# Save the results tables for abundance and prevalence, respectively #  
h_maas.res <- h.maas$results

# Save only the results of with values corresponding to communities found within Site A #
h_comp.res <- c()
for(i in 1:nrow(h_maas.res)){
  if(h_maas.res$value[i] == "HHI" | h_maas.res$value[i] == "HLO"){
    h_comp.res <- rbind(h_comp.res, h_maas.res[i,])
  }
}

# Save only the significant results #
h_sigs.res <- c()
for(i in 1:nrow(h_comp.res)){
  if(h_comp.res$qval[i] < 0.05){
    h_sigs.res <- rbind(h_sigs.res, h_comp.res[i,])
  }
}

# Save the significant results that correspond to taxa found in the mycobloom community #
h_mycs.res <- c()
for(i in 1:nrow(h_sigs.res)){
  if(h_sigs.res$feature[i] %in% myc.name){
    h_mycs.res <- rbind(h_mycs.res, h_sigs.res[i,])
  }
}


# All Samples Across Treatments Test #
all.maas <- Maaslin2(input_data = nomy$otu,
                               input_metadata = nomy$met,
                               output = './maas/all.maas',
                               fixed_effects = "Treats",
                               random_effects = "Sites",
                               normalization = "TSS",
                               transform = "LOG",
                               correction = "BH",
                               reference = c('Treats,Control'),
                               plot_scatter = FALSE,
                               plot_heatmap = FALSE,) 

# Save the results tables for abundance and prevalence, respectively #  
all_maas.res <- all.maas$results

# Change the names such that they match the original ASV names #
all_maas.res$feature <- sub("^(ASV[0-9]+).([^_]+).$", "\\1(\\2)", all_maas.res$feature)

# Save only the significant results #
all_sigs.res <- c()
for(i in 1:nrow(all_maas.res)){
  if(all_maas.res$qval[i] < 0.05){
    all_sigs.res <- rbind(all_sigs.res, all_maas.res[i,])
  }
}

# Save the significant results that correspond to taxa found in the mycobloom community #
all_mycs.res <- c()
for(i in 1:nrow(all_sigs.res)){
  if(all_sigs.res$feature[i] %in% myc.name){
    all_mycs.res <- rbind(all_mycs.res, all_sigs.res[i,])
  }
}

# Make sepearte data.frames for comparisons between high and low communities #
high_maas.res <- filter(all_maas.res, value == "High")
low_maas.res <- filter(all_maas.res, value == "Low")

# take the negative log10 of each qval #
high_maas.res$`-log10(qval)` <- sapply(high_maas.res$qval, function(x) -log(x, base = 10))
low_maas.res$`-log10(qval)` <- sapply(low_maas.res$qval, function(x) -log(x, base = 10))

# Show the fold change for each ASV by Treatment comparison #
for(i in 1:nrow(high_maas.res)){
  if(high_maas.res$coef[i] >= 0){
    high_maas.res$back[i] <- exp(high_maas.res$coef[i])
  } else{
    high_maas.res$back[i] <- -1 * exp(-1 * high_maas.res$coef[i]) 
  }
}

for(i in 1:nrow(low_maas.res)){
  if(low_maas.res$coef[i] >= 0){
    low_maas.res$back[i] <- exp(low_maas.res$coef[i])
  } else{
    low_maas.res$back[i] <- -1 * exp(-1 * low_maas.res$coef[i]) 
  }
}

# Add annotations for colors #
for(i in 1:nrow(high_maas.res)){
  if(high_maas.res$qval[i] < 0.05 & high_maas.res$coef[i] > 0){
    high_maas.res$Direction[i] <- 'High'
  } else if(high_maas.res$qval[i] < 0.05 & high_maas.res$coef[i] < 0){
    high_maas.res$Direction[i] <- 'Control'
  } else{
    high_maas.res$Direction[i] <- 'Neither'
  }
}

for(i in 1:nrow(low_maas.res)){
  if(low_maas.res$qval[i] < 0.05 & low_maas.res$coef[i] > 0){
    low_maas.res$Direction[i] <- 'Low'
  } else if(low_maas.res$qval[i] < 0.05 & low_maas.res$coef[i] < 0){
    low_maas.res$Direction[i] <- 'Control'
  } else{
    low_maas.res$Direction[i] <- 'Neither'
  }
}

# Add colors for samples that only correspond to taxa from the mycobloom community #
for(i in 1:nrow(high_maas.res)){
  if(high_maas.res$feature[i] %in% myc.name){
    high_maas.res$Myc[i] <- 'myc'
  }else{
    high_maas.res$Myc[i] <- 'not'
  }
}

for(i in 1:nrow(low_maas.res)){
  if(low_maas.res$feature[i] %in% myc.name){
    low_maas.res$Myc[i] <- 'myc'
  }else{
    low_maas.res$Myc[i] <- 'not'
  }
}

# Add final annotations so that only significant mycobloom taxa are represented #
for(i in 1:nrow(high_maas.res)){
  if(high_maas.res$Direction[i] == 'High' & high_maas.res$Myc[i] == 'myc'){
    high_maas.res$sigmyc[i] <- 'yes'
    high_maas.res$feats[i] <- high_maas.res$feature[i]
  }else {
    high_maas.res$sigmyc[i] <- 'no'
    high_maas.res$feats[i] <- ''
  }
}

# Add final annotations so that only significant mycobloom taxa are represented #
for(i in 1:nrow(low_maas.res)){
  if(low_maas.res$Direction[i] == 'Low' & low_maas.res$Myc[i] == 'myc'){
    low_maas.res$sigmyc[i] <- 'yes'
    low_maas.res$feats[i] <- low_maas.res$feature[i]
  }else {
    low_maas.res$sigmyc[i] <- 'no'
    low_maas.res$feats[i] <- ''
  }
}

high_maas.res$sigmyc <- factor(high_maas.res$sigmyc, levels = c('yes', 'no'))
low_maas.res$sigmyc <- factor(low_maas.res$sigmyc, levels = c('yes', 'no'))

# Make a volcano plot #
if(!requireNamespace("ggrepel")) install.packages("ggrepel")
library(ggrepel); packageVersion("ggrepel")
if(!requireNamespace('ggtext')) install.packages('ggtext')
library(ggtext); packageVersion('ggtext')

high.volc <- ggplot(high_maas.res, aes(x = coef, y = `-log10(qval)`, color = sigmyc, label = feats)) +
  geom_hline(yintercept = -log(0.05, 10), col = 'gray', linewidth = unit(1, 'mm'), linetype = 'dashed') +
  geom_vline(xintercept = 0, col = 'gray', linewidth = unit(1,'mm'), linetype = 'dashed') +
  geom_point(size = 4) +
  scale_color_manual(labels = c('myc', 'not'), values = c('red', 'gray')) +
  theme_prism() +
  geom_label_repel(size = 5, family = "Liberation Sans", fontface = "bold", min.segment.length = unit(0, 'lines'), box.padding = unit(0.5, 'lines')) +
  scale_x_continuous(name = 'log<sub>2</sub>-fold change', breaks = seq(-4,4,1), limits = c(-4, 4.1)) +
  scale_y_continuous(breaks = seq(0,3.5,0.5), expand = expansion(add = c(0,0)), limits = c(0,3.5)) +
  theme(legend.position = 'none',
        axis.title.x = ggtext::element_markdown(family = 'Liberation Sans', face = 'bold', vjust = 0.1, size = 24),
        axis.title.y = ggtext::element_markdown(family = 'Liberation Sans', face = 'bold', vjust = 0.1, size = 24),
        axis.text.x = element_text(size = 14, family = 'Liberation Sans'),
        axis.text.y = element_text(size = 14, family = 'Liberation Sans'))
                    
low.volc <- ggplot(low_maas.res, aes(x = coef, y = `-log10(qval)`, color = sigmyc, label = feats)) +
  geom_hline(yintercept = -log(0.05, 10), col = 'gray', linewidth = unit(1, 'mm'), linetype = 'dashed') +
  geom_vline(xintercept = 0, col = 'gray', linewidth = unit(1,'mm'), linetype = 'dashed') +
  geom_point(size = 4) +
  scale_color_manual(, values = c('gray', 'red')) +
  theme_prism() +
  geom_label_repel(size = 5, family = "Liberation Sans", fontface = "bold", min.segment.length = unit(0, 'lines'), box.padding = unit(0.5, 'lines')) +
  scale_x_continuous(name = 'log<sub>2</sub>-fold change', breaks = seq(-4,4,1), limits = c(-4.1, 4.1)) +
  scale_y_continuous(breaks = seq(0,2,0.5), expand = expansion(add = c(0,0)), limits = c(0,2.1)) +
  theme(legend.position = 'none',
        axis.title.x = ggtext::element_markdown(family = 'Liberation Sans', face = 'bold', vjust = 0.1, size = 24),
        axis.title.y = ggtext::element_markdown(family = 'Liberation Sans', face = 'bold', vjust = 0.1, size = 24),
        axis.text.x = element_text(size = 14, family = 'Liberation Sans'),
        axis.text.y = element_text(size = 14, family = 'Liberation Sans'))

#### Correlation with Aboveground Diversity ####
# Load in the aboveground diversity data #
raw_above.df <- read.csv2('./KSP24_count.csv', sep = ',')
raw_above.df <- t(raw_above.df)
colnames(raw_above.df) <- raw_above.df[1,]

# Replace the empty entries with 0s #
for(i in 1:nrow(raw_above.df)){
  for(j in 1:ncol(raw_above.df)){
    if(raw_above.df[i,j] == ""){
      raw_above.df[i,j] = 0 
    }
  }
}

filt_above.df <- as.data.frame(raw_above.df[-1,])

# Make the counts numeric #
filt_above.df[,5:ncol(filt_above.df)] <- lapply(filt_above.df[,5:ncol(filt_above.df)], as.numeric)

# Order the data.frame by Site #
filt_above.df <- filt_above.df[sort(rownames(filt_above.df), decreasing = FALSE), ]

# Make an "OTU" and "MET" table of the observations #
above.otu <- filt_above.df[,5:ncol(filt_above.df)]
above.met <- filt_above.df[,1:4]

above.dist <- vegdist(above.otu, method = 'bray')
above.dist <- as.data.frame(as.matrix(above.dist))
below.dist <- vegdist(t(final_ksp$otu), method = 'bray')
below.dist <- as.data.frame(as.matrix(below.dist))
below.dist <- below.dist[-73,-73]

dist_avg <- function(dist, lets, numbs) {
  # Function that takes the within group average of a distance matrix of each group specified by their sample names and creates a new distance matrix #  
  within <- list()
  between <- list()
  groups <- c()
  final <- matrix(nrow = nrow(dist), ncol = ncol(dist))
  for(i in LETTERS[1:lets]){
    for(j in numbs){
      unq <- paste0(i,j)
      groups <- c(groups, paste0(i,j))
      within.dist <- dist[grep(paste0(i,j,'.*'), rownames(dist), value = TRUE),grep(paste0(i,j,'.*'), colnames(dist), value = TRUE)]
      if(length(rownames(within.dist)) > 1){
        within[[paste0(i,j)]] <- mean(as.dist(within.dist))
      } else{
        within[[paste0(i,j)]] <- NA
      }
    }
  }
  final <- matrix(nrow = length(groups), ncol = length(groups))
  final <- as.data.frame(final)
  rownames(final) <- groups; colnames(final) <- groups
  for(i in groups){
    final[i,i] <- within[[i]]
  }
  group_comb <- combn(groups, 2)
  for(i in 1:ncol(group_comb)){
    g1 <- group_comb[1,i]
    g2 <- group_comb[2,i]
    between.dist <- dist[grep(paste0(g1, '.*'), rownames(dist), value = TRUE), grep(paste0(g2, '.*'), colnames(dist), value = TRUE)]
    final[g1,g2] <- mean(as.dist(between.dist))
  }
  final <- t(final)
  for(i in 1:nrow(final)){
    for(j in 1:ncol(final)){
      if(is.na(final[i,j])){
        final[i,j] <- final[j,i]
      }
    }
  }
  return(final)
}

final_below.dist <- dist_avg(below.dist, 8, c("C", "LO", "HI"))
final_above.dist <- dist_avg(above.dist, 8, c(".C", ".low", ".high"))

rownames(final_above.dist) <- rownames(final_below.dist); colnames(final_above.dist) <- colnames(final_below.dist)

beta.man <- mantel(as.matrix(final_above.dist), as.matrix(final_below.dist), method = 'spearman', permutations = 9999)

data.framify <- function(mat_df1, mat_df2){
  # function that converts a data.frame of a distance matrix into a dataframe in which all values are saved as one column of a data.frame #
  joined <- matrix(nrow = 1, ncol = 7)
  joined <- as.data.frame(joined)
  rowcount <- 1
  colnames(joined) <- c('Above', 'Below', 'Rowname', 'Colname', 'Comparison', 'Thorough', 'Di')
  for(i in 1:nrow(mat_df1)){
    for(j in 1:ncol(mat_df1)){
      if(!is.na(mat_df1[i,j])){
        joined[rowcount,1] <- mat_df1[i,j]
        joined[rowcount,2] <- mat_df2[i,j]
        joined[rowcount,3] <- rownames(mat_df1)[i]
        joined[rowcount,4] <- colnames(mat_df1)[j]
        if(rownames(mat_df1)[i] == colnames(mat_df1)[j]){
          joined[rowcount,5] = 'Within'
        } else{
          joined[rowcount,5] = 'Between'
        }
        joined[rowcount,6] <- paste0(substr(joined[rowcount,3],2,2), '-', substr(joined[rowcount,4],2,2))
        rownames(joined)[rowcount] <- paste0(joined[rowcount,3], '-', joined[rowcount,4])
        if(substr(joined[rowcount,3], 1,1) == substr(joined[rowcount,4],1,1) & substr(joined[rowcount,3],2,2) == substr(joined[rowcount,4],2,2)){
          joined[rowcount,7] <- 'WSWT'
        } else if(substr(joined[rowcount,3], 1,1) != substr(joined[rowcount,4],1,1) & substr(joined[rowcount,3],2,2) == substr(joined[rowcount,4],2,2)) {
          joined[rowcount,7] <- 'ASWT'
        } else if(substr(joined[rowcount,3], 1,1) == substr(joined[rowcount,4],1,1) & substr(joined[rowcount,3],2,2) != substr(joined[rowcount,4],2,2)){
          joined[rowcount,7] <- 'WSAT'
        } else{
          joined[rowcount,7] <- 'ASAT'
        }
        rowcount <- rowcount + 1
      }
    }
  }
  return(joined)
}

beta.df <- data.framify(final_above.dist, final_below.dist)

within_site <- function(df){
  # function that only saves the pairwise comparisons within site #
  within <- matrix(nrow = 1, ncol = 7)
  within = as.data.frame(within)
  colnames(within) <- colnames(df)
  rowcount <- 1
  for(i in 1:nrow(df)){
    if(substr(df$Rowname[i],1,1) == substr(df$Colname[i],1,1)){
      within[rowcount,] <- df[i,]
      rownames(within)[rowcount] <- paste0(df$Rowname[i], '-', df$Colname[i])
      rowcount <- rowcount + 1
    }
  }
  return(within)
}

beta_within.df <- within_site(beta.df)

beta.fit <- lm(Above~Below, beta.df)
beta.aov <- Anova(beta.fit, type = 2)

beta.df$Di <- factor(beta.df$Di, levels = c('WSWT', 'WSAT', 'ASWT', 'ASAT'))

beta_man.plot <- ggplot(beta.df, aes(x = Below, y = Above)) +
  geom_abline(intercept = coef(beta.fit)[1],
              slope = coef(beta.fit)[2],
              color = 'black') +
  geom_point(size = 5) +
  theme_prism() +
  scale_y_continuous(name = "Aboveground Biodiversity Dissimilarity (Bray-Curtis)", breaks = seq(0.2,1,0.2), limits = c(0.2,1)) +
  scale_x_continuous(name = "Belowground Biodiversity Dissimilarity (Bray-Curtis)", breaks = seq(0.2,1,0.2), limits = c(0.3,1)) +
  annotate("richtext", x = 0.65, y = 0.2,
           label = paste0("Mantel r = ", round(beta.man$statistic, 3), "; <em>p</em> = ", round(beta.man$signif, 3)),
           size = 14, fontface = 'bold', family = 'Liberation Sans') +
  theme(legend.text = element_text(size = 16, color = 'black', face = 'bold', family = 'Liberation Sans'),
        legend.title = element_text(size = 20, color = 'black', face = 'bold', family = 'Liberation Sans'),
        legend.key.spacing.y = unit(3, 'mm'),
        legend.position = "none",,
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))

save.image("./ksp_amf.RData")

# This just outputs the final time in which the Rscript finishes running #
cat("\n## Script finished at", Sys.time(), "\n")
sink(type = "message")
sink()
