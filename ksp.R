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

library(dada2)

if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
library(phyloseq); packageVersion("phyloseq")

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr); packageVersion('dplyr')
load("ksp_amf.RData")
# Here we actually make the phyloseq object (raw_ksp.ps) that contains the ASV table (nochim_ksp.st), taxonomy table (ksp.taxa), #
# and metadata table (ksp.met) #
nochim_ksp.st <- t(nochim_ksp.st)
raw_ksp.ps <- phyloseq(otu_table(nochim_ksp.st, taxa_are_rows = TRUE),
                       tax_table(ksp.taxa),
                       sample_data(ksp.met))

resave(raw_ksp.ps, file = "./abridged.RData")

# Filter out samples that did not have any matches to the MaarJAM database #
raw_ksp.ps <- tax_table(subset_taxa(ksp.ps, !is.na(Family)))

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
decompose_ps(myc.ps, "myc")

# Each component can be accessed by typing the name of the decomposed phyloseq object ("myc") followed by an accessor ("$") and the three letter abbreviation of the object #
# ("tax" for taxonomy table, "otu" for asv table, "met" for metadata table, "dna" for DNAStringSet, and "fra" for the Frankenstein table) #
myc$fra

resave(myc, file = "./abridged.RData")


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
ksp.hits <- matrix(nrow = nrow(ksp$tax), ncol = 12)
ksp.hits <- as.data.frame(ksp.hits) 
hold <- c()
for(i in 1:length(ksp$dna)){
  hold <- predict(blast.db, ksp$dna[i])
  ksp.hits[i,] <- hold[1,]
  ksp$tax$Best_Hit[i] <- hold[1, 2]
}

# Filter out reads that do not correspond to a NCBI entry #
filt_ksp.tax <- filter(ksp$tax, !is.na(ksp$tax$Best_Hit))

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
decompose_ps(final_ksp.ps, final_ksp)
out.dna <- readDNAStringSet('./reference/outgroup.fasta')
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
ksp.tre$tip.label <- sub("^(ASV[0-9]+)_([^_]+)_$", "\\1(\\2)", ksp_fun.tre$tip.label)

# Save the tree into the phyloseq object #
phy_tree(final_ksp.ps) <- ksp.tre

# Now we can actually do the filtering using the tree by first labelling our outgroups #
out.nam <- c("Outgroup1", "Outgroup2")

# Find the most recent common ancestor node for the two outgroups #
out.mrca <- getMRCA(ksp.tre, out.nam)

# Find all of the descendants of this MRCA #
out.des <- getDescendants(ksp.tre, out.mrca)
out.tip <- phy_tree(ksp.tre$tip.label[out.des[out.des <= length(ksp.tre$tip.label)]])

# Finally, remove the taxa denoted in the tips denoted in out.tip from the phyloseq object #
final_ksp.ps <- subset_taxa(final_ksp.ps, !taxa_names(final_ksp.ps) %in% out.tip)
                    
resave(final_ksp.ps, file = "./abridged.RData") 
save.image("ksp_amf.RData")

#### Alpha Diversity Measurement and Visualization ####
ksp.rich <- estimate_richness(final_ksp.ps) # automatically performs alpha diversity calculations #
ksp.rich <- cbind(ksp.met, ksp.rich)
View(ksp.rich)

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
ksp_rich.mnsd <- filter(ksp_rich.mnsd, Treatment != 'MycoBloom')

# Individual Site  Diversity # 
a.rich <- filter(ksp.rich, Site == 'A')
b.rich <- filter(ksp.rich, Site == 'B')
c.rich <- filter(ksp.rich, Site == 'C')
d.rich <- filter(ksp.rich, Site == 'D')
e.rich <- filter(ksp.rich, Site == 'E')
f.rich <- filter(ksp.rich, Site == 'F')
g.rich <- filter(ksp.rich, Site == 'G')
h.rich <- filter(ksp.rich, Site == 'H')

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
d_sha.hsd # Low vs. control #
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
summary(a.sha)
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

ggplot(plot.rich, aes(x = Site, y = sha.mean, fill = Treatment, color = Treatment)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  theme_prism() +
  ylab('Shannon Diversity') +
  scale_fill_manual(values = c("white", "gray", "#4D4D4D", 'black')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
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
d_sim.hsd # Low vs. control #
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
plot.rich <- cbind(ksp_rich.mnsd, sim.let)

ggplot(plot.rich, aes(x = Site, y = evn.mean, fill = Treatment, color = Treatment)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  theme_prism() +
  ylab("Simpson's Diversity") +
  scale_fill_manual(values = c("white", "gray", "#4D4D4D", 'black')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  ggtitle("Simpson's Diversity") +
  geom_text(aes(label = sim.let, y = evn.mean + evn.sd + 0.1), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 12) +
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
plot.rich <- cbind(ksp_rich.mnsd, cha.let)

ggplot(plot.rich, aes(x = Site, y = cha.mean, fill = Treatment, color = Treatment)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  theme_prism() +
  ylab("ASV Richness") +
  scale_fill_manual(values = c("white", "gray", "#4D4D4D", 'black')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  ggtitle("ksp Sample Chao1 Diversity") +
  geom_text(aes(label = cha.let, y = cha.mean + cha.sd + 0.1), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 12) +
  geom_errorbar(aes(ymin = cha.mean - cha.sd, ymax = cha.mean + cha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2)

#### Beta Diversity ####
ksp_prop.ps <- transform_sample_counts(ksp.ps, function(otu) otu/sum(otu))
ord.pcoa.wuni <- ordinate(ksp_prop.ps, method="PCoA", distance="bray")
ksp.bdist <- phyloseq::distance(ksp.ps, method = "bray", weighted = F)

if(!requireNamespace("vegan")) install.package("vegan")
library(vegan); packageVersion('vegan')
ksp.perm <- adonis2(ksp.bdist~Site*Treatment, data = ksp.met)
ksp.perm

plot_ordination(ksp_prop.ps, ord.pcoa.wuni, color="Site", shape = 'Treatment', title="ksp Sample PCoA") +
  theme_prism() +
  geom_point(size = 6) 

#### Per Sample Analysis ####
# A #
a.ps <- subset_samples(ksp.ps, Site == 'A')
a.met <- as(sample_data(a.ps), 'data.frame')
a_prop.ps <- transform_sample_counts(a.ps, function(otu) otu/sum(otu))
a_ord.pcoa.wuni <- ordinate(a_prop.ps, method="PCoA", distance="bray")
a.bdist <- phyloseq::distance(a.ps, method = "bray", weighted = F)

a.perm <- adonis2(a.bdist~Treatment, data = a.met)
a.perm

plot_ordination(a_prop.ps, a_ord.pcoa.wuni, color="Treatment", title="A Samples PCoA") +
  theme_prism() +
  geom_point(size = 6) +
  annotate(geom = 'text',x = 0.1, y = -0.45, label = 'P-value = 0.527', size = 8)

# B #
b.ps <- subset_samples(ksp.ps, Site == 'B')
b.met <- as(sample_data(b.ps), 'data.frame')
b_prop.ps <- transform_sample_counts(b.ps, function(otu) otu/sum(otu))
b_ord.pcoa.wuni <- ordinate(b_prop.ps, method="PCoA", distance="bray")
b.bdist <- phyloseq::distance(b.ps, method = "bray", weighted = F)

b.perm <- adonis2(b.bdist~Treatment, data = b.met)
b.perm

plot_ordination(b_prop.ps, b_ord.pcoa.wuni, color="Treatment", title="B Samples PCoA") +
  theme_prism() +
  geom_point(size = 6) +
  annotate(geom = 'text',x = -0.2, y = -0.4, label = 'P-value = 0.31', size = 8)

# C #
c.ps <- subset_samples(ksp.ps, Site == 'C')
c.met <- as(sample_data(c.ps), 'data.frame')
c_prop.ps <- transform_sample_counts(c.ps, function(otu) otu/sum(otu))
c_ord.pcoa.wuni <- ordinate(c_prop.ps, method="PCoA", distance="bray")
c.bdist <- phyloseq::distance(c.ps, method = "bray", weighted = F)

c.perm <- adonis2(c.bdist~Treatment, data = c.met)
c.perm

plot_ordination(c_prop.ps, c_ord.pcoa.wuni, color="Treatment", title="C Samples PCoA") +
  theme_prism() +
  geom_point(size = 6) +
  annotate(geom = 'text',x = -0.2, y = -0.4, label = 'P-value = 0.267', size = 8)

# D #
d.ps <- subset_samples(ksp.ps, Site == 'D')
d.met <- as(sample_data(d.ps), 'data.frame')
d_prop.ps <- transform_sample_counts(d.ps, function(otu) otu/sum(otu))
d_ord.pcoa.wuni <- ordinate(d_prop.ps, method="PCoA", distance="bray")
d.bdist <- phyloseq::distance(d.ps, method = "bray", weighted = F)

d.perm <- adonis2(d.bdist~Treatment, data = d.met)
d.perm

plot_ordination(d_prop.ps, d_ord.pcoa.wuni, color="Treatment", title="D Samples PCoA") +
  theme_prism() +
  geom_point(size = 6) +
  annotate(geom = 'text',x = 0.0, y = -0.4, label = 'P-value = 0.228', size = 8)

# E #
e.ps <- subset_samples(ksp.ps, Site == 'E')
e.met <- as(sample_data(e.ps), 'data.frame')
e_prop.ps <- transform_sample_counts(e.ps, function(otu) otu/sum(otu))
e_ord.pcoa.wuni <- ordinate(e_prop.ps, method="PCoA", distance="bray")
e.bdist <- phyloseq::distance(e.ps, method = "bray", weighted = F)

e.perm <- adonis2(e.bdist~Treatment, data = e.met)
e.perm

plot_ordination(e_prop.ps, e_ord.pcoa.wuni, color="Treatment", title="E Samples PCoA") +
  theme_prism() +
  geom_point(size = 6) +
  annotate(geom = 'text',x = 0.0, y = -0.3, label = 'P-value = 0.235', size = 8)

# F #
f.ps <- subset_samples(ksp.ps, Site == 'F')
f.met <- as(sample_data(f.ps), 'data.frame')
f_prop.ps <- transform_sample_counts(f.ps, function(otu) otu/sum(otu))
f_ord.pcoa.wuni <- ordinate(f_prop.ps, method="PCoA", distance="bray")
f.bdist <- phyloseq::distance(f.ps, method = "bray", weighted = F)

f.perm <- adonis2(f.bdist~Treatment, data = f.met)
f.perm

plot_ordination(f_prop.ps, f_ord.pcoa.wuni, color="Treatment", title="F Samples PCoA") +
  theme_prism() +
  geom_point(size = 6) +
  annotate(geom = 'text',x = 0.0, y = -0.3, label = 'P-value = 0.525', size = 8)

# G #
g.ps <- subset_samples(ksp.ps, Site == 'G')
g.met <- as(sample_data(g.ps), 'data.frame')
g_prop.ps <- transform_sample_counts(g.ps, function(otu) otu/sum(otu))
g_ord.pcoa.wuni <- ordinate(g_prop.ps, method="PCoA", distance="bray")
g.bdist <- phyloseq::distance(g.ps, method = "bray", weighted = F)

g.perm <- adonis2(g.bdist~Treatment, data = g.met)
g.perm

plot_ordination(g_prop.ps, g_ord.pcoa.wuni, color="Treatment", title="G Samples PCoA") +
  theme_prism() +
  geom_point(size = 6) +
  annotate(geom = 'text',x = -0.1, y = -0.5, label = 'P-value = 0.536', size = 8)

# H #
h.ps <- subset_samples(ksp.ps, Site == 'H')
h.met <- as(sample_data(h.ps), 'data.frame')
h_prop.ps <- transform_sample_counts(h.ps, function(otu) otu/sum(otu))
h_ord.pcoa.wuni <- ordinate(h_prop.ps, method="PCoA", distance="bray")
h.bdist <- phyloseq::distance(h.ps, method = "bray", weighted = F)

h.perm <- adonis2(h.bdist~Treatment, data = h.met)
h.perm

plot_ordination(h_prop.ps, h_ord.pcoa.wuni, color="Treatment", title="H Samples PCoA") +
  theme_prism() +
  geom_point(size = 6) +
  annotate(geom = 'text',x = -0.1, y = -0.45, label = 'P-value = 0.514', size = 8)

#### Stacked Histograms ####
# First we start by making a color pallette for each unique ASV #
if(!requireNamespace("Polychrome")) install.packages("Polychrome")
library(Polychrome); packageVersion("Polychrome")
ksp.colr <- createPalette(ntaxa(final_ksp.ps),  c("#ff0000", "#00ff00", "#0000ff"))
ksp.colr <- as.data.frame(ksp.colr)
rownames(soil_nod.colr) <- taxa_names(final_ksp.ps)

# Add a gray color for "Other"
soil_nod.colr[118,] <- "#D4D4D4" 
rownames(soil_nod.colr)[118] <- "Other" 

# Save 'ASV' as it's own unique level of taxonomy #
tax_table(final_ksp.ps)$ASV <- taxa_names(final_ksp.ps)

# Construct a phyloseq object and data frame object that contains just the desired samples #
hg_myc.ps <- subset_samples(final_ksp.ps, Treatment == "MycoBloom1")

# Collect only the top 19 ASVs, and group the remaining ASVs into "Other" and save their names in their orders of abundance to make a subsetted palette for plotting #
hg_myc.ps <- aggregate_top_taxa2(hg_myc.ps, 19, "ASV")
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