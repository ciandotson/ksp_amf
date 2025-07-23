# The Effects of Locally-Cultivated Mycorrhizal Inoculant on Kankakee Sands Preserve AMF Community Diversity 

## Introduction

This repository contains a semi-automated workflow that analyzes the Arbuscular Mycorrhizal Fungal (AMF) communities of the restored sand prairies of Kankakee Sands Preserves (KSP). This pipeline was developed to analyze the data found in [this publication](https://www.youtube.com/watch?v=dQw4w9WgXcQ).

## How to Run This Pipeline 
In terms of what is required of the user, most of the work is in the frontend. In summary, the user must:  

- Clone the Repository  
- Download the Raw Reads from Sequence Read Archive (SRA)
- Initialize the Conda Environment
- Run the Pipeline.

Once the user performs these tasks, the R script `ksp_amf.R` will perform the whole analysis and output the reproduced results with the associated publication pertaining to the whole AMF community data. Below, the user may find a step-by-step guide to reproducing these results.

# 1. Cloning the Repository
Aside from the raw reads, which is explained in step 2, all of the data and scripts necessary to run this pipeline is found in this repository. To clone this repository, open your terminal and run the following commands:

```bash
cd
git clone https://github.com/ciandotson/ksp_amf.git
```

These commands will automatically download all of the source code directly into your local network. You can access the contents by running the command:

```bash
cd ~/ksp_amf
```

For the duration of the workflow, please make sure you are in this repository.

# 2. Downloading the Raw Reads from SRA
As the raw reads are too large to be stored in a Github repository, we have to download the reads from NCBI's Sequence Read Archive (SRA). This can be done by running the following code:

> Will update this once the reads are uploaded

Though not required, I would recommend moving the reads into the current directory. You can do this using the commands:

```bash
mkdir ~/ksp_amf/reads
mv -r ~/path/to/your/reads ~/ksp_amf/reads
```

This makes step 4 a bit easier, as you can then just copy and paste the final command directly into your terminal directly as opposed to having to adjust the filepath to where the reads are in your environment.

# 3. Initializing the Conda Environment
Though the vast majority of bioinformatic tools used in this pipeline are automatically downloaded in R if needed, there are a few bash packages (namely cutadapt, mafft, and iqtree), that need to be installed prior to running the pipeline. The easiest way to do this is by creating a conda environment that has these packages installed. If you do not have conda installed, a tutorial is found [here](https://www.anaconda.com/docs/getting-started/miniconda/install). Once installed, you can run the following commands to both make an environment and download the necessary packages to the environment.

```bash
conda create -n PSF
conda activate PSF
conda install cutadapt
conda install mafft
conda install iqtree
```

With this, all of the necessary data and tools are downloaded!
# 4. Running the Pipeline
To use the downloaded tools and data, all you have to do is run the following command from the repository:

```bash
Rscript ksp_amf.R --forward ~/ksp_amf/reads/dir_name | cat > ksp_amf.log
```

If you had not previously moved the files to this repository, change `~/ksp_amf/reads/dir_name` to where the raw forward reads are located in your local environment.

This command calls the main R script `ksp_amf.R` to take as input the file paths of the raw reads and output a series of files that pertain to the results of the initial publication, including:
- `./reads`: a directory that contains all of the read data, including reads with adapter sequences removed, reads processed by `dada2`, and alignment and tree files.
- `./ksp_amf.RData`: an R environment that contains the saved image of the entire working environment of the pipeline.
- `./abridged.RData`: an R environment that is a subset of `ksp_amf.RData`, which saves the most important objects such as `phyloseq` objects, `ggplot2` plots, and outputs from major statistical tests. This is done because the entire R environment is too large for most laptops to process, so if the user wants to rerun portions of the code without HPC, this environment can be used.
