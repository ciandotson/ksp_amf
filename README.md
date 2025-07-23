# ksp_amf

## Introduction

This repository contains a semi-automated workflow that analyzes the Arbuscular Mycorrhizal Fungal (AMF) communities of the restored sand prairies of Kankakee Sands Preserves (KSP). This pipeline was developed to analyze the data found in [this publication](https://www.youtube.com/watch?v=dQw4w9WgXcQ).

## How to Run This Pipeline 
In terms of what is required of the user, most of the work is in the frontend. In summary, the user must:  

- Clone the repository  
- Download the raw reads from Sequence Read Archive (SRA)
- Initialize the Conda Environment
- Run the initialization command.

Once the user performs these tasks, the R script `ksp_amf.R`, will perform the whole analysis and output the reproduced results with the associated publication pertaining to the whole AMF community data. Below, the user may find a step-by-step guide to reproducing these results.

# 1. Cloning the Repository
Aside from the raw reads, which is explained in step 2, all of the data and scripts necessary to run this pipeline is found in this repository. To clone this repository, open your terminal and run the following commands:

  `cd`
  
  `git clone https://github.com/ciandotson/ksp_amf.git`

