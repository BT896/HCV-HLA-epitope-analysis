# HCV-HLA-epitope-analysis
# Overview

This repository contains R code for visualising HLA-associated epitope hits across longitudinally sampled hepatitis C virus (HCV) isolates.

The script illustrates the distribution and persistence of detected HLA-associated epitope regions across serial samples from two hosts (P102 and P110).

The analysis generates a genome-wide visualisation of detected epitope regions and exports observed epitope amino-acid sequences in CSV format.

The example datasets correspond to hosts included in the accompanying manuscript.

# Repository structure

HCV-HLA-epitope-analysis/
├── README.md
├── scripts/
│   └── plot_and_export_epitopes.R
├── data/
│   ├── epitope_hits_with_coordinates_P102.txt
│   └── epitope_hits_with_coordinates_P110.txt
└── results/

# Requirements

The analysis requires R and the following packages:
ggplot2
dplyr
tidyr
patchwork
stringr

Additional packages may be required depending on the final version of the script.

Install the packages using: 
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", "stringr"))

# Input files

The script uses two tab-delimited input files:
epitope_hits_with_coordinates_P102.txt
epitope_hits_with_coordinates_P110.txt

Each file contains detected HLA-associated epitope hits and their corresponding HCV amino-acid coordinates.

The script can be run once the appropriate input files have been placed in the data/ directory.

The required columns are:
Fasta_Header (Sample identifier containing the host and collection date)
Epitope_Hit (Detected epitope amino-acid sequence)
HLA_Association (Associated HLA serotype)
Coordinate1 (Epitope start coordinate)
Coordinate2 (Epitope end coordinate)

Dates are extracted from the sample identifiers in DDMMYY format.

# Running the analysis

Download or clone the repository and open R or RStudio.

Set the working directory to the repository's root folder, then run:
source("scripts/plot_and_export_epitopes.R")

Alternatively, from a terminal:
Rscript scripts/plot_and_export_epitopes.R

# Outputs

The script produces two files in the results/ directory.

1. Genome-wide epitope visualisation

P102_P110_HLA_epitope_map.png

A genome-wide representation of HLA-associated epitope hits across serial HCV isolates. Epitope regions are displayed according to their amino-acid coordinates and sample collection dates.

2. Epitope amino-acid sequence table

HLA_variable_epitopes_P102_P110.csv

A CSV file containing epitope coordinates and observed amino-acid sequences arranged by sample.

The CSV can be used for subsequent alignment-based analyses.

Important: Blank sequence cells indicate that no epitope hit was recorded for that sample and region. They do not necessarily indicate missing sequence data.

# Interpretation

The script evaluates variability in detected epitope-hit presence across longitudinal samples. This should not be interpreted as a direct measure of amino-acid sequence conservation.

For regions without detected epitope hits, complete amino-acid sequences may need to be obtained separately from the corresponding consensus sequence alignments.

# Publication

This code accompanies the (currently) draft manuscript:

["Molecular epidemiology reveals ongoing HCV exposure following acute infection: deep-sequencing evidence from high-risk UK cohorts"]

Licence

[The R code in this repository is released under the MIT Licence. See the LICENSE file for details.]
