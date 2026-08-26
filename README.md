# Nematode Macroevolutionary Mega-Phylogeny Pipeline 🪱🧬

![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-brightgreen.svg)

## Overview
This repository contains the datasets and advanced computational R scripts used to analyze the largest and most densely sampled macroevolutionary framework of the phylum Nematoda to date (3,654 unique species). 

Our pipeline transitions from stringent bioinformatic curation to rigorous hypothesis testing. It is specifically designed to perform **Ancestral State Reconstruction (ASR)**, test the strict irreversibility of **Dollo's Law**, and evaluate deep evolutionary transitions between free-living and parasitic lifestyles while heavily mitigating methodological artifacts like overfitting and taxon-sampling bias.

## Key Features (Reviewer-Proof Validation)
The core analysis script incorporates a multi-layered robustness validation module:
- **Formal Model Selection:** AICc-based comparison among Equal Rates (ER), All-Rates-Different (ARD), and Dollo-like Irreversible models.
- **Hidden State Avoidance Check:** Comparative evaluation with `corHMM` to prevent parameter inflation and overfitting.
- **Parallelized Jackknife Cross-Validation:** 100-iteration stress test (dropping 10% of taxa per run) to ensure the stability of the root state and transition rates (q01/q10).
- **Dynamic Stratified Root Sensitivity:** Automated alteration of the most ancient topological nodes to test the depth of the ancestral signal.
- **Phylogenetic Signal Calculation:** Evaluating trait conservatism utilizing Fritz and Purvis' D-statistic (`caper`).

## Repository Structure
```text
├── data/
│   ├── Pruned_Tree.nwk              # Maximum-Likelihood Mega-phylogeny (3,654 taxa)
│   └── Matched_Trait_Data.csv       # Ecological traits (Lifestyle, Habitat)
├── scripts/
│   └── MacroEvol_Analysis_V11.R     # Main macroevolutionary R pipeline
├── results/                         # Generated automatically after running the script
│   ├── Tables/                      # AICc, Jackknife raw data, Node probabilities
│   └── Figures/                     # Jackknife density, Rate stability, PieChart Tree
└── README.md
Prerequisites & Installation
The pipeline is written in R. To run the analyses, you will need the following dependencies. You can install them directly in your R console:

R
install.packages(c("ape", "phytools", "phangorn", "caper", "corHMM", "parallel", "ggplot2"))
Usage
Clone this repository to your local machine:

Bash
git clone [https://github.com/YourUsername/Nema-MegaPhylo.git](https://github.com/YourUsername/Nema-MegaPhylo.git)
Set your working directory in R to the cloned repository folder.

Source or run the main pipeline script:

R
source("scripts/MacroEvol_Analysis_V11.R")
(Note: The script utilizes mclapply for Unix/Mac and a PSOCK cluster via parLapply for Windows to ensure seamless cross-platform parallel processing during the 100x Jackknife permutations).

Expected Outputs
Running the pipeline will automatically generate a highly detailed Robustness_Validation_Report_v11.txt alongside several publication-ready outputs:

ASR_PieChart_Tree.pdf: High-resolution phylogeny mapping ancestral states across all nodes.

Jackknife_Density_Plot.pdf: Visual distribution of root probabilities across 100 permutations.

Transition_Rates_Stability.pdf: Line plots demonstrating the stability of transition rates (q01 vs q10).

Various CSV files containing trait distributions, AICc scores, and raw node probabilities.

Citation
If you use these scripts or datasets in your research, please cite our upcoming paper:

Noshadi, A., et al. (2026). [Insert Exact Title of the Paper Here]. Journal Name. DOI: [To be added]

License
This project is licensed under the MIT License - see the LICENSE file for details.
