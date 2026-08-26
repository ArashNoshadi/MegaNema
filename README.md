
# Nematoda Megaphylogeny & Macroevolutionary Pipeline 🪱

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)](https://www.r-project.org/)
[![Status](https://img.shields.io/badge/Status-Active-success)]()

> A comprehensive computational pipeline for resolving the deep evolutionary history, ecological transitions, and genomic paradoxes of the phylum Nematoda using a 3,654-taxa ribosomal supermatrix.

## 📖 Overview

This repository contains the data, scripts, and analytical pipelines used in our study investigating the macroevolutionary trajectory of Nematoda. By exhaustively curating ~6 million public sequence records and assembling contiguous ribosomal operons (18S, ITS1, 5.8S, ITS2, 28S), we constructed a high-fidelity 3,654-taxa megaphylogeny. 

Our computational framework transitions from stringent phylogenomic curation to rigorous hypothesis testing, specifically addressing:
1. **Dollo's Law & Parasitism:** Challenging the strict irreversibility of parasitic transitions.
2. **Early Burst Diversification:** Quantifying whether the basal emergence of parasitism catalyzed rapid adaptive radiations.
3. **Genome Size Paradox:** Evaluating the decoupling of structural genomic constraints and evolutionary rates across distinct trophic strategies.

## 🚀 Key Analytical Modules

The analytical pipeline is built primarily in R and encompasses several advanced macroevolutionary models:

* **Formal Model Selection (AICc):** Systematic comparison of Equal Rates (ER), All-Rates-Different (ARD), and Dollo-like Irreversible models.
* **Hidden State Avoidance Check:** Evaluation of `corHMM` models to prevent parameter inflation and overfitting.
* **Robustness & Validation (Reviewer-Proofing):**
  * **Parallel Jackknife Cross-Validation:** 100 iterations with 10% random taxon dropout.
  * **Dynamic Stratified Root Sensitivity:** Stress-testing basal clade topology.
  * **Phylogenetic Signal (D-Statistic):** Measuring trait conservatism (Fritz and Purvis' D).
* **Diversification Modeling:** Pybus & Harvey’s Gamma ($\gamma$), Colless/Sackin indices, and LTT trajectory analysis.
* **Continuous Trait Evolution (PGLS):** Testing Ornstein-Uhlenbeck vs. Brownian motion models for genome architecture.

## 📂 Repository Structure

```text
├── data/
│   ├── Pruned_Tree.nwk                  # The 3,654-taxa ML tree (IQ-TREE output)
│   ├── Matched_Trait_Data.csv           # Ecological & genomic trait metadata
│   └── Genomes_Size_DB.csv              # Annotated genome sizes (Mb)
├── scripts/
│   ├── 01_Operon_Assembly_Curator.py    # De novo assembly & NUMT/Chimera removal
│   ├── 02_Macroevolutionary_ASR_v11.R   # Core ASR, Jackknifing, and Validation pipeline
│   └── 03_PGLS_Diversification.R        # Trait evolution and Gamma statistic tests
├── results/
│   ├── reports/                         # Automated textual reports (e.g., Validation_Report_v11.txt)
│   ├── figures/                         # High-res PDFs (PieChart trees, Density plots)
│   └── tables/                          # CSV outputs (Model selection, Node probabilities)
└── README.md

```

## ⚙️ Prerequisites & Installation

### R Dependencies

Ensure you have R (>= 4.0.0) installed. The macroevolutionary scripts require the following packages:

```R
install.packages(c("ape", "phytools", "phangorn", "caper", "corHMM", "ggplot2"))

```

### Upstream Software

* **IQ-TREE 2:** For maximum-likelihood tree inference and ultrafast bootstrapping.
* **Python 3.8+:** For the initial sequence parsing and curation pipeline.

## 💻 Usage

### 1. Data Preparation

Place your raw Newick tree and CSV metadata into the `data/` directory. Ensure taxonomic names match exactly between the tree tips and the dataframe row names.

### 2. Running the Core Validation Pipeline

Execute the main R script to reproduce the Ancestral State Reconstruction (ASR) and validation metrics:

```bash
Rscript scripts/02_Macroevolutionary_ASR_v11.R

```

**Expected Outputs (generated in the working directory):**

* `Model_Selection.csv`: AICc comparison scores.
* `Jackknife_Results.csv` & `Jackknife_Density_Plot.pdf`: Robustness distribution.
* `Transition_Rates_Stability.pdf`: Fluctuation tracking of $q_{01}$ and $q_{10}$ rates.
* `ASR_PieChart_Tree.pdf`: Full topology with annotated ancestral states.
* `Robustness_Validation_Report_v11.txt`: Comprehensive text summary of the statistical defense.

## 📊 Interpreting the Output

The pipeline automatically generates a `Robustness_Validation_Report_v11.txt` which interprets the statistical outputs. For example, a severe AICc penalty for the Irreversible model (+248 ΔAICc) mathematically refutes Dollo's law for the dataset, supporting the biological reality of secondary reversions.

## 🤝 Contributing & Collaborations

This project is part of an ongoing effort to map the deep evolutionary history of Nematoda. Future extensions (e.g., Sequential Bayesian Dating integrating whole-genome scaffolds) are actively being developed. If you are interested in collaborating or extending this framework, please feel free to open an issue or submit a pull request.

## ✉️ Contact

**Arash Noshadi**

MSc Researcher in Agricultural Biotechnology, University of Tehran

Email: [Your_Email@example.com]

LinkedIn: [Your_LinkedIn_Profile]

## 📜 Citation

If you use this pipeline or data in your research, please cite our upcoming manuscript:

> Noshadi, A., et al. (In Prep). *Dynamic orchestration of macroevolutionary transitions in Nematoda: Re-evaluating Dollo's Law and the Genome Size Paradox.*

---

*Built with rigorous statistical screening to ensure biologically meaningful macroevolutionary inferences.*

```

```
