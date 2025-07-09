# READ ME

<a href="https://creativecommons.org">Untitled</a> © 1999 by <a href="https://creativecommons.org">Jane Doe</a> is licensed under <a href="https://creativecommons.org/licenses/by/4.0/">CC BY 4.0</a><img src="https://mirrors.creativecommons.org/presskit/icons/cc.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/by.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;">


# Ohnologues vs BioMart: A Comparative Analysis of Gene Paralogues

## 📘 Overview

This project compares paralogue gene sets identified via **BioMart** and **ohnologues** (resulting from ancient whole-genome duplications) using data-driven visualizations and gene similarity thresholds. It investigates how these paralogues overlap, how they relate to essentiality (via the **FUSIL** dataset), and how similarity thresholds influence paralogue detection.

## 📁 Data Sources

- **BioMart Paralogues**  
  Retrieved from Ensembl BioMart. Contains pairwise paralogue relationships and percentage identity.

- **Ohnologues (Relaxed 2R Pairs)**  
  Downloaded from [ohnologs.org](http://ohnologs.org). Represents genes derived from ancient whole-genome duplications.

- **Protein-Coding Gene List**  
  Used to filter valid gene symbols. Source: local Ensembl data file.

- **FUSIL Dataset**  
  Provides gene essentiality classifications in mammals.

## 🧪 Methodology

1. **Data Loading**
   - Loads all required datasets: FUSIL, BioMart paralogues, ohnologues, and the protein-coding gene list.

2. **Preprocessing**
   - Cleans column names, filters by gene symbols, handles missing values.
   - Merges paralogue data with gene essentiality (FUSIL).

3. **Analysis**
   - Compares BioMart vs Ohnologues across similarity thresholds (30%, 50%, 70%).
   - Generates Venn diagrams to show gene overlap.
   - Creates scatter plots to compare paralogue counts per gene.

4. **Visualization**
   - Density plots show the distribution of % similarity across paralogue groups:
     - In BioMart only
     - In Ohnologues only
     - In both
     - In neither

## 📊 Output

The notebook generates:
- Venn diagrams of gene overlap per threshold
- Scatter plots of paralogue counts
- Density plots of paralogue similarity distributions

These visualizations help highlight how strict similarity thresholds affect paralogue detection and how ohnologues differ from BioMart's paralogues.

## 💡 Motivation

Understanding differences between duplication annotations (ohnologues vs general paralogues) can help:
- Interpret gene essentiality
- Infer evolutionary history
- Design gene-based screens or knockout studies

## ▶️ How to Run

1. Install required libraries:

```r
install.packages(c("tidyverse", "VennDiagram", "readxl", "writexl", "grid", "UpSetR"))
