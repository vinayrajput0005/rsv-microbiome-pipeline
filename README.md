# Comprehensive Microbial Genomics Analysis Pipeline (RSV Cohort)

[![R Version](https://img.shields.io/badge/R-%E2%89%A4%204.0-blue.svg)](https://www.r-project.org/)
[![Bioinformatics](https://img.shields.io/badge/Domain-Bioinformatics%20%2F%20Genomics-green.svg)]()
[![Microbiome](https://img.shields.io/badge/Library-Phyloseq%20%2F%20DESeq2-orange.svg)]()

An R workflow tailored for deep molecular ecological profiling of Next-Generation Sequencing (NGS) microbiome data. This pipeline was developed to characterize and deconstruct upper respiratory or gut microbiota composition across different clinical presentations of **Respiratory Syncytial Virus (RSV A, RSV B)** compared against a **Recovered** cohort.

The architecture integrates alpha/beta diversity modeling, non-parametric multi-factor PERMANOVA, robust negative binomial dispersion models (DESeq2), and non-metric multidimensional scaling (NMDS) combined with continuous clinical covariate vector fitting (`envfit`) to map viral load dynamics (qPCR Ct values).

---

## 📑 Table of Contents
1. Workflow Architecture
2. Prerequisites & Dependencies
3. Input Data Specification
4. Pipeline Components
5. Key Features & Statistical Rigor
6. Expected Output Deliverables
7. Execution Guide
8. Developer & Citation Information

---

## 🧬 Workflow Architecture

The script processes data sequentially through a modular framework to guarantee reproducibility, prevent factor-level dropouts, and maximize biological insight:

[Raw BIOM + Tax Tables] ──> [Rigorous Metadata Cleaning] ──> [Dual-Barrier Abundance Filtering]
                                                                        │
       ┌─────────────────────────┬──────────────────────────────────────┼────────────────────────┐
       ▼                         ▼                                      ▼                        ▼
[Alpha Diversity]        [Beta Diversity]                       [DESeq2 Profiling]       [Core Microbiome Landscape]
 - Shannon / Simpson      - Compositional Bray-Curtis            - Negative Binomial      - Prevalence (>=70%)
 - Rarefaction Modeling   - Multi-factor PERMANOVA                 Dispersion Models      - Median Abundance (>=0.1%)
 - Kruskal-Wallis Test    - PCoA Space Ordination                - Dynamic Contrasts      - 3-Way Venn Intersections
                                                                        │
                                                                        ▼
                                                             [Clinical Burden Mapping]
                                                              - RSV A / RSV B Subsets
                                                              - NMDS Ordination
                                                              - envfit Vector Mapping (qPCR Ct)

---

## 📦 Prerequisites & Dependencies

The analytics engine relies heavily on the Bioconductor and CRAN ecosystems. Ensure you have R (v4.0 or greater) installed along with the following libraries:

### Bioconductor Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("phyloseq", "DESeq2", "microbiome"))

### CRAN Packages
install.packages(c("ggplot2", "vegan", "ggpubr", "dplyr", "tidyr", 
                   "viridis", "rstatix", "RColorBrewer", "pheatmap", 
                   "ggVennDiagram", "ggrepel", "scales", "reshape2"))

---

## 💾 Input Data Specification

The script looks directly for a formatted BIOM file in your working directory:
* Filename: otu_table_merged.biom
* Taxonomy Mapping: Must include standard tier labeling containing a field mapped to "Genus".
* Sample Metadata Fields Required:
    - Virus: Categorical factor containing "RSV A", "RSV B", and "Recovered".
    - Symptoms_Severity: Categorical factor containing "Mild", "Moderate", "Severe", "Recovered".
    - Status: Categorical factor containing "Positive", "Recovered".
    - Age_group: Categorical covariate used as a confounding control in PERMANOVA blocks.
    - Ct_value: Numeric continuous variable representing diagnostic qPCR cycle threshold values (viral load indicator).

---

## ⚙️ Pipeline Components

### 1. Data Sanitization & Prevalence Filtering
* Enforces explicit experimental factor level orders to fix downslope comparison anchors.
* Executes dual-barrier filtration: purges operational taxonomic units (OTUs) with total absolute counts <= 10, and removes low-prevalence taxa missing from >90% of global sample libraries.

### 2. Alpha and Beta Ecology Dynamics
* Rarefies depth stochastically to the lowest sample sequencing boundary to guarantee even-depth alpha testing.
* Deconstructs spatial ecology using a Compositional Bray-Curtis distance matrix processed through PERMANOVA (adonis2) controlling for age brackets.

### 3. Differential Abundance Modeling (DESeq2)
* Translates taxonomic tables seamlessly into negative binomial distribution objects.
* Deploys custom data-slicing functions to automatically run multiple contrast groups (RSV A vs RSV B, RSV A vs Recovered, etc.) while protecting against low sample size exceptions.

### 4. Continuous Vector Fitting (envfit)
* Subsets diagnostic-positive samples to extract core physiological insights.
* Projects continuous variable vectors (Ct_value) onto a Non-Metric Multidimensional Scaling (NMDS) manifold to visualize exactly how viral load curves shift community geography.

### 5. Multi-Threshold Core Extraction
* Defines the core homeostatic microbiome using strict dual boundaries: a minimum 70% population prevalence inside a cohort group AND a minimum 0.1% median relative abundance. 
* Generates comparative group exports alongside standard 3-way Venn intersection matrices.

---

## 📊 Key Features & Statistical Rigor

* Reproducibility Anchor: set.seed(123) is hardcoded prior to any stochastic operation (Rarefaction, PERMANOVA permutations, NMDS stress optimization loops) to enforce run-to-run consistency.
* Outlier Resiliency: Uses median relative abundance instead of mean relative abundance for core microbiome selection, protecting the target profiles from being biased by high-abundance anomalies in singular samples.
* Automated Factor Dropping: Uses strict droplevels() calls inside sub-population analyses to eliminate dangling factors that cause downstream contrast model failures.

---

## 📁 Expected Output Deliverables

Upon successful run completion, the pipeline populates your workspace with high-resolution graphics (600 DPI publication standard) and comma-separated numeric tables:

### 📈 Generated Visualizations
* Ecology_Alpha_Diversity_Shannon.png: Boxplots enriched with Kruskal-Wallis omnibus metrics and jitter configurations.
* Ecology_Beta_PCoA_BrayCurtis.png: PCoA coordinate space plotting accompanied by R2 and p-value PERMANOVA subtitle keys.
* Overlap_Distribution_Barplot.png & Overlap_Venn_Intersection.png: Distribution visualizations isolating core vs unique shared candidates.
* Clinical_NMDS_Infectious_Burden_Map.png: High-fidelity NMDS layout carrying the mapped envfit continuous vector gradient for viral loads.
* Core_Microbiome_Publication_Quality.pdf / .png: Highly detailed scatter manifold highlighting core taxonomic guilds faceted by virus status.

### 📋 Generated Data Tables
* DESeq2_All_Features_*.csv & DESeq2_Signif_Features_*.csv: Full and statistically significant (padj < 0.05, |log2FC| > 1) differentially abundant genera.
* Overlap_Shared_Taxa_Intersection.csv / Overlap_Unique_Taxa_*.csv: CSV lists detailing structural intersections.
* Core_Microbiome_Comprehensive_Matrix.csv: Integrated file mapping all parsed background data.
* Core_Microbiome_Restricted_Only.csv & Core_Microbiome_Export_*.csv: Explicitly isolated core ecosystem lists split globally and by group.

---

## 🚀 Execution Guide

Make the script executable and run it natively from your terminal interface:
$ chmod +x microbiome_pipeline.R
$ ./microbiome_pipeline.R

Alternatively, invoke it directly within your local R session or RStudio environment:
> source("microbiome_pipeline.R")

---

## 🧑‍💻 Developer & Citation Information

* Lead Architecture Developer: Vinay Rajput, PhD
* Specialization: Next-Generation Sequencing (NGS) Analysis, Bioinformatics, & Microbial Genomics
* Domain Focus: Environmental Surveillance, Public Health Intelligence, One Health Frameworks

If utilizing this codebase or analytics framework in secondary research pipelines or clinical surveillance setups, please ensure correct attribution to the repository structure. For collaborative integrations or structural custom adjustments, please open an issue tracking card in the project repository dashboard.
