# AD_Oligos: Single-cell Analysis of Oligodendrocytes in Alzheimer's Disease

[![Python](https://img.shields.io/badge/Python-3.12-blue.svg)](https://python.org)
[![Scanpy](https://img.shields.io/badge/Scanpy-1.9.0+-green.svg)](https://scanpy.readthedocs.io/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains comprehensive single-cell RNA sequencing analysis of oligodendrocytes from the SEA-AD (Seattle Alzheimer's Disease) dataset, focusing on two brain regions: the dorsolateral prefrontal cortex (DLPFC) and middle temporal gyrus (MTG) across different Braak stages of Alzheimer's disease progression.

## 🧠 Overview

Alzheimer's disease is characterized by progressive neurodegeneration, and oligodendrocytes play crucial roles in maintaining brain homeostasis through myelination and metabolic support. This project analyzes how oligodendrocyte gene expression patterns change across different stages of Alzheimer's disease pathology (Braak stages) in two key brain regions.

## 📁 Repository Structure

```
AD_Oligos/
├── src/                                    # Analysis notebooks
│   ├── DLPFC_SEA-AD_OL.ipynb              # DLPFC oligodendrocyte analysis
│   ├── MTG_SEA-AD_Oligos.ipynb            # MTG oligodendrocyte analysis
├── output/                                 # Analysis results (generated)
│   ├── DLPFC/                             # DLPFC-specific outputs
│   │   ├── processed_DLPFC_SEA-AD_oligos.h5ad (not pushed to github)
│   │   ├── DLPFC_cluster_braak_composition_*.csv
│   │   ├── DLPFC_deg_results_*.csv       # Single-cell DEG results
│   │   ├── DLPFC_pseudobulk_deg_results_*.csv  # Pseudo-bulk DEG results
│   │   └── DLPFC_pseudobulk_sample_metadata.csv
│   └── MTG/                               # MTG-specific outputs
│       ├── processed_MTG_SEA-AD_oligos.h5ad (not pushed to github)
│       ├── MTG_cluster_braak_composition_*.csv
│       ├── MTG_deg_results_*.csv         # Single-cell DEG results
│       ├── MTG_pseudobulk_deg_results_*.csv    # Pseudo-bulk DEG results
│       └── MTG_pseudobulk_sample_metadata.csv
├── Methods_Oligos_SEA-AD_Oligos.md        # Detailed methods section
├── environment.yml                         # Conda environment specification
├── .env                                    # Environment variables (data paths)
├── LICENSE
└── README.md
```

## 🚀 Quick Start

### Prerequisites

- Python 3.12+
- Conda environment with required packages (see below)
- Access to SEA-AD dataset files

### Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/yourusername/AD_Oligos.git
   cd AD_Oligos
   ```

2. **Set up conda environment:**
   ```bash
   # Create environment from YAML file
   conda env create -f environment.yml
   conda activate snRNAseq
   ```

3. **Configure environment:**
   ```bash
   # Create .env file with data path
   echo "DATA_PATH=/path/to/your/SEA-AD/data" > .env
   ```

4. **Run analysis:**
   ```bash
   # Start Jupyter
   jupyter notebook
   
   # Open and run notebooks in src/ directory
   # - DLPFC_SEA-AD_OL.ipynb for dorsolateral prefrontal cortex analysis
   # - MTG_SEA-AD_Oligos.ipynb for middle temporal gyrus analysis
   ```

## 📊 Analysis Pipeline

### Data Processing
- **Quality Control**: Cell and gene filtering based on expression thresholds
- **Normalization**: Log-transformation and scaling for downstream analysis
- **Filtering**: Focus on Braak stages: Reference, Braak III, and Braak VI

### Clustering & Visualization
- **Dimensionality Reduction**: PCA and UMAP embedding
- **Clustering**: Leiden clustering with resolution 0.05
- **Visualization**: UMAP plots colored by Braak stage

### Differential Expression Analysis

**Single-cell DEG:**
- **Method**: Wilcoxon rank-sum test on log-normalized data
- **Level**: Cell-level analysis capturing heterogeneity
- **Output**: Gene-level statistics with adjusted p-values
- **Two analyses**: Standard grouping (Reference vs Braak III/VI) and pooled control (Control=Ref+B0+BII vs Braak III/VI)

**Pseudo-bulk DEG (DESeq2):**
- **Method**: DESeq2 negative binomial generalized linear model
- **Normalization**: DESeq2 median-of-ratios normalization
- **Dispersion**: Empirical Bayes shrinkage for robust variance estimates
- **Advantage**: Gold-standard RNA-seq analysis treating donors as biological replicates
- **Output**: Shrunken fold changes with FDR-adjusted p-values
- **Two analyses**: Standard grouping and pooled control (Reference+Braak0+BraakII) for increased power

### Multiple Testing Correction Strategy

- **Within-comparison control**: Every differential expression test (single-cell and DESeq2) reports Benjamini–Hochberg FDR-adjusted p-values (`pval_adj` / `padj`). The adjustment is applied independently for each comparison (e.g., Reference vs Braak III) across the genes tested in that contrast.
- **Across-comparison control**: No additional correction is currently applied across the broader set of analyses (single-cell vs pseudo-bulk, standard vs pooled, MTG vs DLPFC). These comparisons are partially overlapping and exploratory, so a single global adjustment would be overly conservative and reduce power substantially.
- **Recommended practice**: Interpret within-comparison FDR < 0.05 as evidence specific to that contrast, and prioritize genes that replicate across multiple analyses. If stricter family-wise control is required, consider applying Holm/BH procedures post hoc to the minimum adjusted p-values per gene or restricting to a predefined analysis hierarchy.

**Comparisons for both methods:**
- Reference vs Braak III
- Braak III vs Braak VI

## 📈 Key Results

The analysis generates comprehensive outputs including:

- **Cluster Composition Analysis**: Statistical comparison of cell populations across Braak stages with chi-square testing
- **Single-cell Differential Gene Expression**: Cell-level identification of genes with altered expression in AD progression
- **Pseudo-bulk Differential Gene Expression**: Donor-aggregated analysis providing increased statistical power
- **Volcano Plots**: Visualization of significant gene expression changes for both single-cell and pseudo-bulk results
- **UMAP Embeddings**: Spatial organization of cells by disease stage
- **Sample Metadata**: Donor-level summaries including cell counts per Braak stage

### Output Files Per Region (DLPFC/MTG):

**Clustering & Composition:**
- `*_cluster_braak_composition_counts.csv`: Cell counts per cluster and Braak stage
- `*_cluster_braak_composition_proportions.csv`: Proportions of Braak stages within each cluster

**Single-cell DEG (Standard Grouping):**
- `*_deg_results_Reference_vs_Braak_III.csv`: Genes differentially expressed between Reference and Braak III
- `*_deg_results_Braak_III_vs_Braak_VI.csv`: Genes differentially expressed between Braak III and Braak VI

**Single-cell DEG (Pooled Control):**
- `*_singlecell_pooled_Control_vs_Braak_III.csv`: Single-cell DEG (Control=Ref+B0+BII vs Braak III)
- `*_singlecell_pooled_Control_vs_Braak_VI.csv`: Single-cell DEG (Control vs Braak VI)
- `*_singlecell_pooled_Braak_III_vs_Braak_VI.csv`: Single-cell DEG (Braak III vs VI)

**Pseudo-bulk DEG (DESeq2 - Standard Grouping):**
- `*_deseq2_standard_Reference_vs_Braak_III.csv`: DESeq2 results (Reference: 3 donors)
- `*_deseq2_standard_Braak_III_vs_Braak_VI.csv`: DESeq2 results (well-powered)

**Pseudo-bulk DEG (DESeq2 - Pooled Control):**
- `*_deseq2_pooled_Control_vs_Braak_III.csv`: DESeq2 results (Control=Ref+Braak0+BraakII: 9 donors)
- `*_deseq2_pooled_Control_vs_Braak_VI.csv`: DESeq2 results (increased power)
- `*_deseq2_pooled_Braak_III_vs_Braak_VI.csv`: DESeq2 results

**Metadata & Summaries:**
- `*_pseudobulk_sample_metadata.csv`: Sample-level metadata with donor IDs and cell counts
- `*_deseq2_summary.csv`: Comparison of significant genes across all analyses

## 🛠️ Computational Environment

- **Platform**: Stanford Sherlock high-performance computing cluster
- **Environment**: Conda (snRNAseq) with Python 3.12
- **Key Packages**: Scanpy, pandas, matplotlib, seaborn, scipy
- **Configuration**: High-resolution plotting (150 DPI display, 300 DPI saved figures)

## 📋 Data Requirements

### Input Files
- `OLIGO_DLPFC_SEA-AD.h5ad`: DLPFC oligodendrocyte data
- `OLIGO_MTG_SEA-AD.h5ad`: MTG oligodendrocyte data

### Data Format
- **Format**: AnnData objects (.h5ad)
- **Metadata**: Braak stage annotations required
- **Genes**: Gene symbols in `feature_name` column preferred

## 📝 Methods

Detailed methods are available in [`Methods_Oligos_SEA-AD_Oligos.md`](Methods_Oligos_SEA-AD_Oligos.md), including:

- Quality control parameters
- Clustering algorithms and parameters
- Single-cell statistical testing approaches
- Pseudo-bulk aggregation methodology
- DESeq2 differential expression analysis (standard and pooled groupings)
- Visualization methods
- Computational specifications

**DESeq2 Implementation**: See [`DESeq2_IMPLEMENTATION_INSTRUCTIONS.md`](DESeq2_IMPLEMENTATION_INSTRUCTIONS.md) for detailed instructions on the pseudo-bulk analysis pipeline.

### Key Methodological Notes:

**Why Both Single-cell and Pseudo-bulk DEG?**

The analysis includes both approaches because they complement each other:

1. **Single-cell DEG** (Wilcoxon rank-sum test):
   - Captures cell-level heterogeneity
   - Tests differences between individual cells
   - Can detect genes with variable expression patterns
   - May have inflated significance due to treating cells as independent
   - Performed with both standard and pooled groupings for comprehensive analysis

2. **Pseudo-bulk DEG** (DESeq2 on donor-aggregated counts):
   - Treats donors as biological replicates (proper experimental unit)
   - Uses negative binomial model appropriate for count data
   - Empirical Bayes shrinkage stabilizes fold change estimates
   - Accounts for mean-variance relationship and dispersion
   - Gold standard for RNA-seq differential expression
   - Two complementary analyses:
     * **Standard**: Preserves biological distinctions (Reference n=3, caution advised)
     * **Pooled**: Combines early-stage groups (Reference+Braak0+BraakII) for better power (n=9)

**Interpretation Guidelines:**
- **Most robust genes**: Significant in single-cell AND DESeq2 (either grouping)
- **Standard vs Pooled DESeq2**: Standard preserves distinctions but has limited power for Reference; Pooled increases power but assumes Reference/Braak0/BraakII are similar (pre-clinical stages)
- Genes significant only in single-cell may reflect cell-state transitions or pseudo-replication
- Genes significant only in pseudo-bulk have consistent donor-level changes
- **Recommended**: Focus on genes significant across multiple methods for validation

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **SEA-AD Consortium** for providing the single-cell RNA sequencing data
- **Stanford Sherlock** for computational resources
- **Scanpy Community** for excellent single-cell analysis tools

## 📧 Contact

For questions or collaborations, please open an issue or contact the maintainers.

## 📚 References

Citation data to follow.
 

---

**Note**: This analysis is part of ongoing research into Alzheimer's disease pathology. Please ensure appropriate data access permissions and ethical approvals before reproducing these analyses.