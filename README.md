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

**Pseudo-bulk DEG:**
- **Method**: Aggregation by donor followed by Mann-Whitney U test
- **Normalization**: Counts per million (CPM)
- **Advantage**: Increased statistical power by treating donors as biological replicates
- **Output**: Gene-level statistics with robust fold change estimates

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

**Single-cell DEG:**
- `*_deg_results_Reference_vs_Braak_III.csv`: Genes differentially expressed between Reference and Braak III
- `*_deg_results_Braak_III_vs_Braak_VI.csv`: Genes differentially expressed between Braak III and Braak VI

**Pseudo-bulk DEG:**
- `*_pseudobulk_deg_results_Reference_vs_Braak_III.csv`: Pseudo-bulk DEG results with log2FC and adjusted p-values
- `*_pseudobulk_deg_results_Braak_III_vs_Braak_VI.csv`: Pseudo-bulk DEG results with log2FC and adjusted p-values
- `*_pseudobulk_sample_metadata.csv`: Sample-level metadata including donor IDs and cell counts

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
- Pseudo-bulk aggregation and differential expression methodology
- Visualization methods
- Computational specifications

### Key Methodological Notes:

**Why Both Single-cell and Pseudo-bulk DEG?**

The analysis includes both approaches because they complement each other:

1. **Single-cell DEG** (Wilcoxon rank-sum test):
   - Captures cell-level heterogeneity
   - Tests differences between individual cells
   - Can detect genes with variable expression patterns
   - May have inflated significance due to treating cells as independent

2. **Pseudo-bulk DEG** (Mann-Whitney U test on aggregated counts):
   - Treats donors as biological replicates (proper experimental unit)
   - Reduces technical noise through aggregation
   - More conservative and statistically robust
   - Better accounts for donor-to-donor variability
   - Recommended for identifying reliable disease-associated genes

**Interpretation Guidelines:**
- Genes significant in **both** analyses are most robust and reliable
- Genes significant only in single-cell may reflect cell-state transitions
- Genes significant only in pseudo-bulk may have consistent directional changes masked by cellular heterogeneity in single-cell analysis

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