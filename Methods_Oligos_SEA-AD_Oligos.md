# Methods

## Single-cell RNA Sequencing Analysis of Oligodendrocytes in the Dorsolateral Prefrontal Cortex and Middle Temporal Gyrus

### Data Acquisition and Preprocessing

Oligodendrocyte single-cell RNA sequencing data from two brain regions was obtained from the SEA-AD (Seattle Alzheimer's Disease) dataset: (1) the dorsolateral prefrontal cortex (DLPFC) and (2) the middle temporal gyrus (MTG). Both datasets contained expression profiles for oligodendrocytes across different Braak stages of Alzheimer's disease progression. Data was loaded as AnnData objects and inspected for quality metrics including cell and gene counts, available metadata columns, and data normalization status.

### Quality Control and Data Filtering

Data normalization status was assessed by examining the distribution of expression values and checking for the presence of raw count data. For datasets containing raw counts, standard quality control metrics were calculated including mitochondrial gene expression percentages. Cells with fewer than 200 expressed genes and genes expressed in fewer than 3 cells were filtered out. Cells with greater than 20% mitochondrial gene expression were excluded from downstream analysis.

### Data Normalization and Feature Selection

For raw count data, expression values were normalized to 10,000 transcripts per cell and log-transformed using log1p transformation. Highly variable genes were identified using Scanpy's implementation with parameters: minimum mean expression of 0.0125, maximum mean expression of 3, and minimum dispersion of 0.5. Only highly variable genes were retained for downstream clustering and differential expression analysis.

### Braak Stage Filtering

Cells were filtered to include only those from three specific Braak stages: Reference (no neuropathology), Braak III (moderate neurofibrillary tangle pathology), and Braak VI (severe neurofibrillary tangle pathology). This filtering was performed using exact string matching on the 'Braak stage' metadata column to ensure accurate classification.

### Clustering Analysis

Dimensionality reduction and clustering were performed using standard single-cell analysis pipelines. Expression data was scaled to unit variance with a maximum value of 10. Principal component analysis (PCA) was performed using the ARPACK solver, and the first 40 principal components were used to construct a neighborhood graph with 10 nearest neighbors. Uniform Manifold Approximation and Projection (UMAP) was used for two-dimensional visualization. Leiden clustering was performed with a resolution parameter of 0.05 to identify distinct oligodendrocyte subpopulations.

### Cluster Composition Analysis

The composition of each cluster was analyzed with respect to Braak stage using cross-tabulation analysis. Absolute counts and proportions of cells from each Braak stage within each cluster were calculated and visualized using stacked bar plots. Statistical significance of differences in cluster composition across Braak stages was assessed using chi-square tests of independence.

### Differential Gene Expression Analysis

Differential gene expression analysis was performed between Braak stages using the Wilcoxon rank-sum test implemented in Scanpy, with identical procedures applied to both DLPFC and MTG datasets. To improve interpretability, gene symbols were used as primary identifiers instead of Ensembl IDs. Two pairwise comparisons were performed for each brain region: (1) Reference versus Braak III, and (2) Braak III versus Braak VI. For each comparison, genes with adjusted p-values < 0.05 and absolute log2 fold changes > 0.5 were considered significantly differentially expressed.

### Visualization and Statistical Analysis

Volcano plots were generated to visualize differential gene expression results, with genes colored according to their expression patterns: red points indicating higher expression in the comparison group, blue points indicating higher expression in the reference group. Statistical significance thresholds were indicated by dashed lines. Top differentially expressed genes were annotated on plots for easy identification.

UMAP plots were generated to visualize the spatial distribution of cells colored by Braak stage for both brain regions, providing insight into potential clustering of cells from different disease stages and enabling comparison of cellular organization between DLPFC and MTG.

### Data Processing and Reproducibility

All analyses were performed using Python with the Scanpy package (version 1.9.0 or later) for single-cell analysis, pandas for data manipulation, matplotlib and seaborn for visualization, and scipy for statistical testing. Raw data was preserved throughout the analysis pipeline, with normalized and processed versions stored separately. Results were saved as both processed AnnData objects and individual CSV files for each differential expression comparison, with separate output directories for DLPFC and MTG analyses. All output files were prefixed with their respective brain region identifiers (DLPFC_ or MTG_) to ensure clear organization and reproducibility.

### Computational Environment

Analysis was conducted Stanford's Sherlock cluster in a Jupyter notebook environment on the Stanford Sherlock high-performance computing cluster. The analysis utilized a dedicated conda environment with Python 3.12, configured specifically for single-cell RNA sequencing analysis. Key software packages included Scanpy (version 1.9.0 or later) for single-cell analysis, pandas for data manipulation, matplotlib and seaborn for visualization, and scipy for statistical testing. The conda environment also included dependencies for AnnData handling, UMAP implementation, and statistical analysis libraries.

Statistical analyses were performed using default parameters unless otherwise specified, with multiple testing correction applied using the Benjamini-Hochberg method for differential expression analysis.

The analysis pipeline was designed to be reproducible across different computing environments, with all software versions and dependencies documented in the conda environment specification. Data paths were managed through environment variables to ensure portability between different systems and users.
