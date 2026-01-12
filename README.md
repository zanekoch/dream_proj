# DREAM Complex and Somatic Mutations Research Project

This repository contains the analysis code and documentation for the manuscript **"The DREAM complex links somatic mutation, lifespan, and disease"** by Koch et al.

## Project Overview

This bioinformatics research project investigates the relationship between DREAM (Dp, Rb-like-1, E2f, And MuvB) complex activity and somatic mutation burden across multiple species and datasets. The study demonstrates that DREAM complex activity significantly impacts lifetime somatic mutation burden, with effects linked to altered lifespan and age-related disease pathology.

## Figure Generation Guide

This section maps each figure in the manuscript to the specific notebook(s) where the analyses and visualizations were generated.

### Main Figures

#### Figure 1: DREAM Activity Measurement and Validation
- **Figure 1a-b (DREAM activity concept/heatmap)**: `notebooks/062624_TMS_analysis2.ipynb`
- **Figure 1c-d (DREAM activity validation - harmine/INDY treatment and lin-52 mutation)**: `notebooks/040324_DREAM_activity_validation.ipynb` *(notebook not included in repository)*

#### Figure 2: DREAM Activity and Somatic Mutations in Single Cells (Tabula Muris Senis)
- **All panels (2a-h)**: `notebooks/062624_TMS_analysis2.ipynb`
- Mutation burden by age across tissues (2a)
- DREAM activity vs mutation burden across tissues (2b-e)
- Liver hepatocytes DREAM activity vs mutations by age (2f)
- Cell type and tissue analysis (2g-h)

#### Figure 3: Cross-Species Analysis
- **All panels (3a-f)**: `notebooks/022924_cross_species.ipynb`
- Species phylogeny and DREAM activity (3a-b)
- DREAM activity vs maximum lifespan by tissue (3c-d)
- Cross-tissue DREAM coordination (3e)
- DREAM activity vs species mutation rates (3f)

#### Figure 4: Mouse Lifespan and High-Fat Diet Study
- **All panels (4a-c)**: `notebooks/071524_mouse_lifespan.ipynb`
- Kaplan-Meier survival curves by DREAM activity (4a)
- Cox proportional hazards analysis (4b)
- DREAM activity vs lifespan difference (4c)

#### Figure 5: Human Somatic Mutations and Alzheimer's Disease
- **Panels (5a-e)**: `notebooks/111424_human_SComatic_mutations.ipynb`; **Panels (5f-g)** `notebooks/112024_synapseAD.ipynb`
- Human tissue mutation rates (5b-c)
- DREAM activity vs mutation burden (5d)
- DREAM activity vs neuropathology (5e)
- Alzheimer's disease risk analysis (5f-g)

#### Figure 6: DREAM Knockout Mouse Study
- **All panels (6a-h)**: `notebooks/081025_dream_mice_udseq2.ipynb`
- DREAM knockout experimental design (6a)
- Single-base substitution analysis (6b-d)
- Insertion/deletion mutation analysis (6e-g)
- Summary model (6h)

### Supplementary Figures

#### Supplementary Figure 1: DREAM Activity Validation
- **Panels 1a-b**: - `notebooks/010324_cptac3_analysis.ipynb`
- DREAM activity vs protein abundance and phosphorylation status 
- **Panels 1c-d**: `notebooks/062624_TMS_analysis2.ipynb`
- DREAM complex protein correlations
- DREAM activity vs proliferation markers
- **Panel 1e-j**: - `notebooks/122225_revisions_TMS.ipynb`
- DREAM activity vs proliferation markers  contd.
- DNA repair gene expression in high vs. low DREAM

#### Supplementary Figure 2: Additional TMS Analysis
- **All panels**: `notebooks/062624_TMS_analysis2.ipynb`
- DREAM activity distributions by organ and age
- Additional mutation burden correlations
- Proliferation rate corrections

#### Supplementary Figure 3: Cross-Species Analysis Details
- **All panels**: `notebooks/022924_cross_species.ipynb`
- Body weight corrections
- Additional species analysis

#### Supplementary Figure 4: Species Stress Resistance
- **All panels**: `notebooks/022924_cross_species.ipynb`
- Stress resistance vs DREAM activity across species

#### Supplementary Figure 5: Mouse Lifespan Study Details
- **All panels**: `notebooks/071524_mouse_lifespan.ipynb`
- Strain characteristics and survival curves
- Physiological covariates

#### Supplementary Figure 6: Human Mutation Analysis Details
- **All panels**: `notebooks/111424_human_SComatic_mutations.ipynb` and `notebooks/112024_synapseAD.ipynb`
- Single cell type mutation rates
- Literature comparisons
- Additional AD analysis

#### Supplementary Figure 7: Duplex sequencing 
- **All panels**: `notebooks/081025_dream_mice_udseq2.ipynb`
- Duplex sequencing modelling and outlier detection

#### Supplementary Figure 8: The association of DREAM-related pathways with lifespan and somatic mutation burden
- **All panels**: `notebooks/120925_revisions.ipynb`
- Controlling for parallel transcriptional programs does not abrogrogate DREAM's assoc. with lifespan or mutation rate.


## Source Code Documentation

### Core Data Loading Architecture (`source/`)

#### `read_data.py`
Main entry point for data loading. Contains the `DatasetLoader` factory class that handles loading different types of datasets (expression, mutation, methylation) across multiple data sources.

**Key Features:**
- Unified interface for loading diverse datasets
- Supports Tabula Muris Senis, cross-species, human datasets
- Handles both single-cell and bulk data

#### Dataset Classes
- **`expr_dataset.py`**: Expression data handling for bulk RNA-seq data
- **`sc_expr_dataset.py`**: Single-cell expression data processing (e.g., Tabula Muris Senis, Tabula Sapiens)
- **`mutation_dataset.py`**: Somatic mutation data processing and analysis
- **`methyl_dataset.py`**: DNA methylation data handling
- **`meta_expr_dataset.py`**: Meta-analysis utilities for combining multiple expression datasets

#### Analysis Utilities
- **`utils.py`**: Core utility functions including:
  - `read_dream_files()`: Load DREAM gene sets
  - Statistical plotting functions (exponential fits, major axis regression)
  - Cross-species gene mapping utilities
- **`single_molecule_seq_plotting.py`**: Specialized functions for UDSeq (single-molecule sequencing) data analysis and visualization
- **`ad_risk_models.py`**: Alzheimer's disease risk modeling and regression analysis

### Analysis Scripts

#### Python Scripts (`python_scripts/`)
- **`calc_ssgsea_dream_activity.py`**: Calculate DREAM activity scores using single-sample GSEA
- **`run_sigprofiler_extractor.py`**: Extract mutational signatures using SigProfiler
- **`run_cosmic_fit_with_ffpe.py`**: Fit COSMIC mutational signatures with FFPE correction
- **`random_background_*_corr_*.py`**: Random background analysis for various datasets

#### R Scripts (`r_scripts/`)
- **`SCNormalization.R`**: Single-cell RNA-seq normalization procedures

#### Bash Scripts (`bash_scripts/`)
- **`align_jung_rnaseq.sh`**: RNA-seq alignment pipeline
- **`call_variants_jung_rnaseq_oneChr.sh`**: Variant calling from RNA-seq data
- **`create_mutation_matrix.sh`**: Generate mutation count matrices
- **Various processing scripts**: BAM file processing, VCF handling, genome indexing

#### SLURM Scripts (`slurm_scripts/`)
- **`submit_py_jobs.job`**: Submit Python analysis jobs to HPC cluster
- **`submit_r_jobs.job`**: Submit R analysis jobs to HPC cluster
- **`align_jung_rnaseq.sh`**: HPC-optimized RNA-seq alignment

## Environment Setup

### Using uv (Recommended)

[uv](https://github.com/astral-sh/uv) is a fast Python package manager. The project environment is located at `dream_proj_env3_uv/`.

```bash
# Activate the environment
source dream_proj_env3_uv/bin/activate

# To use in Jupyter notebooks, register the kernel (one-time setup):
python -m ipykernel install --user --name dream_proj_env3_uv --display-name "Python (dream_proj_env3_uv)"

# Then select "Python (dream_proj_env3_uv)" as the kernel in your notebook
```

#### Installing uv (if needed)
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

#### Recreating the environment
```bash
uv venv dream_proj_env3_uv --python 3.10
uv pip install -r requirements_uv.txt --python dream_proj_env3_uv/bin/python
```

### Key Dependencies
- **Python**: pandas, numpy, scipy, matplotlib, seaborn, scanpy, anndata
- **Bioinformatics**: pybiomart, gseapy, sigprofiler, cosmic-fit
- **Statistics**: statsmodels, lifelines, scikit-learn
- **Visualization**: colorcet, upsetplot

## Data Architecture

### Key Datasets
- **Tabula Muris Senis**: Single-cell atlas of aging in 21 mouse tissues
- **Cross-Species**: Expression data from 92 mammalian species
- **Human Data**: Tabula Sapiens and SEA-AD (Seattle Alzheimer's Disease Brain Cell Atlas)
- **Mouse Lifespan**: 50 inbred strains with survival data
- **DREAM Knockout**: UDSeq mutation data from DREAM loss-of-function mice

### Gene Sets
- **DREAM Complex Targets**: 328 genes with validated DREAM binding sites (from Bujarrabal-Dueso et al. 2024, `./data/bujarrabal_dueso/tableS12_dream_promoter_binding.csv`)

## Key Findings Summary

1. **Cellular Level**: DREAM activity associates with increased somatic mutation rates in single cells
2. **Species Level**: Lower DREAM activity predicts longer maximum lifespan across 92 mammalian species
3. **Disease Level**: Reduced DREAM activity in humans predicts later Alzheimer's disease onset
4. **Experimental Validation**: DREAM knockout mice show 4.2% reduction in single-base substitutions and 19.6% reduction in indels

## Citation

If you use this code or data, please cite:

> Koch, Z. et al. (2025). The DREAM complex links somatic mutation, lifespan, and disease. https://www.biorxiv.org/content/10.1101/2025.09.15.676396v1

## Contact

- **Corresponding Author**: Trey Ideker (tideker@health.ucsd.edu)
- **First Author**: Zane Koch (zkoch@ucsd.edu)
