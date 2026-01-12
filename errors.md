# Environment Setup and Import Errors

This document lists issues encountered when setting up the environment and testing notebook imports.

## Environment Setup

The environment was created successfully using:
```bash
uv venv dream_proj_env3_uv --python 3.10
uv pip install -r requirements_uv.txt --python dream_proj_env3_uv/bin/python
```

All 218 packages installed successfully with no errors.

---

## Import Testing Results

### Notebook Import Tests

All 9 notebooks had their import cells tested successfully:

| Notebook | Status |
|----------|--------|
| `010324_cptac3_analysis.ipynb` | Success |
| `022924_cross_species.ipynb` | Success |
| `062624_TMS_analysis2.ipynb` | Success |
| `071524_mouse_lifespan.ipynb` | Success |
| `081025_dream_mice_udseq2.ipynb` | Success |
| `111424_human_SComatic_mutations.ipynb` | Success |
| `112024_synapseAD.ipynb` | Success |
| `120925_revisions.ipynb` | Success |
| `122225_revisions_TMS.ipynb` | Success |

### Source Module Import Tests

All source modules import successfully:

| Module | Status |
|--------|--------|
| `read_data` | Success |
| `utils` | Success |
| `expr_dataset` | Success |
| `sc_expr_dataset` | Success |
| `meta_expr_dataset` | Success |
| `single_molecule_seq_plotting` | Success |
| `mutation_dataset` | Success |
| `methyl_dataset` | Success |
| `ad_risk_models` | Success |

---

## Issues Found

### 1. Missing Notebook Referenced in README

**File:** `notebooks/040324_DREAM_activity_validation.ipynb`

**Impact:** The README references this notebook for Figure 1c-d (DREAM activity validation - harmine/INDY treatment and lin-52 mutation), but the file does not exist in the repository.

**Recommendation:** Either add the notebook to the repository or update the README to note it's not included.

---

### 2. Missing Utility Files

**Directory:** `utilities/` exists but is empty (only contains `.gitkeep`)

The following utility files are referenced in source code but are missing:

| File | Referenced By |
|------|--------------|
| `human_mouse_ensembl_genes.txt.gz` | `expr_dataset.py`, `methyl_dataset.py` |
| `human_rat_ensembl_genes.txt` | `expr_dataset.py`, `methyl_dataset.py` |
| `human_worm_WB_ensembl_genes.tsv` | `expr_dataset.py`, `methyl_dataset.py` |
| `gencode.v46.basic.annotation.gtf.gz` | `expr_dataset.py` |
| `gencode.vM34.basic.annotation.gtf.gz` | `sc_expr_dataset.py` |
| `cell_cycle_genesets/cell_cycle_genes_no_dream.txt` | Python scripts |
| `ensembl_to_geneSymbol.tsv` | Notebooks |

**Impact:** These files are needed for gene mapping, annotation, and analysis functions. The modules will import successfully, but functions that use these files will fail at runtime.

**Recommendation:** Add these utility files to the repository or provide instructions for how users can obtain them.

---

### 3. Missing Figure Subdirectories

**Directory:** `figures/` exists but is missing required subdirectories

The following subdirectories are referenced in notebooks but do not exist:

| Directory | Referenced By |
|-----------|--------------|
| `figures/supplementary/` | Multiple notebooks |
| `figures/fig_knockout_mice/` | `081025_dream_mice_udseq2.ipynb` |

**Impact:** Notebooks that attempt to save figures to these directories will fail with "No such file or directory" errors.

**Recommendation:** Create these subdirectories or update notebooks to create them automatically.

---

### 4. Missing Data Directories

The following data directories are referenced in `read_data.py` but do not exist in `data/`:

| Directory | Description |
|-----------|-------------|
| `data/aon_2020` | Missing dataset |
| `data/barger_2008` | Missing dataset |
| `data/cross_species` | Cross-species expression data |
| `data/eisenberg_2016` | Missing dataset |
| `data/hashimoto_2019` | Missing dataset |
| `data/ma_sc_rat_CR` | Rat caloric restriction single-cell data |
| `data/neff_2013` | Missing dataset |
| `data/otero-garcia_2022` | Missing dataset |
| `data/pearson_2008` | Missing dataset |
| `data/zhang_2023` | Missing dataset |
| `data/alexandrov_collab_2025` | UDSeq data (referenced in `081025_dream_mice_udseq2.ipynb`) |

**Impact:** Loading these datasets via `DatasetLoader` will fail. Most notebooks will not be fully runnable without the corresponding data.

**Recommendation:** These data directories need to be added for full reproducibility. Consider:
- Adding the data files to the repository (if size permits)
- Providing download links or instructions
- Using a data sharing platform (e.g., Zenodo, Figshare)

---

## Summary

| Category | Issue Count | Severity |
|----------|-------------|----------|
| Environment setup | 0 | N/A |
| Import errors | 0 | N/A |
| Missing notebooks | 1 | Medium |
| Missing utility files | 7 | High |
| Missing figure directories | 2 | Low |
| Missing data directories | 11 | High |

**Overall:** The environment installs correctly and all imports work. However, notebooks will not be fully runnable due to missing data files, utility files, and figure directories.
