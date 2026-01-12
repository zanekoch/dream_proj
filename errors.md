# Environment Setup and Import Errors

This document lists issues encountered when setting up the environment and testing notebook imports.

## Environment Setup

The environment was created successfully using:
```bash
uv venv dream_proj_env3_uv --python 3.10
uv pip install -r requirements_uv.txt --python dream_proj_env3_uv/bin/python
```

All 218 packages installed successfully.

---

## Import Testing Results

### Successful Notebooks

The following notebooks had their imports tested successfully:
- `010324_cptac3_analysis.ipynb` - All imports work
- `022924_cross_species.ipynb` - All imports work (after fix)
- `062624_TMS_analysis2.ipynb` - All imports work
- `071524_mouse_lifespan.ipynb` - All imports work
- `081025_dream_mice_udseq2.ipynb` - All imports work
- `111424_human_SComatic_mutations.ipynb` - All imports work
- `112024_synapseAD.ipynb` - All imports work
- `120925_revisions.ipynb` - All imports work
- `122225_revisions_TMS.ipynb` - All imports work

### Source Module Import Results

| Module | Status |
|--------|--------|
| `read_data` | Success |
| `utils` | Success |
| `expr_dataset.ExpressionDataset` | Success |
| `sc_expr_dataset.ScExpressionDataset` | Success |
| `meta_expr_dataset.MetaExpressionDataset` | Success |
| `single_molecule_seq_plotting` | Success |
| `mutation_dataset.MutationDataset` | Success |
| `methyl_dataset.MethylationDataset` | Success |
| `ad_risk_models` | Success |

---

## Issues Found and Fixed

### 1. Missing Notebook Referenced in README - FIXED

**Status:** Fixed (README updated to note notebook not included)

The README references `notebooks/040324_DREAM_activity_validation.ipynb` for Figure 1c-d, but this notebook does not exist in the repository. The README has been updated to indicate this notebook is not included.

---

### 2. Missing Import Cell in `022924_cross_species.ipynb` - FIXED

**Status:** Fixed

Added import cell at the beginning of the notebook with all necessary imports.

---

### 3. Missing `utilities/` Directory - FIXED (empty)

**Status:** Fixed (created empty directory)

The `utilities/` directory has been created but is empty. The following files need to be added for full functionality:

| File | Referenced By |
|------|--------------|
| `human_mouse_ensembl_genes.txt.gz` | `expr_dataset.py`, `methyl_dataset.py` |
| `human_rat_ensembl_genes.txt` | `expr_dataset.py`, `methyl_dataset.py` |
| `human_worm_WB_ensembl_genes.tsv` | `expr_dataset.py`, `methyl_dataset.py` |
| `gencode.v46.basic.annotation.gtf.gz` | `expr_dataset.py` |
| `gencode.vM34.basic.annotation.gtf.gz` | `sc_expr_dataset.py` |
| `cell_cycle_genesets/cell_cycle_genes_no_dream.txt` | `random_background_synapseHarmonization.py` |
| `ensembl_to_geneSymbol.tsv` | Referenced in notebooks |

---

### 4. Missing `figures/` Directory - FIXED (empty)

**Status:** Fixed (created empty directory with subdirectories)

Created:
- `figures/`
- `figures/supplementary/`
- `figures/fig_knockout_mice/`

---

## Remaining Issues (Not Fixed)

### Missing Data Directories

The following data directories referenced in `read_data.py` are missing. These were not fixed as requested:

| Directory | Status |
|-----------|--------|
| `data/aon_2020` | Missing |
| `data/barger_2008` | Missing |
| `data/cross_species` | Missing |
| `data/eisenberg_2016` | Missing |
| `data/hashimoto_2019` | Missing |
| `data/ma_sc_rat_CR` | Missing |
| `data/neff_2013` | Missing |
| `data/otero-garcia_2022` | Missing |
| `data/pearson_2008` | Missing |
| `data/zhang_2023` | Missing |
| `data/alexandrov_collab_2025` | Missing (referenced in `081025_dream_mice_udseq2.ipynb`) |

---

## Summary

| Issue | Status |
|-------|--------|
| Missing notebook in README | Fixed (noted as not included) |
| Missing import cell in cross_species notebook | Fixed |
| Missing utilities/ directory | Fixed (empty) |
| Missing figures/ directory | Fixed (empty with subdirs) |
| Missing data directories | Not fixed (as requested) |
| Missing utility files | Not fixed (directory created empty) |
