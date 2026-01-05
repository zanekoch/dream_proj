# Plan: Fix Hard-Coded Paths

This document outlines the plan to convert all hard-coded absolute paths to relative paths based on the repository root.

## Strategy

All paths will be relative to the repo root (`dream_proj/`). Each file will determine the repo root dynamically:

**For Python source files:**
```python
import os
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
```

**For notebooks (in `notebooks/` directory):**
```python
import os
REPO_ROOT = os.path.dirname(os.getcwd()) if os.path.basename(os.getcwd()) == 'notebooks' else os.getcwd()
# Or more reliably:
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) if '__file__' in dir() else os.path.abspath('..')
```

---

## Files to Edit

### 1. `source/read_data.py` (~45 paths)

| Line | Current Path | New Path |
|------|-------------|----------|
| 177 | `/cellar/users/zkoch/dream/data/petljak_2019` | `os.path.join(REPO_ROOT, 'data/petljak_2019')` |
| 228 | `/cellar/users/zkoch/dream/data/crosby_2022/GSE215974_Processed_data_FE1_hydrazine.xlsx` | `os.path.join(REPO_ROOT, 'data/crosby_2022/GSE215974_Processed_data_FE1_hydrazine.xlsx')` |
| 264 | `/cellar/users/zkoch/dream/data/liu_2021/GSE161789_raw_counts_and_FPKM_for_RNA-seq.xlsx` | `os.path.join(REPO_ROOT, 'data/liu_2021/GSE161789_raw_counts_and_FPKM_for_RNA-seq.xlsx')` |
| 299 | `/cellar/users/zkoch/dream/data/tabula_sapiens` | `os.path.join(REPO_ROOT, 'data/tabula_sapiens')` |
| 524 | `/cellar/users/zkoch/dream/data/SEA-AD/aws_snRNA_seq` | `os.path.join(REPO_ROOT, 'data/SEA-AD/aws_snRNA_seq')` |
| 736-737 | `/cellar/users/zkoch/dream/data/synapse_MIT_ROSMAP_Multiomics/dileep_2023/...` | `os.path.join(REPO_ROOT, 'data/synapse_MIT_ROSMAP_Multiomics/dileep_2023/...')` |
| 866 | `/cellar/users/zkoch/dream/data/synapse_rna_seq_harmonization` | `os.path.join(REPO_ROOT, 'data/synapse_rna_seq_harmonization')` |
| 893 | `/cellar/users/zkoch/dream/data/paine_2024` | `os.path.join(REPO_ROOT, 'data/paine_2024')` |
| 921 | `/cellar/users/zkoch/dream/data/lu_2022` | `os.path.join(REPO_ROOT, 'data/lu_2022')` |
| 1054 | `/cellar/users/zkoch/dream/data/lu_2014` | `os.path.join(REPO_ROOT, 'data/lu_2014')` |
| 1080 | `/cellar/users/zkoch/dream/data/hashimoto_2019/hashimoto_2019.h5ad` | `os.path.join(REPO_ROOT, 'data/hashimoto_2019/hashimoto_2019.h5ad')` |
| 1095 | `/cellar/users/zkoch/dream/data/otero-garcia_2022/...h5ad` | `os.path.join(REPO_ROOT, 'data/otero-garcia_2022/...h5ad')` |
| 1113 | `/cellar/users/zkoch/dream/data/mammalian_methylation_consort/` | `os.path.join(REPO_ROOT, 'data/mammalian_methylation_consort/')` |
| 1162 | `/cellar/users/zkoch/dream/data/gyenis_2023` | `os.path.join(REPO_ROOT, 'data/gyenis_2023')` |
| 1179 | `/cellar/users/zkoch/dream/data/motrpac_2024/transcriptomics/results` | `os.path.join(REPO_ROOT, 'data/motrpac_2024/transcriptomics/results')` |
| 1207 | `/cellar/users/zkoch/dream/data/uxa_2019` | `os.path.join(REPO_ROOT, 'data/uxa_2019')` |
| 1233 | `/cellar/users/zkoch/dream/data/bujarrabal_dueso` | `os.path.join(REPO_ROOT, 'data/bujarrabal_dueso')` |
| 1259 | `/cellar/users/zkoch/dream/data/bujarrabal_dueso` | `os.path.join(REPO_ROOT, 'data/bujarrabal_dueso')` |
| 1285 | `/cellar/users/zkoch/dream/data/cao_et_al_2024` | `os.path.join(REPO_ROOT, 'data/cao_et_al_2024')` |
| 1336 | `/cellar/users/zkoch/dream/data/liu_2023` | `os.path.join(REPO_ROOT, 'data/liu_2023')` |
| 1361-1366 | `/cellar/users/zkoch/dream/data/cross_species/nebulas_sc/...` | `os.path.join(REPO_ROOT, 'data/cross_species/nebulas_sc/...')` |
| 1496 | `/cellar/users/zkoch/dream/data/cross_species/nebulas/...` | `os.path.join(REPO_ROOT, 'data/cross_species/nebulas/...')` |
| 1548 | `/cellar/users/zkoch/dream/data/pearson_2008` | `os.path.join(REPO_ROOT, 'data/pearson_2008')` |
| 1576 | `/cellar/users/zkoch/dream/data/barger_2008` | `os.path.join(REPO_ROOT, 'data/barger_2008')` |
| 1605 | `/cellar/users/zkoch/dream/data/aon_2020` | `os.path.join(REPO_ROOT, 'data/aon_2020')` |
| 1640 | `/cellar/users/zkoch/dream/data/gtex` | `os.path.join(REPO_ROOT, 'data/gtex')` |
| 1683 | `/cellar/users/zkoch/dream/data/eisenberg_2016` | `os.path.join(REPO_ROOT, 'data/eisenberg_2016')` |
| 1711 | `/cellar/users/zkoch/dream/data/neff_2013` | `os.path.join(REPO_ROOT, 'data/neff_2013')` |
| 1739 | `/cellar/users/zkoch/dream/data/zhang_2023` | `os.path.join(REPO_ROOT, 'data/zhang_2023')` |
| 1767 | `/cellar/users/zkoch/dream/data/mercken_2014` | `os.path.join(REPO_ROOT, 'data/mercken_2014')` |
| 1795 | `/cellar/users/zkoch/dream/data/yu_2012` | `os.path.join(REPO_ROOT, 'data/yu_2012')` |
| 1824 | `/cellar/users/zkoch/dream/data/fok_cr_2014` | `os.path.join(REPO_ROOT, 'data/fok_cr_2014')` |
| 1854 | `/cellar/users/zkoch/dream/data/fok_2014` | `os.path.join(REPO_ROOT, 'data/fok_2014')` |
| 1884 | `/cellar/users/zkoch/dream/data/fok_2014` | `os.path.join(REPO_ROOT, 'data/fok_2014')` |
| 1913 | `/cellar/users/zkoch/dream/data/zhou_2012` | `os.path.join(REPO_ROOT, 'data/zhou_2012')` |
| 1963 | `/cellar/users/zkoch/dream/data/ma_sc_rat_CR` | `os.path.join(REPO_ROOT, 'data/ma_sc_rat_CR')` |
| 1989 | `/cellar/users/zkoch/dream/data/palovics_parabiosis/all_tissues.h5ad` | `os.path.join(REPO_ROOT, 'data/palovics_parabiosis/all_tissues.h5ad')` |
| 2054 | `/cellar/users/zkoch/dream/data/martin_montalvo` | `os.path.join(REPO_ROOT, 'data/martin_montalvo')` |
| 2083 | `/cellar/users/zkoch/dream/data/boutant_nestle` | `os.path.join(REPO_ROOT, 'data/boutant_nestle')` |
| 2112 | `/cellar/users/zkoch/dream/data/tyshkovskiy` | `os.path.join(REPO_ROOT, 'data/tyshkovskiy')` |
| 2144 | `/cellar/users/zkoch/dream/data/msalt` | `os.path.join(REPO_ROOT, 'data/msalt')` |
| 2188 | `/cellar/users/zkoch/dream/data/tcga/processed_data` | `os.path.join(REPO_ROOT, 'data/tcga/processed_data')` |

**Action:** Add `REPO_ROOT` definition at top of file, then find/replace all paths.

---

### 2. `source/expr_dataset.py` (~15 paths)

| Line | Current Path | New Path |
|------|-------------|----------|
| 438 | `/cellar/users/zkoch/dream/data/liu_2023/life_history_traits.xlsx` | `os.path.join(REPO_ROOT, 'data/liu_2023/life_history_traits.xlsx')` |
| 465 | `/cellar/users/zkoch/dream/data/cao_et_al_2024/lesion_count_greater_than_one_support.txt` | `os.path.join(REPO_ROOT, 'data/cao_et_al_2024/lesion_count_greater_than_one_support.txt')` |
| 469 | `/cellar/users/zkoch/dream/data/cao_et_al_2024/lesion_count_greater_than_5_support.txt` | `os.path.join(REPO_ROOT, 'data/cao_et_al_2024/lesion_count_greater_than_5_support.txt')` |
| 608 | `/cellar/users/zkoch/dream/data/uxa_2019/supp_table_s2.xlsx` | `os.path.join(REPO_ROOT, 'data/uxa_2019/supp_table_s2.xlsx')` |
| 629 | `/cellar/users/zkoch/dream/data/motrpac_2024/proteomics_untargeted/...` | `os.path.join(REPO_ROOT, 'data/motrpac_2024/proteomics_untargeted/...')` |
| 663 | `/cellar/users/zkoch/dream/data/gyenis_2023` | `os.path.join(REPO_ROOT, 'data/gyenis_2023')` |
| 891 | `/cellar/users/zkoch/dream/utilities/human_mouse_ensembl_genes.txt.gz` | `os.path.join(REPO_ROOT, 'utilities/human_mouse_ensembl_genes.txt.gz')` |
| 903 | `/cellar/users/zkoch/dream/utilities/human_mouse_ensembl_genes.txt.gz` | `os.path.join(REPO_ROOT, 'utilities/human_mouse_ensembl_genes.txt.gz')` |
| 919 | `/cellar/users/zkoch/dream/utilities/human_rat_ensembl_genes.txt` | `os.path.join(REPO_ROOT, 'utilities/human_rat_ensembl_genes.txt')` |
| 939 | `/cellar/users/zkoch/dream/utilities/human_worm_WB_ensembl_genes.tsv` | `os.path.join(REPO_ROOT, 'utilities/human_worm_WB_ensembl_genes.tsv')` |
| 964 | `/cellar/users/zkoch/dream/utilities/human_mouse_ensembl_genes.txt.gz` | `os.path.join(REPO_ROOT, 'utilities/human_mouse_ensembl_genes.txt.gz')` |
| 981 | `/cellar/users/zkoch/dream/utilities/human_mouse_ensembl_genes.txt.gz` | `os.path.join(REPO_ROOT, 'utilities/human_mouse_ensembl_genes.txt.gz')` |
| 1004 | `/cellar/users/zkoch/dream/utilities/human_rat_ensembl_genes.txt` | `os.path.join(REPO_ROOT, 'utilities/human_rat_ensembl_genes.txt')` |
| 1031 | `/cellar/users/zkoch/dream/utilities/human_worm_WB_ensembl_genes.tsv` | `os.path.join(REPO_ROOT, 'utilities/human_worm_WB_ensembl_genes.tsv')` |
| 1607 | `/cellar/users/zkoch/dream/utilities/gencode.v46.basic.annotation.gtf.gz` | `os.path.join(REPO_ROOT, 'utilities/gencode.v46.basic.annotation.gtf.gz')` |

**Action:** Add `REPO_ROOT` definition at top of file, then find/replace all paths.

---

### 3. `source/methyl_dataset.py` (~5 paths)

| Line | Current Path | New Path |
|------|-------------|----------|
| 50 | `/cellar/users/zkoch/dream/data/mammalian_methylation_consort/li_2024/...xlsx` | `os.path.join(REPO_ROOT, 'data/mammalian_methylation_consort/li_2024/...xlsx')` |
| 78 | `/cellar/users/zkoch/dream/utilities/human_mouse_ensembl_genes.txt.gz` | `os.path.join(REPO_ROOT, 'utilities/human_mouse_ensembl_genes.txt.gz')` |
| 89 | `/cellar/users/zkoch/dream/utilities/human_rat_ensembl_genes.txt` | `os.path.join(REPO_ROOT, 'utilities/human_rat_ensembl_genes.txt')` |
| 106 | `/cellar/users/zkoch/dream/utilities/human_worm_WB_ensembl_genes.tsv` | `os.path.join(REPO_ROOT, 'utilities/human_worm_WB_ensembl_genes.tsv')` |

**Action:** Add `REPO_ROOT` definition at top of file, then find/replace all paths.

---

### 4. `source/random_background_synapseHarmonization.py` (~6 paths)

| Line | Current Path | New Path |
|------|-------------|----------|
| 14 | `source_path = "/cellar/users/zkoch/dream"` | `REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))` |
| 150 | `/cellar/users/zkoch/dream/utilities/cell_cycle_genesets/cell_cycle_genes_no_dream.txt` | `os.path.join(REPO_ROOT, 'utilities/cell_cycle_genesets/cell_cycle_genes_no_dream.txt')` |
| 172-179 | `/cellar/users/zkoch/dream/data/synapse_rna_seq_harmonization/...` | `os.path.join(REPO_ROOT, 'data/synapse_rna_seq_harmonization/...')` |

**Action:** Replace `source_path` with `REPO_ROOT`, update all paths.

---

### 5. Notebooks

For each notebook, add this cell near the top (after imports):

```python
import os
# set repo root for relative paths
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname('__file__'), '..'))
```

Then replace all hard-coded paths.

#### 5.1 `notebooks/081025_dream_mice_udseq2.ipynb`

| Location | Current | New |
|----------|---------|-----|
| Cell with fig_dir | `fig_dir = "/cellar/users/zkoch/dream/figures/fig_knockout_mice"` | `fig_dir = os.path.join(REPO_ROOT, 'figures/fig_knockout_mice')` |
| id_activity_fn | `/cellar/users/zkoch/dream/data/alexandrov_collab_2025/...` | `os.path.join(REPO_ROOT, 'data/alexandrov_collab_2025/...')` |
| sbs_activity_fn | `/cellar/users/zkoch/dream/data/alexandrov_collab_2025/...` | `os.path.join(REPO_ROOT, 'data/alexandrov_collab_2025/...')` |
| savefig calls (3) | `/cellar/users/zkoch/dream/figures/supplementary/...` | `os.path.join(REPO_ROOT, 'figures/supplementary/...')` |

#### 5.2 `notebooks/111424_human_SComatic_mutations.ipynb`

| Location | Current | New |
|----------|---------|-----|
| 8 savefig() calls | `/cellar/users/zkoch/dream/figures/...` | `os.path.join(REPO_ROOT, 'figures/...')` |

#### 5.3 `notebooks/112024_synapseAD.ipynb`

| Location | Current | New |
|----------|---------|-----|
| 8 savefig() calls | `/cellar/users/zkoch/dream/figures/...` | `os.path.join(REPO_ROOT, 'figures/...')` |

#### 5.4 `notebooks/071524_mouse_lifespan.ipynb`

| Location | Current | New |
|----------|---------|-----|
| 12 savefig() calls | `/cellar/users/zkoch/dream/figures/...` | `os.path.join(REPO_ROOT, 'figures/...')` |

#### 5.5 `notebooks/022924_cross_species.ipynb`

| Location | Current | New |
|----------|---------|-----|
| ~20 savefig() calls | `/cellar/users/zkoch/dream/figures/...` | `os.path.join(REPO_ROOT, 'figures/...')` |

#### 5.6 `notebooks/062624_TMS_analysis2.ipynb`

| Location | Current | New |
|----------|---------|-----|
| ~15 savefig() calls | `/cellar/users/zkoch/dream/figures/...` | `os.path.join(REPO_ROOT, 'figures/...')` |

#### 5.7 `notebooks/120925_revisions.ipynb`

| Location | Current | New |
|----------|---------|-----|
| ensembl_to_geneSymbol path | `/cellar/users/zkoch/dream/utilities/ensembl_to_geneSymbol.tsv` | `os.path.join(REPO_ROOT, 'utilities/ensembl_to_geneSymbol.tsv')` |
| chip_df read_excel | `/cellar/users/zkoch/dream/data/litovchick_2007/mmc2.xls` | `os.path.join(REPO_ROOT, 'data/litovchick_2007/mmc2.xls')` |
| to_csv calls | `/cellar/users/zkoch/dream/data/upsetplot_reviewer_sets/...` | `os.path.join(REPO_ROOT, 'data/upsetplot_reviewer_sets/...')` |
| savefig() calls | `/cellar/users/zkoch/dream/figures/...` | `os.path.join(REPO_ROOT, 'figures/...')` |

#### 5.8 `notebooks/122225_revisions_TMS.ipynb`

| Location | Current | New |
|----------|---------|-----|
| ~4 savefig() calls | `/cellar/users/zkoch/dream/figures/...` | `os.path.join(REPO_ROOT, 'figures/...')` |

---

## Execution Order

1. **Source files first** (they may be imported by notebooks):
   - [x] `source/read_data.py`
   - [x] `source/expr_dataset.py`
   - [x] `source/methyl_dataset.py`
   - [x] `source/random_background_synapseHarmonization.py`
   - [x] `source/utils.py` (additional file not in original plan)
   - [x] `source/sc_expr_dataset.py` (additional file not in original plan)

2. **Notebooks second**:
   - [x] `notebooks/081025_dream_mice_udseq2.ipynb`
   - [x] `notebooks/111424_human_SComatic_mutations.ipynb`
   - [x] `notebooks/112024_synapseAD.ipynb`
   - [x] `notebooks/071524_mouse_lifespan.ipynb`
   - [x] `notebooks/022924_cross_species.ipynb`
   - [x] `notebooks/062624_TMS_analysis2.ipynb`
   - [x] `notebooks/120925_revisions.ipynb`
   - [x] `notebooks/122225_revisions_TMS.ipynb`
   - [x] `notebooks/010324_cptac3_analysis.ipynb` (additional file not in original plan)

3. **Python scripts** (additional files not in original plan):
   - [x] `python_scripts/random_background_dreamVmut_corr_SEAAD.py`
   - [x] `python_scripts/run_sigprofiler_extractor.py`
   - [x] `python_scripts/random_background_synapseHarmonization.py`
   - [x] `python_scripts/run_cosmic_fit_with_ffpe.py`
   - [x] `python_scripts/random_background_dreamVmut_corr_TMS.py`
   - [x] `python_scripts/calc_ssgsea_dream_activity.py`

---

## Notes

- The stderr output in notebooks showing `/cellar/users/zkoch/miniconda3/...` paths are from Python warnings and don't need to be edited - they will change automatically when run in a different environment.
- Some paths may be in commented-out code - these should still be updated for consistency.
- After edits, run each notebook to verify paths resolve correctly.
