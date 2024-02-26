Each number name file is on .h5ad downloaded from the figshare (https://figshare.com/articles/dataset/count-data-brain-annotated_h5ad/17111009)
- these were then combined, the TMS samples removed, and written to all_tissues_only_parabiosis_no_tms.h5ad

# Overview

* Data is  available in `h5ad` and `rds` formats
* Cell quality control thresholds:
   * cells with less than 500 genes or 5k reads are discarded
   * mito/ribo/ercc read thresholds are set for cell QC, 10%,10%,30% respectively
* Data is Gencode vM19 aligned and merged with cells from TMS male samples
* Count data includes the raw counts without any transformation
* Cell attributes:
  * `facs_selection`
  * `plate`: plate id
  * `pair`: pair id for parabionts, `TMS` indicated for mice from the aging cohort
  * `mouse_id`
  * `data`: `pb` or `tms` indicating the data source
  * `condition` : `Y`,`A`,`IY`,`HY`,`IA`,`HA` indicating the 6 different conditions within the data
  * `tissue`
  * `subtissue`: indicates the different regions of the brain. 'nan' in case of other tissues
  * `age`
  * `cell_type`: cell type annotation