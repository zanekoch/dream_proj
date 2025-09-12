## CPTAC genomic Data downloading
### Step 1: Data using TCGABiolinks
Run TCGA_data_downloader.R (Use conda environment `r_go`) to get the desired project and datatypes, populating a `GDCdata` folder with the data. 
### Step 2: Metadata/samplesheet
Select the same samples on the TCGA GDC website. Add them to cart. Go to cart. Download clinical (aka metdata) and sample sheet. Upload using sftp to this cluster. Unpack the clinical and use clinical.csv. This creates `/dream/data/tcga/manual_download/CPTAC-3/clinical.tsv` and `/dream/data/tcga/manual_download/CPTAC-3/sample_sheet.tsv`. Remove rows of the sample sheet that are indistinguishable but map to different files. Currently unknown what causes this but my choice was to keep the first row and remove the rest. This creates `/dream/data/tcga/manual_download/CPTAC-3/sample_sheetsample_sheet_no_duplicate_rows.tsv`.
### Step 3: Data processing
Run `~/dream/source/process_tcgabiolinks.py` (Use conda environment `dream_proj`) with the unduplicated samplesheet, metadata, and path to each project-datatype pair to create data matrices with synchronized sample IDs and such. Outputs to `/dream/data/tcga/processed/CPTAC-3_mutation.parquet` and `/dream/data/tcga/processed/CPTAC-3_expression.parquet` etc.
### Step 4: Manual BS
Do they steps in `/cellar/users/zkoch/dream/notebooks/010224_process_tcga_data.ipynb` to clean the data and form a nice metadata table. This creates `/dream/data/tcga/processed/CPTAC-3_metadata.parquet`.
### Step 5: Protoemtics
- Download the pdc client `wget https://pdc-download-clients.s3.amazonaws.com/pdc-client_v1.0.7_Ubuntu_x64.zip`
- log2 intensity values for kidney cancers were downloaded from `https://proteomic.datacommons.cancer.gov/pdc/analysis/dd0a228f-1fb3-11e9-b7f8-0a80fada099c?StudyName=CPTAC%20CCRCC%20Discovery%20Study%20-%20Phosphoproteome` as cptac_phospo_ccrcc.gct
- mapping from log2 intensity value aliquot_id to sample_ids was downloaded from `https://proteomic.datacommons.cancer.gov/pdc/study/PDC000128` the biospecimen page 


### Log
#### CPTAC-3: 
Added methylation, somatic mutaiton, and expression files to cart. Downloaded sample sheet and clinical (1,235 samples and 5,646 files) to `tcga/manual_download/CPTAC-3`.
##### Cell cycle processing
Done using scanpy in `/cellar/users/zkoch/dream/notebooks/010824_cellCycle_deconv.ipynb` and the cell cycle assignments are in `/cellar/users/zkoch/dream/data/tcga/processed_data/CPTAC-3_cell_cycle_scores.parquet` (using cellXgene2 environment)
Done using scanpy in `/cellar/users/zkoch/dream/notebooks/010824_cellCycle_deconv.ipynb` and the cell cycle assignments are in `/cellar/users/zkoch/dream/data/tcga/processed_data/CPTAC-3_cell_cycle_scores.parquet` (using cellXgene2 environment)