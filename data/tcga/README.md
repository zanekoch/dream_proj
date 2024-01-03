## Data downloading
### Step 1: Data using TCGABiolinks
Run TCGA_data_downloader.R to get the desired project and datatypes, populating a `GDCdata` folder with the data.
### Step 2: Metadata/samplesheet
Select the same samples on the TCGA GDC website. Add them to cart. Go to cart. Download clinical (aka metdata) and sample sheet. Upload using sftp to this cluster. Unpack the clinical and use clinical.csv 
### Step 3: Data processing
Run `/cellar/users/zkoch/dream/source/process_tcgabiolinks.py` with the samplesheet, metdata, and path to each project-datatype pair to create data matrices with synchronized sample IDs and such. Outputs to `data/tcga/processed/`.


### Log
#### CPTAC-3: Added methylation, somatic mutaiton, and expression files to cart. Downloaded sample sheet and clinical (1,235 samples and 5646 files) to `tcga/manual_download/CPTAC-3`