Install synpase client into dream_proj_env
`pip install synapseclient`

foloow these (https://python-docs.synapse.org/tutorials/authentication/) instruction to creaet access token and login

Use download.py to download what is in your cart

Then move it out of subdirs with:
`find . -mindepth 2 -type f -exec mv {} . \;`

File descriptions
- all_brain_regions_filt_preprocessed_scanpy_norm.final_noMB.cell_labels.tsv: cells labels and umap coordinates
- all_brain_regions_filt_preprocessed_scanpy_fullmatrix.h5ad: full matrix of sc data
- MIT_ROSMAP_Multiomics_assay_snRNAseq_metadata.csv: mapping of each sample (cell?) to its sequencing method
- 

sc_rosmap_preprocessed.pkl: preprocessed with ScExpressionDataset.preprocess