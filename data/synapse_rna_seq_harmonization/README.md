Install synpase client into dream_proj_env
`pip install synapseclient`

foloow these (https://python-docs.synapse.org/tutorials/authentication/) instruction to creaet access token and login

Use download.py to download what is in your cart

Then move it out of subdirs with:
`find . -mindepth 2 -type f -exec mv {} . \;`

NOTE:
seems like the main diff between *_counts_matrix.txt and *_counts_matrix_clean.txt is that non clean has syn IDs and clean has other types of IDs


ROSMAP
- metadata in RNAseq_Harmonization_ROSMAP_combined_metadata.csv 
- RNA-seq which has sample names matching metadata in the ROSMAP*_clean.txt