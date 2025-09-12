Install synpase client into dream_proj_env
`pip install synapseclient`

foloow these (https://python-docs.synapse.org/tutorials/authentication/) instruction to creaet access token and login

Use download.py to download what is in your cart

Then move it out of subdirs with:
`find . -mindepth 2 -type f -exec mv {} . \;`

Files are provided per chromosome

File descriptions
- File type Description
*.recalibrated_variants.annotated.clinical.txt	All low frequency HIGH/MODERATE annotated variants with possible clinical impact (from ClinVar ) in text file format
*.recalibratedvariants.annotated.codingrare.txt	All HIGH/MODERATE annotated variants with less than 5% allele frequency in 1000genomes and ExAC in text file format
*.recalibrated_variants.annotated.coding.txt	All annotated variants with HIGH/MODERATE impact in text file format
*.recalibrated_variants.annotated.txt	All variants with annotations in text file format
*.recalibrated_variants.annotated.vcf.gz	All variants with annotations in variant call format (VCF)
*.recalibrated_variants.annotated.vcf.gz.tbi	Index file for annotated VCF file
*.recalibrated_variants.vcf.gz	All variants in variant call format
*.recalibrated_variants.vcf.gz.tbi	Index file for VCF file with all variants