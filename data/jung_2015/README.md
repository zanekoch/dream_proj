From jung et al. 2015 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4551908/)
2 samples ~10 years apart from 10 older humans. RNA-seq + DNAm + histone marks. 

# ./raw_RNA_seq 
fastq files. Downloaded from https://www.ebi.ac.uk/ena/browser/view/PRJNA223350 with `download.sh`

# ./RNA_seq
RNA-seq files, to be put in a matrix. Downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51518

# methyl_histone_GSE51517_series_matrix.txt.gz
DNAm and histone mark data, downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi

# sample_name_mapping.tsv
mapping between GEO, SRR, and sample name, downloaded from https://www.ebi.ac.uk/ena/browser/view/PRJNA223350

# STAR_index
where index files from running star on Jung 2015 data go.
Created by running /cellar/users/zkoch/dream/bash_scripts/create_star_index.sh

# STAR_results
Where the aligned RNA-seq is.
Created by running /cellar/users/zkoch/dream/bash_scripts/align_jung_rnaseq.sh using sbatch (/cellar/users/zkoch/dream/slurm_scripts/submit_bash_jobs.job)

# RNA-seq STAR alignment
Used this tutorial https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html
