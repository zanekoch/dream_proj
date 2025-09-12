#!/bin/sh

# use STAR to align .fastq.gz files to reference genome
raw_data_dir=/cellar/users/zkoch/dream/data/jung_2015/raw_RNA_seq
ref_dir=/cellar/users/zkoch/dream/data/reference_files
star_idx_dir=/cellar/users/zkoch/dream/data/jung_2015/STAR_index
aligned_results_dir=/cellar/users/zkoch/dream/data/jung_2015/STAR_results


# Step 1: create genome indexes
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $star_idx_dir  --genomeFastaFiles $ref_dir/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile $ref_dir/Homo_sapiens.GRCh38.111.gtf --sjdbOverhang 100