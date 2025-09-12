#!/bin/sh

# use STAR to align .fastq.gz files to reference genome
raw_data_dir=/cellar/users/zkoch/dream/data/jung_2015/raw_RNA_seq
ref_dir=/cellar/users/zkoch/dream/data/reference_files
star_idx_dir=/cellar/users/zkoch/dream/data/jung_2015/STAR_index
aligned_results_dir=/cellar/users/zkoch/dream/data/jung_2015/STAR_results_two_pass
fastq_fn=$1

# get basename of fastq_fn
base_fn=$(basename $fastq_fn)
# remove .fastq.gz
base_fn="${base_fn%.fastq}"
echo $base_fn

STAR --genomeDir $star_idx_dir \
    --runThreadN 2 \
    --readFilesIn $fastq_fn \
    --outFileNamePrefix $aligned_results_dir/$base_fn \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --twopassMode Basic

        
