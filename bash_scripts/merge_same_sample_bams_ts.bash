#!/bin/bash

# Usage: ./merge_bam_files.sh /path/to/bam/files

# Check if the directory path is provided
if [ -z "$1" ]; then
    echo "Please provide the path to the directory containing BAM files."
    echo "Usage: $0 /path/to/bam/files"
    exit 1
fi

# Directory containing BAM files
bam_dir="$1"
cd "$bam_dir" || { echo "Directory not found: $bam_dir"; exit 1; }
# sample of choice to merge
sample_to_merge="$2"

# set output bam file path
output_bam="/cellar/users/zkoch/dream/data/tabula_sapiens/merged_bams/${sample_to_merge}_outs_possorted_genome_bam.bam"

# Ensure samtools is installed
module load samtools
if ! command -v samtools &> /dev/null; then
    echo "samtools could not be found. Please install samtools to use this script."
    exit 1
fi

# get list of bam files for the sample of choice
bam_files=$(ls -1 ${sample_to_merge}_*.bam)

# Perform the merge using samtools
echo "Merging BAM files for sample ${sample_to_merge}:"
echo "$bam_files"
samtools merge -@ 16 --write-index -f "$output_bam" $bam_files

# check if the merge was successful
if [ $? -eq 0 ]; then
    echo "Successfully created $output_bam"
else
    echo "Error merging files for sample $sample_to_merge"
fi
