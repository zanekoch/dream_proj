#!/bin/bash

# VCF Liftover Script: mm39 to mm10
# This script lifts over VCF coordinates from mouse genome mm39 to mm10 using CrossMap
# Processes all VCF files matching the pattern *.vcf

# Set variables
INPUT_PATTERN="./*.vcf"  # Pattern to match VCF files
CHAIN_FILE="/cellar/users/zkoch/dream/utilities/mm39ToMm10.over.chain.gz"
REFERENCE_FASTA="/cellar/users/zkoch/dream/utilities/mm10.fa"

# Create output directory
mkdir -p liftover_output

# Download required files if they don't exist
echo "Checking for required files..."

if [ ! -f "$CHAIN_FILE" ]; then
    echo "Downloading chain file..."
    wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/liftOver/mm39ToMm10.over.chain.gz
fi

if [ ! -f "$REFERENCE_FASTA" ]; then
    echo "Downloading mm10 reference genome..."
    echo "Note: You may want to download individual chromosomes or use a pre-existing reference"
    echo "Example: wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"
    echo "Then: gunzip mm10.fa.gz"
    echo "Please ensure you have the mm10 reference FASTA file before running CrossMap"
fi

# Check if CrossMap is installed
if ! command -v CrossMap &> /dev/null; then
    echo "CrossMap is not installed. Please install it using:"
    echo "conda install -c bioconda crossmap"
    echo "or"
    echo "pip install CrossMap"
    exit 1
fi

# Check if any VCF files exist
vcf_files=($INPUT_PATTERN)
if [ ! -f "${vcf_files[0]}" ]; then
    echo "No VCF files found matching pattern: $INPUT_PATTERN"
    exit 1
fi

echo "Found ${#vcf_files[@]} VCF file(s) to process"

# Process each VCF file
for vcf_file in $INPUT_PATTERN; do
    if [ -f "$vcf_file" ]; then
        echo "Processing: $vcf_file"
        
        # Generate output filenames based on input filename
        base_name=$(basename "$vcf_file" .vcf)
        output_vcf="${base_name}_mm10.vcf"
        unmapped_file="${base_name}_unmapped.vcf"
        
        # Run CrossMap VCF liftover
        echo "Running CrossMap VCF liftover from mm39 to mm10 for $vcf_file..."
        CrossMap vcf \
                "$CHAIN_FILE" \
                "$vcf_file" \
                "$REFERENCE_FASTA" \
                "liftover_output/$output_vcf"
        
        # Check if liftover was successful
        if [ $? -eq 0 ]; then
            echo "Liftover completed successfully for $vcf_file!"
            echo "Output VCF: liftover_output/$output_vcf"
            echo "Unmapped variants: liftover_output/$unmapped_file"
            
            # Optional: compress output VCF
            echo "Compressing output VCF..."
            bgzip "liftover_output/$output_vcf"
            tabix -p vcf "liftover_output/${output_vcf}.gz"
            
            echo "Compressed VCF: liftover_output/${output_vcf}.gz"
            echo "Index file: liftover_output/${output_vcf}.gz.tbi"
            echo "---"
        else
            echo "Liftover failed for $vcf_file. Continuing with next file..."
            echo "---"
        fi
    fi
done

echo "All liftover processes complete!"
