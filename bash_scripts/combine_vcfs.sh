#!/bin/bash

# Create output directory
OUTPUT_DIR=~/dream/data/alexandrov_collab_2025/dupcaller_output/vcfs/combined_vcfs
mkdir -p "$OUTPUT_DIR"

# Change to VCF directory
cd ~/dream/data/alexandrov_collab_2025/dupcaller_output/vcfs/

# Get unique sample names (excluding mm10 variants and bcfnormed files)
SAMPLES=$(ls *.vcf | grep -E "^S[0-9]+_(snv|indel)\.vcf$" | sed 's/_\(snv\|indel\)\.vcf//' | sort -u)

for sample in $SAMPLES; do
    echo "Processing sample: $sample"
    
    snv_file="${sample}_snv.vcf"
    indel_file="${sample}_indel.vcf"
    output_file="$OUTPUT_DIR/${sample}_combined.vcf"
    
    if [[ -f "$snv_file" && -f "$indel_file" ]]; then
        # Extract header from SNV file
        grep "^#" "$snv_file" > "$output_file"
        
        # Add non-header lines from SNV file
        grep -v "^#" "$snv_file" >> "$output_file"
        
        # Add non-header lines from INDEL file
        grep -v "^#" "$indel_file" >> "$output_file"
        
        echo "  Created: $output_file"
    else
        echo "  Warning: Missing files for $sample"
        if [[ ! -f "$snv_file" ]]; then
            echo "    Missing: $snv_file"
        fi
        if [[ ! -f "$indel_file" ]]; then
            echo "    Missing: $indel_file"
        fi
    fi
done

echo "Done! Combined VCF files are in: $OUTPUT_DIR"