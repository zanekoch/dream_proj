#!/usr/bin/env python
"""
Fix malformed deletion entries in VCF files after CrossMap liftover
Addresses the issue where deletions have empty ALT alleles and single-base REF alleles
"""

import sys
import argparse
import pysam
from pathlib import Path

def fix_deletion_format(vcf_in_path, vcf_out_path, reference_fasta):
    """
    Fix malformed deletion entries in VCF file
    
    Parameters:
    -----------
    vcf_in_path : str
        Path to input VCF file with malformed deletions
    vcf_out_path : str
        Path to output fixed VCF file
    reference_fasta : str
        Path to reference genome FASTA file (mm10)
    """
    
    # Open reference genome
    ref = pysam.FastaFile(reference_fasta)
    
    # Open input VCF
    vcf_in = pysam.VariantFile(vcf_in_path)
    
    # Create output VCF with same header
    vcf_out = pysam.VariantFile(vcf_out_path, 'w', header=vcf_in.header)
    
    fixed_count = 0
    total_deletions = 0
    
    for record in vcf_in:
        needs_fix = False
        
        # Check if this is a malformed deletion
        # Symptoms: empty ALT allele or ALT is "." and REF is single base
        if record.alts:
            for alt in record.alts:
                if alt == "" or (alt == "." and len(record.ref) == 1):
                    needs_fix = True
                    total_deletions += 1
                    break
        
        if needs_fix:
            try:
                # For deletions, we need to:
                # 1. Get the deleted sequence from INFO field if available
                # 2. Or reconstruct from original coordinates
                
                # Try to get deletion length from INFO field
                del_len = None
                if 'SVLEN' in record.info:
                    del_len = abs(record.info['SVLEN'])
                elif 'END' in record.info:
                    del_len = record.info['END'] - record.pos
                
                if del_len and del_len > 0:
                    # Reconstruct proper deletion format
                    # Get reference sequence including the deletion
                    chrom = record.chrom
                    # Position in VCF is 1-based, pysam expects 0-based
                    start = record.pos - 1
                    end = start + del_len + 1  # +1 to include anchor base
                    
                    # Fetch reference sequence
                    ref_seq = ref.fetch(chrom, start, end).upper()
                    
                    if len(ref_seq) > 1:
                        # Create new record with proper format
                        # REF should be the full deleted sequence plus anchor
                        # ALT should be just the anchor base
                        record.ref = ref_seq
                        record.alts = (ref_seq[0],)  # First base only for deletion
                        fixed_count += 1
                        print(f"Fixed deletion at {chrom}:{record.pos} - REF: {ref_seq}, ALT: {ref_seq[0]}")
                
            except Exception as e:
                print(f"Warning: Could not fix deletion at {record.chrom}:{record.pos} - {e}")
        
        # Write record (fixed or original)
        vcf_out.write(record)
    
    vcf_in.close()
    vcf_out.close()
    ref.close()
    
    print(f"\nFixed {fixed_count} out of {total_deletions} malformed deletions")
    return fixed_count, total_deletions


def validate_vcf(vcf_path):
    """
    Validate that VCF deletions are properly formatted
    """
    vcf = pysam.VariantFile(vcf_path)
    
    issues = []
    for record in vcf:
        if record.alts:
            for alt in record.alts:
                # Check for empty or malformed deletion
                if alt == "" or alt == ".":
                    issues.append(f"{record.chrom}:{record.pos} has empty ALT allele")
                elif len(record.ref) > len(alt) and len(record.ref) == 1:
                    issues.append(f"{record.chrom}:{record.pos} deletion has single-base REF")
    
    vcf.close()
    return issues


def main():
    parser = argparse.ArgumentParser(description='Fix malformed deletions in VCF after CrossMap liftover')
    parser.add_argument('input_vcf', help='Input VCF file with malformed deletions')
    parser.add_argument('output_vcf', help='Output fixed VCF file')
    parser.add_argument('reference', help='Reference genome FASTA file (mm10)')
    parser.add_argument('--validate', action='store_true', help='Validate output VCF after fixing')
    
    args = parser.parse_args()
    
    # Fix the VCF
    print(f"Processing {args.input_vcf}...")
    fixed, total = fix_deletion_format(args.input_vcf, args.output_vcf, args.reference)
    
    # Optionally validate
    if args.validate:
        print(f"\nValidating {args.output_vcf}...")
        issues = validate_vcf(args.output_vcf)
        if issues:
            print(f"Found {len(issues)} remaining issues:")
            for issue in issues[:10]:  # Show first 10
                print(f"  - {issue}")
        else:
            print("Validation passed - no malformed deletions found!")
    
    print(f"\nFixed VCF written to: {args.output_vcf}")


if __name__ == "__main__":
    main()