
# repo root for relative paths
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#!/usr/bin/env python3
"""
Script to run SigProfilerAssignment cosmic_fit with custom FFPE + COSMIC signatures
"""

import os
import sys
from SigProfilerAssignment import Analyzer as Analyze

def main():
    # Input files
    signature_database = os.path.join(REPO_ROOT, "data/alexandrov_collab_2025/output_SigProfilerAssignment/cosmic_and_ffpe_signatures_with_SBS87.tsv")
    input_matrix = os.path.join(REPO_ROOT, "data/alexandrov_collab_2025/dupcaller_output/vcfs/mouse_udseq.SBS96.all")
    
    # Output directory
    output_dir = os.path.join(REPO_ROOT, "data/alexandrov_collab_2025/output_SigProfilerAssignment/cosmic_fit_results_only_extracted87")
    
    # Check if input files exist
    if not os.path.exists(signature_database):
        print(f"Error: Signature database file not found: {signature_database}")
        sys.exit(1)
    
    if not os.path.exists(input_matrix):
        print(f"Error: Input matrix file not found: {input_matrix}")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    print("Running SigProfilerAssignment cosmic_fit...")
    print(f"Signature database: {signature_database}")
    print(f"Input matrix: {input_matrix}")
    print(f"Output directory: {output_dir}")
    
    try:
        # Run cosmic_fit with custom signature database
        # Note: No explicit CPU parameter found in documentation
        # The tool may use internal parallelization
        Analyze.cosmic_fit(
            samples=input_matrix,
            output=output_dir,
            input_type="matrix",
            signature_database=signature_database,
            context_type="96",
            verbose=True
        )
        
        print("SigProfilerAssignment completed successfully!")
        print(f"Results saved to: {output_dir}")
        
    except Exception as e:
        print(f"Error running SigProfilerAssignment: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()