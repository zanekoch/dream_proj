
# repo root for relative paths
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#!/usr/bin/env python
"""
SigProfilerExtractor analysis for UDSeq mice mutation data
Extracts mutational signatures from SBS, DBS, and ID matrices
"""

import os
import glob
import sys
from SigProfilerExtractor import sigpro as sig
from SigProfilerMatrixGenerator import install as genInstall

def run_signature_extraction(matrix_path, output_name, context_type_desc):
    """Run SigProfilerExtractor on a single matrix file"""
    print(f"\nRunning SigProfilerExtractor on {context_type_desc}...")
    print(f"Input path: {matrix_path}")
    print(f"Output name: {output_name}")
    
    try:
        sig.sigProfilerExtractor(
            "matrix",                               # Input type: matrix files
            output_name,                           # Output folder name
            matrix_path,                           # Path to matrix file
            reference_genome="mm10",               # Reference genome (mm10 for mouse)
            minimum_signatures=1,                  # Minimum number of signatures
            maximum_signatures=15,                 # Maximum number of signatures  
            nmf_replicates=100,                   # Number of NMF iterations for stability
            resample=True,                        # Bootstrap resampling
            cpu=8,                               # Use available CPUs
            gpu=False,                           # Don't use GPU
            nmf_init="random",                   # Random initialization
            precision="single",                  # Single precision for speed
            matrix_normalization="gmm",          # Gaussian mixture model normalization
            seeds="random",                      # Random seeds
            min_nmf_iterations=10000,            # Minimum NMF iterations
            max_nmf_iterations=1000000,          # Maximum NMF iterations
            nmf_test_conv=10000,                 # Convergence test iterations
            nmf_tolerance=1e-15,                 # Convergence tolerance
            get_all_signature_matrices=False     # Don't save all intermediate matrices
        )
        
        print(f"SigProfilerExtractor completed successfully for {context_type_desc}!")
        print(f"Results saved in: {output_name}/")
        return True
        
    except Exception as e:
        print(f"ERROR running SigProfilerExtractor on {context_type_desc}: {e}")
        return False

def main():
    print("Starting SigProfilerExtractor analysis for SBS, DBS, and ID matrices...")
    
    # Install mm10 reference genome if not already installed
    print("\nInstalling/checking mm10 reference genome...")
    try:
        genInstall.install('mm10')
        print("mm10 reference genome ready")
    except Exception as e:
        print(f"Error installing mm10 genome: {e}")
        sys.exit(1)
    
    # Define matrix paths and output names
    """
    # not ffpe corrected mutations
    base_path = os.path.join(REPO_ROOT, "data/alexandrov_collab_2025/dupcaller_output/vcfs/combined_vcfs/output")
    matrices = [
        {
            "path": f"{base_path}/SBS/mouse_udseq.SBS96.all",
            "output": "sbs_signatures",
            "desc": "SBS96 (Single Base Substitutions)"
        },
        {
            "path": f"{base_path}/DBS/mouse_udseq.DBS78.all", 
            "output": "dbs_signatures",
            "desc": "DBS78 (Dinucleotide Substitutions)"
        },
        {
            "path": f"{base_path}/ID/mouse_udseq.ID83.all",
            "output": "id_signatures", 
            "desc": "ID83 (Insertions/Deletions)"
        }
    ]"""
    
    # ffpe corrected mutations
    base_path = os.path.join(REPO_ROOT, "data/alexandrov_collab_2025/FFPEsig/FFPEsig_OUTPUT/")
    matrices = [
        {
            "path": f"{base_path}/corrected_sbs_df_for_sigprofiler.csv",
            "output": "sbs_signatures",
            "desc": "SBS96 (Single Base Substitutions) - FFPE corrected"
        }
    ]
    
    success_count = 0
    
    # Run signature extraction for each matrix type
    for matrix in matrices:
        if os.path.exists(matrix["path"]):
            if run_signature_extraction(matrix["path"], matrix["output"], matrix["desc"]):
                success_count += 1
        else:
            print(f"WARNING: Matrix file not found: {matrix['path']}")
    
    print(f"\n=== SUMMARY ===")
    print(f"Successfully processed {success_count}/{len(matrices)} matrix types")
    
    if success_count > 0:
        print(f"Signature extraction results saved in respective output directories:")
        for matrix in matrices:
            if os.path.exists(matrix["path"]):
                print(f"  - {matrix['desc']}: {matrix['output']}/")
    else:
        print("No matrices were successfully processed")
        sys.exit(1)

if __name__ == "__main__":
    main()