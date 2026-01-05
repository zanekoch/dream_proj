# use cellXgene2 conda env
import pandas as pd
import scanpy as sc
import os
from matplotlib import pyplot as plt
import numpy as np
import statsmodels.formula.api as smf
import random
import gseapy
import sys
import argparse
import pickle

# repo root for relative paths
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# add source directory to path
if os.path.join(REPO_ROOT, 'source') not in sys.path:
    sys.path.append(os.path.join(REPO_ROOT, 'source'))
# read source files
import read_data

def get_ssgsea(
    gene_set: list, expr_object: sc.AnnData, regress: bool, ssgsea : bool, n_cores: int
    ):
    genes_w_expr = expr_object.columns[expr_object.columns.str.startswith('ENSG')]
    gene_set_w_expr = list(set(gene_set).intersection(genes_w_expr))
    # get the mean expression of these genes across all cells
    if ssgsea == False:
        expr_object['expr'] = expr_object[gene_set_w_expr].mean(axis=1)
    # or do ssgsea
    else:
        gene_set_dict = {
            'gene_set': gene_set
            }
        ssgsea = gseapy.ssgsea(
            data= expr_object[genes_w_expr].T, # make sure we don't include meta cols
            gene_sets=gene_set_dict,
            outdir=None, no_plot=True, verbose = False,
            threads = n_cores
            )
        results_df = ssgsea.res2d
        results_df.set_index('Name', inplace=True)
        results_df.drop(columns = ['Term','ES'], inplace=True)
        results_df.rename(
            columns = {'NES':'expr'},
            inplace=True
            )
        # convert cols to float
        results_df = results_df.astype(float)
        # merge with expr_object
        expr_object['expr'] = results_df['expr']
    if regress:
        expr_ols = smf.ols(
            formula='expr ~ n_genes_expressed * total_seq_depth', data=expr_object
            ).fit()
        expr_object['expr_resid'] = expr_ols.resid
        expr_object['expr_resid_scaled'] = -1 * expr_object['expr_resid']
        expr_object['expr_resid_scaled'] = expr_object['expr_resid_scaled'] - min(expr_object['expr_resid_scaled'])
        to_return = expr_object['expr_resid_scaled'].copy(deep=True)
        # drop
        expr_object.drop(['expr_resid', 'expr_resid_scaled'], axis=1, inplace=True)
        return to_return
    else:
        expr_object['expr_scaled'] = -1 * expr_object['expr']
        expr_object['expr_scaled'] = expr_object['expr_scaled'] - min(expr_object['expr_scaled'])
        to_return = expr_object['expr_scaled'].copy(deep=True)
        # drop
        expr_object.drop(['expr','expr_scaled'], axis=1, inplace=True)
        return to_return
    

def create_random_background(
    num_iter,
    expr_object,
    regress,
    ssgsea,
    n_dream_genes_w_expr,
    choose_from_background = [],
    n_cores = 1
    ):
    
    # randomly choose sets of n_dream_genes_w_expr from expressed genes
    random.seed(42)
    if choose_from_background == []:
        genes_to_choose_from = set(expr_object.columns[expr_object.columns.str.startswith('ENSG')])
    else:
        genes_to_choose_from = set(choose_from_background)
        print(f"using {len(genes_to_choose_from)} provided genes to choose from")
        
    # do ssgsea on each random set of genes
    random_results = []
    for i in range(num_iter):
        # choose random genes
        random_genes = random.sample(list(genes_to_choose_from), n_dream_genes_w_expr)
        # do ssgsea
        random_results.append(get_ssgsea(random_genes, expr_object, regress, ssgsea, n_cores))
        if i % 10 == 0:
            print(f"finished {i} of {num_iter}", flush=True)
        
    random_background_results_df = pd.concat(random_results, axis = 1)
    random_background_results_df.columns = [f'random_{n}' for n in range(len(random_results))]
    return random_background_results_df

def main():
    # create argument parser
    parser = argparse.ArgumentParser(
        description='Create random background for DREAM vs mutation burden correlation'
        )
    parser.add_argument(
        '--start_sample_num', type=int,
        help='Start sample number for random background generation'
        )
    parser.add_argument(
        '--end_sample_num', type=int,
        help='End sample number for random background generation'
        )
    parser.add_argument(
        '--regress', type=str, default='True',
        help='Whether to regress out n_counts and n_genes'
        )
    parser.add_argument(
        '--use-cell-cycle-genes', type=str, default='False',
        help='Whether to use cell cycle genes for random background generation'
        )
    parser.add_argument(
        '--n-cores', type=int, default=1,
        help='Number of cores to use for random background generation'
        )
    # parge
    args = parser.parse_args()
    start_sample_num = args.start_sample_num
    end_sample_num = args.end_sample_num
    regress = args.regress
    use_cell_cycle_genes = args.use_cell_cycle_genes
    n_cores = args.n_cores
    if regress == 'True':
        regress = True
    else:
        regress = False
    print(f"Using {n_cores} cores")
    
    # load expression data
    print("loading expression data", flush=True)
    loader = read_data.DatasetLoader("synapse_rna_seq_harmonization")
    harmonization = loader.load_dataset()
    harmonization.calc_total_seq_depth()
    harmonization.get_dream_gene_expression()
    
    choose_from_background = []
    if use_cell_cycle_genes == 'True':
        cell_cycle_genes = pd.read_csv(
            os.path.join(REPO_ROOT, "utilities/cell_cycle_genesets/cell_cycle_genes_no_dream.txt"),
            header=None
            )[0].values.tolist()
        choose_from_background = cell_cycle_genes
            
    print("generating random backgrounds", flush=True)
    # 1000 cells 500 times took 1 hour
    # so 10k cells 100 times should take 2
    
    background_df = create_random_background(
        num_iter = 500,
        expr_object = harmonization.expression_df.iloc[start_cell_num:end_cell_num],
        regress = regress,
        ssgsea = True,
        n_dream_genes_w_expr = len(harmonization.dream_regulated_genes_w_expression),
        choose_from_background = choose_from_background,
        n_cores = n_cores
        )
    
    # save to parquet
    if regress:
        background_df.to_parquet(
            os.path.join(REPO_ROOT, f'data/synapse_rna_seq_harmonization/random_background/sea-ad_random_background500iter_{start_cell_num}-{end_cell_num}cells.parquet')
            )
        print(f"wrote to {os.path.join(REPO_ROOT, f'data/synapse_rna_seq_harmonization/random_background/sea-ad_random_background500iter_{start_cell_num}-{end_cell_num}cells.parquet')}")
    else:
        background_df.to_parquet(
            os.path.join(REPO_ROOT, f'data/synapse_rna_seq_harmonization/random_background/synapse_random_background500iter_{start_cell_num}-{end_cell_num}cells_noregress.parquet')
            )
        print(f"wrote to {os.path.join(REPO_ROOT, f'data/synapse_rna_seq_harmonization/random_background/synapse_random_background500iter_{start_cell_num}-{end_cell_num}cells_noregress.parquet')}")

if __name__ == "__main__":
    main()