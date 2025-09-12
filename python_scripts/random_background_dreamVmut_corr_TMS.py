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

# add source directory to path
source_path = "/cellar/users/zkoch/dream"#os.path.abspath(os.path.join('..'))
if source_path not in sys.path:
    sys.path.append(os.path.join(source_path, 'source'))
# read source files
import read_data

def get_ssgsea(
    gene_set: list, expr_object: sc.AnnData, regress: bool, ssgsea : bool
    ):
    genes_w_expr = expr_object.var_names
    gene_set_w_expr = list(set(gene_set).intersection(genes_w_expr))
    # get the mean expression of these genes across all cells
    if ssgsea == False:
        expr_object.obs['expr'] = expr_object[:, gene_set_w_expr].X.mean(axis=1)
    else:
        # get dream gene to check for enrichment in expression
        gene_set_dict = {
            'gene_set': gene_set
            }
        # run ssgsea
        
        # get number of system threads
        threads = os.cpu_count()
        df = expr_object.to_df().T
        ssgsea = gseapy.ssgsea(
            data=df, # get a df version of the expr
            gene_sets=gene_set_dict,
            outdir=None, no_plot=True, verbose = False,
            threads = threads
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
        expr_object.obs['expr'] = results_df['expr']
    if regress:
        expr_ols = smf.ols(
            formula='expr ~ n_counts * n_genes ', data=expr_object.obs
            ).fit()
        expr_object.obs['expr_resid'] = expr_ols.resid
        expr_object.obs['expr_resid_scaled'] = -1 * expr_object.obs['expr_resid']
        expr_object.obs['expr_resid_scaled'] = expr_object.obs['expr_resid_scaled'] - min(expr_object.obs['expr_resid_scaled'])
        to_return = expr_object.obs['expr_resid_scaled'].copy(deep=True)
        # drop
        expr_object.obs.drop(['expr_resid', 'expr_resid_scaled'], axis=1, inplace=True)
        return to_return
    else:
        expr_object.obs['expr_scaled'] = -1 * expr_object.obs['expr']
        expr_object.obs['expr_scaled'] = expr_object.obs['expr_scaled'] - min(expr_object.obs['expr_scaled'])
        to_return = expr_object.obs['expr_scaled'].copy(deep=True)
        # drop
        expr_object.obs.drop(['expr','expr_scaled'], axis=1, inplace=True)
        return to_return
    

def create_random_background(num_iter, expr_object, regress, ssgsea, choose_from_background = []):
    n_dream_genes_w_expr = 439 #len(dream_genes_w_expr) 
    # randomly choose sets of n_dream_genes_w_expr from facs genes
    random.seed(42)
    if choose_from_background == []:
        genes_to_choose_from = set(expr_object.var_names)
    else:
        genes_to_choose_from = set(choose_from_background)
        print(f"using {len(genes_to_choose_from)} provided genes to choose from")
        
    random_expressions = []
    for i in range(num_iter):
        # choose random genes
        random_genes = random.sample(list(genes_to_choose_from), n_dream_genes_w_expr)
        random_expressions.append(get_ssgsea(random_genes, expr_object, regress, ssgsea))
        if i % 50 == 0:
            print(f"finished {i} of {num_iter}", flush=True)
        
    random_background_mean_expression_df = pd.concat(random_expressions, axis = 1)
    random_background_mean_expression_df.columns = [f'random_{n}' for n in range(len(random_expressions))]
    # add necessary cols
    #random_background_mean_expression_df['mutation_count_resid'] = expr_object.obs['mutation_count_resid']
    #random_background_mean_expression_df['mutation_count_resid_scaled'] = expr_object.obs['mutation_count_resid_scaled']
    random_background_mean_expression_df['age_months'] = expr_object.obs['age_months']
    random_background_mean_expression_df['tissue'] = expr_object.obs['tissue']
    random_background_mean_expression_df['subtissue'] = expr_object.obs['subtissue']
    random_background_mean_expression_df['cell_type'] = expr_object.obs['cell_type']
    random_background_mean_expression_df['n_genes'] = expr_object.obs['n_genes']
    random_background_mean_expression_df['n_counts'] = expr_object.obs['n_counts']
    # comment/uncomment for facs_raw/facs

    #random_background_mean_expression_df['mean_ddr_gene_expr_resid'] = expr_object.obs['mean_ddr_gene_expr_resid']

    return random_background_mean_expression_df

def main():
    # create argument parser
    parser = argparse.ArgumentParser(description='Create random background for DREAM vs mutation burden correlation')
    parser.add_argument(
        '--start_cell_num', type=int,
        help='Start cell number for random background generation'
        )
    parser.add_argument(
        '--end_cell_num', type=int,
        help='End cell number for random background generation'
        )
    parser.add_argument(
        '--regress', type=str, default='True',
        help='Whether to regress out n_counts and n_genes'
        )
    # parge
    args = parser.parse_args()
    start_cell_num = args.start_cell_num
    end_cell_num = args.end_cell_num
    regress = args.regress
    if regress == 'True':
        regress = True
    else:
        regress = False
    
    # load data
    loader = read_data.DatasetLoader("TMS")
    tms = loader.load_dataset()
    """All done in loader now
    print("calculating mutation burden", flush=True)
    tms.calculate_mutation_burden('mutation_count_per_kb_top50expr', top_perc_expr=.5, max_vaf = 0.6)
    tms.get_dream_gene_expression()
    # get normalized enrichment score (but actually already calcd)
    tms.dream_enrichment_ssgsea()"""
    
    dream_col = 'DREAM_normalized_enrichment_score_resid'
    mut_col = 'mutation_count_per_kb_top50expr'
    
    print("generating random backgrounds", flush=True)
    # 1000 cells 500 times took 1 hour
    background_df = create_random_background(
        100, tms.adata[start_cell_num:end_cell_num, :],
        regress = regress, ssgsea = True
        )
    # add mutation col and dream col
    background_df[mut_col] = tms.dream_expression[start_cell_num:end_cell_num, :].obs[mut_col]
    background_df[dream_col] = tms.dream_expression[start_cell_num:end_cell_num, :].obs[dream_col]
    if regress:
        background_df.to_parquet(
            f'/cellar/users/zkoch/dream/data/tabula_muris_senis/random_background_new_DREAM/tms_random_background500iter_{start_cell_num}-{end_cell_num}cells.parquet'
            )
        print(f"wrote to /cellar/users/zkoch/dream/data/tabula_muris_senis/random_background_new_DREAM/tms_random_background500iter_{start_cell_num}-{end_cell_num}cells.parquet")
    else:
        background_df.to_parquet(
            f'/cellar/users/zkoch/dream/data/tabula_muris_senis/random_background_new_DREAM/tms_random_background500iter_{start_cell_num}-{end_cell_num}cells_noregress.parquet'
            )
        print(f"wrote to /cellar/users/zkoch/dream/data/tabula_muris_senis/random_background_new_DREAM/tms_random_background500iter_{start_cell_num}-{end_cell_num}cells_noregress.parquet")

if __name__ == "__main__":
    main()