# use dream_proj_env conda environment
import pandas as pd
import numpy as np
import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr, ttest_ind, mannwhitneyu,linregress,ttest_ind_from_stats
import glob
import gseapy
import scanpy as sc
from pybiomart import Dataset
import random
from pylr2 import regress2
import statsmodels.formula.api as smf
import colorcet as cc
import psynlig

source_path = "/cellar/users/zkoch/dream"
if source_path not in sys.path:
    sys.path.append(os.path.join(source_path, 'source'))
# read source files
import read_data
from expr_dataset import ExpressionDataset
from meta_expr_dataset import MetaExpressionDataset
import utils

loader = read_data.DatasetLoader("tabula_sapiens")
ts = loader.load_dataset()

# read in first command line arg
cell_num_start = int(sys.argv[1])
# subset sc_rosmap to only include the cells in the range cell_num_start to cell_num_start + 100000
ts.adata = ts.adata[cell_num_start: cell_num_start + 10000]
# do the ssgsea
ts.get_n_genes_and_counts()
ts.get_dream_gene_expression()
ts.dream_enrichment_ssgsea()
# save sc_rosmap to pickle 
import pickle
with open(f"/cellar/users/zkoch/dream/data/tabula_sapiens/ssgsea/sc_tabula_sapiens_preprocessed_ssgsea_cellnumstart{cell_num_start}.pkl", "wb") as f:
    pickle.dump(ts, f)
