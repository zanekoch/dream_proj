import pandas as pd
import os
import pandas as pd
from collections import defaultdict
import glob
import dask.dataframe as dd
import numpy as np
from pybiomart import Dataset
import gseapy
import scanpy as sc 
from typing import Union, List
import statsmodels.formula.api as smf

"""
To add a new expression dataset
1. Add to DatasetLoader.load_dataset
2. Write a load_{dataset_name} function
3. Edit DatasetLoader.read_soft if need be
4. Add the dataset to ExpressionDataset.__init__, making the metadata mergable with the expression df
"""

class MetaExpressionDataset:
    """ Class to represent multiple expression datasets """
    def __init__(
        self,
        dataset_names : List[str]
        ) -> None:
        self.dataset_names = dataset_names
        self.datasets = {}
        # load each dataset
        for dataset_name in self.dataset_names:
            loader = DatasetLoader(dataset_name)
            expression_dset = loader.load_dataset()           
            self.datasets[dataset_name] = expression_dset
        
    def scale_by_total_seq_depth(self):
        for dataset_name in self.dataset_names:
            self.datasets[dataset_name].scale_by_total_seq_depth()
            
    def log_scale_expr(self):
        for dataset_name in self.dataset_names:
            self.datasets[dataset_name].log_scale_expr()
            
    def get_dream_gene_expression(self):
        for dataset_name in self.dataset_names:
            self.datasets[dataset_name].get_dream_gene_expression()
            
    def dream_enrichment_ssgsea(self):
        for dataset_name in self.dataset_names:
            self.datasets[dataset_name].dream_enrichment_ssgsea()
            
    def run_GSEA(self):
        all_gs_results_dfs = []
        for dataset_name in self.dataset_names:
            if dataset_name == 'tyshkovskiy':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query(
                        "condition != 'Control'"
                        )['condition'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'age'
                    )
            elif dataset_name == 'boutant_nestle':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'tissue'
                    )
            elif dataset_name == 'martin_montalvo':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'tissue'
                    )
            elif dataset_name == 'zhou_2012':
                self.datasets[dataset_name].expression_df['is_high_fat'] = self.datasets[dataset_name].expression_df['condition'].str.contains('HF')
                self.datasets[dataset_name].expression_df['is_high_fat'] = self.datasets[dataset_name].expression_df['is_high_fat'].map({True:'HF', False:'LF'})
                self.datasets[dataset_name].metadata_df['is_high_fat'] = self.datasets[dataset_name].metadata_df['condition'].str.contains('HF')
                self.datasets[dataset_name].metadata_df['is_high_fat'] = self.datasets[dataset_name].metadata_df['is_high_fat'].map({True:'HF', False:'LF'})
                self.datasets[dataset_name].expression_df['condition2'] = self.datasets[dataset_name].expression_df['condition'].map({
                    'LF':'Control', 'LF+Exercise':'Exercise','LF+CR':'CR', 'HF+Exercise':'Exercise', 'HF':'Control', 'HF+CR':'CR'
                })
                self.datasets[dataset_name].metadata_df['condition2'] = self.datasets[dataset_name].metadata_df['condition'].map({
                    'LF':'Control', 'LF+Exercise':'Exercise','LF+CR':'CR', 'HF+Exercise':'Exercise', 'HF':'Control', 'HF+CR':'CR'
                })
                self.datasets[dataset_name].meta_cols.append('condition2')
                self.datasets[dataset_name].meta_cols.append('is_high_fat')

                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition2',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition2 != 'Control'")['condition2'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'is_high_fat'
                    )
            elif dataset_name == 'fok_chronic_2014' or dataset_name == 'fok_short_term_2014':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'sex'
                    )
            elif dataset_name == 'fok_cr_2014':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control'
                )
            else:
                raise NotImplementedError(f"run_GSEA not implemented for {dataset_name}")
            self.datasets[dataset_name].all_gs_results_df['dataset'] = dataset_name
            all_gs_results_dfs.append(self.datasets[dataset_name].all_gs_results_df)
        self.all_gs_results_df = pd.concat(all_gs_results_dfs)
        self.all_gs_results_df.reset_index(inplace=True, drop = True)
        
class ScExpressionDataset:
    """ Class to represent a single cell expression dataset """
    def __init__(
        self,
        # load anndata.AnnData object
        adata: sc.AnnData,
        gene_species: str, # species of gene names
        dataset: str,
        expr_counts_for_mutation: pd.DataFrame = None,
        mutation_counts: pd.DataFrame = None
        ) -> None:
        self.adata = adata
        self.gene_species = gene_species
        self.dataset = dataset
        self.expr_counts_for_mutation = expr_counts_for_mutation
        self.mutation_counts = mutation_counts
        self.already_converted = False
        if self.dataset == 'palovics_parabiosis':
            self.adata.obs['condition'] = self.adata.obs['condition'].map({
                'IY': 'Isochronic Young', 'HY': 'Heterochronic Young',
                'IA': 'Isochronic Aged', 'HA': 'Heterochronic Aged',
                'Y': 'Young-TMS', 'A': 'Aged-TMS'
                })
        """elif self.dataset == 'ma_sc_rat_CR':
            self.adata.obs['tissue'] = self.adata.obs.index.str.split('_')[1]
            self.adata.obs['condition'] = self.adata.obs.index.str.split('_')[-1].str.split('-')[0]"""
            
    def pre_process(self, log_scale: bool = True) -> None:
        """Preprocess the adata object
        ### Returns:
        None
        """
        if self.dataset == 'tms' or self.dataset == 'palovics_parabiosis':
            raise NotImplementedError("Preprocessing for TMS and palovics_parabiosis not recommended as they are saved as processed")
        # change all var names to upper case
        self.adata.var_names = self.adata.var_names.str.upper()
        # filter out cells with less than 200 genes or genes expressed in less than 3 cells
        sc.pp.filter_cells(self.adata, min_genes=200)
        sc.pp.filter_genes(self.adata, min_cells=3)
        # calculate QC metrics
        sc.pp.calculate_qc_metrics(
            self.adata, percent_top=None, log1p=False, inplace=True
            )
        # filter cells with more than 10,000,000 total counts
        self.adata.obs['outlier_total_counts'] = self.adata.obs.total_counts > 10000000
        print(f"filtering out {sum(self.adata.obs['outlier_total_counts'])} cells with more than 10,000,000 total counts out of {self.adata.shape[0]} cells")
        self.adata = self.adata[~self.adata.obs['outlier_total_counts'], :]
        if log_scale:
            # log scale the expression
            sc.pp.log1p(self.adata)
        # scale by total sequence depth
        sc.pp.scale(self.adata, max_value=10, zero_center=False)
        
    def read_dream_files(self) -> None:
        """Read in DREAM files
        ### Returns:
        None
        """
        self.dream_regulated_genes = pd.read_csv(
            "/cellar/users/zkoch/dream/data/bujarrabal_dueso/tableS12_dream_promoter_binding.csv", index_col=1
            )
        # replace spaces with underscores and make lowercase
        self.dream_regulated_genes.columns = [
            x.lower().replace(" ", "_") for x in self.dream_regulated_genes.columns
            ]
        # TODO: may need to add gene converter read in

    def get_dream_gene_expression(
        self,
        row_limiting_query: str = None,
        convert = True
        ) -> pd.DataFrame:
        """Get expression of DREAM genes
        ### Parameters:
        row_limiting_query : str
            Query to limit rows of dream_regulated_genes dataframe
        convert : bool
            Wether to convert the gene names of adata to match the dream gene names
        ### Returns:
        None
        """
        # check if dream files have been read in
        if not hasattr(self, "dream_regulated_genes"):
            self.read_dream_files()
        # get the DREAM regulated genes we want, index is gene names
        if row_limiting_query != None:
            dream_regulated_genes_names = self.dream_regulated_genes.query(row_limiting_query).index
        else:
            dream_regulated_genes_names = self.dream_regulated_genes.index
        # TODO: may need to add gene convertersion for non-TMS data eventually
        if convert and self.gene_species != "human" and not self.already_converted:
            # do gene conversion using 
            self._convert_genes_to_human_ensembl()
            
        self.dream_regulated_genes_w_expression = list(
            set(self.adata.var_names).intersection(
                set(dream_regulated_genes_names)
                )
            )
        print(f"Found {len(self.dream_regulated_genes_w_expression)} DREAM genes with expression")
        # create a subset of adata that is just the dream genes
        self.dream_expression = self.adata[
            :, self.dream_regulated_genes_w_expression
            ].copy()
        # calculate mean
        self.dream_expression.obs['mean_dream_activity'] = np.mean(
            self.dream_expression.X, axis=1
            )
        # get resid of mean dream activity
        self.scale_by_seq_depth(col_name = 'mean_dream_activity')
        self.adata.obs['mean_dream_activity_resid'] = self.dream_expression.obs['mean_dream_activity_resid']

    def scale_by_seq_depth(self, col_name: str):
        mut_ols = smf.ols(
            formula=f'{col_name} ~ total_counts * n_genes',
            data=self.dream_expression.obs
            ).fit()
        self.dream_expression.obs[f'{col_name}_resid'] = mut_ols.resid
        # invert residuals making highest values low and lowest values high
        self.dream_expression.obs[f'{col_name}_resid'] = -1 * self.dream_expression.obs[f'{col_name}_resid']
        # then add the min value to make all values positive
        self.dream_expression.obs[f'{col_name}_resid'] = self.dream_expression.obs[f'{col_name}_resid'] - min(self.dream_expression.obs[f'{col_name}_resid'])
        print(f"scaled {col_name} by sequence depth and created {col_name}_resid")
    
    def _convert_genes_to_human_ensembl(self) -> None:
        """Convert the gene names of adata to human ensembl gene names"""
        self.already_converted = True
        if self.gene_species == "mouse":
            species_name = 'mmusculus_gene_ensembl'
            dataset = Dataset(name='mmusculus_gene_ensembl', host='http://www.ensembl.org')
            raise NotImplementedError("_convert_genes_to_human_ensembl not implemented for mouse")
        elif self.gene_species == "rat":
            species_name = 'rnorvegicus_gene_ensembl'
            dataset = Dataset(name='rnorvegicus_gene_ensembl', host='http://www.ensembl.org')
            # get a mapping from ensembl rat gene names to human gene names
            ensembl_gene_mapper = dataset.query(
                attributes=['ensembl_gene_id', 'external_gene_name',
                            'hsapiens_homolog_orthology_type',
                            'hsapiens_homolog_orthology_confidence',
                            'hsapiens_homolog_associated_gene_name']
                )
            ensembl_gene_mapper.dropna(subset = 'Human homology type', inplace = True)
            ensembl_gene_mapper.drop_duplicates(
                subset = 'Gene stable ID', inplace = True
                )
            self.adata.var['human_gene_name'] = self.adata.var[
                'gene_ids'
                ].map(ensembl_gene_mapper.set_index('Gene stable ID')['Human gene name'])
            self.adata.var.dropna(subset = ['human_gene_name'], inplace = True)
            self.adata.var.drop_duplicates(subset = 'human_gene_name', inplace = True)
            self.adata = self.adata[:, self.adata.var.index]
            # remove genes dropped from var
            self.adata.var.set_index('human_gene_name', inplace = True, drop = False)
            
            self.adata = self.adata[:, self.adata.var.index]
        
            self.adata.var_names = self.adata.var.index
            
        else:
            raise NotImplementedError(f"Species {self.gene_species} not implemented")
            
    def map_illm_probe_to_ensembl(self):
        dataset = Dataset(name='mmusculus_gene_ensembl', host='http://www.ensembl.org')
        ensembl_gene_mapper = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name','illumina_mouseref_8'])
        self.expression_df.columns = self.expression_df.columns.map(
                ensembl_gene_mapper.set_index('ILLUMINA MouseRef 8 probe')['Gene stable ID'].to_dict()
                )
            # drop columns with nan
        self.expression_df = self.expression_df.loc[:, ~self.expression_df.columns.isna()]
            # combine columns with duplicate names by summing
        self.expression_df = self.expression_df.groupby(
                self.expression_df.columns, axis=1
                ).sum()
        
    def dream_enrichment_ssgsea(self) -> None:
        """Run ssgsea on the DREAM genes"""
        if self.dream_regulated_genes_w_expression is None:
            # run get_dream_gene_expression
            self.get_dream_gene_expression()
        # get dream gene to check for enrichment in expression
        dream_gene_set_dict = {
            'dream_reg_genes': self.dream_regulated_genes_w_expression
            }
        # run ssgsea
        ssgsea = gseapy.ssgsea(
            data=self.adata.to_df().T, # get a df version of the expr
            gene_sets=dream_gene_set_dict,
            outdir=None, no_plot=True
            )
        results_df = ssgsea.res2d
        results_df.set_index('Name', inplace=True)
        results_df.drop(columns = ['Term'], inplace=True)
        results_df.rename(
            columns = {'NES':'DREAM_normalized_enrichment_score', 'ES': 'DREAM_enrichment_score'},
            inplace=True
            )
        # convert cols to float
        results_df = results_df.astype(float)
        # merge to dream expression on index
        self.dream_expression.obs = self.dream_expression.obs.merge(
            results_df, left_index=True, right_index=True
            )

    def _get_gene_length(
        self, 
        gene_name:str,
        gencode_m34:pd.DataFrame,
        type: str) -> int:
        """Get the length of a gene from gencode_m34"""
        gene = gencode_m34.query("gene_name == @gene_name and type == @type")
        if gene.shape[0] > 1:
            #print(gene)
            # use just first row
            gene = gene.iloc[0]
        elif gene.shape[0] == 0:
            return -1
        try:
            length = gene['end'].values[0] - gene['start'].values[0]
        except:
            length = gene['end'] - gene['start']
        return length

    def _read_gencode_m34(self) -> pd.DataFrame:
        """Read in gencode m34 and process it"""
        gencode_m34 = pd.read_csv(
            "/cellar/users/zkoch/dream/utilities/gencode.vM34.basic.annotation.gtf.gz",
            sep='\t', skiprows=5, header=None
            )
        gencode_m34.columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
        # split attributes column
        gencode_m34['gene_id'] = gencode_m34['attributes'].str.split('gene_id "').str[1].str.split('"').str[0]
        gencode_m34['gene_name'] = gencode_m34['attributes'].str.split('gene_name "').str[1].str.split('"').str[0]
        gencode_m34['gene_name'] = gencode_m34['gene_name'].str.upper()
        return gencode_m34

    def calculate_mutation_burden(
        self, 
        col_name: str,
        top_perc_expr: float = 0.5,
        max_vaf: float = 0.6,
        ) -> None:
        """Calculate the mutation burden for each cell per kb. Only considers mutations to genes which were the top_perc_expr expressed across all ages
        ### Parameters:
        col_name : str
            Column name in adata.obs to store this mutation burden
        top_perc_expr : float
            The top percentage of genes to consider for mutation burden. Genes with a lower expression than the top_perc_expr in any age group are filtered.
        max_vaf : float
            The maximum variant allele frequency to consider a mutation, below is filtered
        ### Returns:
        None
        """
        # process self.expr_counts_for_mutation and self.mutation_counts
        # drop cells without a cell ontology class
        self.expr_counts_for_mutation = self.expr_counts_for_mutation[
            self.expr_counts_for_mutation.obs['cell_ontology_class']!='nan'
            ]
        self.mutation_counts = self.mutation_counts[
            self.mutation_counts.obs.index.isin(self.expr_counts_for_mutation.obs.index)
            ]
        # keep only genes with at least 1 mutation
        sc.pp.filter_genes(self.mutation_counts, min_cells=1)
        self.expr_counts_for_mutation = self.expr_counts_for_mutation[
            :,self.expr_counts_for_mutation.var_names.isin(self.mutation_counts.var_names)
            ]
        # filter genes that are expressed in fewer than 3 cells
        sc.pp.filter_genes(self.expr_counts_for_mutation, min_cells=3)
        self.mutation_counts = self.mutation_counts[
            :,self.mutation_counts.var_names.isin(self.expr_counts_for_mutation.var_names)
            ]
        # find mutations present in more than max_vaf of cells
        countsgenemutations = pd.DataFrame(
            index = self.mutation_counts.var_names,columns = ['3m','18m','21m','24m']
            )
        # get the count of mutations in each gene in each age
        aux = self.mutation_counts[self.mutation_counts.obs['age']=='3m']
        countsgenemutations['3m'] = (np.sum(aux.X>0,axis=0)/aux.shape[0]).transpose()
        aux = self.mutation_counts[self.mutation_counts.obs['age']=='18m']
        countsgenemutations['18m'] = (np.sum(aux.X>0,axis=0)/aux.shape[0]).transpose()
        aux = self.mutation_counts[self.mutation_counts.obs['age']=='21m']
        countsgenemutations['21m'] = (np.sum(aux.X>0,axis=0)/aux.shape[0]).transpose()
        aux = self.mutation_counts[self.mutation_counts.obs['age']=='24m']
        countsgenemutations['24m'] = (np.sum(aux.X>0,axis=0)/aux.shape[0]).transpose()
        # get the average mutation count across all ages
        countsgenemutations['avg'] = countsgenemutations.sum(axis=1)/4
        countsgenemutations.sort_values(by='avg',ascending=False)
        # filter out genes with a mutation in more than max_vaf of cells across all ages
        germline_mutations = countsgenemutations.loc[
            countsgenemutations['avg'] > max_vaf
            ].index
        genestoremove = list(set(germline_mutations))
        print("removing", len(genestoremove), "genes with germline mutations")
        sc.pp.filter_cells(self.expr_counts_for_mutation,min_genes=1)
        sc.pp.filter_cells(self.mutation_counts,min_genes=0)
        # remove these from both mutation counts and expression counts
        self.mutation_counts = self.mutation_counts[
            :,~self.mutation_counts.var_names.isin(genestoremove)
            ]
        self.expr_counts_for_mutation = self.expr_counts_for_mutation[
            :,~self.expr_counts_for_mutation.var_names.isin(genestoremove)
            ]

        """sc.pp.filter_cells(self.expr_counts_for_mutation,min_genes=0)
        sc.pp.filter_cells(self.mutation_counts,min_genes=0)"""
        # convert expression values to a dataframe
        counts = self.expr_counts_for_mutation.X.todense()
        counts = pd.DataFrame(
            counts, index = self.expr_counts_for_mutation.obs_names,
            columns=self.expr_counts_for_mutation.var_names
            )
        # iterate across each tissue and age to find the genes that are expressed in the top_perc_expr of cells
        genes_tissue_dict = {}
        for t in list(set(self.expr_counts_for_mutation.obs['tissue'])):
            # select this tissue
            adata_aux = self.expr_counts_for_mutation[
                self.expr_counts_for_mutation.obs['tissue']==t
                ]
            if t == 'Mammary_Gland':
                # select this age
                adata_aux3 = adata_aux[adata_aux.obs['age']=='3m']
                counts_aux3 = counts.loc[adata_aux3.obs_names]
                # count the number of cells with expression for each gene
                counts_aux3.loc['sum'] = counts_aux3.astype(bool).sum(axis=0)
                # filter out genes that are expressed in fewer than top_perc_expr of cells
                minimum_expression_count = counts_aux3.shape[0] * (1 - top_perc_expr)
                a3 = counts_aux3.loc['sum'] > minimum_expression_count
                a3 = list(a3[a3==True].index)
                # 18 months
                adata_aux18 = adata_aux[adata_aux.obs['age']=='18m']
                counts_aux18 = counts.loc[adata_aux18.obs_names]
                counts_aux18.loc['sum'] = counts_aux18.astype(bool).sum(axis=0)
                minimum_expression_count = counts_aux18.shape[0] * (1 - top_perc_expr)
                a18 = counts_aux18.loc['sum']>minimum_expression_count
                a18 = list(a18[a18==True].index)
                # 21 months
                adata_aux21 = adata_aux[adata_aux.obs['age']=='21m']
                counts_aux21 = counts.loc[adata_aux21.obs_names]
                counts_aux21.loc['sum'] = counts_aux21.astype(bool).sum(axis=0)
                minimum_expression_count = counts_aux21.shape[0] * (1 - top_perc_expr)
                a21 = counts_aux21.loc['sum'] > minimum_expression_count
                a21 = list(a21[a21==True].index)
                # get the intersection of genes expressed in the top_perc_expr of cells in all ages
                genes_tissue_dict[t] = list(set(a3) & set(a18) & set(a21))
            else:
                # select this age
                adata_aux3 = adata_aux[adata_aux.obs['age']=='3m']
                counts_aux3 = counts.loc[adata_aux3.obs_names]
                # count the number of cells with expression for each gene
                counts_aux3.loc['sum'] = counts_aux3.astype(bool).sum(axis=0)
                # filter out genes that are expressed in fewer than top_perc_expr of cells
                minimum_expression_count = counts_aux3.shape[0] * (1 - top_perc_expr)
                a3 = counts_aux3.loc['sum'] > minimum_expression_count
                a3 = list(a3[a3==True].index)
                # 18 months
                adata_aux18 = adata_aux[adata_aux.obs['age']=='18m']
                counts_aux18 = counts.loc[adata_aux18.obs_names]
                counts_aux18.loc['sum'] = counts_aux18.astype(bool).sum(axis=0)
                minimum_expression_count = counts_aux18.shape[0] * (1 - top_perc_expr)
                a18 = counts_aux18.loc['sum']>minimum_expression_count
                a18 = list(a18[a18==True].index)
                # 24 months
                adata_aux24 = adata_aux[adata_aux.obs['age']=='24m']
                counts_aux24 = counts.loc[adata_aux24.obs_names]
                counts_aux24.loc['sum'] = counts_aux24.astype(bool).sum(axis=0)
                minimum_expression_count = counts_aux24.shape[0] * (1 - top_perc_expr)
                a24 = counts_aux24.loc['sum'] > minimum_expression_count
                a24 = list(a24[a24==True].index)
                # get the intersection of genes expressed in the top_perc_expr of cells in all ages
                genes_tissue_dict[t] = list(set(a3) & set(a18) & set(a24))
        # get one list of all genes we are calling muts from
        self.gene_w_called_muts = []
        for t in genes_tissue_dict.keys():
            self.gene_w_called_muts.extend(genes_tissue_dict[t])
        self.gene_w_called_muts = list(set(self.gene_w_called_muts))
        
        # now get the gene lengths of the genes we are calling muts from
        gencode_m34 = self._read_gencode_m34()
        # calculate gene lengths
        gene_w_called_muts_df = pd.DataFrame(self.gene_w_called_muts, columns=['gene_name'])
        gene_w_called_muts_df['transcript_length'] = gene_w_called_muts_df.apply(
            lambda x: self._get_gene_length(x['gene_name'], gencode_m34, 'transcript'), axis=1
            )
        # for each tissue, get the total length of genes with called mutations and mut burden
        self.mutation_counts.obs[col_name] = -1
        total_mut_called_gene_length = {}
        for t in genes_tissue_dict.keys():
            # get gene legnths
            genes_called_from = genes_tissue_dict[t]
            genes_called_from_w_lengths = set(
                gene_w_called_muts_df['gene_name'].to_list()
                ).intersection(set(genes_called_from))
            if len(genes_called_from_w_lengths) != len(genes_called_from):
                print(t, len(genes_called_from), len(genes_called_from_w_lengths))
            total_mut_called_gene_length[t] = sum(
                gene_w_called_muts_df.query("gene_name in @genes_called_from_w_lengths")['transcript_length']
                )
            # get mutation count sum
            mut_count_this_tissue = self.mutation_counts[
                self.mutation_counts.obs['tissue'] == t, genes_called_from
                ].X.sum(axis=1).A1
            # scale by the total bps
            scaled_counts = mut_count_this_tissue / (total_mut_called_gene_length[t] / 1000)
            self.mutation_counts.obs.loc[
                self.mutation_counts.obs['tissue'] == t, col_name
                ] = scaled_counts
        # replace _ in key with space
        new_keys = [x.replace('_', ' ') for x in total_mut_called_gene_length.keys()]
        total_mut_called_gene_length = dict(zip(new_keys, total_mut_called_gene_length.values()))
        # add to adata.obs
        self.mutation_counts.obs['short_cell_name'] = self.mutation_counts.obs[
            'cell'
            ].str.split('_').str[:2].str.join('_')
        self.adata.obs = self.adata.obs.merge(
            self.mutation_counts.obs[['short_cell_name', col_name]],
            on='short_cell_name', how='left'
            )
              
class ExpressionDataset:
    """ Class to represent an expression dataset """
    def __init__(
        self, 
        expression_df: pd.DataFrame, 
        species: str,
        metadata_df: pd.DataFrame,
        dataset: str
        ) -> None:
        """Constructor for ExpressionDataset
        ### Parameters:
        expression_df : pd.DataFrame
            Expression dataframe, samples x genes
        species : str
            String identifyin the species used to define genes in expression_df
        metadata_df : pd.DataFrame
            Metadata dataframe
        dataset : str
            Dataset name
        ### Returns:
        None
        """
        self.expression_df = expression_df
        self.species = species
        self.metadata_df = metadata_df
        self.dataset = dataset
        # remove the trailing .N or .NN from genes if present
        self.expression_df.columns = self.expression_df.columns.str.replace(
            r'\.\d+$', '', regex=True
            )
        # do again incase multiple .N or .NN
        self.expression_df.columns = self.expression_df.columns.str.split(".").str[0]
        self.dream_regulated_genes_w_expression = None
        
        if self.dataset == 'tcga':
            # add ages, tissue, cancer type, tissue, sample_type, and case_id 
            self.meta_cols = ['age_at_index', 'primary_diagnosis', 'tissue_or_organ_of_origin', 'case_id', 'sample_type']
            self.expression_df = self.expression_df.merge(
                self.metadata_df[self.meta_cols],
                left_index=True, right_index=True
                )
            # add a column specifying the number of samples a case has
            self.get_num_samples_per_case()
            # simplify tissue and cancer names
            self.simplify_tissue_and_cancer_names()
        elif self.dataset == 'mSalt':
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_name'
                )
        elif self.dataset == 'tyshkovskiy':
            self.metadata_df.rename(columns = {'intervention':'condition'}, inplace=True)
            self.metadata_df.loc[
                self.metadata_df['condition'].str.contains('Control'), 'condition'
                ] = 'Control'
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            # remove sample_title from metacols
            self.meta_cols.remove('sample_title')
            self.expression_df.set_index('sample_title', inplace=True)
            # turn age into int
            self.expression_df['age'] = self.expression_df['age'].str.strip('months').astype(str)
        elif self.dataset == 'martin_montalvo':
            self.map_illm_probe_to_ensembl()
            # combine w metadata
            self.metadata_df['condition'] = self.metadata_df['sample_title'].apply(lambda x: x.split('-')[-2])
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'SD':'Control', 'MET':'Metformin'}
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
        elif self.dataset == 'boutant_nestle':
            # map the probe names to ensembl
            self.map_illm_probe_to_ensembl()
            # combine with metadata
            self.metadata_df.rename({'intervention':'condition'}, axis=1, inplace=True)
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'WT mice':'Control', 'EX mice':'Exercise', 'CR mice':'CR', 'SIRT1Tg/Tg mice':'SIRT1Tg/Tg'}
                )
            self.metadata_df['sample_title'] = self.metadata_df['sample_title'].str.split(' ').str[-1]
            # make tissue column
            self.metadata_df['tissue'] = self.metadata_df['sample_title'].str.split('_').str[0]
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            self.expression_df.set_index('sample_title', inplace=True)
            # remove sample_title from metacols
            self.meta_cols.remove('sample_title')
        elif self.dataset == 'zhou_2012':
            # convert age to int
            self.metadata_df['age'] = self.metadata_df['age'].str.strip('months')
            self.metadata_df['age'] = self.metadata_df['age'].str.strip('weeks')
            self.metadata_df['age'] = self.metadata_df['age'].astype(int)
            # last 4 rows are in months
            self.metadata_df.rename(columns = {'age':'age_months'}, inplace=True)
            self.metadata_df['condition'] = self.metadata_df['sample_title'].apply(
                lambda x: x.split('_')[0]
            )
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'LF+Ex':'LF+Exercise', 'HF+Ex':'HF+Exercise'}
            )
            # map the probe names to ensembl
            self.map_affy_probe_to_ensembl()
            # combine with metadata
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
        elif self.dataset == 'fok_chronic_2014':
            self.map_illm_probe_to_ensembl()
            self.metadata_df['condition'] = self.metadata_df['sample_title'].apply(lambda x: x.split(" ")[0])
            self.metadata_df['condition'] = self.metadata_df['condition'].str[:4]
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'Cont':'Control', 'Rapa':'Rapamycin'}
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
        elif self.dataset == 'fok_short_term_2014':
            self.map_illm_probe_to_ensembl()
            self.metadata_df['condition'] = self.metadata_df['sample_title'].apply(lambda x: x.split(" ")[0])
            self.metadata_df['condition'] = self.metadata_df['condition'].str[:4]
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'Cont':'Control', 'Rapa':'Rapamycin'}
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
        elif self.dataset == 'fok_cr_2014':
            self.map_illm_probe_to_ensembl()
            # keep everything else stuff after last space
            self.metadata_df['condition'] = self.metadata_df['sample_title'].apply(
                lambda x: ' '.join(x.split(" ")[:-1])
            )
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'Ad Libitum':'Control', 'Dietary Restriction':'CR', 'Rapamycin + Dietary Restriction':'Rapamycin+CR'}
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
        else:
            raise NotImplementedError(f"Dataset {self.dataset} not implemented")

    def map_illm_probe_to_ensembl(self):
        """Map illumina probe names to ensembl gene names"""
        dataset = Dataset(name='mmusculus_gene_ensembl', host='http://www.ensembl.org')
        ensembl_gene_mapper = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name','illumina_mouseref_8'])
        self.expression_df.columns = self.expression_df.columns.map(
                ensembl_gene_mapper.set_index('ILLUMINA MouseRef 8 probe')['Gene stable ID'].to_dict()
                )
            # drop columns with nan
        self.expression_df = self.expression_df.loc[:, ~self.expression_df.columns.isna()]
            # combine columns with duplicate names by summing
        self.expression_df = self.expression_df.groupby(
                self.expression_df.columns, axis=1
                ).sum()
        
    def map_affy_probe_to_ensembl(self):
        """Map affy probe names to ensembl gene names"""
        dataset = Dataset(name='mmusculus_gene_ensembl', host='http://www.ensembl.org')
        ensembl_gene_mapper = dataset.query(
            attributes=['ensembl_gene_id', 'external_gene_name','affy_mouse430_2']
            )
        self.expression_df.columns = self.expression_df.columns.map(
                ensembl_gene_mapper.set_index('AFFY Mouse430 2 probe')['Gene stable ID'].to_dict()
                )
            # drop columns with nan
        self.expression_df = self.expression_df.loc[:, ~self.expression_df.columns.isna()]
            # combine columns with duplicate names by summing
        self.expression_df = self.expression_df.groupby(
                self.expression_df.columns, axis=1
                ).sum()
                
    def get_num_samples_per_case(self) -> None:
        """Add a column to self.expression specifying if that row's sample_id is unique to a case_id
        ### Returns:
        None
        """
        # get the count of sample_ids per case_id
        sample_id_count = self.metadata_df.reset_index().groupby(
            "case_id"
            )['sample_id'].nunique()
        # add a column to self.expression_df specifying if that row's sample_id is unique to a case_id
        self.expression_df['num_samples_this_case'] = self.expression_df[
            'case_id'
            ].map(sample_id_count.to_dict())
        # also add to metadata_df
        self.metadata_df['num_samples_this_case'] = self.metadata_df[
            'case_id'
            ].map(sample_id_count.to_dict())
        self.meta_cols.append('num_samples_this_case')
    
    def read_dream_files(self) -> None:
        """Read in DREAM files
        If the expression df is not human genes, convert the dream regulation genes to apprioriate species (self.species)
        ### Returns:
        None
        """
        self.dream_regulated_genes = pd.read_csv(
            "/cellar/users/zkoch/dream/data/bujarrabal_dueso/tableS12_dream_promoter_binding.csv", index_col=0
            )
        # replace spaces with underscores and make lowercase
        self.dream_regulated_genes.columns = [
            x.lower().replace(" ", "_") for x in self.dream_regulated_genes.columns
            ]
        # dream genes are human genes
        if self.species == "human":
            pass 
        # can convert human genes to mouse
        elif self.species == "mouse":
            self.gene_converter = pd.read_csv(
                "/cellar/users/zkoch/dream/utilities/human_mouse_ensembl_genes.txt.gz",
                sep="\t", index_col=0, 
                )
        else :
            raise NotImplementedError(f"Species {self.species} not implemented")
    
    def get_dream_gene_expression(
        self,
        row_limiting_query: str = None
        ) -> pd.DataFrame:
        """Get expression of DREAM genes
        ### Parameters:
        row_limiting_query : str
            Query to limit rows of dream_regulated_genes dataframe
        ### Returns:
        None
        """
        # check if dream files have been read in
        if not hasattr(self, "dream_regulated_genes"):
            self.read_dream_files()
        # get the DREAM regulated genes we wnat
        if row_limiting_query != None:
            dream_regulated_genes_names = self.dream_regulated_genes.query(row_limiting_query).index
        else:
            dream_regulated_genes_names = self.dream_regulated_genes.index
        # convert them if need be
        if self.species != "human":
            need_to_convert_genes = True
        else:
            need_to_convert_genes = False
        # if we need to convert genes
        if need_to_convert_genes:
            if self.species == "mouse":
                # convert to mouse genes
                dream_regulated_genes_names_converted = self.gene_converter.loc[
                    dream_regulated_genes_names, "Mouse gene stable ID"
                    ]
                dream_regulated_genes_names_converted.dropna(inplace=True)
                # get expression of converted genes
                self.dream_regulated_genes_w_expression = list(
                    set(dream_regulated_genes_names_converted).intersection(
                        set(self.expression_df.columns)
                        )
                    )
                print("Converted DREAM genes to mouse genes")
        else:
            self.dream_regulated_genes_w_expression = list(
                set(dream_regulated_genes_names).intersection(
                    set(self.expression_df.columns)
                    )
                )
            print("Did not need to convert DREAM genes")
        print(f"Found {len(self.dream_regulated_genes_w_expression)} DREAM genes with expression")
        # select the dream regulated genes from the expression df
        dream_expression = self.expression_df[
            self.dream_regulated_genes_w_expression + self.meta_cols
            ].copy(deep = True)
        self.dream_expression = dream_expression
        self.dream_expression['mean_dream_reg_expr'] = self.dream_expression[
            self.dream_regulated_genes_w_expression
            ].mean(axis=1)
        self.scale_dream_by_seq_depth(col_name = 'mean_dream_reg_expr')

    def scale_dream_by_seq_depth(self, col_name: str):
        mut_ols = smf.ols(
            formula=f'{col_name} ~ total_seq_depth * n_genes_expressed',
            data=self.dream_expression
            ).fit()
        self.dream_expression[f'{col_name}_resid'] = mut_ols.resid
        # invert residuals making highest values low and lowest values high
        self.dream_expression[f'{col_name}_resid'] = -1 * self.dream_expression[f'{col_name}_resid']
        # then add the min value to make all values positive
        self.dream_expression[f'{col_name}_resid'] = self.dream_expression[f'{col_name}_resid'] - min(self.dream_expression[f'{col_name}_resid'])
        print(f"scaled {col_name} by sequence depth and created {col_name}_resid")
        
    def dream_enrichment_ssgsea(self) -> None:
        """Run ssgsea on the DREAM genes"""
        if self.dream_regulated_genes_w_expression is None:
            # run get_dream_gene_expression
            self.get_dream_gene_expression()
        # run ssgsea
        dream_gene_set_dict = {'dream_reg_genes':self.dream_regulated_genes_w_expression}
        # select non-meta columns
        expr_df = self.expression_df.loc[:, ~self.expression_df.columns.isin(self.meta_cols)]
        ssgsea = gseapy.ssgsea(
            data=expr_df.T, gene_sets=dream_gene_set_dict,
            outdir=None, no_plot=True
            )
        results_df = ssgsea.res2d
        results_df.set_index('Name', inplace=True)
        results_df.drop(columns = ['Term'], inplace=True)
        results_df.rename(
            columns = {'NES':'DREAM_normalized_enrichment_score', 'ES': 'DREAM_enrichment_score'},
            inplace=True
            )
        # convert cols to float
        results_df = results_df.astype(float)
        # merge to dream expression on index
        self.dream_expression = self.dream_expression.merge(
            results_df, left_index=True, right_index=True
            )
        self.scale_dream_by_seq_depth(col_name = 'DREAM_normalized_enrichment_score')
        self.scale_dream_by_seq_depth(col_name = 'DREAM_enrichment_score')
        
    def test_differential_dream_expression(
        self,
        class_col: str,
        treatment_classes: list,
        control_class: str,
        secondary_grouping_col: str = ''
        ) -> None:
        """
        Run GSEA comparing treatment_classes to control_class. Optionally, group by secondary_grouping like tissue or sex
        ### Parameters:
        class_col : str
            Column name in metadata_df to get classes from
        treatment_classes : list
            List of classes to compare to control_class
        control_class : str
            Class to compare treatment_classes to
        secondary_grouping_col : str
            Column name in metadata_df to group by
        ### Returns:
            None
        """
        dream_gene_set_dict = {'dream_reg_genes':self.dream_regulated_genes_w_expression}
        gs_results_dfs = []
        if secondary_grouping_col == '':
            for treatment_class in treatment_classes:
                # select this class's rows
                expr_df = self.expression_df.query(
                    f"{class_col} == '{treatment_class}' or {class_col} == '{control_class}'"
                    )
                # get class col
                classes = expr_df[class_col].to_list()
                # select non-meta columns
                expr_df = expr_df.loc[:, ~expr_df.columns.isin(self.meta_cols)]
                expr_df = expr_df.T
                # run
                gs_res = gseapy.GSEA(
                    data=expr_df, gene_sets=dream_gene_set_dict,
                    classes=classes, permutation_type = 'gene_set', method='s2n'
                    )
                gs_res.pheno_pos = treatment_class
                gs_res.pheno_neg = control_class
                gs_res.run()
                gs_results_df = gs_res.res2d
                gs_results_df['treatment_class'] = treatment_class
                gs_results_df['control_class'] = control_class
                gs_results_df['secondary_grouping'] = ''
                gs_results_dfs.append(gs_results_df)
                print("finished", treatment_class, "vs. ", control_class)
            # check if class has attribute all_gs_results_df
            if not hasattr(self, 'all_gs_results_df'):
                self.all_gs_results_df = pd.concat(gs_results_dfs)
                # make treatment_class and control_class the first columns
                self.all_gs_results_df = self.all_gs_results_df[
                    ['treatment_class', 'control_class', 'secondary_grouping'] + self.all_gs_results_df.columns.to_list()[:-3]
                ]
            else:
                all_gs_results_df = pd.concat(gs_results_dfs)
                all_gs_results_df = all_gs_results_df[
                    ['treatment_class', 'control_class', 'secondary_grouping'] + all_gs_results_df.columns.to_list()[:-3]
                ]
                # add to all_gs_results_df
                self.all_gs_results_df = pd.concat(
                    [self.all_gs_results_df, all_gs_results_df]
                    )
        else:
            for secondary_grouping_val in self.expression_df[secondary_grouping_col].unique():
                for treatment_class in treatment_classes:
                    # select this class's rows
                    expr_df = self.expression_df.query(
                        f"{secondary_grouping_col} == '{secondary_grouping_val}' and ({class_col} == '{treatment_class}' or {class_col} == '{control_class}')"
                        )
                    # get class col
                    classes = expr_df[class_col].to_list()
                    if len(set(classes)) != 2:
                        print("missing 2 classes")
                        print(secondary_grouping_val)
                        print(treatment_class)
                        print(control_class)
                        print(expr_df)
                        print(classes)
                        continue
                    # select non-meta columns
                    expr_df = expr_df.loc[:, ~expr_df.columns.isin(self.meta_cols)]
                    expr_df = expr_df.T
                    # run
                    try:
                        gs_res = gseapy.GSEA(
                            data=expr_df, gene_sets=dream_gene_set_dict,
                            classes=classes, permutation_type = 'gene_set', method='s2n'
                            )
                        gs_res.pheno_pos = treatment_class
                        gs_res.pheno_neg = control_class
                        gs_res.run()
                    except:
                        print(secondary_grouping_val)
                        print(treatment_class)
                        print(control_class)
                        print(expr_df)
                        print(classes)
                        raise ValueError
                    
                    gs_results_df = gs_res.res2d
                    gs_results_df['treatment_class'] = treatment_class
                    gs_results_df['control_class'] = control_class
                    gs_results_df['secondary_grouping'] = secondary_grouping_val
                    gs_results_dfs.append(gs_results_df)
                    print("finished", treatment_class, "vs. ", control_class, "in", secondary_grouping_val)
            # check if class has attribute all_gs_results_df
            if not hasattr(self, 'all_gs_results_df'):
                self.all_gs_results_df = pd.concat(gs_results_dfs)
                # make treatment_class and control_class the first columns
                self.all_gs_results_df = self.all_gs_results_df[
                    ['treatment_class', 'control_class', 'secondary_grouping'] + self.all_gs_results_df.columns.to_list()[:-3]
                ]
            else:
                all_gs_results_df = pd.concat(gs_results_dfs)
                all_gs_results_df = all_gs_results_df[
                    ['treatment_class', 'control_class', 'secondary_grouping'] + all_gs_results_df.columns.to_list()[:-3]
                ]
                # add to all_gs_results_df
                self.all_gs_results_df = pd.concat(
                    [self.all_gs_results_df, all_gs_results_df]
                    )
    
    def simplify_tissue_and_cancer_names(self) -> None:
        """Simplify tissue names
        ### Returns:
        None
        """
        # map to simple tissue names
        simple_tissue_map = {
            # kidney
            'Kidney, NOS': 'Kidney',
            # lung
            'Upper lobe, lung': 'Lung',
            'Lower lobe, lung' : 'Lung',
            'Lung, NOS' : 'Lung',
            'Middle lobe, lung' : 'Lung',
            # brain
            'Brain, NOS' : 'Brain',
            'Parietal lobe' : 'Brain',
            'Temporal lobe' : 'Brain',
            'Frontal lobe' : 'Brain',
            'Occipital lobe' : 'Brain',
            # uterus
            'Corpus uteri' : 'Uterus',
            'Endometrium' : 'Uterus',
            # pancreas 
            'Pancreas, NOS' : 'Pancreas',
            'Head of pancreas' : 'Pancreas',
            'Tail of pancreas' : 'Pancreas',
            'Body of pancreas' : 'Pancreas',
            # head and neck
            'Larynx, NOS': 'Mouth',
            'Tongue, NOS' : 'Mouth',
            'Gum, NOS' : 'Mouth',
            'Overlapping lesion of lip, oral cavity and pharynx' : 'Mouth',
            'Floor of mouth, NOS' : 'Mouth',
            'Cheek mucosa' : 'Mouth',
            'Tonsil, NOS' : 'Mouth',
            'Head, face or neck, NOS' : 'Mouth',
            'Lip, NOS' : 'Mouth',
            'Oropharynx, NOS' : 'Mouth',
            'Base of tongue, NOS' : 'Mouth'
            }    
        self.expression_df['tissue'] = self.expression_df['tissue_or_organ_of_origin'].map(simple_tissue_map)
        # drop uknown, none, or nan
        self.expression_df = self.expression_df.query(
            "tissue != 'Unknown' and tissue != 'None'"
            )
        self.expression_df.dropna(subset=['tissue'], inplace=True)
        simple_cancer_map = {
            'Renal cell carcinoma, NOS' : 'Renal cell carcinoma',
            'Endometrioid adenocarcinoma, NOS' : 'Endometrioid adenocarcinoma', # uteran
            'Adenocarcinoma, NOS' : 'Adenocarcinoma',
            'Glioblastoma' : 'Glioblastoma',
            'Squamous cell carcinoma, NOS' : 'Squamous cell carcinoma', # skin
            'Infiltrating duct carcinoma, NOS' : 'Infiltrating duct carcinoma', # breast
            'Oligodendroglioma, anaplastic' : 'Oligodendroglioma',
            'Oligodendroglioma, NOS' : 'Oligodendroglioma',
        }
        self.expression_df['cancer_type'] = self.expression_df['primary_diagnosis'].map(simple_cancer_map)
        # drop uknown, none, or nan
        self.expression_df = self.expression_df.query(
            "cancer_type != 'Unknown' and cancer_type != 'None'"
            )
        self.expression_df.dropna(subset=['cancer_type'], inplace=True)
        self.meta_cols.extend(['tissue', 'cancer_type'])
    
    def scale_by_total_seq_depth(self) -> None:
        """Scale the expression df by the total sequence depth of each sample
        ### Returns:
        None
        """
        all_columns = self.expression_df.columns
        # remove meta cols
        non_meta_cols = [x for x in all_columns if x not in self.meta_cols]
        # get sum of all expression columns
        self.expression_df['total_seq_depth'] = self.expression_df.loc[:, non_meta_cols].sum(axis=1)
        # claculate the number of genes with nonzeor expression
        self.expression_df['n_genes_expressed'] = self.expression_df[non_meta_cols].gt(0).sum(axis=1)
        self.meta_cols = self.meta_cols + ['total_seq_depth', 'n_genes_expressed']
        
        # scale the numerical columns by the total sequence depth
        scaled_expression = self.expression_df[non_meta_cols].div(
            self.expression_df['total_seq_depth'], axis=0
            )
        # set the scaled expression columns to the scaled values
        self.expression_df = pd.concat(
            [scaled_expression, self.expression_df[self.meta_cols]], axis=1
            )
            
    def log_scale_expr(self) -> None:
        """Log scale the expression df
        ### Returns:
        None
        """
        all_columns = self.expression_df.columns
        # remove meta cols
        non_meta_cols = [x for x in all_columns if x not in self.meta_cols]
        # log scale the numerical columns
        scaled_expression_cols = np.log2(
            self.expression_df[non_meta_cols] + 1
            )
        # min-max scale to be in range 0 to 10
        scaled_expression_cols = scaled_expression_cols / scaled_expression_cols.max()
        scaled_expression_cols = scaled_expression_cols * 10
        # reconstruct the df
        self.expression_df = pd.concat(
            [scaled_expression_cols, self.expression_df[self.meta_cols]], axis=1
        )
    def batch_correct(
        self, 
        batch_col : str
        ) -> None:
        """Batch correct the expression df"""
        # get the expression columns
        all_columns = self.expression_df.columns
        # remove meta cols
        non_meta_cols = [x for x in all_columns if x not in self.meta_cols]
        # get the expression columns
        expression = self.expression_df[non_meta_cols]
        # batch correct
        expression = sc.pp.combat(
            expression.T, batch=self.expression_df[batch_col]
            ).T

class MutationDataset:
    """ Class to represent a mutation dataset """
    def __init__(
        self, 
        mutation_df: pd.DataFrame, 
        metadata_df: pd.DataFrame
        ) -> None:
        """Constructor for MutationDataset
        ### Parameters:
        mutation_df : pd.DataFrame
            Mutation dataframe, samples x genes
        metadata_df : pd.DataFrame
            Metadata dataframe
        ### Returns:
        None
        """
        self.mutation_df = mutation_df
        self.metadata_df = metadata_df
        # keep columns we want and create self.mutation_burden df
        self.process_mutations()
        # add a column specifying the number of samples a case has
        self.get_num_samples_per_case()    
        # simplify tissue and cancer names
        self.simplify_tissue_and_cancer_names()
    
    def process_mutations(self) -> None:
        """Select columns we want and calculate the mutation burden for each sample
        ### Returns:
        None
        """
        columns_we_want = [
            'first_sample_id', 'first_case_id', 'Chromosome', 
            'Start_Position', 'End_Position', 'Variant_Type',
            'Reference_Allele',
            'Tumor_Seq_Allele1', # genotype of tumor on first allele
            'Tumor_Seq_Allele2', # genotype of tumor on second allele
            'dbSNP_RS', # wether the mutation is in dbSNP
            'HGVSc', 'HGVSp_Short', # coding/protein variant in HGVS format
            # HGVSc mutation is 99.99% the same as Reference_Allele>Allele or its complement
            't_depth', # depth in tumor bam
            't_ref_count', # depth supporting reference
            't_alt_count', # depth supporting variant allele
            'Allele', # variant allele (same as tumor_seq_allele2, why?)
            'Gene', # gene name, Ensembl ID
            'One_Consequence', # consequence of variant (e.g., missense_variant)
            'Protein_position', # position of variant in protein/total length of protein
            ]
        self.mutation_df = self.mutation_df[columns_we_want].copy(deep=True)
        # rename columns
        self.mutation_df.rename(
            {'first_sample_id': 'sample_id', 'first_case_id': 'case_id'}, axis=1, inplace=True
            )
        # add columns of interest
        self.mutation_df['maf'] = self.mutation_df['t_alt_count'] / self.mutation_df['t_depth']
        # add ref>var column
        self.mutation_df['ref>var'] = self.mutation_df['Reference_Allele'] + '>' + self.mutation_df['Allele']
        # add a chr:start
        self.mutation_df['chr:start'] = self.mutation_df['Chromosome'].astype(str) + ':' + self.mutation_df['Start_Position'].astype(str)
        
        # create mutation burden df
        self.mutation_burden = self.mutation_df.groupby('sample_id').size().to_frame('mutation_burden')
        # add ages, tissue, cancer type, and case_id (for uniqueness)
        self.mutation_burden = self.mutation_burden.merge(
            self.metadata_df[['age_at_index', 'primary_diagnosis', 'tissue_or_organ_of_origin', 'case_id']], left_index=True, right_index=True
            )
        
    def get_num_samples_per_case(self) -> None:
        """Add a column to self.mutation_df specifying if that row's sample_id is unique to a case_id
        ### Returns:
        None
        """
        # get the count of sample_ids per case_id
        sample_id_count = self.mutation_df.groupby(
            'case_id'
            )['sample_id'].nunique()
        # add a column to self.mutation_df specifying if that row's sample_id is unique to a case_id
        self.mutation_df['num_samples_this_case'] = self.mutation_df[
            'case_id'
            ].map(sample_id_count.to_dict())
        # and to mutation_burden
        self.mutation_burden['num_samples_this_case'] = self.mutation_burden[
            'case_id'
            ].map(sample_id_count.to_dict())
        # also add to metadata_df
        self.metadata_df['num_samples_this_case'] = self.metadata_df[
            'case_id'
            ].map(sample_id_count.to_dict())

    def simplify_tissue_and_cancer_names(self) -> None:
        """Simplify tissue names
        ### Returns:
        None
        """
        # map to simple tissue names
        simple_tissue_map = {
            # kidney
            'Kidney, NOS': 'Kidney',
            # lung
            'Upper lobe, lung': 'Lung',
            'Lower lobe, lung' : 'Lung',
            'Lung, NOS' : 'Lung',
            'Middle lobe, lung' : 'Lung',
            # brain
            'Brain, NOS' : 'Brain',
            'Parietal lobe' : 'Brain',
            'Temporal lobe' : 'Brain',
            'Frontal lobe' : 'Brain',
            'Occipital lobe' : 'Brain',
            # uterus
            'Corpus uteri' : 'Uterus',
            'Endometrium' : 'Uterus',
            # pancreas 
            'Pancreas, NOS' : 'Pancreas',
            'Head of pancreas' : 'Pancreas',
            'Tail of pancreas' : 'Pancreas',
            'Body of pancreas' : 'Pancreas',
            # head and neck
            'Larynx, NOS': 'Mouth',
            'Tongue, NOS' : 'Mouth',
            'Gum, NOS' : 'Mouth',
            'Overlapping lesion of lip, oral cavity and pharynx' : 'Mouth',
            'Floor of mouth, NOS' : 'Mouth',
            'Cheek mucosa' : 'Mouth',
            'Tonsil, NOS' : 'Mouth',
            'Head, face or neck, NOS' : 'Mouth',
            'Lip, NOS' : 'Mouth',
            'Oropharynx, NOS' : 'Mouth',
            'Base of tongue, NOS' : 'Mouth'
            }    
        self.mutation_burden['tissue'] = self.mutation_burden['tissue_or_organ_of_origin'].map(simple_tissue_map)
        # drop uknown, none, or nan
        self.mutation_burden = self.mutation_burden.query(
            "tissue != 'Unknown' and tissue != 'None'"
            )
        self.mutation_burden.dropna(subset=['tissue'], inplace=True)
        simple_cancer_map = {
            'Renal cell carcinoma, NOS' : 'Renal cell carcinoma',
            'Endometrioid adenocarcinoma, NOS' : 'Endometrioid adenocarcinoma', # uteran
            'Adenocarcinoma, NOS' : 'Adenocarcinoma',
            'Glioblastoma' : 'Glioblastoma',
            'Squamous cell carcinoma, NOS' : 'Squamous cell carcinoma', # skin
            'Infiltrating duct carcinoma, NOS' : 'Infiltrating duct carcinoma', # breast
            'Oligodendroglioma, anaplastic' : 'Oligodendroglioma',
            'Oligodendroglioma, NOS' : 'Oligodendroglioma',
        }
        self.mutation_burden['cancer_type'] = self.mutation_burden['primary_diagnosis'].map(simple_cancer_map)
        # drop uknown, none, or nan
        self.mutation_burden = self.mutation_burden.query(
            "cancer_type != 'Unknown' and cancer_type != 'None'"
            )
        self.mutation_burden.dropna(subset=['cancer_type'], inplace=True)
        
class MethylationDataset:
    """ Class to represent a methylation dataset """
    def __init__(
        self, 
        methylation_df: pd.DataFrame, 
        metadata_df: pd.DataFrame
        ) -> None:
        """Constructor for MethylationDataset
        ### Parameters:
        methylation_df : pd.DataFrame
            Methylation dataframe, samples x genes
        metadata_df : pd.DataFrame
            Metadata dataframe
        ### Returns:
        None
        """
        self.methylation_df = methylation_df
        self.metadata_df = metadata_df

class DatasetLoader:
    """Class to load a dataset"""
    def __init__(
        self,
        dataset_name:str,
        load_expression: bool = True,
        load_mutation: bool = True,
        load_methylation: bool = True
        ) -> None:
        """Constructor for DatasetLoader:
        Optionally choose which datasets (in case you don't want to wait for methylation for example)
        ### Parameters
        dataset_name : str
            Name of the dataset to load
        load_expression : bool
            Whether to load expression data (if available)
        load_mutation : bool
            Whether to load mutation data (if available)
        load_methylation : bool
            Whether to load methylation data (if available)
        ### Returns
        None
        """
        self.dataset_name = dataset_name
        self.load_expression = load_expression
        self.load_mutation = load_mutation
        self.load_methylation = load_methylation

    def load_dataset(self) -> Union[
        List[Union[ScExpressionDataset, ExpressionDataset, MutationDataset, MethylationDataset]],
        ScExpressionDataset
        ]:
        """Load the dataset based on the dataset_name
        ### Returns
        A list of dataset objects or a single dataset object : list | ScExpressionDataset | ExpressionDataset | MutationDataset | MethylationDataset
            The contents of this list depends on the dataset. mSalt it is a list of two ExpressionDataset objects, across_species and treated_mice. For CPTAC-3 it is an [ExpressionDataset, MutationDataset, MethylationDataset]
        """
        print(f"Loading dataset: {self.dataset_name}")
        if self.dataset_name == "mSalt":
            dataset_list = self.load_mSalt()
        elif self.dataset_name == 'CPTAC-3':
            dataset_list = self.load_tcga()
        elif self.dataset_name == 'tyshkovskiy':
            dataset_list = self.load_tyshkovskiy()
        elif self.dataset_name == 'boutant_nestle':
            dataset_list = self.load_boutant_nestle()
        elif self.dataset_name == 'martin_montalvo':
            dataset_list = self.load_martin_montalvo()
        elif self.dataset_name == 'TMS':
            dataset_list = self.load_tms()
        elif self.dataset_name == 'palovics_parabiosis':
            dataset_list = self.load_palovics_parabiosis()
        elif self.dataset_name == 'ma_sc_rat_CR':
            dataset_list = self.load_ma_sc_rat_CR()
        elif self.dataset_name == 'zhou_2012':
            dataset_list = self.load_zhou_2012()
        elif self.dataset_name == 'fok_chronic_2014':
            dataset_list = self.load_fok_chronic_2014()
        elif self.dataset_name == 'fok_short_term_2014':
            dataset_list = self.load_fok_short_term_2014()
        elif self.dataset_name == 'fok_cr_2014':
            dataset_list = self.load_fok_cr_2014()
        else:
            raise NotImplementedError(f"Dataset {self.dataset_name} not implemented")
        return dataset_list
    
    def load_fok_cr_2014(self) -> ExpressionDataset:
        """Load the Fok CR 2014 dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        # paths
        dataset_path = "/cellar/users/zkoch/dream/data/fok_cr_2014"
        treated_mice_fn = "GSE40977_expression.csv"
        treated_mice_metadata_fn = "GSE40977_family.soft"

        if not self.load_expression:
            raise ValueError("load_expression must be True to load Fok cr 2014")
        # load data
        treated_mice_metadata = self.read_soft(
            os.path.join(dataset_path, treated_mice_metadata_fn)
            )
        treated_mice_expr = self.read_soft_expression(
            os.path.join(dataset_path, treated_mice_fn), metadata = treated_mice_metadata,
            array_type='illumina'
            )
        # create expression dataset object
        treated_mice = ExpressionDataset(
            expression_df=treated_mice_expr,
            species="mouse",
            metadata_df=treated_mice_metadata,
            dataset="fok_cr_2014"
            )
        return treated_mice
    
    def load_fok_chronic_2014(self) -> ExpressionDataset:
        """Load the Fok 2014 dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        # paths
        dataset_path = "/cellar/users/zkoch/dream/data/fok_2014"
        treated_mice_fn2 = "GSE48333_expression.csv"
        treated_mice_metadata_fn2 = "GSE48333_family.soft"

        if not self.load_expression:
            raise ValueError("load_expression must be True to load Fok 2014")
        
        treated_mice_metadata = self.read_soft(
            os.path.join(dataset_path, treated_mice_metadata_fn2)
            )
        treated_mice_expr = self.read_soft_expression(
            os.path.join(dataset_path, treated_mice_fn2), metadata = treated_mice_metadata,
            array_type='illumina'
            )
        # create expression dataset object
        treated_mice = ExpressionDataset(
            expression_df=treated_mice_expr,
            species="mouse",
            metadata_df=treated_mice_metadata,
            dataset="fok_chronic_2014"
            )
        return treated_mice
    
    def load_fok_short_term_2014(self) -> ExpressionDataset:
        """Load the Fok 2014 dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        # paths
        dataset_path = "/cellar/users/zkoch/dream/data/fok_2014"
        treated_mice_fn1 = "GSE48331_expression.csv"
        treated_mice_metadata_fn1 = "GSE48331_family.soft"

        if not self.load_expression:
            raise ValueError("load_expression must be True to load Fok 2014")
        # load data
        treated_mice_metadata = self.read_soft(
            os.path.join(dataset_path, treated_mice_metadata_fn1)
            )
        treated_mice_expr = self.read_soft_expression(
            os.path.join(dataset_path, treated_mice_fn1), metadata = treated_mice_metadata,
            array_type='illumina'
            )
        # create expression dataset object
        treated_mice = ExpressionDataset(
            expression_df=treated_mice_expr,
            species="mouse",
            metadata_df=treated_mice_metadata,
            dataset="fok_short_term_2014"
            )
        return treated_mice
    
    def load_zhou_2012(self) -> ExpressionDataset:
        """Load the Zhou 2012 dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        dataset_path = "/cellar/users/zkoch/dream/data/zhou_2012"
        
        treated_mice_fn1 = "GSE36838_expression.csv" # 18 mice
        treated_mice_metadata_fn1 = "GSE36838_family.soft"
        
        treated_mice_fn2 = "GSE36836_expression.csv" # 4 pooled bois
        treated_mice_metadata_fn2 = "GSE36836_family.soft"
        
        if not self.load_expression:
            raise ValueError("load_expression must be True to load Zhou 2012")
        # load data
        treated_mice_metadata1 = self.read_soft(
            os.path.join(dataset_path, treated_mice_metadata_fn1)
            )
        treated_mice_expr1 = self.read_soft_expression(
            os.path.join(dataset_path, treated_mice_fn1), metadata = treated_mice_metadata1,
            array_type='affymetrix'
            )
        # decided not to use the other 4 mice bc they have big batch effect that I don't know how to correct
        """treated_mice_metadata2 = self.read_soft(
            os.path.join(dataset_path, treated_mice_metadata_fn2)
            )
        treated_mice_expr2 = self.read_soft_expression(
            os.path.join(dataset_path, treated_mice_fn2), metadata = treated_mice_metadata2,
            array_type='affymetrix'
            )
        # concat the metadata
        treated_mice_metadata = pd.concat([treated_mice_metadata1, treated_mice_metadata2])
        treated_mice_metadata.reset_index(drop=True, inplace=True)
        # merge expression data, matching columns
        treated_mice_expr = pd.merge(
            treated_mice_expr1.T, treated_mice_expr2.T, left_index=True, right_index=True
            ).T"""
        # create expression dataset object
        treated_mice = ExpressionDataset(
            expression_df=treated_mice_expr1,
            species="mouse",
            metadata_df=treated_mice_metadata1,
            dataset="zhou_2012"
            )
        return treated_mice
    
    def load_ma_sc_rat_CR(self) -> ScExpressionDataset:
        """Load the Ma et al. 2020 CR rat dataset
        ### Returns
        adata : ScExpressionDataset
            Expression dataset for CR rat
        """
        from glob import glob
        # directory filled with barcodes.tsv.gz, genes.tsv.gz, and matrix.mtx.gz files
        data_dir = "/cellar/users/zkoch/dream/data/ma_sc_rat_CR"
        adata_fn = 'all_tissues_with_ssgsea.h5ad'
        metadata_fn = 'GSE137869_family.soft'
        adata = sc.read_h5ad(os.path.join(data_dir, adata_fn))
        # read metadata
        # already combined in adata
        """metadata = self.read_soft(os.path.join(data_dir, metadata_fn))
        metadata['sample'] = metadata['sample_geo_accession'] + '_' + metadata['sample_title']
        # join with adata.obs
        adata.obs = adata.obs.merge(
            metadata, left_on = 'sample',
            right_on='sample', how='left'
            )"""
        ma_sc_rat_CR = ScExpressionDataset(
            adata=adata,
            gene_species="rat",
            dataset="ma_sc_rat_CR"
            )
        return ma_sc_rat_CR
    
    def load_palovics_parabiosis(self) -> ScExpressionDataset:
        """Load the Palovics Parabiosis dataset
        ### Returns
        adata : ScExpressionDataset
            Expression dataset for parabiosis
        """
        fn = "/cellar/users/zkoch/dream/data/palovics_parabiosis/all_tissues.h5ad"
        adata = sc.read_h5ad(fn)
        palovics_parabiosis = ScExpressionDataset(
            adata=adata,
            gene_species="mouse",
            dataset="palovics_parabiosis"
            )
        return palovics_parabiosis
        
    def load_tms(self) -> ScExpressionDataset:
        """Load the Tabula Muris Senis dataset
        ### Returns
        facs : ScExpressionDataset
            SCExpression dataset for FACS
        """
        tms_data_dir = "/cellar/users/zkoch/dream/data/tabula_muris_senis"
        facs_fn = "tabula-muris-senis-zane-processed-facs_with_mutation_counts_ssgsea.h5ad"
        expr_counts_for_mutation_fn = "adata_with_ercc_genecode_counts_for_gatk_with_metadata.h5ad"
        mutation_counts_fn = "adata_with_ercc_gatk_all_data_with_metadata.h5ad"
        # read annData
        facs = sc.read_h5ad(os.path.join(tms_data_dir, facs_fn))
        # read mutations
        expr_counts_for_mutation = sc.read_h5ad(
            os.path.join(tms_data_dir, expr_counts_for_mutation_fn)
            )
        expr_counts_for_mutation.var_names = [x.upper() for x in expr_counts_for_mutation.var_names]

        mutation_counts = sc.read_h5ad(
            os.path.join(tms_data_dir, mutation_counts_fn)
            )
        mutation_counts.var_names = [x.upper() for x in mutation_counts.var_names]

        # create ScExpressionDataset object
        tms = ScExpressionDataset(
            adata=facs,
            gene_species="human",
            dataset="tms",
            expr_counts_for_mutation = expr_counts_for_mutation,
            mutation_counts = mutation_counts
            )
        return tms

    def load_martin_montalvo(self) -> ExpressionDataset:
        """Load the Martin Montalvo dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        # paths
        dataset_path = "/cellar/users/zkoch/dream/data/martin_montalvo"
        treated_mice_fn = "GSE40936_expression.soft"
        treated_mice_metadata_fn = "GSE40936_family.soft"

        if not self.load_expression:
            raise ValueError("load_expression must be True to load Martin Montalvo")
        # load data
        treated_mice_metadata = self.read_soft(
            os.path.join(dataset_path, treated_mice_metadata_fn)
            )
        treated_mice_expr = self.read_soft_expression(
            os.path.join(dataset_path, treated_mice_fn), metadata = treated_mice_metadata,
            array_type='illumina'
            )
        # create expression dataset object
        treated_mice = ExpressionDataset(
            expression_df=treated_mice_expr,
            species="mouse",
            metadata_df=treated_mice_metadata,
            dataset="martin_montalvo"
            )
        return treated_mice
    
    def load_boutant_nestle(self) -> ExpressionDataset:
        """Load the Boutant Nestle dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        dataset_path = "/cellar/users/zkoch/dream/data/boutant_nestle"
        treated_mice_fn = "GSE70857_non-normalized.txt"
        treated_mice_metadata_fn = "GSE70857_family.soft"
        if not self.load_expression:
            raise ValueError("load_expression must be True to load Boutant")
        # load data
        treated_mice_expr = pd.read_csv(
            os.path.join(dataset_path, treated_mice_fn),
            index_col=0, sep = "\t"
            )
        treated_mice_metadata = self.read_soft(
            os.path.join(dataset_path, treated_mice_metadata_fn)
            )
        
        # create expression dataset objects
        treated_mice = ExpressionDataset(
            expression_df=treated_mice_expr.T,
            species="mouse",
            metadata_df=treated_mice_metadata,
            dataset="boutant_nestle"
            )
        return treated_mice
        
    def load_tyshkovskiy(self) -> ExpressionDataset:
        """Load the Tyshkovskiy dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        dataset_path = "/cellar/users/zkoch/dream/data/tyshkovskiy"
        treated_mice_fn = "GSE131754_Interventions_assigned_reads.txt.gz"
        treated_mice_metadata_fn = "GSE131754_family.soft"
        if not self.load_expression:
            raise ValueError("load_expression must be True to load Tyshkovskiy")
        # load data
        treated_mice_expr = pd.read_csv(
            os.path.join(dataset_path, treated_mice_fn),
            index_col=0, sep = "\t"
            )
        treated_mice_metadata = self.read_soft(
            os.path.join(dataset_path, treated_mice_metadata_fn)
            )
        # create expression dataset objects
        treated_mice = ExpressionDataset(
            expression_df=treated_mice_expr.T,
            species="mouse",
            metadata_df=treated_mice_metadata,
            dataset="tyshkovskiy"
            )
        return treated_mice

    def load_mSalt(self) -> List[ExpressionDataset]:
        """Load the mSalt dataset
        ### Returns
        across_species : ExpressionDataset
            Expression dataset for across species
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        # paths
        dataset_path = "/cellar/users/zkoch/dream/data/msalt"
        across_species_fn = "across_species.csv"
        across_species_metadata_fn = "across_species_meta.soft"
        treated_mice_fn = "treated_mice.csv"
        treated_mice_metadata_fn = "treated_mice_meta.soft"
        if not self.load_expression:
            raise ValueError("load_expression must be True to load mSalt")
        # load data
        across_species_expr = pd.read_csv(
            os.path.join(dataset_path, across_species_fn),
            index_col=0
            )
        treated_mice_expr = pd.read_csv(
            os.path.join(dataset_path, treated_mice_fn),
            index_col=0
            )
        # read metadata from .soft files
        across_species_metadata = self.read_soft(
            os.path.join(dataset_path, across_species_metadata_fn)
            )
        treated_mice_metadata = self.read_soft(
            os.path.join(dataset_path, treated_mice_metadata_fn)
            )
        # create expression dataset objects
        across_species = ExpressionDataset(
            expression_df=across_species_expr.T,
            species="mouse", # df is indxd by mouse genes
            metadata_df=across_species_metadata,
            dataset="mSalt"
            )
        treated_mice = ExpressionDataset(
            expression_df=treated_mice_expr.T,
            species="mouse",
            metadata_df=treated_mice_metadata,
            dataset="mSalt"
            )
        return across_species, treated_mice
        
    def load_tcga(self) -> List[Union[ExpressionDataset, MutationDataset, MethylationDataset]]:
        """Load a TCGA dataset based on self.dataset_name
        ### Returns
        A list of a subset of ExpressionDataset, MutationDataset, MethylationDataset
        """
        # find the files for this dataset
        processed_data_dir = "/cellar/users/zkoch/dream/data/tcga/processed_data"
        dataset_files = glob.glob(os.path.join(processed_data_dir, self.dataset_name +  "*"))
        # initialize empty fns
        metadata_fn, expression_fn, mutation_fn, methylation_dir, methylation_fn = "","","","",""
        for dataset_fn in dataset_files:
            # if the fn contains metadata
            if "metadata" in dataset_fn:
                metadata_fn = dataset_fn
            # only read in parquet expression
            elif "expression" in dataset_fn and ".tsv" not in dataset_fn:
                expression_fn = dataset_fn
            elif "mutation" in dataset_fn:
                mutation_fn = dataset_fn
            elif "methylation_dask" in dataset_fn:
                methylation_dir = dataset_fn
            elif "methylation.parquet" in dataset_fn:
                methylation_fn = dataset_fn
            else:
                print(f"Unknown file type: {dataset_fn}")
        # create datasets based on what files we have
        expression_dataset, mutation_dataset, methylation_dataset = None, None, None
        if metadata_fn == "":
            raise ValueError(f"Metadata file not found for {self.dataset_name}")
        else:
            metadata_df = pd.read_parquet(metadata_fn)
            print(f"Loaded metadata for {self.dataset_name}")
        if expression_fn != "" and self.load_expression:
            expression_df = pd.read_parquet(expression_fn)
            # get metadata rows related to expression
            expr_metadata_df = metadata_df.query(
                "data_type == 'expression'"
                ).reset_index(drop=True)
            expr_metadata_df.set_index("sample_id", inplace=True)
            expression_dataset = ExpressionDataset(
                expression_df=expression_df.T,
                species="human",
                metadata_df=expr_metadata_df,
                dataset="tcga"
                )
            print(f"Created ExpressionDataset for {self.dataset_name}")
        if mutation_fn != "" and self.load_mutation:
            mutation_df = pd.read_parquet(mutation_fn)
            mutation_metadata_df = metadata_df.query(
                "data_type == 'mutation'"
                ).reset_index(drop=True)
            mutation_metadata_df.set_index("sample_id", inplace=True)
            mutation_dataset = MutationDataset(
                mutation_df=mutation_df,
                metadata_df=mutation_metadata_df
                )
            print(f"Created MutationDataset for {self.dataset_name}")
        if methylation_dir != "" and self.load_methylation:
            print("reading dask methylation df")
            methylation_dd = dd.read_parquet(methylation_fn)
            # convert to pandas df
            methylation_df = methylation_dd.compute()
            # get metadata rows related to methylation
            methylation_metadata_df = metadata_df.query(
                "data_type == 'methylation'"
                ).reset_index(drop=True)
            methylation_metadata_df.set_index("sample_id", inplace=True)
            methylation_dataset = MethylationDataset(
                methylation_df=methylation_df.T,
                metadata_df=methylation_metadata_df
                )
            print(f"Created MethylationDataset for {self.dataset_name}")
        # load from methylation fn if we don't have a methylation dir
        elif methylation_fn != "" and self.load_methylation:
            methylation_df = pd.read_parquet(methylation_fn)
            # get metadata rows related to methylation
            methylation_metadata_df = metadata_df.query(
                "data_type == 'methylation'"
                ).reset_index(drop=True)
            methylation_metadata_df.set_index("sample_id", inplace=True)
            methylation_dataset = MethylationDataset(
                methylation_df=methylation_df,
                metadata_df=methylation_metadata_df
                )
            print(f"Created MethylationDataset for {self.dataset_name}")
            
        return [expression_dataset, mutation_dataset, methylation_dataset]
        
    def read_soft(
        self,
        fn : str
        ) -> pd.DataFrame:
        """Read metadata from .soft file
        ### Parameters
        fn : str
            Path to .soft file
        ### Returns
        metadata : pd.DataFrame
            Metadata dataframe
        """
        # create defaultDict filled with empty lists
        metadata = defaultdict(list)
        
        with open(fn, "r") as f:
            # iterate over lines
            lines = f.readlines()
            for i in range(len(lines)):
                # until we find a new sample, keep going
                if lines[i].startswith("^SAMPLE"):
                    # once we find a new sample, add it as a key
                    sample_name = lines[i].strip("^SAMPLE = ").strip()
                    metadata['sample'].append(sample_name)
                    # then continue looping over subsequent lines
                    not_next_sample = True
                    while not_next_sample:
                        i += 1
                        # if we find a new sample, break
                        if lines[i].startswith("^SAMPLE") or i == len(lines)-1:
                            i -= 1
                            not_next_sample = False
                            break
                        # otherwise, add things we want
                        else:
                            """elif lines[i].startswith("!Sample_description"):
                                    metadata['sample_name'].append(
                                        lines[i].strip("!Sample_description = ").strip()
                                        )"""
                            if lines[i].startswith("!Sample_geo_accession"):
                                metadata['sample_geo_accession'].append(
                                    lines[i].strip("!Sample_geo_accession = ").strip()
                                    )
                            if lines[i].startswith("!Sample_title"):
                                metadata['sample_title'].append(
                                    lines[i].strip("!Sample_title = ").strip()
                                    )
                            elif lines[i].startswith("!Sample_characteristics_ch1 = Sex:"):
                                # why have to do this way?
                                to_remove = "!Sample_characteristics_ch1 = Sex: " 
                                metadata['sex'].append(
                                    lines[i][len(to_remove):].strip()
                                    )
                            elif lines[i].startswith("!Sample_characteristics_ch1 = gender:"):
                                # why have to do this way?
                                to_remove = "!Sample_characteristics_ch1 = gender: " 
                                metadata['sex'].append(
                                    lines[i][len(to_remove):].strip()
                                    )
                            elif lines[i].startswith("!Sample_characteristics_ch1 = genotype:"):
                                # why have to do this way?
                                to_remove = "!Sample_characteristics_ch1 = genotype: " 
                                metadata['intervention'].append(
                                    lines[i][len(to_remove):].strip()
                                    )
                            elif lines[i].startswith("!Sample_characteristics_ch1 = age:"):
                                to_remove = "!Sample_characteristics_ch1 = age: " 
                                metadata['age'].append(
                                    lines[i][len(to_remove):].strip()
                                    )
                            elif lines[i].startswith("!Sample_characteristics_ch1 = treatment"):
                                to_remove = "!Sample_characteristics_ch1 = treatment: "
                                metadata['treatment'].append(
                                    lines[i][len(to_remove):].strip()
                                    )
                            elif lines[i].startswith("!Sample_characteristics_ch1 = intervention duration"):
                                to_remove = "!Sample_characteristics_ch1 = intervention duration: "
                                metadata['intervention_duration'].append(
                                    lines[i][len(to_remove):].strip()
                                    )
                            elif lines[i].startswith("!Sample_characteristics_ch1 = intervention"):
                                to_remove = "!Sample_characteristics_ch1 = intervention: "
                                metadata['intervention'].append(
                                    lines[i][len(to_remove):].strip()
                                    )
                            elif lines[i].startswith("!Sample_characteristics_ch1 = dose"):
                                to_remove = "!Sample_characteristics_ch1 = dose: "
                                metadata['dose'].append(
                                    lines[i][len(to_remove):].strip()
                                    )
                            elif lines[i].startswith("!Sample_characteristics_ch1 = Sex"):
                                to_remove = "!Sample_characteristics_ch1 = Sex: "
                                metadata['sex'].append(
                                    lines[i][len(to_remove):].strip()
                                    )      
                            elif lines[i].startswith("!Sample_characteristics_ch1 = tissue"):
                                to_remove = "!Sample_characteristics_ch1 = tissue: "
                                metadata['tissue'].append(
                                    lines[i][len(to_remove):].strip()
                                    )  
                            elif lines[i].startswith("!Sample_characteristics_ch1 = Stage"):
                                to_remove = "!Sample_characteristics_ch1 = Stage: "
                                metadata['condition'].append(
                                    lines[i][len(to_remove):].strip()
                                    )  
                                      
        # turn into dataframe
        # if the lists are not the same lenght pad with nans
        max_len = max([len(x) for x in metadata.values()])
        for key in metadata.keys():
            if len(metadata[key]) < max_len:
                metadata[key] = metadata[key] + [np.nan] * (max_len - len(metadata[key]))
        metadata = pd.DataFrame(metadata)
        # drop rows that are all nan
        metadata.dropna(how='any', inplace=True, axis = 0)
        return metadata

    def read_soft_expression(
        self,
        fn : str,
        metadata : pd.DataFrame,
        array_type : str = "illumina"
        ) -> pd.DataFrame:
        """Read expresson from .soft file
        ### Parameters
        fn : str
            Path to exxpression.soft file
        metadata : pd.DataFrame
            Metadata dataframe to get sample names and number frmo
        array_type : str
            Type of array used to generate the expression data [illumina | affymetrix]
        ### Returns
        expression_df : pd.DataFrame
            expression_df dataframe
        """
        #samples = treated_mice_metadata['sample_title'] 
        expression_df = pd.read_csv(fn, sep = '\t', index_col = None, header=None)
        if array_type == "illumina":
            if len(expression_df.columns) == 3:
                expression_df.columns = ['probe_name', 'raw_expression', 'z_log_expr']
            elif len(expression_df.columns) == 2:
                expression_df.columns = ['probe_name', 'raw_expression']
            else:
                raise ValueError(f"Number of columns in expression_df not recognized: {len(expression_df.columns)}")
        elif array_type == "affymetrix":
            expression_df.columns = ['probe_name', 'raw_expression']
        else:
            raise ValueError(f"array_type {array_type} not recognized")
        length = expression_df.shape[0]
        num_samples = metadata.shape[0]
        genes_per_sample = length / num_samples
        # check if number of genes per sample is integer
        assert genes_per_sample.is_integer(), f"Number of genes per sample is not integer: {genes_per_sample}"
        # if it is reshape to samples x genes
        expression_df = pd.DataFrame(
            expression_df['raw_expression'].values.reshape((num_samples, int(genes_per_sample))),
            index = metadata['sample_title'],
            columns = expression_df['probe_name'].values[:int(genes_per_sample)]
        )
        expression_df.fillna(0, inplace = True)
        return expression_df

