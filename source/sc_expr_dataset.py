import pandas as pd
import scanpy as sc 
from pybiomart import Dataset
import statsmodels.formula.api as smf
import gseapy
import numpy as np
import utils
import colorcet as cc
import matplotlib.pyplot as plt
import seaborn as sns

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
        
    def get_n_genes_and_counts(self) -> None:
        """Get the number of genes detected in each cell and total counts"""
        self.adata.obs['n_genes'] = (self.adata.X > 0).sum(axis=1)
        self.adata.obs['n_counts'] = self.adata.X.sum(axis=1)

    def get_dream_gene_expression(
        self,
        convert = True,
        use_new: bool = False
        ) -> pd.DataFrame:
        """Get expression of DREAM genes
        ### Parameters:
        convert : bool
            Wether to convert the gene names of adata to match the dream gene names
        use_new : bool
            Wether to use the new dream file
        ### Returns:
        None
        """
        self._read_and_convert_dream(convert, use_new)
        # create a subset of adata that is just the dream genes
        self.dream_expression = self.adata[
            :, self.dream_regulated_genes_w_expression
            ].copy()
        # calculate mean
        self.dream_expression.obs['mean_dream_activity'] = np.mean(
            self.dream_expression.X, axis=1
            )
        # get resid of mean dream activity
        self.scale_dream_by_seq_depth(col_name = 'mean_dream_activity')
        self.adata.obs['mean_dream_activity_resid'] = self.dream_expression.obs['mean_dream_activity_resid']

    def _read_and_convert_dream(
        self, 
        convert: bool = True,
        use_new: bool = False
        ):
        print("using new dream file: ", use_new)
        self.dream_regulated_genes = utils.read_dream_files(use_new)
        # get the DREAM regulated genes we want, index is gene names

        if convert and self.gene_species != "human" and not self.already_converted:
            # do gene conversion using 
            self._convert_genes_to_human_ensembl()
        dream_regulated_genes_names = self.dream_regulated_genes.index
        
        # intersection of symbols
        self.dream_regulated_genes_w_expression = list(
            set(self.adata.var_names).intersection(
                set(dream_regulated_genes_names)
                )
            )
        print(f"Found {len(self.dream_regulated_genes_w_expression)} DREAM genes with expression")
        
    def scale_dream_by_seq_depth(self, col_name: str, eq: str = None, invert: bool = True) -> None:
        if eq is None:
            # use n_genes for tms and n_genes_by_counts for rosmap
            eq = f'{col_name} ~ n_counts * n_genes'
            # total_counts * n_genes_by_counts'  
            # n_genes is same as n_genes_by_counts (calc'd by CZI but is not changed)
            # total_counts is the total counts of the cell after my QC
            # n_counts is total before any preproc or QC
        mut_ols = smf.ols(
            formula=eq,
            data=self.dream_expression.obs
            ).fit()
        self.dream_expression.obs[f'{col_name}_resid'] = mut_ols.resid
        # invert residuals making highest values low and lowest values high
        if invert:
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
            # get a mapping from ensembl mouse gene names to human gene names
            ensembl_gene_mapper = dataset.query(
                attributes=['ensembl_gene_id', 'external_gene_name',
                            'hsapiens_homolog_orthology_type',
                            'hsapiens_homolog_orthology_confidence',
                            'hsapiens_homolog_associated_gene_name']
                )
            # convert both to uppercase
            ensembl_gene_mapper['Gene name'] = ensembl_gene_mapper['Gene name'].str.upper()
            self.adata.var.index = self.adata.var.index.str.upper()
            # only keep genes with a human homolog and unduplicate
            ensembl_gene_mapper.dropna(subset = 'Human homology type', inplace = True)
            ensembl_gene_mapper.drop_duplicates(
                subset = 'Gene name', inplace = True
                )
            # map the gene names through mouse gene names to human gene names
            self.adata.var['human_gene_name'] = self.adata.var.index.map(
                ensembl_gene_mapper.set_index('Gene name')['Human gene name']
                )
            self.adata.var.dropna(subset = ['human_gene_name'], inplace = True)
            self.adata.var.drop_duplicates(subset = 'human_gene_name', inplace = True)
            self.adata = self.adata[:, self.adata.var.index]
            # remove genes dropped from var
            self.adata.var.set_index('human_gene_name', inplace = True, drop = False)
            self.adata = self.adata[:, self.adata.var.index]
            self.adata.var_names = self.adata.var.index
            self.adata = self.adata[:, self.adata.var_names]
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
        elif self.gene_species == 'symbol':
            # use gene symbols instead of ensembl
            self.dream_regulated_genes.set_index('gene_name', inplace=True, drop = False)
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
        
    def dream_enrichment_ssgsea(
        self,
        use_new: bool = False,
        threads: int = 1
        ) -> None:
        """Run ssgsea on the DREAM genes
        ### Parameters:
        use_new : bool
            Wether to use the new dream file
        ### Returns:
        None
        """
        
        print("using new dream file: ", use_new)
        self.dream_regulated_genes = utils.read_dream_files(use_new)
        # get dream gene to check for enrichment in expression
        dream_gene_set_dict = {
            'dream_reg_genes': self.dream_regulated_genes_w_expression
            }
        # run ssgsea
        ssgsea = gseapy.ssgsea(
            data=self.adata.to_df().T, # get a df version of the expr
            gene_sets=dream_gene_set_dict,
            outdir=None, no_plot=True, verbose=True,
            threads=threads
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
        self.scale_dream_by_seq_depth(col_name = 'DREAM_normalized_enrichment_score')
        self.scale_dream_by_seq_depth(col_name = 'DREAM_enrichment_score')
        
    

    def _get_gene_length(
        self, 
        gene_name:str,
        gencode_m34:pd.DataFrame,
        type: str
        ) -> int:
        """Get the length of a gene from gencode_m34
        ### Parameters:
        gene_name : str
            Gene name to get length of
        gencode_m34 : pd.DataFrame
            Gencode m34 dataframe
        type : str
            Feature of gene to get length of, 'gene' or 'transcript'
        
        """
        gene = gencode_m34.query("gene_name == @gene_name and type == @type")
        if gene.shape[0] > 1:
            #print(gene)
            # use just first row
            gene = gene.iloc[0]
        elif gene.shape[0] == 0:
            return -1
        # 
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
    
    def calc_expression_change_by_length_bin(
        self,
        length_col = 'gene_length_bp',
        bin_type = 'quantile',
        dream_col = 'DREAM_normalized_enrichment_score_resid',
        n_bins = 5,
        skip_binning = False
        ):
        """Assign genes to length bins and calculate change in expression of each gene across age
        ### Parameters
        length_col : str
            Column name in adata.var that contains gene lengths ('gene_length_bp' or 'transcript_length_bp')
        bin_type : str
            Type of binning to use, 'quantile' or 'equal'
        dream_col : str
            Column name in adata.obs that contains DREAM expression
        n_bins : int
            Number of bins to use
        skip_binning : bool
            Skip binning and use existing bins
        ### Returns
        None
        """
        if not skip_binning:
            expression_bin_col = length_col + '_bin'
            # split genes into bins based on gene length
            if bin_type == 'quantile':
                self.adata.var[expression_bin_col] = pd.qcut(self.adata.var[length_col], q=n_bins)
            elif bin_type == 'equal':
                self.adata.var[expression_bin_col] = pd.cut(self.adata.var[length_col], bins=n_bins)
            # get list of genes in each group
            gene_groups = self.adata.var.groupby(expression_bin_col).groups
            # get the mean expression of each group of genes in each cell
            for group, genes in gene_groups.items():
                self.adata.obs[f'mean_expr_{group}'] = self.adata[:, genes].X.mean(axis=1)
        # for each gene, 
        
    
    def calc_expression_by_length_bin(
        self,
        length_col = 'gene_length_bp', 
        bin_type = 'quantile',
        bin_edges = [],
        dream_col = 'DREAM_normalized_enrichment_score_resid',
        n_bins = 5,
        skip_binning = False
        ):
        """Assign genes to length bins and calculate mean expression of each bin in each cell.
        ### Parameters
        length_col : str
            Column name in adata.var that contains gene lengths ('gene_length_bp' or 'transcript_length_bp')
        bin_type : str
            Type of binning to use, 'quantile' or 'equal'
        dream_col : str
            Column name in adata.obs that contains DREAM expression
        n_bins : int
            Number of bins to use
        skip_binning : bool
            Skip binning and use existing bins
        None
        """
        
        if not skip_binning:
            expression_bin_col = length_col + '_bin'
            # split genes into bins based on gene length
            if bin_type == 'quantile':
                self.adata.var[expression_bin_col] = pd.qcut(self.adata.var[length_col], q=n_bins)
            elif bin_type == 'equal':
                if len(bin_edges) > 0:
                    self.adata.var[expression_bin_col] = pd.cut(self.adata.var[length_col], bins=bin_edges)
                else:
                    self.adata.var[expression_bin_col] = pd.cut(self.adata.var[length_col], bins=n_bins)
            # get list of genes in each group
            gene_groups = self.adata.var.groupby(expression_bin_col).groups
            # get the mean expression of each group of genes in each cell
            for group, genes in gene_groups.items():
                self.adata.obs[f'mean_expr_{group}'] = self.adata[:, genes].X.mean(axis=1)
                
        # create long df of mean expression by length bin
        length_cols = [x for x in self.adata.obs.columns if 'mean_expr' in x]
        cols_to_select = ['age_months'] + length_cols + ['cell']
        expression_by_length_and_age = self.adata.obs[
            cols_to_select
            ].melt(
                id_vars=['age_months', 'cell']
                ).rename(
                    columns={'variable': 'length_bin', 'value': 'length_mean_expression'}
                    )
        expression_by_length_and_age.query("age_months != 21", inplace=True)
        
        # add dream expression
        expression_by_length_and_age['DREAM_normalized_enrichment_score_resid'] = expression_by_length_and_age['cell'].map(
            self.dream_expression.obs.set_index('cell')['DREAM_normalized_enrichment_score_resid']
            )
                   
        return expression_by_length_and_age
    
    def calc_mutation_rate(
        self,
        mutation_column: str,
        min_callable_sites: int = 1e5,
        callable_sites_col: str = 'SitesPerCell',
        min_cell_per_subclass: int = 100,
        use_lme: bool = False,
        age_column: str = 'Age at Death',
        cell_type_column: str = 'Subclass',
        donor_column: str = 'Donor ID',
        dream_activity_bin: str = None
    ) -> pd.DataFrame:
        """Calculate mutation rates per year for each cell subclass.
        
        Args:
            mutation_column: Column name containing mutation counts
            min_callable_sites: Minimum number of callable sites required per cell
            callable_sites_col: Column name containing number of callable sites
            min_cell_per_subclass: Minimum number of cells required per subclass
            use_lme: Whether to use linear mixed effects model
            age_column: Column name containing age information
            cell_type_column: Column name containing cell subclass information
            donor_column: Column name containing donor ID information
            
        Returns:
            DataFrame containing mutation rates and statistics per subclass
        """
        slopes = []
        slope_se = []
        r2s = []
        pvals = []
        
        # Subset to cells with minimum callable sites
        regr_df = self.dream_expression.obs.query(f"{callable_sites_col} >= @min_callable_sites")
        regr_df.dropna(subset=[mutation_column, age_column], inplace=True)
        
        # TODO: EXPERIMENTAL
        regr_df = regr_df.query("mutation_burden_per_genome_nanfilled < 1e5")
        
        # groupby 'DREAM_normalized_enrichment_score_resid' and keep top 20% of cells by tissue
        if dream_activity_bin is not None:
            # Group by tissue and donor_id and only create bins where there is variation
            regr_df['dream_activity_bin'] = regr_df.groupby(['tissue']).apply(
                lambda x: pd.qcut(x['DREAM_normalized_enrichment_score_resid'], 5, labels=['Low', 'Low-Medium', 'Medium','Medium-High', 'High']) 
                if len(x['DREAM_normalized_enrichment_score_resid'].unique()) > 1 
                else pd.Series(['Medium'] * len(x), index=x.index)
            ).reset_index(level=[0,1], drop=True)
            regr_df = regr_df.query("dream_activity_bin == @dream_activity_bin")
        
        # Filter celltypes with too few cells
        ct_counts = regr_df.query(f"{mutation_column} > 0")[cell_type_column].value_counts()
        celltypes_to_keep = ct_counts[ct_counts >= min_cell_per_subclass].index
        regr_df = regr_df.loc[regr_df[cell_type_column].isin(celltypes_to_keep)]
        
        # Convert celltype to categorical and clean unused categories
        regr_df[cell_type_column] = regr_df[cell_type_column].astype('category')
        regr_df[cell_type_column] = regr_df[cell_type_column].cat.remove_unused_categories()
        
        # Calculate rates for each celltype
        for celltype in regr_df[cell_type_column].unique():
            celltype_data = regr_df.query(f"{cell_type_column} == @celltype").reset_index(drop=True)
            if use_lme:
                try:
                    # Linear mixed effects model
                    celltype_data = celltype_data.dropna(subset=[mutation_column, donor_column, age_column])
                    lme = smf.mixedlm(
                        f'{mutation_column} ~ {age_column} + 0',
                        data=celltype_data,
                        groups=celltype_data[donor_column]
                    ).fit()
                    std_err = float(lme.summary().tables[1].loc[age_column, 'Std.Err.'])
                    slope = lme.params[age_column]
                    slopes.append(slope)
                    slope_se.append(std_err)
                    r2s.append(None)  # No R2 for LME
                    pvals.append(lme.pvalues[age_column])
                except:
                    print(celltype)
                    print(lme.summary())
                    slopes.append(np.nan)
                    slope_se.append(np.nan)
                    r2s.append(np.nan)
                    pvals.append(np.nan)
                    # pop cell type
                    continue
                    
            else:
                # Ordinary least squares
                ols = smf.ols(
                    f'{mutation_column} ~ {age_column} + 0',
                    data=celltype_data
                ).fit()
                slopes.append(ols.params[age_column])
                std_err = float(ols.bse[age_column])
                slope_se.append(std_err)
                r2s.append(ols.rsquared)
                pvals.append(ols.pvalues[age_column])
        
        # Create results dataframe
        results = pd.DataFrame({
            cell_type_column: regr_df[cell_type_column].unique(),
            f'{mutation_column}_per_year': slopes,
            f'{mutation_column}_per_year_se': slope_se,
            f'{mutation_column}_per_year_pval': pvals,
        })
        
        if not use_lme:
            results[f'{mutation_column}_per_year_r2'] = r2s
            
        return results
    
    def plot_mutation_rate(
        self,
        mutation_rate_df: pd.DataFrame,
        mutation_column: str,
        cell_type_column: str = 'Subclass',
        lineage_column: str = 'lineage',
        subset_to_lineage: str = None,
        axes: plt.Axes = None,
        ):
        """Plot mutation rate per year for each cell type.
        
        Args:
            mutation_rate_df: DataFrame containing mutation rates and statistics per subclass (output of calc_mutation_rate)
            cell_type_column: Column name containing cell subclass information
        """
        if subset_to_lineage is not None:
            # check if lineage_column is in mutation_rate_df
            if lineage_column not in mutation_rate_df.columns:
                mapper = self.dream_expression.obs[[cell_type_column, lineage_column]].drop_duplicates().reset_index(drop = True).set_index(cell_type_column)
                mutation_rate_df['lineage'] = mutation_rate_df[cell_type_column].map(mapper['lineage'])
            mutation_rate_df = mutation_rate_df.query(f"{lineage_column} == @subset_to_lineage")
            mutation_rate_df[cell_type_column] = mutation_rate_df[cell_type_column].astype('category')
            mutation_rate_df[cell_type_column] = mutation_rate_df[cell_type_column].cat.remove_unused_categories()
            
        # sort by mutation rate
        mutation_rate_df = mutation_rate_df.sort_values(by=f'{mutation_column}_per_year', ascending=False)
        mutation_rate_df.reset_index(drop=True, inplace=True)
        
        # create categorical color map from cc.glasbey_cool that is len(mutation_rate_df[cell_type_column])
        cmap = cc.glasbey[:len(mutation_rate_df[cell_type_column])]
        # plot
        if axes is None:
            fig, ax = plt.subplots(figsize=(10, 4))
        else:
            ax = axes
        ax.bar(
            mutation_rate_df[cell_type_column], mutation_rate_df[f'{mutation_column}_per_year'],
            yerr = mutation_rate_df[f'{mutation_column}_per_year_se'],
            capsize = 3.5, color = cmap
        )
        # angle x-axis labels
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        ax.set_ylabel(f'Somatic mutations per year')
        sns.despine(ax=ax)
        return ax