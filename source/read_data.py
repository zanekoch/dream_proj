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
import anndata as ad
from typing import Union, List
import statsmodels.formula.api as smf
import pickle

# my library code 
from expr_dataset import ExpressionDataset
from sc_expr_dataset import ScExpressionDataset
from mutation_dataset import MutationDataset
from methyl_dataset import MethylationDataset
import utils

"""
To add a new expression dataset
1. Add to DatasetLoader.load_dataset
2. Write a load_{dataset_name} function
3. Edit DatasetLoader.read_soft if need be
4. Add the dataset to ExpressionDataset.__init__, making the metadata mergable with the expression df
5. Add dataset to MetaExpressionDataset.run_GSEA with proper control etc.
"""

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

    def load_dataset(
        self,
        kws: dict = {
            'read_whole_dataset': False,
            'start_cell_num': -1,
            'end_cell_num': -1
            }
        ) -> Union[
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
        elif self.dataset_name == 'yu_2012':
            dataset_list = self.load_yu_2012()
        elif self.dataset_name == 'mercken_2014':
            dataset_list = self.load_mercken_2014()
        elif self.dataset_name == 'zhang_2023':
            dataset_list = self.load_zhang_2023()
        elif self.dataset_name == 'neff_2013':
            dataset_list = self.load_neff_2013()
        elif self.dataset_name == 'eisenberg_2016':
            dataset_list = self.load_eisenberg_2016()
        elif self.dataset_name == 'gtex':
            dataset_list = self.load_gtex(kws)
        elif self.dataset_name == 'aon_2020':
            dataset_list = self.load_aon_2020()
        elif self.dataset_name == 'barger_2008':
            dataset_list = self.load_barger_2008()
        elif self.dataset_name == 'pearson_2008':
            dataset_list = self.load_pearson_2008()
        elif self.dataset_name == 'nebulas':
            dataset_list = self.load_nebulas()
        elif self.dataset_name == 'nebulas_sc':
            dataset_list = self.load_nebulas_sc()
        elif self.dataset_name == 'liu_2023':
            dataset_list = self.load_liu_2023()
        elif self.dataset_name == 'cao_2024':
            dataset_list = self.load_cao_2024()
        elif self.dataset_name == 'bujarrabal_dueso_2023_u2osRNA':
            dataset_list = self.load_bujarrabal_dueso_2023()
        elif self.dataset_name == 'bujarrabal_dueso_2023_wormRNA':
            dataset_list = self.load_bujarrabal_dueso_2023_wormRNA()
        elif self.dataset_name == 'uxa_2019':
            dataset_list = self.load_uxa_2019()
        elif self.dataset_name == 'motrpac_2024':
            dataset_list = self.load_motrpac_2024()
        elif self.dataset_name == 'gyenis_2023':
            dataset_list = self.load_gyenis_2023(kws)
        elif self.dataset_name == 'mammalian_methylation_consort':
            dataset_list = self.load_mammalian_methylation_consort()
        elif self.dataset_name == 'otero-garcia_2022':
            dataset_list = self.load_otero_garcia_2022()
        elif self.dataset_name == 'hashimoto_2019':
            dataset_list = self.load_hashimoto_2019()
        elif self.dataset_name == 'lu_2014':
            dataset_list = self.load_lu_2014()
        elif self.dataset_name == 'williams_2022':
            dataset_list = self.load_williams_2022()
        elif self.dataset_name == 'lu_2022':
            dataset_list = self.load_lu_2022()
        elif self.dataset_name == 'paine_2024':
            dataset_list = self.load_paine_2024()
        elif self.dataset_name == 'synapse_rna_seq_harmonization':
            dataset_list = self.load_synapse_rna_seq_harmonization()
        elif self.dataset_name == 'synapse_MIT_ROSMAP_Multiomics_10x':
            dataset_list = self.load_synapse_MIT_ROSMAP_Multiomics_10x()
        elif self.dataset_name == 'synapse_MIT_ROSMAP_Multiomics_smartseq':
            dataset_list = self.load_synapse_MIT_ROSMAP_Multiomics_smartseq()
        elif self.dataset_name == 'SEA-AD':
            dataset_list = self.load_SEA_AD(read_whole_dataset = kws['read_whole_dataset'], start_cell_num = kws['start_cell_num'], end_cell_num = kws['end_cell_num'])
        elif self.dataset_name == 'tabula_sapiens':
            dataset_list = self.load_tabula_sapiens()
        elif self.dataset_name == 'liu_2021':
            dataset_list = self.load_liu_2021()
        elif self.dataset_name == 'crosby_2022':
            dataset_list = self.load_crosby_2022()
        elif self.dataset_name == 'petljak_2019':
            dataset_list = self.load_petljak_2019()
        else:
            raise NotImplementedError(f"Dataset {self.dataset_name} not implemented")
        return dataset_list

    def load_petljak_2019(self) -> ExpressionDataset:
        """Load the Petljak 2019 cell line dataset
        ### Returns
        cell_lines : ExpressionDataset
            Expression dataset for cell lines
        """
        data_dir = '/cellar/users/zkoch/dream/data/petljak_2019'
        
        # read in mapping of COSMIC id to cell line name
        id_to_cell_line = pd.read_excel(os.path.join(data_dir, 'TableS1E.xlsx'))
        # set second row as column names
        id_to_cell_line.columns = id_to_cell_line.iloc[1]
        id_to_cell_line = id_to_cell_line.iloc[3:]
        # drop first column
        id_to_cell_line = id_to_cell_line.iloc[:, 1:]
        # set COSMIC id as index
        id_to_cell_line.set_index('COSMIC identifier', inplace=True, drop=True)
        # set as string
        id_to_cell_line.index = id_to_cell_line.index.astype(str)
        mapper = id_to_cell_line['Sample Name'].to_dict()
        # set index to be Sample Name
        id_to_cell_line.set_index('Sample Name', inplace=True, drop=True)
        # drop nan columns
        id_to_cell_line = id_to_cell_line.loc[:, id_to_cell_line.columns.notna()]
        
        # read expression data
        expr_fn = os.path.join(data_dir, 'Cell_line_RMA_proc_basalExp.txt')
        expr_df = pd.read_csv(expr_fn, sep='\t').T
        # set first row to be column names
        expr_df.columns = expr_df.iloc[0]
        expr_df = expr_df.iloc[2:]
        # strip 'DATA.' from start of cell line names index
        expr_df.index = expr_df.index.str.strip('DATA.')
        expr_df.index = expr_df.index.astype(str)
        # map COSMIC IDs to cell line names
        expr_df.index = expr_df.index.map(mapper)
        # drop rows with missing values in index
        expr_df = expr_df.loc[expr_df.index.notna()]
        # drop columns with missing values
        expr_df = expr_df.loc[:, expr_df.columns.notna()]

        # create ExpressionDataset object
        cell_lines = ExpressionDataset(
            expression_df=expr_df,
            species="symbol",
            metadata_df=id_to_cell_line,
            dataset="petljak_2019"
        )
        return cell_lines

    def load_crosby_2022(self) -> ExpressionDataset:
        """Load the Crosby 2022 dataset from GSE215974_Processed_data_FE1_hydrazine.xlsx
        
        Returns:
            ExpressionDataset: An ExpressionDataset object containing the Crosby 2022 data
        """
        # path to dataset
        fn = "/cellar/users/zkoch/dream/data/crosby_2022/GSE215974_Processed_data_FE1_hydrazine.xlsx"
        
        # read excel file
        expression_df = pd.read_excel(fn)
        
        # set gene names as index
        expression_df.set_index(expression_df.columns[0], inplace=True)
        
        # transpose so samples are rows
        expression_df = expression_df.T
        
        # create metadata from sample names
        metadata_df = pd.DataFrame(index=expression_df.index)
        # treatment is the number after the last underscore
        metadata_df['treatment'] = metadata_df.index.str.split('145_').str[-1].str[:-1]
        metadata_df['washout'] = metadata_df.index.str.contains('plus20hr')
        
    
        # create ExpressionDataset object
        dataset = ExpressionDataset(
            expression_df=expression_df,
            species="mouse", 
            metadata_df=metadata_df,
            dataset="crosby_2022"
        )
        
        return dataset

    def load_liu_2021(self) -> ExpressionDataset:
        """
        Load the Liu 2021 dataset from GSE161789_raw_counts_and_FPKM_for_RNA-seq.xlsx
        
        Returns:
            ExpressionDataset: An ExpressionDataset object containing the Liu 2021 data
        """
        # path to the dataset
        fn = "/cellar/users/zkoch/dream/data/liu_2021/GSE161789_raw_counts_and_FPKM_for_RNA-seq.xlsx"
        
        # read the excel file
        expression_df = pd.read_excel(fn)
        
        # create metadata dataframe with sample information
        metadata_df = pd.DataFrame(index=expression_df.columns[1:17])  # assuming first column is gene IDs
        metadata_df['treatment'] = metadata_df.index.str.split('_').str[1]
        metadata_df['value_type'] = metadata_df.index.str.split('_').str[-1]
        metadata_df['replicate'] = metadata_df.index.str.split('_').str[2]

        # transpore and drop last 9 columns
        expression_df.set_index(expression_df.columns[0], inplace=True)
        expression_df = expression_df.iloc[:, :-9]
        expression_df = expression_df.T
 
        
        # create the ExpressionDataset object
        # assuming human species, modify if different
        dataset = ExpressionDataset(
            expression_df=expression_df,  # set first column as index
            species="human",
            metadata_df=metadata_df,
            dataset="liu_2021"
        )
        
        return dataset


    def load_tabula_sapiens(self) -> ScExpressionDataset:
        """Load the Tabula Sapiens dataset
        ### Returns
        adata : ScExpressionDataset
            Expression dataset for Tabula Sapiens
        """
        data_dir = '/cellar/users/zkoch/dream/data/tabula_sapiens'
        adata = sc.read_h5ad(os.path.join(data_dir, 'tabula_sapiens.h5ad'))

        # subset to 10x 3'
        # find rows where assay is 10x 3' v3
        assay_mask = adata.obs['assay'] == '10x 3\' v3'
        adata = adata[assay_mask, :]

        # tissue (was tissue_in_publication) is a more general tissue name than 'tissue_specific' (was tissue), for example 'Eye' includes tissue_specifc: 'eye', 'conjunctiva', 'lacrimal gland', 'cornea', 'retinal neural layer', 'sclera'
        adata.obs.rename(columns = {'tissue': 'tissue_specific', 'tissue_in_publication': 'tissue'}, inplace = True)

        # convert from category to string
        adata.obs['tissue'] = adata.obs['tissue'].astype(str)
        adata.obs['tissue_specific'] = adata.obs['tissue_specific'].astype(str)

        # get age as float
        adata.obs['age'] = adata.obs['development_stage'].str.split('-').str[0].astype(float)

        # broader cell types (cell_type) will allow more mutations to be detected than more specific cell types (cell_ontology_class) because there are less branches for the mutation to happen "before"
        # based on this, what do we want?
        # we want to distinguish trends within tissues we want as many cell types per tissue as possible, but more cell types will discard more mutations
        # but this is probably better than having less resolution to compare within tissues
        # plus scRNA-seq mutations are already generally too high
        # => use granular cell types ('cell_ontology_class')
        adata.obs.rename({'cell_type': 'broader_cell_type', 'cell_ontology_class': 'cell_type'}, inplace = True)
        adata.obs['cell_type'] = adata.obs['cell_type'].astype(str)

        # set gene names as index
        """adata.var.reset_index(inplace = True)
        adata.var.set_index('feature_name', inplace = True)"""
        
        # replace any ' ' or '/' with '_' in the cell_type column to make work with SComatic
        adata.obs['cell_type'] = adata.obs['cell_type'].str.replace(' ', '_')
        adata.obs['cell_type'] = adata.obs['cell_type'].str.replace('/', '_')
        
        ts = ScExpressionDataset(
            adata=adata,
            # ensemble gene ids
            gene_species="human",
            dataset="tabula_sapiens"
            )
        # calculate mean dream expression to populate dream_expression
        ts.get_n_genes_and_counts()
        ts.get_dream_gene_expression()
        # read in ssgsea results
        ssgsea_dream_expression = pd.read_parquet(
            '/cellar/users/zkoch/dream/data/tabula_sapiens/ssgsea_dream_expression.parquet'
            )
        ts.dream_expression.obs = ssgsea_dream_expression
        
        ############# Read in pseudo-bulk mutations ############   
        print("Reading in mutations")
        # find all mutation files
        mutation_fns = glob.glob('/cellar/users/zkoch/dream/data/tabula_sapiens/output_dir/*/Step4_VariantCalling/*.calling.step2.pass.tsv')
        # read in each file
        mutation_dfs = []
        for fn in mutation_fns:
            try:
                df = pd.read_csv(fn, sep='\t', skiprows=27, header = 1)
            except:
                # skip if file is empty
                continue
            
            # get sample id
            sample_id = fn.split('/')[-1].split('.')[0]
            # add to df
            df['mut_sample_id'] = sample_id
            # append to list
            mutation_dfs.append(df)
        mutation_df = pd.concat(mutation_dfs)
        del mutation_dfs

        # get count of mutations in each cell type in each donor
        count_by_type_by_donor = mutation_df.groupby(['mut_sample_id', 'Cell_types']).size().sort_values(ascending = False).reset_index().rename(columns = {0: 'mutation_burden', 'Cell_types': 'cell_type'})
        # to not lose cell ids when we merge 
        ts.dream_expression.obs.reset_index(inplace = True)
        ts.dream_expression.obs.rename(columns = {'index': 'cell_id'}, inplace = True)
        # merge with ts.dream_expression.obs
        ts.dream_expression.obs = ts.dream_expression.obs.merge(
                count_by_type_by_donor, left_on = ['donor_id', 'cell_type'],
                right_on = ['mut_sample_id', 'cell_type'], how = 'left'
                )   
        # drop mut_sample_id
        ts.dream_expression.obs.drop(columns = ['mut_sample_id'], inplace = True)
        ts.dream_expression.obs.set_index('cell_id', inplace = True)

        ########## Read in sc genotypes ############
        print("Reading in sc genotypes")
        sc_genotype_fns = glob.glob('/cellar/users/zkoch/dream/data/tabula_sapiens/output_dir/*/SingleCellAlleles/*single_cell_genotype.tsv')
        # read in each
        sc_genotype_dfs = []
        for fn in sc_genotype_fns:
            df = pd.read_csv(fn, sep='\t')
            # get sample id
            sample_id = fn.split('/')[-3]
            # add to df
            df['mut_sample_id'] = sample_id
            # append to list
            sc_genotype_dfs.append(df)
        sc_genotype_df = pd.concat(sc_genotype_dfs)
        del sc_genotype_dfs
        # remove chrY and chrX
        sc_genotype_df = sc_genotype_df.loc[(sc_genotype_df['#CHROM'] != 'chrY') & (sc_genotype_df['#CHROM'] != 'chrX')]
        # groupby the cell barcode to get the mutation burden per cell
        sc_burden = sc_genotype_df.groupby(['mut_sample_id', 'CB']).size().reset_index().rename(columns = {0: 'sc_mutation_burden'})
        
        # get CB of dream expression
        ts.dream_expression.obs['CB'] = ts.dream_expression.obs.index.str.split('_').str[0]
        # merge with dream expression on CB
        ts.dream_expression.obs = ts.dream_expression.obs.merge(
            sc_burden, left_on = ['CB', 'donor_id'], right_on = ['CB', 'mut_sample_id'], how = 'left'
            )
        ts.dream_expression.obs.drop(columns = ['mut_sample_id'], inplace = True)

        ############ SC mutation burden corroberated by bulk ############
        # subset to mutations found in sc that were also seen in the bulk level for that cell type
        dropped_sc_genotype_df = sc_genotype_df.query(
            "Cell_type_expected == Cell_type_observed and ALT_expected == Base_observed"
            )
        sc_burden = dropped_sc_genotype_df.groupby(['mut_sample_id', 'CB']).size().reset_index().rename(columns = {0: 'sc_mutation_burden_in_bulk'})
        # merge with dream expression on CB
        ts.dream_expression.obs = ts.dream_expression.obs.merge(
            sc_burden, left_on = ['CB', 'donor_id'], right_on = ['CB', 'mut_sample_id'], how = 'left'
            )
        # drop mut_sample_id
        ts.dream_expression.obs.drop(columns = ['mut_sample_id'], inplace = True)

        ############ Read in sc callable sites ################
        print("Reading in sc callable sites (SitesPerCell)")
        sc_callable_sites_fns = glob.glob('/cellar/users/zkoch/dream/data/tabula_sapiens/output_dir/*/UniqueCellCallableSites/*.SitesPerCell.tsv')
        # read in each file
        sc_callable_sites_dfs = []
        for fn in sc_callable_sites_fns:
            df = pd.read_csv(fn, sep=',')
            # get sample id
            sample_id = fn.split('/')[-1].split('.')[0]
            # add to df
            df['mut_sample_id'] = sample_id
            # append to list
            sc_callable_sites_dfs.append(df)
        sc_callable_sites_df = pd.concat(sc_callable_sites_dfs)
        del sc_callable_sites_dfs
        sc_callable_sites_df.reset_index(inplace=True, drop = True)
        # drop rows with duplicated Donor ID and CB
        sc_callable_sites_df = sc_callable_sites_df.drop_duplicates(subset = ['mut_sample_id', 'CB'])
        # merge with dream expression
        ts.dream_expression.obs = ts.dream_expression.obs.merge(
            sc_callable_sites_df[['CB', 'mut_sample_id', 'SitesPerCell']], left_on = ['CB', 'donor_id'], right_on = ['CB', 'mut_sample_id'], how = 'left'
            )
        # drop mut_sample_id
        ts.dream_expression.obs.drop(columns = ['mut_sample_id'], inplace = True)

        ########## Read in callable sites ############
        print("Reading in callable sites")

        callable_sites_fns = glob.glob('/cellar/users/zkoch/dream/data/tabula_sapiens/output_dir/*/CellTypeCallableSites/*coverage_cell_count.report.tsv')
        # read in each file
        callable_site_dfs = []
        for fn in callable_sites_fns:
            df = pd.read_csv(fn, sep='\t')
            # get sample id
            sample_id = fn.split('/')[-1].split('.')[0]
            # add to df
            df['mut_sample_id'] = sample_id
            # append to list
            callable_site_dfs.append(df)
        callable_site_df_ungrouped = pd.concat(callable_site_dfs)
        del callable_site_dfs

        # condense this into the number of callable sites per cell type per sample
        ##INFO=Dp,Description=Depth of coverage (reads) in the cell type supporting the variant
        ##INFO=Nc,Description=Number of distinct cells found in the cell type with the mutation
        callable_site_df = callable_site_df_ungrouped.groupby(['Cell_types', 'mut_sample_id'])['DP'].sum().reset_index()
        # rename to callable_sites
        callable_site_df.rename(columns = {'DP': 'callable_sites', 'Cell_types': 'cell_type', 'mut_sample_id': 'donor_id'}, inplace=True)
        # drop rows with duplicated Donor ID and Subclass
        callable_site_df = callable_site_df[~callable_site_df.duplicated(subset = ['donor_id', 'cell_type'])]
        # merge with dream expression
        ts.dream_expression.obs = ts.dream_expression.obs.merge(
            callable_site_df, left_on = ['donor_id', 'cell_type'],
            right_on = ['donor_id', 'cell_type'], how = 'left'
            )

        # scale mutation rates 
        ts.dream_expression.obs['mutation_burden_per_kb'] = ts.dream_expression.obs['mutation_burden'] / (ts.dream_expression.obs['callable_sites']/1000)
        # sc mutations scaled by sc callable sites
        ts.dream_expression.obs['sc_mutation_burden_per_kb'] = ts.dream_expression.obs['sc_mutation_burden'] / (ts.dream_expression.obs['SitesPerCell']/1000)
        # sc mutations scaled by bulk callable sites
        ts.dream_expression.obs['sc_mutation_burden_per_bulk_kb'] = ts.dream_expression.obs['sc_mutation_burden'] / (ts.dream_expression.obs['callable_sites']/1000)
        # sc mutations also found in bulk, scaled by sc callable sites
        ts.dream_expression.obs['sc_mutation_burden_in_bulk_per_kb'] = ts.dream_expression.obs['sc_mutation_burden_in_bulk'] / (ts.dream_expression.obs['SitesPerCell']/1000)
        # sc mutations also found in bulk, scaled by bulk callable sites
        ts.dream_expression.obs['sc_mutation_burden_in_bulk_per_bulk_kb'] = ts.dream_expression.obs['sc_mutation_burden_in_bulk'] / (ts.dream_expression.obs['callable_sites']/1000)
        # sc mutations also found in bulk, scaled to per genome by sc callable sites
        ts.dream_expression.obs['sc_mutation_burden_in_bulk_per_genome'] = 3e9 * (ts.dream_expression.obs['sc_mutation_burden_in_bulk'] / (ts.dream_expression.obs['SitesPerCell']))
        # sc mutations scaled to per genome by sc callable sites
        ts.dream_expression.obs['sc_mutation_burden_per_genome'] = 3e9 * (ts.dream_expression.obs['sc_mutation_burden'] / (ts.dream_expression.obs['SitesPerCell']))
        # mutation burden scaled to per genome
        ts.dream_expression.obs['mutation_burden_per_genome'] = 3e9 * (ts.dream_expression.obs['mutation_burden'] / (ts.dream_expression.obs['callable_sites']))


        ts.dream_expression.obs['sc_mutation_burden_in_bulk_per_genome_per_year'] = ts.dream_expression.obs['sc_mutation_burden_in_bulk_per_genome'] / ts.dream_expression.obs['age']

        
        return ts


    def load_SEA_AD(
        self,
        read_whole_dataset: bool = False,
        start_cell_num: int = -1,
        end_cell_num: int = -1
        ) -> ScExpressionDataset:
        """Load the SEA-AD dataset
        ### Returns
        adata : ScExpressionDataset
            Expression dataset for SEA-AD
        ### Parameters
        read_whole_dataset : bool
            Whether to read the whole dataset
        start_cell_num : int
            The number of cells to start at
        end_cell_num : int
            The number of cells to end at
        """
        data_dir = '/cellar/users/zkoch/dream/data/SEA-AD/aws_snRNA_seq'
        if read_whole_dataset:
            print("Loading whole dataset")
            adata = sc.read_h5ad(
                os.path.join(data_dir, 'sc_sea_ad_preprocessed_ssgsea_adata.h5ad')
                )
            if start_cell_num != -1 and end_cell_num != -1:
                adata = utils.read_subset_h5ad(
                    os.path.join(data_dir, 'sc_sea_ad_preprocessed_ssgsea_adata.h5ad'),
                    start_cell_num,
                    end_cell_num
                    )
            sea_ad = ScExpressionDataset(
                adata=adata,
                gene_species="symbol",
                dataset="SEA-AD"
                )
        else:
            print("WARNING: loading only DREAM and pretending its both to save time")
            dream = sc.read_h5ad(os.path.join(data_dir, 'sc_sea_ad_preprocessed_ssgsea_dream_expression.h5ad'))
            sea_ad = ScExpressionDataset(
                adata=dream,
                gene_species="symbol",
                dataset="SEA-AD"
                )
        
        # preporcessed by running on each cell and aggregating the results
        #sea_ad.adata = sea_ad.adata[cell_num_start: cell_num_start + 10000]
        # do the ssgsea
        #sea_ad.get_n_genes_and_counts()
        #sea_ad.get_dream_gene_expression()
        #sea_ad.dream_enrichment_ssgsea()
        # scaled by eq = f'{col_name} ~ n_counts * n_genes'
        
        sea_ad = ScExpressionDataset(
            adata=dream,
            gene_species="symbol",
            dataset="SEA-AD"
            )
        sea_ad._read_and_convert_dream()
        sea_ad.dream_expression = dream
        
        # rename subclasses
        sea_ad.adata.obs['Subclass'] = sea_ad.adata.obs['Subclass'].str.replace(' ', '_')
        sea_ad.adata.obs['Subclass'] = sea_ad.adata.obs['Subclass'].str.replace('/', '_')
        sea_ad.dream_expression.obs['Subclass'] = sea_ad.dream_expression.obs['Subclass'].str.replace(' ', '_')
        sea_ad.dream_expression.obs['Subclass'] = sea_ad.dream_expression.obs['Subclass'].str.replace('/', '_')
        
        
        ############ Read in mutations ############
        print("Reading in mutations")
        # find all mutation files
        mutation_fns = glob.glob('/cellar/users/zkoch/dream/data/SEA-AD/synapse/bams/dlpfc/output_dir/*/Step4_VariantCalling/*.calling.step2.pass.tsv')
        # read in each file
        mutation_dfs = []
        for fn in mutation_fns:
            try:
                df = pd.read_csv(fn, sep='\t', skiprows=27, header = 1)
            except:
                # skip if file is empty
                continue
            
            # get sample id
            sample_id = fn.split('/')[-1].split('.')[0]
            # add to df
            df['mut_sample_id'] = sample_id
            # append to list
            mutation_dfs.append(df)
        mutation_df = pd.concat(mutation_dfs)
        del mutation_dfs

        # read in all the manifest files that synapse output, to match sample ids to donor ids
        manifest_fns = glob.glob('/cellar/users/zkoch/dream/data/SEA-AD/synapse/bams/dlpfc/manifest*')
        manifest_dfs = []
        for fn in manifest_fns:
            df = pd.read_csv(fn, sep=',')
            manifest_dfs.append(df)
        manifest_df = pd.concat(manifest_dfs)
        # convert sample id to donor id in mutation df
        manifest_df['mut_sample_id'] = manifest_df['name'].str.split('-').str[0]
        manifest_df['mut_sample_id'] = manifest_df['mut_sample_id'].str.split('_').str[0]
        manifest_df.set_index('mut_sample_id', inplace=True)
        mut_sample_id_to_donor_id = manifest_df['individualID']
        # drop dullicated index values (these are cases where )
        mut_sample_id_to_donor_id = mut_sample_id_to_donor_id[~mut_sample_id_to_donor_id.index.duplicated()]
        # map to donor id
        mutation_df['Donor ID'] = mutation_df['mut_sample_id'].map(mut_sample_id_to_donor_id)
        # get count of mutations in each cell type in each donor
        count_by_type_by_donor = mutation_df.groupby(['Donor ID', 'Cell_types']).size().reset_index().rename(columns = {0: 'mutation_burden', 'Cell_types': 'Subclass'})
        # merge with dream expression
        sea_ad.dream_expression.obs = sea_ad.dream_expression.obs.merge(
            count_by_type_by_donor, left_on = ['Donor ID', 'Subclass'],
            right_on = ['Donor ID', 'Subclass'], how = 'left'
            )
        
        ########## Read in callable sites ############
        print("Reading in callable sites")
        callable_sites_fns = glob.glob('/cellar/users/zkoch/dream/data/SEA-AD/synapse/bams/dlpfc/output_dir/*/CellTypeCallableSites/*coverage_cell_count.report.tsv')
        # read in each file
        callable_site_dfs = []
        for fn in callable_sites_fns:
            df = pd.read_csv(fn, sep='\t')
            # get sample id
            sample_id = fn.split('/')[-1].split('.')[0]
            # add to df
            df['mut_sample_id'] = sample_id
            # append to list
            callable_site_dfs.append(df)
        callable_site_df_ungrouped = pd.concat(callable_site_dfs)
        del callable_site_dfs

        # condense this into the number of callable sites per cell type per sample
        ##INFO=Dp,Description=Depth of coverage (reads) in the cell type supporting the variant
        ##INFO=Nc,Description=Number of distinct cells found in the cell type with the mutation
        callable_site_df = callable_site_df_ungrouped.groupby(['Cell_types', 'mut_sample_id'])['DP'].sum().reset_index()
        # rename to callable_sites
        callable_site_df.rename(columns = {'DP': 'callable_sites', 'Cell_types': 'Subclass'}, inplace=True)
        # map to donor id
        callable_site_df['Donor ID'] = callable_site_df['mut_sample_id'].map(mut_sample_id_to_donor_id)
        # drop rows with duplicated Donor ID and Subclass
        callable_site_df = callable_site_df[~callable_site_df.duplicated(subset = ['Donor ID', 'Subclass'])]
        # merge with dream expression
        sea_ad.dream_expression.obs = sea_ad.dream_expression.obs.merge(
            callable_site_df, left_on = ['Donor ID', 'Subclass'],
            right_on = ['Donor ID', 'Subclass'], how = 'left'
            )
        
        ########## Read in sc genotypes ############
        print("Reading in sc genotypes")
        sc_genotype_fns = glob.glob('/cellar/users/zkoch/dream/data/SEA-AD/synapse/bams/dlpfc/output_dir/*/SingleCellAlleles/*single_cell_genotype.tsv')
        # read in each
        sc_genotype_dfs = []
        for fn in sc_genotype_fns:
            df = pd.read_csv(fn, sep='\t')
            # get sample id
            sample_id = fn.split('/')[-3]
            # add to df
            df['mut_sample_id'] = sample_id
            # append to list
            sc_genotype_dfs.append(df)
        sc_genotype_df = pd.concat(sc_genotype_dfs)
        del sc_genotype_dfs
        # map to donor id
        sc_genotype_df['Donor ID'] = sc_genotype_df['mut_sample_id'].map(mut_sample_id_to_donor_id)


        # groupby the cell barcode to get the mutation burden per cell
        sc_burden = sc_genotype_df.groupby(['Donor ID', 'CB']).size().reset_index().rename(columns = {0: 'sc_mutation_burden'})
        # merge with dream expression on CB
        sea_ad.dream_expression.obs['CB'] = sea_ad.dream_expression.obs['sample_id'].str.split('-').str[0]
        sea_ad.dream_expression.obs = sea_ad.dream_expression.obs.merge(
            sc_burden, left_on = ['CB', 'Donor ID'], right_on = ['CB', 'Donor ID'], how = 'left'
            )

        # this decides if we only want to consider mutations that are seen in a cell type at the pseudo-bulk level
        ########################
        # select rows where we see the expected mutation call in the expected cell type
        dropped_sc_genotype_df = sc_genotype_df.query("Cell_type_expected == Cell_type_observed and ALT_expected == Base_observed")
        sc_burden = dropped_sc_genotype_df.groupby(['Donor ID', 'CB']).size().reset_index().rename(columns = {0: 'sc_mutation_burden_in_bulk'})
        # merge with dream expression on CB
        sea_ad.dream_expression.obs['CB'] = sea_ad.dream_expression.obs['sample_id'].str.split('-').str[0]
        sea_ad.dream_expression.obs = sea_ad.dream_expression.obs.merge(
            sc_burden, left_on = ['CB', 'Donor ID'], right_on = ['CB', 'Donor ID'], how = 'left'
            )
        
        ############ Read in sc callable sites ################
        print("Reading in sc callable sites")
        sc_callable_sites_fns = glob.glob('/cellar/users/zkoch/dream/data/SEA-AD/synapse/bams/dlpfc/output_dir/*/UniqueCellCallableSites/*.SitesPerCell.tsv')
        # read in each file
        sc_callable_sites_dfs = []
        for fn in sc_callable_sites_fns:
            df = pd.read_csv(fn, sep=',')
            # get sample id
            sample_id = fn.split('/')[-1].split('.')[0]
            # add to df
            df['mut_sample_id'] = sample_id
            # append to list
            sc_callable_sites_dfs.append(df)
        sc_callable_sites_df = pd.concat(sc_callable_sites_dfs)
        del sc_callable_sites_dfs
        sc_callable_sites_df.reset_index(inplace=True, drop = True)
        # map to donor id
        sc_callable_sites_df['Donor ID'] = sc_callable_sites_df['mut_sample_id'].map(mut_sample_id_to_donor_id)
        # drop rows with duplicated Donor ID and CB
        sc_callable_sites_df = sc_callable_sites_df.drop_duplicates(subset = ['Donor ID', 'CB'])
        # merge with dream expression
        sea_ad.dream_expression.obs = sea_ad.dream_expression.obs.merge(
            sc_callable_sites_df[['CB', 'Donor ID', 'SitesPerCell']], left_on = ['CB', 'Donor ID'], right_on = ['CB', 'Donor ID'], how = 'left'
            )
        
        # scale to mutation burden per kb
        sea_ad.dream_expression.obs['mutation_burden_per_kb'] = sea_ad.dream_expression.obs['mutation_burden'] / (sea_ad.dream_expression.obs['callable_sites']/1000)
        # sc mutations scaled by sc callable sites
        sea_ad.dream_expression.obs['sc_mutation_burden_per_kb'] = sea_ad.dream_expression.obs['sc_mutation_burden'] / (sea_ad.dream_expression.obs['SitesPerCell']/1000)
        # sc mutations scaled by bulk callable sites
        sea_ad.dream_expression.obs['sc_mutation_burden_per_bulk_kb'] = sea_ad.dream_expression.obs['sc_mutation_burden'] / (sea_ad.dream_expression.obs['callable_sites']/1000)
        # sc mutations also found in bulk, scaled by sc callable sites
        sea_ad.dream_expression.obs['sc_mutation_burden_in_bulk_per_kb'] = sea_ad.dream_expression.obs['sc_mutation_burden_in_bulk'] / (sea_ad.dream_expression.obs['SitesPerCell']/1000)
        # sc mutations also found in bulk, scaled by bulk callable sites
        sea_ad.dream_expression.obs['sc_mutation_burden_in_bulk_per_bulk_kb'] = sea_ad.dream_expression.obs['sc_mutation_burden_in_bulk'] / (sea_ad.dream_expression.obs['callable_sites']/1000)
        # sc mutations also found in bulk, scaled to per genome by sc callable sites
        sea_ad.dream_expression.obs['sc_mutation_burden_in_bulk_per_genome'] = 3e9 * (sea_ad.dream_expression.obs['sc_mutation_burden_in_bulk'] / (sea_ad.dream_expression.obs['SitesPerCell']))
        return sea_ad
        

    def load_synapse_MIT_ROSMAP_Multiomics_smartseq(self) -> ScExpressionDataset:
        """Load the synapse MIT ROSMAP Multiomics smartseq dataset with the fusion data (this is the expression data from mathys_2023 and fusion data from dileep_2023)
        ### Returns
        adata : ScExpressionDataset
            Expression dataset for synapse MIT ROSMAP Multiomics
        """
        """fusions_smart = pd.read_csv("/cellar/users/zkoch/dream/data/synapse_MIT_ROSMAP_Multiomics/dileep_2023/human_smart-seq2_fusions.maxsensitivity_prediction_numbers.meta.counts.tsv", sep="\t", header=0)
        human_smart_expr_colnames_fn = "/cellar/users/zkoch/dream/data/synapse_MIT_ROSMAP_Multiomics/dileep_2023/human_smart-seq2.expr_colnames.tsv"
        human_smart_expr_symbol_rownames_fn = "/cellar/users/zkoch/dream/data/synapse_MIT_ROSMAP_Multiomics/dileep_2023/human_smart-seq2.expr_symbol_rownames.tsv"
        human_smart_expr_colnames = pd.read_csv(human_smart_expr_colnames_fn, sep="\t", header=None)
        human_smart_expr_symbol_rownames = pd.read_csv(human_smart_expr_symbol_rownames_fn, sep="\t", header=None)
        expr_smart = sc.read_mtx("/cellar/users/zkoch/dream/data/synapse_MIT_ROSMAP_Multiomics/dileep_2023/human_smart-seq2.expr.mtx").T

        # set gene names
        expr_smart.var = human_smart_expr_symbol_rownames#.values
        expr_smart.var.rename(columns = {0:'gene_symbol'}, inplace = True)
        expr_smart.var.set_index('gene_symbol', inplace = True, drop = True)
        # set cell names
        expr_smart.obs = human_smart_expr_colnames#.values
        expr_smart.obs.rename(columns = {0:'cell_name'}, inplace = True)
        expr_smart.obs.set_index('cell_name', inplace = True, drop = True)

        # merge and drop rows that are duplicated
        merged_obs = expr_smart.obs.merge(fusions_smart, left_index = True, right_on = 'sample', how = 'inner').drop_duplicates(subset = ['subject', 'sample', 'count', 'celltype'])
        # subset X to only include the cells that are in the merged obs
        expr_smart = expr_smart[merged_obs.index]
        # set the obs to the merged obs
        expr_smart.obs = merged_obs
        expr_smart.obs.reset_index(inplace = True, drop = True)

        # make genes unique by dropping duplicates
        duplicated_vars = expr_smart.var.index.duplicated()
        expr_smart = expr_smart[:, ~duplicated_vars]

        # rename count to fusion_count
        expr_smart.obs.rename(columns = {'count':'fusion_count'}, inplace = True)

        # map mit ids to rosmap ids
        mit_to_rosmap_mapping = pd.read_csv("/cellar/users/zkoch/dream/data/synapse_MIT_ROSMAP_Multiomics/MIT_ROSMAP_Multiomics_individual_metadata.csv", sep=",", header=0)
        mit_to_rosmap_mapping.set_index("subject", inplace=True)
        # drop duplicate index 
        mit_to_rosmap_mapping = mit_to_rosmap_mapping[~mit_to_rosmap_mapping.index.duplicated(keep='first')]
        expr_smart.obs['individualID'] = ''
        expr_smart.obs['individualID'] = expr_smart.obs['subject'].map(mit_to_rosmap_mapping['individualID'])

        # merge with rosmap metadata
        harmonized_metadata = pd.read_csv("/cellar/users/zkoch/dream/data/synapse_rna_seq_harmonization/RNAseq_Harmonization_ROSMAP_combined_metadata.csv")
        # TODO: might be losing longitudinal data here
        harmonized_metadata.drop_duplicates(subset="individualID", keep="first", inplace=True)
        expr_smart.obs = expr_smart.obs.merge(
            harmonized_metadata, left_on="individualID", right_on="individualID", how="left"
            )
        # drop _x columns and rename _y columns
        expr_smart.obs = expr_smart.obs.loc[:,~expr_smart.obs.columns.str.contains('_x')]
        expr_smart.obs.rename(columns = {col: col.replace('_y','') for col in expr_smart.obs.columns}, inplace = True)
            
        synapse = ScExpressionDataset(
            adata=expr_smart,
            gene_species="symbol",
            dataset="synapse_MIT_ROSMAP_Multiomics_smartseq"
            )
        # then I ran 
        # sc_rosmap_smart.pre_process()
        # sc_rosmap_smart.get_dream_gene_expression()
        # sc_rosmap_smart.dream_enrichment_ssgsea()
        """
        with open("/cellar/users/zkoch/dream/data/synapse_MIT_ROSMAP_Multiomics/sc_rosmap_smart_preprocessed_ssgsea.pkl", "rb") as f:
            synapse = pickle.load(f)
        synapse.dream_expression.obs['age_first_ad_dx'] = synapse.dream_expression.obs['age_first_ad_dx'].replace('90+',90).astype(float)
        synapse.dream_expression.obs['age_death'] = synapse.dream_expression.obs['age_death'].replace('90+',90).astype(float)
        synapse.dream_expression.obs['age_at_visit_max'] = synapse.dream_expression.obs['age_at_visit_max'].replace('90+',90).astype(float)
        return synapse
    
    def load_synapse_MIT_ROSMAP_Multiomics_10x(self) -> ScExpressionDataset:
        """Load the synapse MIT ROSMAP Multiomics 10x dataset with the fusion data (this is the expression data from mathys_2023 and fusion data from dileep_2023)
        ### Returns
        adata : ScExpressionDataset
            Expression dataset for synapse MIT ROSMAP Multiomics
        """
        # read in 10x expression and fusions data
        """print("reading in 10x expression and fusions data, takes ~10 minutes")
        fusions_10x = pd.read_csv("/cellar/users/zkoch/dream/data/synapse_MIT_ROSMAP_Multiomics/dileep_2023/human_10x.aggregated.fusion_predictions.processed_counts.tsv", sep="\t", header=0)
        expr_10x = sc.read_h5ad("/cellar/users/zkoch/dream/data/synapse_MIT_ROSMAP_Multiomics/mathys_2023/PFC427_raw_data.h5ad")
        
        # make barcodes match
        expr_10x.obs['num'] = expr_10x.obs['bc'].str.split('-').str[-1]
        expr_10x.obs['bc_no_num'] = expr_10x.obs['bc'].str.split('-').str[0]
        expr_10x.obs['barcode'] = expr_10x.obs['batch'].astype(str) + '-' + expr_10x.obs['num'] + ':' + expr_10x.obs['bc_no_num']
        # merge fusions with 10x
        expr_10x.obs = expr_10x.obs.merge(
            fusions_10x[['barcode', 'subject', 'celltype', 'cell_type_high_res', 'count', 'age_death', 'msex', 'pmi', 'cogdx.ad', 'braak56' ]],
            left_on="barcode", right_on="barcode", how="left"
            )
        expr_10x.obs.rename(columns={"count": "fusion_count"}, inplace=True)

        # map mit ids to rosmap ids
        mit_to_rosmap_mapping = pd.read_csv("/cellar/users/zkoch/dream/data/synapse_MIT_ROSMAP_Multiomics/MIT_ROSMAP_Multiomics_individual_metadata.csv", sep=",", header=0)
        mit_to_rosmap_mapping.set_index("subject", inplace=True)
        # drop duplicate index 
        mit_to_rosmap_mapping = mit_to_rosmap_mapping[~mit_to_rosmap_mapping.index.duplicated(keep='first')]
        expr_10x.obs['individualID'] = expr_10x.obs['subject'].map(mit_to_rosmap_mapping['individualID'])

        # merge with rosmap metadata
        harmonized_metadata = pd.read_csv("/cellar/users/zkoch/dream/data/synapse_rna_seq_harmonization/RNAseq_Harmonization_ROSMAP_combined_metadata.csv")
        # drop rows with duplicate individualID
        # TODO: might be losing longitudinal data here
        harmonized_metadata.drop_duplicates(subset="individualID", keep="first", inplace=True)
        expr_10x.obs = expr_10x.obs.merge(
            harmonized_metadata, left_on="individualID", right_on="individualID", how="left"
            )
        
        synapse = ScExpressionDataset(
            adata=expr_10x,
            gene_species="symbol",
            dataset="synapse_MIT_ROSMAP_Multiomics_10x"
            )"""
        print("loading synapse_MIT_ROSMAP_Multiomics_10x from pickle")
        with open("/cellar/users/zkoch/dream/data/synapse_MIT_ROSMAP_Multiomics/sc_rosmap_preprocessed_ssgsea.pkl", "rb") as f:
            synapse = pickle.load(f)
        # find all columns that end in _x and drop them
        to_drop = synapse.dream_expression.obs.columns[synapse.dream_expression.obs.columns.str.endswith('_x')]
        synapse.dream_expression.obs.drop(columns = to_drop, inplace = True)
        # rename columns that end in _y to remove the _y
        synapse.dream_expression.obs.rename(columns = lambda x: x.replace('_y',''), inplace = True)
        # replace 90+ with 90 and convert to numeric
        synapse.dream_expression.obs['age_first_ad_dx'] = synapse.dream_expression.obs['age_first_ad_dx'].replace('90+',90).astype(float)
        synapse.dream_expression.obs['age_death'] = synapse.dream_expression.obs['age_death'].replace('90+',90).astype(float)
        synapse.dream_expression.obs['age_at_visit_max'] = synapse.dream_expression.obs['age_at_visit_max'].replace('90+',90).astype(float)
        return synapse


    def load_synapse_rna_seq_harmonization(self) -> ExpressionDataset:
        """Load the synapse rna seq harmonization dataset
        ### Returns
        expression : ExpressionDataset
            Expression dataset for synapse rna seq harmonization
        """
        data_dir = '/cellar/users/zkoch/dream/data/synapse_rna_seq_harmonization'
        expr_fns = glob.glob(os.path.join(data_dir, 'ROSMAP_*_gene_all_counts_matrix_clean.txt'))
        # read in data
        expr_datasets = []
        for fn in expr_fns:
            expr_datasets.append(pd.read_csv(fn, sep='\t', index_col=0))
        expr_datasets = pd.concat(expr_datasets, axis=1)
        expr = expr_datasets.T
        # read metadata 
        metadata = pd.read_csv(
            os.path.join(data_dir, 'RNAseq_Harmonization_ROSMAP_combined_metadata.csv'),
            sep=',', index_col=0
            )
        expression = ExpressionDataset(
            expression_df=expr,
            species="human",
            metadata_df=metadata,
            dataset="synapse_rna_seq_harmonization"
            )
        return expression

    def load_paine_2024(self) -> ExpressionDataset:
        """Load the Paine 2024 dataset
        ### Returns
        expression : ExpressionDataset
            Expression dataset for Paine 2024
        """
        data_dir = "/cellar/users/zkoch/dream/data/paine_2024"
        fn = "GSE231853_family.soft"
        # read metadata from soft file
        metadata_df = self.read_soft(
            os.path.join(data_dir, fn)
            )
        # read expression data
        fn = "GSE231853_counts_TGFB.tsv"
        expr = pd.read_csv(
            os.path.join(data_dir, fn),
            sep = "\t", index_col=0
            ) 
        
        expression = ExpressionDataset(
            expression_df=expr.T,
            species="human",
            metadata_df=metadata_df,
            dataset="paine_2024"
            )
        return expression
    
    
    def load_lu_2022(self) -> ExpressionDataset:
        """Load the Lu 2022 dataset
        ### Returns
        expression : ExpressionDataset
            Expression dataset for Lu 2022
        """
        data_dir = "/cellar/users/zkoch/dream/data/lu_2022"
        # match files ending in .reads
        fns = glob.glob(os.path.join(data_dir, "*.reads"))
        expr_dfs = []
        for fn in fns:
            # read each in as tsv
            df = pd.read_csv(fn, sep = "\t", index_col = 0)
            # append
            expr_dfs.append(df.T)
        # concat one on top of the other, matching on columns
        expr_df = pd.concat(expr_dfs)
        # split index on the first number to get species name
        species = expr_df.index.str.split('.').str[2]
        expr_df['species'] = species
        # map species
        species_dict = {
            'ASM': 'Spiny mouse',
            'beaver': 'American beaver',
            'BMR': 'Blind mole rat',
            'BTR': 'Bushy-tailed woodrat',
            'Cap': 'Capybara',
            'CHH': 'Mouse',#'Chinese hamster',
            'Chin': 'Long-tailed chinchilla',
            'CHIPMUNK': 'Eastern chipmunk',
            'DEGUS': 'Degu',
            'DM': 'Deer mouse',
            'DMR': 'Damaraland mole-rat',
            'Ellobius': 'Transcaucasian mole vole',
            'GH': 'Golden hamster',
            'GP': 'Guinea pig',
            'MUSKR': 'Muskrat',
            'WM': 'Mouse', #'Woodchuck',
            'NUTRIA': 'Nutria',
            'NMR': 'Naked mole rat',
            'paca': 'Lowland paca',
            'rat': 'Rat',
            'RED': 'American red squirrel',
            'GRAY': 'Eastern gray squirrel',
            'star': 'Star-nosed mole',
            'WC': 'Woodchuck', #'Northern short-tailed shrew',
            'EMR': 'Eastern mole',
            'ht2': 'Striped dwarf hamster'
        }
        expr_df['species'] = expr_df['species'].map(species_dict)
        tissue = expr_df.index.str.split('(\d+)', expand=True).get_level_values(2).str.split('.').get_level_values(0).str[0].str.split('_').str[0]
        expr_df['tissue'] = tissue
        # map tissue
        tissue_map = {
            'Lung': 'Lung',
            'Liver':'Liver',
            'LiverR':'Liver',
            'Brain':'Brain',
            'Kidney':'Kidney',
            'Kidny':'Kidney',
            'Skin':'Skin',
            'Heart':'Heart',
        }
        expr_df['tissue'] = expr_df['tissue'].map(tissue_map)
        expr_df.dropna(subset = ['species', 'tissue'], inplace = True)
        metadata_df = expr_df[
            ['species', 'tissue', 'WEIGHT',	'LOG2ratio', 'Intercept', 'MLS']
            ].copy(deep = True)
        # create expression dataset object
        expression = ExpressionDataset(
            expression_df=expr_df,
            species="symbol",
            metadata_df=metadata_df,
            dataset="lu_2022"
            )
        return expression
        
    
    def load_williams_2022(self) -> ExpressionDataset:
        """Load the Williams 2022 dataset
        ### Returns
        mice : ExpressionDataset
            Expression dataset for mice of diff strains
        """
        # Read the processed omics data
        all_data = pd.read_csv(
            "/cellar/users/zkoch/dream/data/williams_2022/aData_S1_AllOmicsandPhenotypeData.csv",
            header=0, index_col = None
            )
        # for phenotype and metadata, 3rd column is index
        all_data.set_index(all_data.columns[2], inplace = True)
        # first 8 rows are metadata
        metadata = all_data.iloc[:8, 2:].copy(deep = True).T
        # reanme IDFemaleAgingColony to OmicsEarTag
        metadata.rename(columns = {'IDFemaleAgingColony': 'OmicsEarTag'}, inplace = True)
        # set the index to OmicsEarTag
        metadata.set_index('OmicsEarTag', inplace = True)
        # drop Order3 column
        metadata.drop(columns = 'Order3', inplace = True)
        # phenotypic data like bodyweight, HDL and amalylase are rows 14 to 53, but also select rows 1 and 7 since that is ear tags
        rows_to_select = [1, 7] + list(range(12, 52))
        pheno_data = all_data.iloc[rows_to_select, 2:].copy(deep = True).T
        pheno_data.rename(columns = {'IDFemaleAgingColony': 'OmicsEarTag'}, inplace = True)
        # set the index to OmicsEarTag and Pheno_AnimalEarTag
        pheno_data.set_index(['OmicsEarTag'], inplace = True)

        # merge the metadata and phenotypic data
        metadata = metadata.merge(pheno_data, left_index = True, right_index = True)
        # drop EarTagCurrent_y and rename EarTagCurrent_x to EarTagCurrent
        metadata.drop(columns = 'EarTagCurrent_y', inplace = True)
        metadata.rename(columns = {'EarTagCurrent_x': 'EarTagCurrent'}, inplace = True)
        
        # select RNA-seq data
        # column +1 column has ensembl gene ids
        all_data.set_index('InjectionOrder2', inplace = True, drop = False)
        rows_to_select = [1, 7] + list(range(32497, 75817))
        expr_df = all_data.iloc[rows_to_select, 2:].copy(deep = True).T
        expr_df.rename(columns = {'Orig_ID': 'OmicsEarTag', 'Name1': 'EarTagCurrent'}, inplace = True)
        expr_df.set_index('OmicsEarTag', inplace = True, drop = True)

        # split on  _ID_ at the beginning of the column names and remove it
        new_cols = expr_df.columns.str.split('_ID_', expand = True).droplevel(0).to_list()
        # reset first one to EarTagCurrent
        new_cols[0] = 'EarTagCurrent'
        expr_df.columns = new_cols
        expression = ExpressionDataset(
            expression_df=expr_df,
            species="mouse",
            metadata_df=metadata,
            dataset="williams_2022",
            )
        return expression
    
    def load_lu_2014(self) -> ExpressionDataset:
        """Load the Lu 2014 dataset
        ### Returns
        expression : ExpressionDataset
            Expression dataset for Lu 2014
        """
        data_dir = "/cellar/users/zkoch/dream/data/lu_2014"
        fn = "GSE53890_family.soft"
        # read metadata from soft file
        metadata_df = self.read_soft(
            os.path.join(data_dir, fn)
            )
        # read expression data from soft
        fn = "GSE53890_series_matrix_ensembl_mapped.txt.gz"
        expr = pd.read_csv(
            os.path.join(data_dir, fn),
            sep = ",", index_col=0
            ) 
        expression = ExpressionDataset(
            expression_df=expr.T,
            species="human",
            metadata_df=metadata_df,
            dataset="lu_2014"
            )
        return expression
    
    def load_hashimoto_2019(self) -> ScExpressionDataset:
        """Load the Hashimoto 2019 dataset
        ### Returns
        adata : ScExpressionDataset
            Expression dataset for AD
        """
        fn = "/cellar/users/zkoch/dream/data/hashimoto_2019/hashimoto_2019.h5ad"
        adata = sc.read_h5ad(fn)
        hashimoto = ScExpressionDataset(
            adata=adata,
            gene_species="human",
            dataset="hashimoto_2019"
            )
        return hashimoto
    
    def load_otero_garcia_2022(self) -> ScExpressionDataset:
        """Load the otero_garcia_2022 dataset
        ### Returns
        adata : ScExpressionDataset
            Expression dataset for AD
        """
        fn = "/cellar/users/zkoch/dream/data/otero-garcia_2022/70170717-45c4-4891-9b14-fb795ecc3d94.h5ad"
        """adata = sc.read_h5ad(fn)
        otero_garcia = ScExpressionDataset(
            adata=adata,
            gene_species="human",
            dataset="otero-garcia_2022"
            )"""
        # pickle the object and save it
        with open("/cellar/users/zkoch/dream/data/otero-garcia_2022/otero_garcia_preproc_w_ssgsea.pkl", "rb") as f:
            otero_garcia = pickle.load(f)
        return otero_garcia
      
    def load_mammalian_methylation_consort(self) -> MethylationDataset:
        """Load the mammalian methylation consortium dataset
        ### Returns
        methylation : MethylationDataset
            Methylation dataset for mammalian methylation consortium
        """
        data_dir = "/cellar/users/zkoch/dream/data/mammalian_methylation_consort/"
        
        # read in cpg manifest
        manifest = pd.read_csv(
            os.path.join(
                data_dir,
                'shorvath/shorvath-MammalianMethylationConsortium-67ea35d/annotations/manifest_HorvathMammalMethylChip40.csv'
                ),
            low_memory=False
            )
        manifest.drop(columns='Unnamed: 0', inplace=True)
        # create mapping from cpg to gene
        cpg_to_gene_map = manifest[
            ['IlmnID', 'Human.Hg38_SYMBOL']
            ]
        cpg_to_gene_map.set_index('IlmnID', inplace=True)
        # read in methylation
        mmc_methyl_df = pd.read_parquet(
            os.path.join(
                data_dir, 'GSE223748_datBetaNormalized.parquet'
                )
            )
        # read in sampleSheet.csv
        sample_sheet = pd.read_csv(
            os.path.join(
                data_dir, 'sampleSheet.csv'
                )
            )
        methylation_dset = MethylationDataset(
            methyl_df=mmc_methyl_df,
            metadata_df=sample_sheet,
            probe_map=cpg_to_gene_map,
            dataset="mammalian_methylation_consort",
            species = 'symbol'
            )
        return methylation_dset
  
    def load_gyenis_2023(
        self,
        kws : dict
        ) -> ExpressionDataset:
        """Load the Gyenis 2023 dataset
        ### Parameters
        Keyword arguments
            - strand: 'forward' or 'reverse'
        ### Returns
        nascent_mice : ExpressionDataset
            Expression dataset for nascent rna-seq
        """ 
        data_dir = "/cellar/users/zkoch/dream/data/gyenis_2023"
        eu_seq_fn = "Eu_20bins.xls"
        eu_seq = pd.read_excel(os.path.join(data_dir, eu_seq_fn), sheet_name=kws['strand'])
        mice = ExpressionDataset(
            expression_df=eu_seq,
            species="symbol",
            metadata_df=None,
            dataset="gyenis_2023"
            )
        return mice
    
    def load_motrpac_2024(self) -> ExpressionDataset:
        """Load the Motrpac 2024 dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        rna_dir = "/cellar/users/zkoch/dream/data/motrpac_2024/transcriptomics/results"
        rna_fns = glob.glob(os.path.join(rna_dir, "t*","transcript-rna-seq","*count.txt"))
        # read in each
        rna_dfs = []
        for fn in rna_fns:
            df = pd.read_csv(fn, sep='\t')
            df.set_index('gene_id', inplace=True)
            df = df.T
            # get tissue
            tissue_name = os.path.basename(fn).split('-')[2].split('_')[0]
            df['tissue'] = tissue_name
            rna_dfs.append(df)
        rna_df = pd.concat(rna_dfs)

        mice = ExpressionDataset(
            expression_df=rna_df,
            species="rat",
            metadata_df=None,
            dataset="motrpac_2024"
            )
        return mice
        
    def load_uxa_2019(self) -> ExpressionDataset:
        """Load the Uxa 2019 dataset
        ### Returns
        treated_cells : ExpressionDataset
            Expression dataset for treated mice
        """
        dataset_path = "/cellar/users/zkoch/dream/data/uxa_2019"
        treated_cells_fn = "supp_table_s2.xlsx"
        if not self.load_expression:
            raise ValueError("load_expression must be True to load Uxa 2019")
        # load data
        uxa_df = pd.read_excel(os.path.join(dataset_path, treated_cells_fn), sheet_name = None)
        # create df of samples and their gene expression data
        expr_by_sample_df = uxa_df['allgenes'].iloc[:, -20:]
        # and ENSG  
        expr_by_sample_df['ENSG'] = uxa_df['allgenes']['ENSG']
        expr_by_sample_df.set_index('ENSG', inplace = True)
        # create expression dataset object
        treated_mice = ExpressionDataset(
            expression_df=expr_by_sample_df.T,
            species="human",
            metadata_df=None,
            dataset="uxa_2019"
            )
        return treated_mice
    
    def load_bujarrabal_dueso_2023_wormRNA(self) -> ExpressionDataset:
        """Load the Bujarrabal Dueso 2023 wormRNA dataset
        ### Returns
        treated_cells : ExpressionDataset
            Expression dataset for treated cells
        """
        dataset_path = "/cellar/users/zkoch/dream/data/bujarrabal_dueso"
        treated_cells_fn = "GSE152235_FPKM_lin-52_worms.csv.gz"
        if not self.load_expression:
            raise ValueError("load_expression must be True to load Bujarrabal Dueso 2023")
        
        treated_cells_expr = pd.read_csv(
            os.path.join(dataset_path, treated_cells_fn),
            index_col=0, sep = ","
            )
        # drop last 11 columns
        treated_cells_expr = treated_cells_expr.iloc[:, :-11]
        # create expression dataset object
        treated_cells = ExpressionDataset(
            expression_df=treated_cells_expr.T,
            species="worm",
            metadata_df=pd.DataFrame(),
            dataset="bujarrabal_dueso_2023_wormRNA"
            )
        return treated_cells
    
    def load_bujarrabal_dueso_2023(self) -> ExpressionDataset:
        """Load the Bujarrabal Dueso 2023 treated cells dataset
        ### Returns
        treated_cells : ExpressionDataset
            Expression dataset for treated cells
        """
        dataset_path = "/cellar/users/zkoch/dream/data/bujarrabal_dueso"
        treated_cells_fn = "GSE168401_Readcounts_Harmine_INDY.csv.gz"
        if not self.load_expression:
            raise ValueError("load_expression must be True to load Bujarrabal Dueso 2023")
        
        treated_cells_expr = pd.read_csv(
            os.path.join(dataset_path, treated_cells_fn),
            index_col=17, sep = ","
            )
        # drop index column
        treated_cells_expr.drop(columns = 'index', inplace = True)
        # create expression dataset object
        treated_cells = ExpressionDataset(
            expression_df=treated_cells_expr.T,
            species="human",
            metadata_df=pd.DataFrame(),
            dataset="bujarrabal_dueso_2023"
            )
        return treated_cells
    
    def load_cao_2024(self) -> ExpressionDataset:
        """Loads the Cao 2024 dataset
        ### Returns
        expr_dataset : ExpressionDataset
            Expression dataset for mice who also had SSB and AP measured
        """
        dataset_path = "/cellar/users/zkoch/dream/data/cao_et_al_2024"
        cao_expr = pd.read_csv(
            os.path.join(dataset_path, "GSE190955_merge.71samples.TPM.genes.txt.gz"),
            sep = "\t", index_col = 0, encoding_errors='replace', low_memory=False
            ).T
        cao_expr =  cao_expr.iloc[1:]
        cao_expr24 = pd.read_csv(
            os.path.join(dataset_path, "GSE239751_mouse.TPM.genes.24M.txt.gz"),
            sep = "\t", index_col = 0, encoding_errors='replace', low_memory=False
            ).T
        cao_expr24 =  cao_expr24.iloc[1:]
        cao_expr = pd.concat([cao_expr, cao_expr24])
        
        # first metadata
        cao_metadata = pd.read_csv(
            "/cellar/users/zkoch/dream/data/cao_et_al_2024/GSE190955-GPL24247_series_matrix.txt.gz",
            index_col=0, sep = "\t",
            skiprows=36, header=None
            )
        cao_metadata = cao_metadata.T
        cao_metadata = cao_metadata.iloc[:,[0,1,7,11,15,17]].copy(deep = True)
        cao_metadata.columns = ['sample_title', 'sample_geo_accession', 'tissue', 'age', 'dna/rna', 'protocol']
        # 24m metadata
        cao_metadata24 = pd.read_csv(
            "/cellar/users/zkoch/dream/data/cao_et_al_2024/GSE239751_series_matrix.txt.gz",
            index_col=0, sep = "\t",
            skiprows=36, header=None
            )
        cao_metadata24 = cao_metadata24.T
        cao_metadata24 = cao_metadata24.iloc[:,[-2,2,6,7,11]].copy(deep = True)
        cao_metadata24.columns = [ 'sample_geo_accession','tissue', 'age', 'dna/rna','protocol']
        sample_titles = cao_metadata24.iloc[0:12]['protocol'].apply(lambda x: x.split(',')[1].strip()).to_list()
        cao_metadata24['sample_title'] = [x + ' SSB' for x in sample_titles] + sample_titles + [x + ' RNA' for x in sample_titles]
        cao_metadata24.reset_index(inplace = True, drop = True)
        # combine metadata
        cao_metadata = pd.concat([cao_metadata, cao_metadata24])
            
        cao_expr = ExpressionDataset(
            expression_df=cao_expr,
            species="symbol",
            metadata_df=cao_metadata,
            dataset="cao_2024"
            )
        return cao_expr
    
    def load_liu_2023(self) -> ExpressionDataset:
        """Load the expression data of 100+ species
        ### Returns
        across_species : ExpressionDataset
            Expression dataset for across species
        """
        dataset_path = "/cellar/users/zkoch/dream/data/liu_2023"
        liu_expr = pd.read_csv(
            os.path.join(dataset_path, "bc_rpkm_tmm_log2_data.tsv"),
            #os.path.join(dataset_path, "rawdata.tsv"),
            sep = "\t", index_col = 0
            )
        liu_meta = pd.read_csv(
            os.path.join(dataset_path, "metainfo.tsv"), 
            sep = "\t", index_col = 0
            )
        # create expression dataset object
        across_species = ExpressionDataset(
            expression_df=liu_expr.T,
            species="human",
            metadata_df=liu_meta,
            dataset="liu_2023"
            )
        return across_species
    
    def load_nebulas_sc(self) -> ScExpressionDataset:
        """Load the many species expression datasets in nebulas as one
        ### Returns
        across_species : ScExpressionDataset
            ScExpression dataset for across species
        """
        dataset_path = "/cellar/users/zkoch/dream/data/cross_species/nebulas_sc/gsapub/ftp/pub/gen"
        # within dataset_path there are numerous subdirectories
        species_dirs = glob.glob(os.path.join(dataset_path, "*"))
        try:
            # read from h5ad
            processed_fn = "/cellar/users/zkoch/dream/data/cross_species/nebulas_sc/all_species_processed_w_ssgsea.h5ad"
            adata = sc.read_h5ad(processed_fn)
        except:
            # iterate across species dirs, getting the barcodes, features, and matrix
            adata_list = []
            for species_dir in species_dirs:
                meta_fn = glob.glob(os.path.join(species_dir, "*txt"))[0]
                meta_df = pd.read_csv(meta_fn, sep = '\t')
                # inside some species there are multiple subdirectories (corresponding to tissues)
                tissue_subdirs = glob.glob(os.path.join(species_dir, "*", "*"))
                for tissue_subdir in tissue_subdirs:
                    # check if the directory contains a .h5ad file
                    h5ad_fns = glob.glob(os.path.join(tissue_subdir, "*.h5ad"))
                    if len(h5ad_fns) > 0:
                        # read in   
                        adata = sc.read_h5ad(h5ad_fns[0])
                        # all cols uppercase
                        adata.var.index = adata.var.index.str.upper()
                        adata.columns = adata.var_names
                        adata_list.append(adata)
                    else: # no cached h5ad so read from mtx
                        # create adata from 10x
                        adata = sc.read_10x_mtx(
                            tissue_subdir, var_names='gene_symbols', cache=True
                            )
                        # find the appropriate row in the metadata file to set obs with
                        sample_id = tissue_subdir.split("/")[-1]
                        # find row with GEN Sample ID == sample_id
                        meta_row = meta_df.query("`GEN Sample ID` == @sample_id")
                        # get number of cells
                        n_cells = adata.shape[0]
                        # repeat the row n_cells times
                        meta_row = pd.concat([meta_row]*n_cells, ignore_index=True)
                        # set this as the obs for every cell in adata
                        adata.obs = meta_row
                        # save as h5ad 
                        adata.write_h5ad(os.path.join(tissue_subdir, f"{sample_id}.h5ad"))
                        adata_list.append(adata)
            # concatenate adata_list, merging on var_names
            adata = ad.concat(adata_list, join = 'outer', index_unique = None)
            # drop any var starting with LOC from adata
            print("!\nWARNING: Dropping any gene names containing LOC\n!")
            adata = adata[:, ~adata.var.index.str.contains("LOC")]
            # change all names to uppercase
            adata.var.index = adata.var.index.str.upper()
            # update the adata columns
            adata.columns = adata.var_names
            # make index unique
            adata.obs.reset_index(inplace=True, drop=True)
        
        # create expression object
        across_species = ScExpressionDataset(
            adata=adata,
            gene_species="symbol",
            dataset="nebulas_sc"
            )
        return across_species

    def map_species_symbol_to_human_ensembl(
        self,
        species: str,
        expr_df: pd.DataFrame
        ) -> pd.DataFrame:
        """Map species gene names containing symbols to human ensembl gene names
        ### Parameters
        species : str
            The species of the expression dataset
        expr_df : pd.DataFrame
            The expression dataframe
        ### Returns
        expr_df : pd.DataFrame
            The expression dataframe with human ensembl gene names as index
        """
        if species == 'Gallus gallus':
            symbols = expr_df.index.str.split("_").str[1].str.split('-').str[0]
        elif species == 'Canis lupus familiaris' or species == 'Ovis aries' or species == 'Bos taurus' or species == 'Macaca mulatta' or species == 'Callithrix jacchus' or species == 'Oryctolagus cuniculus' or species == 'Capra hircus' or species == 'Sus scrofa':
            symbols = expr_df.index.str.split("_").str[1]
        elif species == 'Mesocricetus auratus':
            # symbols are uppercase version of index
            symbols = expr_df.index.str.upper()
        elif species == 'Homo sapiens':
            expr_df.index = expr_df.index.str.split("_").str[0]
            return expr_df
        elif species == 'Mus musculus' or species == 'Rattus norvegicus':
            symbols = expr_df.index.str.split("_").str[1].str.upper()
        elif species == 'Mustela putorius furo':
            # already symbols
            symbols = expr_df.index
        elif species == 'Macaca fascicularis':
            # skip for now because need to do a species specific thang 
            df = pd.DataFrame()
            return df
        else:
            print(expr_df)
            raise NotImplementedError(f"Mapping species {species} to human ensembl not implemented")
        # set symbosl to index and drop nan and duplicates
        expr_df.index = symbols
        # drop rows with nan index
        expr_df = expr_df.loc[~expr_df.index.isna()]
        expr_df = expr_df.groupby(expr_df.index).sum()
        # get the human ensembl gene mapper
        dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
        ensembl_gene_mapper = dataset.query(attributes=[
            'external_gene_name',
            'ensembl_gene_id',
        ])
        # turn into series with unique index for mapping
        ensembl_gene_mapper = ensembl_gene_mapper.set_index('Gene name')['Gene stable ID']
        ensembl_gene_mapper = ensembl_gene_mapper[~ensembl_gene_mapper.index.duplicated(keep='first')]
        # do mapping
        expr_df.index = expr_df.index.map(
            ensembl_gene_mapper
            )
        # drop index with nan
        expr_df = expr_df.loc[~expr_df.index.isna()]
        # combine index with duplicate names by summing
        expr_df = expr_df.groupby(
            expr_df.index, axis=0
            ).sum()
        # raise error if df has less than 5000 rows
        if expr_df.shape[0] < 5000:
            raise ValueError(f"Expression dataframe has less than 5000 rows after mapping {species} symbols to human ensembl")
        return expr_df
            
    def load_nebulas(self) -> ExpressionDataset:
        """Load the many species expression datasets in nebulas as one
        ### Returns
        across_species : ExpressionDataset
            Expression dataset for across species
        """
        dataset_path = "/cellar/users/zkoch/dream/data/cross_species/nebulas/gsapub/ftp/pub/gen"
        # within dataset_path there are numerous subdirectories
        species_dirs = glob.glob(os.path.join(dataset_path, "*"))
        # get metadata and expression for each species
        metadata_dfs = []
        expression_dfs = []
        for species_dir in species_dirs:
            # get files
            metadata_fn = glob.glob(os.path.join(species_dir, "*meta.txt"))[0]
            expression_fn = glob.glob(os.path.join(species_dir,"*", "*GeneMatrix_rawCounts.txt"))[0]
            # read
            metadata_df = pd.read_csv(metadata_fn, index_col=None, sep = "\t")
            expression_df = pd.read_csv(expression_fn, index_col=0, sep="\t")
            # from metadata get the species name
            species = metadata_df.loc[0, "Species"]
            # convert gene names to human ensembl
            expression_df = self.map_species_symbol_to_human_ensembl(
                species, expression_df
                )
            if expression_df.shape[0] ==0:
                print(f"Skipping {species} because no gene names could be mapped to human ensembl or didn't try")
                continue
            print(f"Converted {species} gene names to human ensembl, resulting in {expression_df.shape[0]} genes and {expression_df.shape[1]} samples")
            # drop any column containing 'Unnamed' in name
            expression_df = expression_df.loc[:, ~expression_df.columns.str.contains('Unnamed')]
            # add to dict
            metadata_dfs.append(metadata_df)
            expression_dfs.append(expression_df)
        # merge all expression dataframes together, jooining on index
        expression_df = pd.concat(
            expression_dfs, axis=1
            )
        # same for metadata, but by columns
        metadata_df = pd.concat(
            metadata_dfs, axis=0
            )
        metadata_df = metadata_df[['GEN Sample ID', 'Species' , 'Age','Age unit','Gender', 'Tissue', 'Disease State','BTO Category', 'Development Stage', 'Case/Control','Case Detail', 'Control Detail','#Reads', 'Healthy Condition']]
        # create expression object
        across_species = ExpressionDataset(
            expression_df=expression_df.T,
            species="human",
            metadata_df=metadata_df,
            dataset="nebulas"
            )
        return across_species
        
    def load_pearson_2008(self) -> ExpressionDataset:
        """Load the Pearson 2008 dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        dataset_path = "/cellar/users/zkoch/dream/data/pearson_2008"
        treated_mice_fn = "GSE11845_series_matrix.txt"
        treated_mice_metadata_fn = "GSE11845_family.soft"
        if not self.load_expression:
            raise ValueError("load_expression must be True to load Mercken 2014")
        # load data
        treated_mice_metadata = self.read_soft(
            os.path.join(dataset_path, treated_mice_metadata_fn)
            )
        treated_mice_expr = pd.read_csv(
            os.path.join(dataset_path, treated_mice_fn),
            index_col=0, sep = "\t"
            )
        # create expression dataset object
        treated_mice = ExpressionDataset(
            expression_df=treated_mice_expr.T,
            species="mouse",
            metadata_df=treated_mice_metadata,
            dataset="pearson_2008"
            )
        return treated_mice
    
    def load_barger_2008(self) -> ExpressionDataset:
        """Load the barger 2008 dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        dataset_path = "/cellar/users/zkoch/dream/data/barger_2008"
        treated_mice_fn = "GSE11291_expression.csv" 
        treated_mice_metadata_fn = "GSE11291_family.soft"
        if not self.load_expression:
            raise ValueError("load_expression must be True to load Barger 2012")
        # load data
        treated_mice_metadata = self.read_soft(
            os.path.join(dataset_path, treated_mice_metadata_fn)
            )
        treated_mice_expr = self.read_soft_expression(
            os.path.join(dataset_path, treated_mice_fn), metadata = treated_mice_metadata,
            array_type='affymetrix'
            )
        
        # create expression dataset object
        treated_mice = ExpressionDataset(
            expression_df=treated_mice_expr,
            species="mouse",
            metadata_df=treated_mice_metadata,
            dataset="barger_2008"
            )
        return treated_mice
    
    def load_aon_2020(self) -> ExpressionDataset:
        """Load the Aon 2020 dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        dataset_path = "/cellar/users/zkoch/dream/data/aon_2020"
        treated_mice_fn = "GSE124294_expression.csv"
        treated_mice_metadata_fn = "GSE124294_family.soft"
        if not self.load_expression:
            raise ValueError("load_expression must be True to load Mercken 2014")
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
            dataset="aon_2020"
            )
        return treated_mice
        
    
    def load_gtex(self, kws: dict) -> List[Union[ExpressionDataset, MutationDataset]]:
        """Load the GTEx dataset
        ### Parameters
        kws : dict
            Keyword arguments
            - mutation_dataset: 'nieto' or 'yizhak'
        ### Returns
        expression : ExpressionDataset
            Expression dataset for GTEx
        mutation : MutationDataset
            Mutation dataset for GTEx
        """
        dataset_path = "/cellar/users/zkoch/dream/data/gtex"
        expression_fn = "gtex_expr_reads.parquet"
        metadata_fn = "gtex_metadata.parquet"
        # metadata, identifiers are SUBJID and SAMPID
        metadata_df = pd.read_parquet(os.path.join(dataset_path, metadata_fn))
        if self.load_expression:
            # expression df, columns = samples (but first col is gene names), rows = genes in ENSG format
            expression_df = pd.read_parquet(os.path.join(dataset_path, expression_fn))
            expression_df.set_index("Name", inplace=True) # set ENSEMBL gene names as index
            # create expression dataset object
            expression = ExpressionDataset(
                expression_df=expression_df.T,
                species="human",
                metadata_df=metadata_df,
                dataset="gtex"
                )
        else:
            expression = None
        if self.load_mutation:
            if kws['mutation_dataset'] == 'nieto':
                mutation_fn = "gtex_mutations.parquet"
                # rows are mutations, identifiers are sample_id and subject_id but idk how they matchup
                mutation_df = pd.read_parquet(os.path.join(dataset_path, mutation_fn))
            elif kws['mutation_dataset'] == 'yizhak':
                mutation_fn = "yizhak_et_al_2019_mutation_calls.xlsx"
                mutation_df = pd.read_excel(os.path.join(dataset_path, mutation_fn), skiprows=1)
            # create mutation dataset object
            mutation = MutationDataset(
                mutation_df=mutation_df,
                metadata_df=metadata_df,
                dataset="gtex"
                )
        else:
            mutation = None
            
        return [expression, mutation]
    
    def load_eisenberg_2016(self) -> ExpressionDataset:
        """Load the Eisenberg 2016 dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        dataset_path = "/cellar/users/zkoch/dream/data/eisenberg_2016"
        treated_mice_fn = "GSE86882_expression.csv"
        treated_mice_metadata_fn = "GSE86882_family.soft"
        if not self.load_expression:
            raise ValueError("load_expression must be True to load Eisenberg 2016")
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
            dataset="eisenberg_2016"
            )
        return treated_mice
    
    def load_neff_2013(self) -> ExpressionDataset:
        """Load the Neff 2013 dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        dataset_path = "/cellar/users/zkoch/dream/data/neff_2013"
        treated_mice_fn = "GSE41018_expression.csv"
        treated_mice_metadata_fn = "GSE41018_family.soft"
        if not self.load_expression:
            raise ValueError("load_expression must be True to load Neff 2013")
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
            dataset="neff_2013"
            )
        return treated_mice
    
    def load_zhang_2023(self) -> ExpressionDataset:
        """Load the Zhang 2023 dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        dataset_path = "/cellar/users/zkoch/dream/data/zhang_2023"
        treated_mice_fn = "GSE234563_Matrix.txt"
        treated_mice_metadata_fn = "GSE234563_family.soft"
        if not self.load_expression:
            raise ValueError("load_expression must be True to load Zhang 2023")
        # load data
        treated_mice_metadata = self.read_soft(
            os.path.join(dataset_path, treated_mice_metadata_fn)
            )
        treated_mice_expr = pd.read_csv(
            os.path.join(dataset_path, treated_mice_fn),
            index_col=0, sep = "\t"
            )
        # create expression dataset object
        treated_mice = ExpressionDataset(
            expression_df=treated_mice_expr.T,
            species="symbol",
            metadata_df=treated_mice_metadata,
            dataset="zhang_2023"
            )
        return treated_mice
    
    def load_mercken_2014(self) -> ExpressionDataset:
        """Load the mercken 2014 dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        dataset_path = "/cellar/users/zkoch/dream/data/mercken_2014"
        treated_mice_fn = "GSE49000_expression.csv"
        treated_mice_metadata_fn = "GSE49000_family.soft"
        if not self.load_expression:
            raise ValueError("load_expression must be True to load Mercken 2014")
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
            dataset="mercken_2014"
            )
        return treated_mice
    
    def load_yu_2012(self) -> ExpressionDataset:
        """Load the Yu 2012 dataset
        ### Returns
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        dataset_path = "/cellar/users/zkoch/dream/data/yu_2012"
        treated_mice_fn = "GSE39313_expression.csv"
        treated_mice_metadata_fn = "GSE39313_family.soft"
        if not self.load_expression:
            raise ValueError("load_expression must be True to load Yu 2012")
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
            dataset="yu_2012"
            )
        return treated_mice
    
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
        """tms_data_dir = "/cellar/users/zkoch/dream/data/tabula_muris_senis"
        # already has DREAM_normalized_enrichment_score calculated
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
        
        # from tms_old_geneset.pkl
        tms.get_dream_gene_expression()
        tms.scale_dream_by_seq_depth(
            col_name = 'DREAM_normalized_enrichment_score',
            eq = 'DREAM_normalized_enrichment_score ~ n_counts * n_genes'
        )
        tms.calculate_mutation_burden(
            'mutation_count_per_kb_top50expr', 
            top_perc_expr=.5, max_vaf = 0.6
            )
        tms.dream_expression.obs[mut_col] = tms.adata.obs[mut_col].values
        """
        # the above is saved into this
        with open('/cellar/users/zkoch/dream/data/tabula_muris_senis/tms_old_geneset.pkl', 'rb') as f:
            tms = pickle.load(f)
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
        print("WARNING: not loading treated_mice")
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
        """treated_mice_expr = pd.read_csv(
            os.path.join(dataset_path, treated_mice_fn),
            index_col=0
            )"""
        # read metadata from .soft files
        across_species_metadata = self.read_soft(
            os.path.join(dataset_path, across_species_metadata_fn)
            )
        """treated_mice_metadata = self.read_soft(
            os.path.join(dataset_path, treated_mice_metadata_fn)
            )"""
        # create expression dataset objects
        across_species = ExpressionDataset(
            expression_df=across_species_expr.T,
            species="mouse", # df is indxd by mouse genes
            metadata_df=across_species_metadata,
            dataset="mSalt"
            )
        """treated_mice = ExpressionDataset(
            expression_df=treated_mice_expr.T,
            species="mouse",
            metadata_df=treated_mice_metadata,
            dataset="mSalt"
            )"""
        return across_species#, treated_mice
        
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
                metadata_df=mutation_metadata_df,
                dataset = 'TCGA'
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
                            elif lines[i].startswith("!Sample_data_processing = "):
                                to_remove = "!Sample_data_processing = "
                                metadata['data_processing'].append(
                                    lines[i][len(to_remove):].strip()
                                    )
                            elif lines[i].startswith("!Sample_extract_protocol_ch1 = "):
                                to_remove = "!Sample_extract_protocol_ch1 = "
                                if "12 male mice" in lines[i]:
                                    continue
                                metadata['sample_extraction_protocol'].append(
                                    lines[i][len(to_remove):].strip()
                                    )
        # turn into dataframe
        # if the lists are not the same lenght pad with nans
        max_len = max([len(x) for x in metadata.values()])
        for key in metadata.keys():
            if len(metadata[key]) < max_len:
                metadata[key] = metadata[key] + [np.nan] * (max_len - len(metadata[key]))
        metadata = pd.DataFrame(metadata)
        
        if self.dataset_name == 'cao_2024':
            # drop rows that are nan in sample column
            metadata.dropna(subset=['sample'], inplace=True)
        else:
            # drop rows that have nan
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




"""

def map_species_gene_to_human_ensembl(
        self,
        species: str,
        expr_df: pd.DataFrame
        ) -> pd.DataFrame:
        Map species gene names to human ensembl gene names
        ### Parameters
        species : str
            The species of the expression dataset
        expr_df : pd.DataFrame
            The expression dataframe
        ### Returns
        expr_df : pd.DataFrame
            The expression dataframe with human ensembl gene names as index
        species_to_ensembl_name_mapping = {
            "Canis lupus familiaris": "clfamiliaris_gene_ensembl",
            "Bos taurus": "btaurus_gene_ensembl",
            "Mesocricetus auratus": "mauratus_gene_ensembl",
            "Gallus gallus": "ggallus_gene_ensembl",
        }
        if species == "Mesocricetus auratus":
            ensembl_query_name = 'external_gene_name'
            index_set_name = 'Gene name'
        else:
            ensembl_query_name = 'ensembl_gene_id'
            index_set_name = 'Gene stable ID'
        # convert species name to ensembl name
        ensembl_species_name = species_to_ensembl_name_mapping[species]
        # get the ensembl gene mapper
        dataset = Dataset(name=ensembl_species_name, host='http://www.ensembl.org')
        ensembl_gene_mapper = dataset.query(attributes=[
            ensembl_query_name,
            'hsapiens_homolog_ensembl_gene',
        ])
        # turn into series with unique index for mapping
        ensembl_gene_mapper = ensembl_gene_mapper.set_index(index_set_name)['Human gene stable ID']
        ensembl_gene_mapper = ensembl_gene_mapper[~ensembl_gene_mapper.index.duplicated(keep='first')]
        
        # remove trailing gene names
        expr_df.index = expr_df.index.str.split("_").str[0]
        if species == "Canis lupus familiaris":
            # replace middle part of gene names with "845"
            # idk why
            idx = expr_df.index.to_list()
            idx = [x[:9] + '845' + x[12:] for x in idx]
            expr_df.index = idx
        elif species == "Gallus gallus":
            # replace tenth string index with 1 and T with G
            idx = expr_df.index.to_list()
            idx = [x[:6] + 'G' + x[7:] for x in idx]
            idx = [x[:10] + '1' + x[11:] for x in idx]
            expr_df.index = idx
        print(expr_df)
        print(ensembl_gene_mapper)
        # do mapping
        expr_df.index = expr_df.index.map(
            ensembl_gene_mapper
            )
        # drop index with nan
        expr_df = expr_df.loc[~expr_df.index.isna()]
        # combine index with duplicate names by summing
        expr_df = expr_df.groupby(
            expr_df.index, axis=0
            ).sum()
        return expr_df
        
"""