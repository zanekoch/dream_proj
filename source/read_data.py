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

# my library code 
from expr_dataset import ExpressionDataset
from meta_expr_dataset import MetaExpressionDataset
from sc_expr_dataset import ScExpressionDataset
from mutation_dataset import MutationDataset
from methyl_dataset import MethylationDataset

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
        elif self.dataset_name == 'yu_2012':
            dataset_list = self.load_yu_2012()
        elif self.dataset_name == 'mercken_2014':
            dataset_list = self.load_mercken_2014()
        elif self.dataset_name == 'zhang_2023':
            dataset_list = self.load_zhang_2023()
        else:
            raise NotImplementedError(f"Dataset {self.dataset_name} not implemented")
        return dataset_list
    
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

