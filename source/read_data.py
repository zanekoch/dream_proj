import pandas as pd
import os
import pandas as pd
from collections import defaultdict
import glob
import dask.dataframe as dd

class ExpressionDataset:
    """ Class to represent an expression dataset """
    def __init__(
        self, 
        expression_df: pd.DataFrame, 
        species: str,
        metadata_df: pd.DataFrame
        ) -> None:
        """Constructor for ExpressionDataset
        ### Parameters:
        expression_df : pd.DataFrame
            Expression dataframe, samples x genes
        species : str
            String identifyin the species used to define genes in expression_df
        metadata_df : pd.DataFrame
            Metadata dataframe
        ### Returns:
        None
        """
        self.expression_df = expression_df
        self.species = species
        self.metadata_df = metadata_df
        # remove the trailing .N or .NN from genes if present
        self.expression_df.columns = self.expression_df.columns.str.replace(
            r'\.\d+$', '', regex=True
            )
        # do again incase multiple .N or .NN
        self.expression_df.columns = self.expression_df.columns.str.split(".").str[0]
        # add ages, tissue, cancer type, tissue, sample_type, and case_id 
        self.meta_cols = ['age_at_index', 'primary_diagnosis', 'tissue_or_organ_of_origin', 'case_id', 'sample_type']
        self.expression_df = self.expression_df.merge(
            self.metadata_df[self.meta_cols],
            left_index=True, right_index=True
            )
        # add a column specifying the number of samples a case has
        self.get_num_samples_per_case()
        
    
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

    def load_dataset(self) -> list:
        """Load the dataset based on the dataset_name
        ### Returns
        A list of dataset objects : list
            The contents of this list depends on the dataset. mSalt it is a list of two ExpressionDataset objects, across_species and treated_mice. For CPTAC-3 it is an [ExpressionDataset, MutationDataset, MethylationDataset]
        """
        print(f"Loading dataset: {self.dataset_name}")
        if self.dataset_name == "mSalt":
            dataset_list = self.load_mSalt()
        elif self.dataset_name == 'CPTAC-3':
            dataset_list = self.load_tcga()
        else:
            raise NotImplementedError(f"Dataset {self.dataset_name} not implemented")
        return dataset_list
        

    def load_mSalt(self) -> list[ExpressionDataset]:
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
            metadata_df=across_species_metadata
            )
        treated_mice = ExpressionDataset(
            expression_df=treated_mice_expr.T,
            species="mouse",
            metadata_df=treated_mice_metadata
            )
        return across_species, treated_mice
        
    def load_tcga(self) -> list[ExpressionDataset, MutationDataset, MethylationDataset]:
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
                metadata_df=expr_metadata_df
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
                            elif lines[i].startswith("!Sample_characteristics_ch1 = age:"):
                                to_remove = "!Sample_characteristics_ch1 = age: " 
                                metadata['age'].append(
                                    lines[i][len(to_remove):].strip()
                                    )
                            elif lines[i].startswith("!Sample_description"):
                                metadata['sample_name'].append(
                                    lines[i].strip("!Sample_description = ").strip()
                                    )
                            elif lines[i].startswith("!Sample_characteristics_ch1 = treatment"):
                                to_remove = "!Sample_characteristics_ch1 = treatment: "
                                metadata['treatment'].append(
                                    lines[i][len(to_remove):].strip()
                                    )
                            elif lines[i].startswith("!Sample_characteristics_ch1 = dose"):
                                to_remove = "!Sample_characteristics_ch1 = dose: "
                                metadata['dose'].append(
                                    lines[i][len(to_remove):].strip()
                                    )                   
        # turn into dataframe
        metadata = pd.DataFrame(metadata)
        return metadata

    