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
        
        
    def read_dream_files(self) -> None:
        """Read in DREAM files
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
        self.gene_conversion = pd.read_csv(
            "/cellar/users/zkoch/dream/utilities/human_mouse_ensembl_genes.txt.gz",
            sep="\t", index_col=0, 
            )
    
    def get_dream_gene_expression(
        self,
        row_limiting_query: str = None
        ) -> pd.DataFrame:
        """Get expression of DREAM genes
        ### Returns:
        dream_expression : pd.DataFrame
            Expression values of DREAM genes
        row_limiting_query : str
            Query to limit rows of dream_regulated_genes dataframe
        """
        # check if dream files have been read in
        if not hasattr(self, "dream_regulated_genes"):
            self.read_dream_files()
        # get human genes
        if row_limiting_query:
            dream_regulated_genes_names = self.dream_regulated_genes.query(row_limiting_query).index
        else:
            dream_regulated_genes_names = self.dream_regulated_genes.index
        # convert to mouse genes
        dream_regulated_genes_names_converted = self.gene_conversion.loc[
            dream_regulated_genes_names, "Mouse gene stable ID"
            ]
        dream_regulated_genes_names_converted.dropna(inplace=True)
        # get expression of converted genes
        dream_regulated_genes_w_expression = list(
            set(dream_regulated_genes_names_converted).intersection(set(self.expression_df.columns))
            )
        dream_expression = self.expression_df[dream_regulated_genes_w_expression].copy(deep = True)
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
        ) -> None:
        """Constructor for DatasetLoader
        ### Parameters
        dataset_name : str
            Name of the dataset to load
        
        ### Returns
        None
        """
        self.dataset_name = dataset_name

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
            elif "expression" in dataset_fn:
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
        if expression_fn != "":
            expression_df = pd.read_parquet(expression_fn)
            expression_dataset = ExpressionDataset(
                expression_df=expression_df,
                species="human",
                metadata_df=metadata_df
                )
            print(f"Created ExpressionDataset for {self.dataset_name}")
        if mutation_fn != "":
            mutation_df = pd.read_parquet(mutation_fn)
            mutation_dataset = MutationDataset(
                mutation_df=mutation_df,
                metadata_df=metadata_df
                )
            print(f"Created MutationDataset for {self.dataset_name}")
        if methylation_dir != "":
            print("reading dask methylation df")
            methylation_dd = dd.read_parquet(methylation_fn)
            # convert to pandas df
            methylation_df = methylation_dd.compute()
            methylation_dataset = MethylationDataset(
                methylation_df=methylation_df,
                metadata_df=metadata_df
                )
            print(f"Created MethylationDataset for {self.dataset_name}")
        # load from methylation fn if we don't have a methylation dir
        elif methylation_fn != "":
            methylation_df = pd.read_parquet(methylation_fn)
            methylation_dataset = MethylationDataset(
                methylation_df=methylation_df,
                metadata_df=metadata_df
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

    