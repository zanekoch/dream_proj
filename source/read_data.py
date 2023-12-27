import pandas as pd
import os
import pandas as pd
from collections import defaultdict

class ExpressionDataset:
    """ Class to represent an expression dataset """
    def __init__(
        self, 
        expression_df: pd.DataFrame, 
        expression_species: str,
        metadata_df: pd.DataFrame
        ) -> None:
        """Constructor for ExpressionDataset
        ### Parameters:
        expression_df : pd.DataFrame
            Expression dataframe, samples x genes
        expression_species : str
            Species of expression dataframe
        metadata_df : pd.DataFrame
            Metadata dataframe
        ### Returns:
        None
        """
        self.expression_df = expression_df
        self.expression_species = expression_species
        self.metadata_df = metadata_df
        
        
    def read_dream_files(self) -> None:
        """Read in DREAM files
        ### Returns:
        None
        """
        self.dream_regulated_genes = pd.read_csv(
            "/cellar/users/zkoch/dream/bujarrabal_dueso/tableS12_dream_promoter_binding.csv", index_col=0
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
        return dream_expression

class DatasetLoader:
    """Class to load a dataset"""
    def __init__(
        self,
        dataset_name:str
        ) -> None:
        """Constructor for DatasetLoader
        ### Parameters
        dataset_name : str
            Name of the dataset to load
        ### Returns
        None
        """
        self.dataset_name = dataset_name

    def load_dataset(self) -> list[ExpressionDataset]:
        """Load the dataset based on the dataset_name
        ### Returns
        A list of ExpressionDataset objects : list[ExpressionDataset]
        """
        print(f"Loading dataset: {self.dataset_name}")
        if self.dataset_name == "mSalt":
            across_species, treated_mice = self.load_mSalt()
            return across_species, treated_mice
        else:
            raise NotImplementedError(f"Dataset {self.dataset_name} not implemented")

    def load_mSalt(self) -> list[ExpressionDataset]:
        """Load the mSalt dataset
        ### Returns
        across_species : ExpressionDataset
            Expression dataset for across species
        treated_mice : ExpressionDataset
            Expression dataset for treated mice
        """
        # paths
        dataset_path = "/cellar/users/zkoch/dream/msalt"
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
            expression_species="mouse", # df is indxd by mouse genes
            metadata_df=across_species_metadata
            )
        treated_mice = ExpressionDataset(
            expression_df=treated_mice_expr.T,
            expression_species="mouse",
            metadata_df=treated_mice_metadata
            )
        return across_species, treated_mice
        
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


