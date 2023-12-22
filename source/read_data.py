import pandas as pd
import os
import pandas as pd
from collections import defaultdict

class ExpressionDataset:
    """ Class to represent an expression dataset """
    def __init__(
        self, 
        expression_df: pd.DataFrame, 
        metadata_df: pd.DataFrame
        ) -> None:
        """Constructor for ExpressionDataset
        ### Parameters:
        expression_df : pd.DataFrame
            Expression dataframe
        metadata_df : pd.DataFrame
            Metadata dataframe
        ### Returns:
        None
        """
        self.expression_df = expression_df
        self.metadata_df = metadata_df


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
            expression_df=across_species_expr,
            metadata_df=across_species_metadata
            )
        treated_mice = ExpressionDataset(
            expression_df=treated_mice_expr,
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


