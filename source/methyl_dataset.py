import pandas as pd
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
