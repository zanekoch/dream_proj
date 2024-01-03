import os
import pandas as pd
import argparse
import glob
import argparse

def combine_tcgabiolinks_files(
    path: str,
    datatype: str,
    sample_sheet: pd.DataFrame
    ) -> None:
    """Combine the TCGAbiolinks files into a single file
    ### Parameters
    path : str
        Path to data
    datatype : str
        Data type
    sample_sheet : pd.DataFrame
        Sample sheet, maps file names to case/sample IDs
    """
    if datatype == 'mutation':
        all_files = glob.glob(os.path.join(path, "*/*.maf.gz"))
    elif datatype == 'expression':
        all_files = glob.glob(os.path.join(path, "*/b*.tsv"))
    else:
        raise NotImplementedError(f"Datatype {datatype} not implemented")
    # read each file into a list of dataframes
    dfs = []
    for fn in all_files:
        if datatype == 'expression':
            dfs.append(read_expression_file(fn, sample_sheet))
        elif datatype == 'mutation':
            dfs.append(read_mutation_file(fn, sample_sheet))
    print(f"Read {len(dfs)} files")
    # concatenate the dataframes
    if datatype == 'mutation':
        df = pd.concat(dfs, axis=0).reset_index(drop = True)
    elif datatype == 'expression':
        df = pd.concat(dfs, axis=1)
    return df
    
def read_expression_file(
    fn: str,
    sample_sheet : pd.DataFrame
    ) -> pd.DataFrame:
    """Read a .tsv file of expression data
    ### Parameters
    fn : str
        Path to .tsv file
    sample_sheet : pd.DataFrame
        Dataframe with id matched to case/sample IDs
    ### Returns
    pivoted_df : pd.DataFrame
        Dataframe with genes as index and one column for the TPM of this file's 
    """
    # set row 2 as header and first column as index
    df = pd.read_csv(fn, sep='\t', index_col=None, low_memory=False, header=1)
    # drop first 4 rows
    df.drop(df.index[:4], inplace=True)
    # get the case id etc. from the filename
    fn_no_path = fn.split("/")[-1]
    this_file_sample_sheet = sample_sheet.query(
            "`File Name` == @fn_no_path"
            )
    # based on the file name, get the case id etc. from the sample sheet
    # `first` because of the processing done in /cellar/users/zkoch/dream/notebooks/010224_process_tcga_data.ipynb to remove lists of case/sample IDs
    df['case_id'] = this_file_sample_sheet['first_case_id'].values[0]    
    df['sample_id'] = this_file_sample_sheet['first_sample_id'].values[0]
    df['sample_type'] = this_file_sample_sheet['first_sample_type'].values[0]
    pivoted_df = df[['gene_id', 'tpm_unstranded']].set_index('gene_id')
    pivoted_df.columns = [df['sample_id'].values[0]]
    return pivoted_df

def read_mutation_file(
    fn: str,
    sample_sheet: pd.DataFrame
    ) -> pd.DataFrame:
    """Read a .maf file
    ### Parameters
    fn : str
        Path to .maf file
    sample_sheet : pd.DataFrame
        Dataframe with id matched to case/sample IDs
    ### Returns
    df : pd.DataFrame
        Dataframe with mutation data and case/sample IDs (all col names are relative derived from the sample sheet, e.g., case_id refers to the case ID from the sample sheet)
    """
    # no index column
    df = pd.read_csv(fn, sep='\t', skiprows=7, index_col=False, low_memory=False)
    # if present, drop columns (because they cause pyarrow to crash)
    to_drop = ['SOMATIC','PUBMED', 'PHENO']
    df.drop(to_drop, axis=1, inplace=True, errors='ignore')
    # get the case id from the filename
    fn_no_path = fn.split("/")[-1]
    this_file_sample_sheet = sample_sheet.query(
            "`File Name` == @fn_no_path"
            )
    # based on the file name, get the case id etc. from the sample sheet
    df['case_id'] = this_file_sample_sheet['first_case_id'].values[0]    
    df['sample_id'] = this_file_sample_sheet['Sample ID'].values[0]
    df['sample_type'] = this_file_sample_sheet['Sample Type'].values[0]
    return df

def map_fn_to_caseid(
    fn_no_path : str,
    sample_sheet : pd.DataFrame,
    df : pd.DataFrame
    ) -> pd.DataFrame:
    """Map the filename to the case ID
    ### Parameters
    fn_no_path : str
        Filename without path
    sample_sheet : pd.DataFrame
        Sample sheet
    df : pd.DataFrame
        Dataframe with mutation/expression data
    ### Returns 
    df : pd.DataFrame
        Dataframe with mutation/expression data and case/sample IDs added
    """ 
    this_file_sample_sheet = sample_sheet.query(
            "`File Name` == @fn_no_path"
            )
    # based on the file name, get the case id etc. from the sample sheet
    # `first` because of the processing done in /cellar/users/zkoch/dream/notebooks/010224_process_tcga_data.ipynb to remove lists of case/sample IDs
    df['case_id'] = this_file_sample_sheet['first_case_id'].values[0]    
    df['sample_id'] = this_file_sample_sheet['first_sample_id'].values[0]
    df['sample_type'] = this_file_sample_sheet['first_sample_type'].values[0]
    return df
    
def main(
    path: str = None,
    dataset_name: str = None,
    output_path: str = None,
    datatype: str = None,
    sample_sheet_fn: str = None,
    ) -> pd.DataFrame:
    """Main function
    ### Parameters
    path : str
        Path to data
    dataset_name : str
        Dataset name
    output_path : str
        Output path
    datatype : str
        Data type
    sample_sheet_fn : str
        Path to sample sheet
    ### Returns
    None
    """
    if path is None or dataset_name is None or output_path is None or datatype is None or sample_sheet_fn is None:
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--path",
            type=str,
            required=True,
            help="Path to data"
        )
        parser.add_argument(
            "--dataset",
            type=str,
            required=True,
            help="Dataset name"
        )
        parser.add_argument(
            "--output",
            type=str,
            required=True,
            help="Output path"
        )
        parser.add_argument(
            "--datatype",
            type=str,
            required=True,
            help="Data type"
        )
        parser.add_argument(
            "--sample_sheet",
            type=str,
            required=True,
            help="Path to sample sheet"
        )
        args = parser.parse_args()
        path = args.path
        dataset_name = args.dataset
        output_path = args.output
        datatype = args.datatype
        sample_sheet_fn = args.sample_sheet
        
    sample_sheet = pd.read_csv(
        sample_sheet_fn, sep='\t',
        header  = 0, index_col=False
        )
    # combine files
    df = combine_tcgabiolinks_files(
        path=path,
        datatype=datatype,
        sample_sheet=sample_sheet
    )
    # save as parquet
    #df.to_parquet(os.path.join(output_path, f"{dataset_name}_{datatype}.parquet"))
    #print(f"Saved to {os.path.join(output_path, f'{dataset_name}_{datatype}.parquet')}")
    return df 

if __name__ == "__main__":
    main()
