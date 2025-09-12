import pandas as pd

class MutationDataset:
    """ Class to represent a mutation dataset """
    def __init__(
        self, 
        mutation_df: pd.DataFrame, 
        metadata_df: pd.DataFrame,
        dataset: str
        ) -> None:
        """Constructor for MutationDataset
        ### Parameters:
        mutation_df : pd.DataFrame
            Mutation dataframe, samples x genes
        metadata_df : pd.DataFrame
            Metadata dataframe
        dataset : str
            Dataset name
        ### Returns:
        None
        """
        self.mutation_df = mutation_df
        self.metadata_df = metadata_df
        self.dataset = dataset
        if self.dataset == 'TCGA':
            # keep columns we want and create self.mutation_burden df
            self.process_mutations()
            # add a column specifying the number of samples a case has
            self.get_num_samples_per_case()    
            # simplify tissue and cancer names
            self.simplify_tissue_and_cancer_names()
        elif self.dataset == 'gtex':
            pass
        
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
        