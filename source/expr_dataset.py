import pandas as pd
import scanpy as sc 
from pybiomart import Dataset
import statsmodels.formula.api as smf
import gseapy
import numpy as np
import glob
import os
import utils

class ExpressionDataset:
    """ Class to represent an expression dataset """
    def __init__(
        self, 
        expression_df: pd.DataFrame, 
        species: str,
        metadata_df: pd.DataFrame,
        dataset: str,
        drop_gene_period: bool = True
        ) -> None:
        """Constructor for ExpressionDataset
        ### Parameters:
        expression_df : pd.DataFrame
            Expression dataframe, samples x genes
        species : str
            String identifyin the species used to define genes in expression_df
        metadata_df : pd.DataFrame
            Metadata dataframe
        dataset : str
            Dataset name
        drop_gene_period : bool [default = True]
            Whether to drop the trailing .N or .NN from gene names
        ### Returns:
        None
        """
        self.expression_df = expression_df
        self.species = species
        self.metadata_df = metadata_df
        self.dataset = dataset
        if drop_gene_period:
            # remove the trailing .N or .NN from genes if present
            self.expression_df.columns = self.expression_df.columns.str.replace(
                r'\.\d+$', '', regex=True
                )
            # do again incase multiple .N or .NN
            self.expression_df.columns = self.expression_df.columns.str.split(".").str[0]
        self.dream_regulated_genes_w_expression = None
        self.already_converted_genes = False
        
        if self.dataset == 'tcga':
            # add ages, tissue, cancer type, tissue, sample_type, and case_id 
            self.meta_cols = ['age_at_index', 'primary_diagnosis', 'tissue_or_organ_of_origin', 'case_id', 'sample_type']
            self.expression_df = self.expression_df.merge(
                self.metadata_df[self.meta_cols],
                left_index=True, right_index=True
                )
            # add a column specifying the number of samples a case has
            self.get_num_samples_per_case()
            # simplify tissue and cancer names
            #self.simplify_tissue_and_cancer_names()
        elif self.dataset == 'mSalt':
            # check if any expressio index contains Baboon
            if any(self.expression_df.index.str.contains('Baboon')):
                # first convert mSalt_across_species.expression_df genes names (columns) from mouse ensemble to human ensemble
                ensembl_species_name = 'mmusculus_gene_ensembl'
                dataset = Dataset(name=ensembl_species_name, host='http://www.ensembl.org')
                ensembl_gene_mapper = dataset.query(attributes=[
                    'ensembl_gene_id',
                    'hsapiens_homolog_ensembl_gene',
                ])
                ensembl_gene_mapper.dropna(inplace=True)
                ensembl_gene_mapper = ensembl_gene_mapper.set_index('Gene stable ID')['Human gene stable ID']
                # drop duplicate indices
                ensembl_gene_mapper = ensembl_gene_mapper[~ensembl_gene_mapper.index.duplicated(keep='first')]
                self.expression_df.columns = self.expression_df.columns.map(ensembl_gene_mapper)
                self.expression_df = self.expression_df.loc[:, ~self.expression_df.columns.isna()]
                self.expression_df = self.expression_df.loc[:, ~self.expression_df.columns.duplicated()]

                # do this for the across species
                self.metadata_df['species'] = self.metadata_df['sample_title'].apply(
                    lambda x: x.split(", ")[1][0].upper() + x.split(", ")[1][1:]
                    )
                self.metadata_df['name'] = self.metadata_df['sample_title'].apply(
                    lambda x: x.split(", ")[1] + '_' + x.split(", ")[-1][-1] + '_' + x.split(", ")[0]
                    )
                # capitalize first letter of name
                self.metadata_df['name'] = self.metadata_df['name'].apply(
                    lambda x: x[0].upper() + x[1:]
                    )
                self.metadata_df.rename(
                    {'sex': 'Gender', 'tissue':'Tissue', 'age':'Age', 'species':'Species'},
                    inplace = True, axis = 1
                    )
                self.metadata_df.set_index("name", inplace=True)
                self.meta_cols = self.metadata_df.columns.to_list()
                self.expression_df = self.expression_df.merge(
                    self.metadata_df,
                    left_index=True, right_index=True
                    )
                #self.expression_df.set_index('name', inplace=True)
                #self.meta_cols.remove('name')
            else:
                # and this for the treated mice
                self.meta_cols = self.metadata_df.columns.to_list()
                """print(self.expression_df)
                print(self.metadata_df)
                self.expression_df = self.expression_df.merge(
                    self.metadata_df,
                    left_index=True, right_on = 'sample_name'
                    )"""
        elif self.dataset == 'tyshkovskiy':
            self.metadata_df.rename(columns = {'intervention':'condition'}, inplace=True)
            self.metadata_df.loc[
                self.metadata_df['condition'].str.contains('Control'), 'condition'
                ] = 'Control'
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            # remove sample_title from metacols
            self.meta_cols.remove('sample_title')
            self.expression_df.set_index('sample_title', inplace=True)
            # turn age into int
            self.expression_df['age'] = self.expression_df['age'].str.strip('months').astype(str)
        elif self.dataset == 'martin_montalvo':
            self.map_illm_probe_to_ensembl()
            # combine w metadata
            self.metadata_df['condition'] = self.metadata_df['sample_title'].apply(lambda x: x.split('-')[-2])
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'SD':'Control', 'MET':'Metformin'}
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
        elif self.dataset == 'boutant_nestle':
            # map the probe names to ensembl
            self.map_illm_probe_to_ensembl()
            # combine with metadata
            self.metadata_df.rename({'intervention':'condition'}, axis=1, inplace=True)
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'WT mice':'Control', 'EX mice':'Exercise', 'CR mice':'CR', 'SIRT1Tg/Tg mice':'SIRT1Tg/Tg'}
                )
            self.metadata_df['sample_title'] = self.metadata_df['sample_title'].str.split(' ').str[-1]
            # make tissue column
            self.metadata_df['tissue'] = self.metadata_df['sample_title'].str.split('_').str[0]
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            self.expression_df.set_index('sample_title', inplace=True)
            # remove sample_title from metacols
            self.meta_cols.remove('sample_title')
        elif self.dataset == 'zhou_2012':
            # convert age to int
            self.metadata_df['age'] = self.metadata_df['age'].str.strip('months')
            self.metadata_df['age'] = self.metadata_df['age'].str.strip('weeks')
            self.metadata_df['age'] = self.metadata_df['age'].astype(int)
            # last 4 rows are in months
            self.metadata_df.rename(columns = {'age':'age_months'}, inplace=True)
            self.metadata_df['condition'] = self.metadata_df['sample_title'].apply(
                lambda x: x.split('_')[0]
            )
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'LF+Ex':'LF+Exercise', 'HF+Ex':'HF+Exercise'}
            )
            # map the probe names to ensembl
            self.map_affy_probe_to_ensembl()
            # combine with metadata
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
        elif self.dataset == 'fok_chronic_2014':
            self.map_illm_probe_to_ensembl()
            self.metadata_df['condition'] = self.metadata_df['sample_title'].apply(lambda x: x.split(" ")[0])
            self.metadata_df['condition'] = self.metadata_df['condition'].str[:4]
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'Cont':'Control', 'Rapa':'Rapamycin'}
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
        elif self.dataset == 'fok_short_term_2014':
            self.map_illm_probe_to_ensembl()
            self.metadata_df['condition'] = self.metadata_df['sample_title'].apply(lambda x: x.split(" ")[0])
            self.metadata_df['condition'] = self.metadata_df['condition'].str[:4]
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'Cont':'Control', 'Rapa':'Rapamycin'}
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
        elif self.dataset == 'fok_cr_2014':
            self.map_illm_probe_to_ensembl()
            # keep everything else stuff after last space
            self.metadata_df['condition'] = self.metadata_df['sample_title'].apply(
                lambda x: ' '.join(x.split(" ")[:-1])
            )
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'Ad Libitum':'Control', 'Dietary Restriction':'CR', 'Rapamycin + Dietary Restriction':'Rapamycin+CR'}
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
        elif self.dataset == 'yu_2012':
            self.map_illm_probe_to_ensembl()
            self.metadata_df.rename(columns = {'treatment':'condition'}, inplace=True)
            # for rows with intervention == 'ApoE-FGF21-Tg', set condition to 'ApoE-FGF21-Tg'
            self.metadata_df.loc[
                self.metadata_df['intervention'] == 'ApoE-FGF21-Tg', 'condition'
                ] = 'ApoE-FGF21-Tg'
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'fed ad libitum':'Control', 'ApoE-FGF21-Tg':'ApoE-FGF21-Tg', 
                 '40% calorie restriction for 2 wks':'CR'}
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
            # drop 24 hour fasted animals bc this isn't control or treatment
            self.expression_df = self.expression_df[
                ~self.expression_df['condition'].str.contains('24')
                ]
            self.metadata_df = self.metadata_df[
                ~self.metadata_df['condition'].str.contains('24')
                ]
        elif self.dataset == 'mercken_2014':
            self.map_illm_probe_to_ensembl()
            self.metadata_df.rename(columns = {'treatment':'condition'}, inplace=True)
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'standard AIN-93G diet (SD)':'Control', 'Standard Diet plus SRT2104':'SRT2104', 
                 'Standard AIN-93G diet, 40% Caloric restriction (CR)':'CR'}
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_title'
                )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
        elif self.dataset == 'zhang_2023':
            self.metadata_df.rename(columns = {'intervention':'condition'}, inplace=True)
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'CreER':'Control', }
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df.columns = self.expression_df.columns.str.upper()
            self.expression_df = self.expression_df.merge(
                            self.metadata_df,
                            left_index=True, right_on = 'sample_title'
                            )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
        elif self.dataset == 'neff_2013':
            self.map_illm_probe_to_ensembl()
            # split on _
            self.metadata_df['condition'] = self.metadata_df['sample_title'].apply(
                lambda x: x.split("_")[1]
                )
            # for couple weird ones split again on space
            self.metadata_df['condition'] = self.metadata_df['condition'].apply(
                lambda x: x.split(" ")[-1]
                )
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                {'control':'Control', 'treated':'Rapamycin'}
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df.columns = self.expression_df.columns.str.upper()
            self.expression_df = self.expression_df.merge(
                            self.metadata_df,
                            left_index=True, right_on = 'sample_title'
                            )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
        elif self.dataset == 'eisenberg_2016':
            self.map_illm_probe_to_ensembl()
            self.metadata_df.rename(columns = {'treatment':'condition'}, inplace=True)
            self.metadata_df['condition'] = self.metadata_df['condition'].replace(
                    {'untreated':'Control', 'treated with 3 mM spermidine':'Spermidine'}
                )
            self.metadata_df['tissue'] = 'Heart'
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df.columns = self.expression_df.columns.str.upper()
            self.expression_df = self.expression_df.merge(
                            self.metadata_df,
                            left_index=True, right_on = 'sample_title'
                            )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
            # remove young rows
            self.expression_df = self.expression_df[~self.expression_df.index.str.contains('Young')]
            self.metadata_df = self.metadata_df[~self.metadata_df['sample_title'].str.contains('Young')]
        elif self.dataset == 'gtex':
            # drop first row of expression
            self.expression_df = self.expression_df.iloc[1:]
            # combine duplicate columns by summing
            self.expression_df = self.expression_df.groupby(
                self.expression_df.columns, axis=1
                ).sum()     
            # keep only columns of metadata we care about
            self.metadata_df.set_index('SAMPID', inplace=True)
            self.metadata_df = self.metadata_df[
                ['SUBJID', 'SEX','AGE','DTHHRDY', 'SMTS', 'SMTSD', 'SMRDTTL','SMMNCPB']
                ]
            self.metadata_df.rename(
                columns = {'SMTS':'tissue', 'SMTSD':'tissue_detail', 'SMRDTTL':'total_reads', 'SMMNCPB':'mean_coverage_per_base', 'DTHHRDY':'death_code'}
                , inplace=True
                )
            self.metadata_df['SEX'].replace(
                {1: 'Male', 2: 'Female'}, inplace=True
                )
            self.metadata_df['death_code'].replace(
                {1: 'fast_accident', 2: 'fast_illness', 3: 'medium_illness', 4: 'long_illness', 0: 'ventilator'}, inplace=True
                )
            # split on - to get age range and then take the mean
            self.metadata_df['AGE'] = self.metadata_df['AGE'].str.split('-').apply(
                lambda x: np.mean([int(i) for i in x])
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            # join with expression on index
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_index=True
                )
            
        elif self.dataset == 'aon_2020':
            self.map_illm_probe_to_ensembl()
            self.metadata_df['condition'] = self.metadata_df['sample_title'].apply(lambda x: x.split("_")[0])
            self.metadata_df['condition'].replace(
                {'Ad lib':'Control'}, inplace=True
                )
            self.metadata_df['diet_type'] = self.metadata_df['sample_title'].apply(lambda x: x.split("_")[1])
            self.metadata_df['diet_type'].replace(
                {'NIA':'low_sugar_low_fat', 'WIS':'high_sugar_high_fat'}, inplace=True
                )
            self.metadata_df['tissue'] = 'Liver'
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df.columns = self.expression_df.columns.str.upper()
            self.expression_df = self.expression_df.merge(
                            self.metadata_df,
                            left_index=True, right_on = 'sample_title'
                            )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
        elif self.dataset == 'barger_2008':    
            self.map_affy_probe_to_ensembl()
            self.metadata_df['tissue'] = self.metadata_df['sample_title'].apply(lambda x: x.split("-")[0])
            self.metadata_df['age'] = self.metadata_df['sample_title'].apply(lambda x: x.split("of")[0].split("-")[1].strip())
            self.metadata_df['condition'] = self.metadata_df['sample_title'].apply(lambda x: x.split("age")[-1].split("diet")[0].strip()[1:])
            self.metadata_df['condition'].replace(
                {'control':'Control', 'resveratrol':'Resveratrol',' resveratrol':'Resveratrol'}, inplace=True
                )
            # capitalize each first letter of tissue
            self.metadata_df['tissue'] = self.metadata_df['tissue'].apply(lambda x: ' '.join([i.capitalize() for i in x.split()]))
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df.columns = self.expression_df.columns.str.upper()
            self.expression_df = self.expression_df.merge(
                            self.metadata_df,
                            left_index=True, right_on = 'sample_title'
                            )
            self.expression_df.set_index('sample_title', inplace=True)
            self.meta_cols.remove('sample_title')
            # remove samples with age == '5 months'
            self.expression_df.query("age != '5 months'", inplace=True)
            self.metadata_df.query("age != '5 months'", inplace=True)
        elif self.dataset == 'pearson_2008':    
            self.map_illm_probe_to_ensembl()
            self.metadata_df['tissue'] = self.metadata_df['sample_title'].apply(lambda x: x.split(" ")[0])
            self.metadata_df['tissue'].replace(
                {'H': 'Heart', 'L': 'Liver', 'M': 'Muscle', 'F': 'Fat'}, inplace=True
                )
            self.metadata_df['condition'] = self.metadata_df['sample_title'].apply(lambda x: x.split(" ")[2])
            self.metadata_df['condition'].replace(
                {'SD':'Control', 'SDHR':'Resveratrol', 'SDLR':'Resveratrol',
                 'EOD': 'CR', 'EODHR':'CR + Resveratrol', 'EODLR':'CR + Resveratrol',
                 'HC': 'High fat', 'HCHR':'High fat + Resveratrol', 'HCLR':'High fat + Resveratrol'
                }, inplace=True
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df.columns = self.expression_df.columns.str.upper()
            self.expression_df = self.expression_df.merge(
                            self.metadata_df,
                            left_index=True, right_on = 'sample'
                            )
            self.expression_df.set_index('sample', inplace=True)
            self.meta_cols.remove('sample')
            # remove high fat samples
            self.expression_df.query("condition != 'High fat' and condition != 'High fat + Resveratrol'", inplace=True)
            self.metadata_df.query("condition != 'High fat' and condition != 'High fat + Resveratrol'", inplace=True)
            # remove rows where condition == 'CR' and tissue == 'Fat'
            self.expression_df.query("not (condition == 'CR' and tissue == 'Fat')", inplace=True)
            self.metadata_df.query("not (condition == 'CR' and tissue == 'Fat')", inplace=True)
        elif self.dataset == 'nebulas':
            # fix tissue names
            self.metadata_df['Tissue'].replace({
                'mixed (corpus callosum, frontal cortex, cerebellum, thalamus, hippocampus, caudate, medulla, pituitary, heart, liver, lung, kidney, skeletal muscle, small intestine, testis and ovary)':'Mixed', 'pineal gland':'Brain', 'retina':'Retina', 'brain ': 'Brain', 'oral pharynx':'Head/neck', 'prefrontal cortex': 'Brain', 'nasal pharynx':'Head/neck', 'olfactory bulb':'Head/neck', 'salivary gland':'Tongue', 'tongue':'Tongue','heart':'Heart', 'liver':'Liver', 'lung':'Lung', 'kidney':'Kidney', 'diaphragm':'Diaphram','epididymis':'Testis', 'fat pad': 'Adipose', 'gastric': 'Stomach', 'kidney facia':'Kidney', 'lymph node': 'Lymph Node', 'mesentery': 'Stomach','pancreas': 'Pancreas', 'prostate':'Prostate', 'small intestine': 'Intestine', 'Cerebrum':'Brain', 'Skeletal Muscle':'Muscle'
                }, inplace=True)
            # join
            self.metadata_df.set_index('GEN Sample ID', inplace=True, drop=True)
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_index=True
                )
            # drop samples with '-' or nan tissue
            self.expression_df = self.expression_df.query("Tissue != '-' ")
            self.expression_df.dropna(subset=['Tissue'], inplace=True)
            # also metadata
            self.metadata_df = self.metadata_df.query("Tissue != '-' ")
            self.metadata_df.dropna(subset=['Tissue'], inplace=True)
        elif self.dataset == 'liu_2023':
            # read in lifehistory
            self.lifehistory_df = pd.read_excel(
                "/cellar/users/zkoch/dream/data/liu_2023/life_history_traits.xlsx",
                skiprows = 2
                )
            self.lifehistory_df.columns = self.lifehistory_df.columns[:8].to_list() + ['AW/g'] + self.lifehistory_df.columns[9:].to_list()
            # merge lifehistory with expression
            self.expression_df.reset_index(inplace=True)
            self.expression_df.rename(columns = {"index": "Sample"}, inplace=True)
            self.expression_df['Species'] = self.expression_df['Sample'].apply(lambda x: x.split(".")[0])
            self.lifehistory_df['Species'] = self.lifehistory_df['Scientific_name'].str.replace(" ", "_")
            self.expression_df = self.expression_df.merge(
                self.lifehistory_df, left_on="Species", right_on="Species", how = 'left'
                )
            self.expression_df.set_index("Sample", inplace=True)
            # also merge in metadata
            self.expression_df = self.expression_df.merge(
                self.metadata_df, left_index=True, right_index=True
                )
            self.meta_cols = self.metadata_df.columns.to_list() + self.lifehistory_df.columns.to_list()
            # remove Scientific_name from meta_cols (twice, once for each df)
            self.meta_cols.remove('Scientific_name')
            self.meta_cols.remove('Scientific_name')
            self.expression_df.drop(columns =[ 'Scientific_name_x','Scientific_name_y'], inplace=True)
            #self.meta_cols.remove('Scientific_name_y')
        elif self.dataset == 'cao_2024':
            # read lesions
            """lesion_df = pd.read_csv("/cellar/users/zkoch/dream/data/cao_et_al_2024/lesion_counts_all.txt", sep="\t", header=None)"""
            """lesion_df = pd.read_csv(
                "/cellar/users/zkoch/dream/data/cao_et_al_2024/lesion_count_greater_than_one_support.txt",
                sep="\t", header=None
                )"""
            lesion_df = pd.read_csv(
                "/cellar/users/zkoch/dream/data/cao_et_al_2024/lesion_count_greater_than_5_support.txt",
                sep="\t", header=None
                )
            # make every odd numbered row its own column and keep even numbered rows as 
            names = lesion_df.loc[::2,0].reset_index(drop=True)
            counts = lesion_df.loc[1::2,0].reset_index(drop=True)
            lesion_df = pd.DataFrame({"name": names, "num_lesions": counts})
            lesion_df['name'] = lesion_df['name'].str.strip('./')
            lesion_df['sample'] = lesion_df['name'].str.split("_").str[0]
            lesion_df['num_lesions'] = lesion_df['num_lesions'].apply(
                lambda x: x.split(' ')[0]
            )
            lesion_df['num_lesions'] = lesion_df['num_lesions'].astype(int)
            
            # remove './ from name
            lesion_df['is_SSB'] = lesion_df['name'].str.contains("SSB")
            # for SSB set tissue to SSB
            lesion_df.loc[lesion_df['is_SSB'], 'tissue'] = lesion_df.loc[lesion_df['is_SSB'], 'name'].apply(
                lambda x: x.split(".")[1]
                )
            lesion_df.loc[lesion_df['is_SSB'], 'age'] = lesion_df.loc[lesion_df['is_SSB'], 'name'].apply(
                lambda x: x.split(".")[2]
                )
            lesion_df.loc[lesion_df['is_SSB'], 'age_int'] = lesion_df.loc[lesion_df['is_SSB'], 'age'].apply(
                lambda x: int(x[:-1])
                )
            lesion_df = lesion_df.merge(self.metadata_df, left_on="sample",right_on="sample_geo_accession", how="left")
            
            # combone age_x and age_y, using age_x if it exists and age_y if it does not
            lesion_df['age'] = lesion_df['age_x'].fillna(lesion_df['age_y'])
            # drop age_x, age_y, and tissue_x
            lesion_df = lesion_df.drop(columns=['age_x', 'age_y', 'tissue_x', 'age_int'])
            # rename tissue_y to tissue
            lesion_df = lesion_df.rename(columns={"tissue_y": "tissue"})
            lesion_df.dropna(axis = 0, inplace=True)
            lesion_df['age_int'] = lesion_df['age'].map({
                'age: 3 months' : 3,
                'age: 19 months' : 19,
                'age: 12 months' : 12,
                'age: 22 months' : 22,
                '3M' : 3,
                '19M' : 19,
                '12M' : 12,
                '22M' : 22,
                '24M' : 24,
                'development stage: 3 months': 3,
                'development stage: 12 months': 12,
                'development stage: 19 months': 19,
                'development stage: 22 months': 22,
                'development stage: 24 months': 24,
            })
            
            lesion_df['base_sample_title'] = lesion_df['sample_title'].str.strip(' SSB').apply(
                lambda x: x.split(',')[0]).apply(
                    lambda x: x.split('-rep1')[0]).apply(
                        lambda x: x.split('E-')[-1])
            # pivot so that we have one row per base_sample_title, with mutliple columns for the num_lesions where is_SSB is True or False  
            lesion_count_by_type = lesion_df.pivot_table(index='base_sample_title', columns='is_SSB', values='num_lesions', aggfunc='median')
            lesion_count_by_type.columns = ['AP_count', 'SSB_count']
            lesion_count_by_type = lesion_count_by_type.merge(lesion_df.drop_duplicates('base_sample_title'), left_index=True, right_on='base_sample_title', how='left').reset_index(drop = True)
            lesion_count_by_type.drop(columns = ['sample_title', 'is_SSB', 'num_lesions'], inplace = True)
            self.metadata_df = lesion_count_by_type
            
            
            # fix names to match between expression and metadata
            self.expression_df.reset_index(drop = False, inplace = True)
            self.expression_df['sample_title'] = self.expression_df['index'].apply(
                lambda x: x.split("RNA")[1][1:].split(".")[0]
            )
            # check if sample title has 2 dashes
            self.expression_df['dash_num'] = self.expression_df['sample_title'].apply(
                lambda x: len(x.split("-"))
            )
            # for rows where dash_num is 2, add a '-' to sample_title before the first letter (non number)
            self.expression_df.loc[self.expression_df['dash_num'] == 2, 'sample_title'] = self.expression_df.loc[
                self.expression_df['dash_num'] == 2, 'sample_title'].apply(
                lambda x: x[:2] + '-' + x[2:]
            )
            self.expression_df['sample_title'].replace({
                '3B--1': '3-B-1',
                '3B--2': '3-B-2',
                '3B--3': '3-B-3',
                '3M--1': '3-M-1',
                '3M--2': '3-M-2',
                '3M--3': '3-M-3',
                }, inplace = True)
            self.expression_df['sample_title'] = self.expression_df['sample_title'].astype(str)
            # sum duplicate columns
            self.expression_df = self.expression_df.groupby(level=0, axis=1).sum()
            # drop last 3 columns
            self.expression_df = self.expression_df.iloc[:, :-3]
            sample_title_save = self.expression_df['sample_title'].copy(deep=True)
            # remove weird non float non ascii characters
            self.expression_df = self.expression_df.apply(
                lambda x: pd.to_numeric(x, errors='coerce')
                )
            self.expression_df.dropna(axis = 1, inplace = True, how = 'any')
            self.expression_df.columns = self.expression_df.columns.str.upper()
            self.expression_df['sample_title'] = sample_title_save
            
            # merge
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                            self.metadata_df,
                            left_on = 'sample_title', right_on = 'base_sample_title',
                            how = 'left'
                            )
            self.expression_df.set_index('sample_title', inplace=True)
        elif self.dataset == 'bujarrabal_dueso_2023':
            # creat treatment column, if index contains harmine then 'Harmine', if contains control then control, if contains INDY then indy
            self.expression_df['treatment'] = 'control'
            self.expression_df['treatment'] = self.expression_df['treatment'].where(
                ~self.expression_df.index.str.contains('Harmine'), 'Harmine'
                )
            self.expression_df['treatment'] = self.expression_df['treatment'].where(
                ~self.expression_df.index.str.contains('Harmine_Control'), 'Harmine control'
                )
            self.expression_df['treatment'] = self.expression_df['treatment'].where(
                ~self.expression_df.index.str.contains('INDY'), 'INDY'
                )
            self.expression_df['treatment'] = self.expression_df['treatment'].where(
                ~self.expression_df.index.str.contains('Control_INDY'), 'INDY control'
                )
            # for rows where treatment contains 'INDY' remove last 5 characters
            self.metadata_df = self.expression_df['treatment'].copy(deep=True).to_frame()
            self.meta_cols = self.metadata_df.columns.to_list()
        elif self.dataset == 'bujarrabal_dueso_2023_wormRNA':
            # these values are from the series matrix file in https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152235
            # create treatment column with first 6 values being no treatment, next 3 being 6h post mock-UV-treatment, next 3 being 6h post UV-treatment
            values = ['6h post mock-UV-treatment']*6 + ['6h post UV-treatment']*6
            self.expression_df['treatment'] = values
            # create genotype column with first 3 being WT, next 3 lin-52(n771), next 3 WT, next 3 lin-52(n771)
            values = ['WT']*3 + ['lin-52(n771)']*3 + ['WT']*3 + ['lin-52(n771)']*3
            self.expression_df['genotype'] = values
            # create metadata df
            self.metadata_df = self.expression_df[['treatment', 'genotype']]
            self.meta_cols = self.metadata_df.columns.to_list()
        elif self.dataset == 'uxa_2019':
            # also read in chip binding data
            uxa_df = pd.read_excel("/cellar/users/zkoch/dream/data/uxa_2019/supp_table_s2.xlsx", sheet_name = None)
            # create df of genes with chip-seq data
            chip_binding = uxa_df['allgenes'][['ENSG', 'Gene Name','TP53 ChIP sum (0..15)','DREAM components binding (0..9)', 'RB/E2F binding (0..5)']]
            # replace '.' with np.nan
            chip_binding = chip_binding.replace('.', np.nan)
            chip_binding.dropna(inplace = True)
            # convert to float
            chip_binding['TP53 ChIP sum (0..15)'] = chip_binding['TP53 ChIP sum (0..15)'].astype(float)
            chip_binding['DREAM components binding (0..9)'] = chip_binding['DREAM components binding (0..9)'].astype(float)
            chip_binding['RB/E2F binding (0..5)'] = chip_binding['RB/E2F binding (0..5)'].astype(float)
            chip_binding.set_index('ENSG', inplace=True)
            self.chip_binding = chip_binding
            
            # split expression df index /
            self.expression_df['condition1'] = self.expression_df.index.str.split('/').str[0]
            self.expression_df['condition2'] = self.expression_df.index.str.split('/').str[1].str.split('_').str[0]
            
            # create metadata df
            self.metadata_df = self.expression_df[['condition1', 'condition2']]
            self.meta_cols = self.metadata_df.columns.to_list()
        elif self.dataset == 'motrpac_2024':
            phospho_fns = glob.glob('/cellar/users/zkoch/dream/data/motrpac_2024/proteomics_untargeted/analysis/prot-ph/normalized-data/*normalized-protein-corrected-logratio.txt')
            # read in each
            phospho_dfs = []
            for fn in phospho_fns:
                df = pd.read_csv(fn, sep='\t')
                df.set_index('viallabel', inplace=True)
                df = df.T
                # get tissue
                tissue_name = os.path.basename(fn).split('-')[2].split('_')[0]
                df['tissue'] = tissue_name
                phospho_dfs.append(df)
            self.phospho_df = pd.concat(phospho_dfs)
            # refseq gene names from uniprot
            self.rbl2_rat_refseq = ['NP_112356.1', 'tissue']
            self.lin52_rat_refseq = ['XP_006225911.1', 'XP_006240420.1', 'XP_001056671.2','XP_343090.4','tissue']
            # create metadata df
            self.metadata_df = self.expression_df[['tissue']]
            self.meta_cols = self.metadata_df.columns.to_list()
        elif self.dataset == 'gyenis_2023':
            def process_binned(df):
                df.drop(columns = 'strand', inplace=True)
                new_cols =[]
                for _ in range(6):
                    new_cols.extend(df.columns[:21].to_list())
                df.columns = new_cols
                # stack every 21 columnds of eu_seq_forward on top of eachother
                df_stacked = pd.concat([df.iloc[:, i:i+21] for i in range(0, df.shape[1], 21)], axis=0)
                # assign mice their ids'
                # TODO: make sure this is right
                num_times = int(df_stacked.shape[0]/6)
                new_col =   ['old_1'] * num_times + ['old_2'] * num_times + ['old_3'] * num_times + ['adult_1'] * num_times + ['adult_2'] * num_times + ['adult_3'] * num_times
                df_stacked['mouse_id'] = new_col
                return df_stacked

            data_dir = "/cellar/users/zkoch/dream/data/gyenis_2023"
            polii = "Total_pol2_20bin.xls"
            # for length
            polii_df = pd.read_excel(os.path.join(data_dir, polii), sheet_name="pol2_log2FC")
            gene_to_length_map = polii_df.set_index('Gene')['Length']
            gene_to_length_map.drop_duplicates(inplace=True)
            # stack nascent counts 
            eu_seq_stacked = process_binned(self.expression_df)
            # drop rows where id has a / in it, is length 1, or is nan
            eu_seq_stacked = eu_seq_stacked[
                ~(eu_seq_stacked['id'].str.contains('/') 
                  | (eu_seq_stacked['id'].str.len() == 1)
                  | eu_seq_stacked['id'].isna())
                ]
            # addlengths
            eu_seq_stacked['length_bp'] = eu_seq_stacked['id'].map(gene_to_length_map)
            # make id uppercase
            eu_seq_stacked['id'] = eu_seq_stacked['id'].str.upper()
            # sum columns that start with bin
            bin_cols = [col for col in eu_seq_stacked.columns if col.startswith('bin')]
            eu_seq_stacked['total_count'] = eu_seq_stacked[bin_cols].sum(axis=1)
            # create a normal expr df
            eu_seq_reg_expr = eu_seq_stacked.pivot_table(index='id', columns='mouse_id', values='total_count').T
            eu_seq_reg_expr['age'] = eu_seq_reg_expr.index.str.split('_').str[0]
            # get columsn that start with bin
            bin_cols = [col for col in eu_seq_stacked.columns if col.startswith('bin')]
            # divide the bin cols by 'bin2'
            """eu_seq_stacked[bin_cols] = eu_seq_stacked[bin_cols].div(eu_seq_stacked['bin1'], axis=0)"""
            self.eq_seq  = eu_seq_stacked
            # for plotting by bin: stack again, making index cols id and length
            eu_seq_stacked2 = eu_seq_stacked.set_index(
                ['id', 'length_bp', 'mouse_id']
                ).stack().reset_index().rename({'level_3': 'bin', 0: 'count'}, axis=1)
            eu_seq_stacked2['age'] = eu_seq_stacked2['mouse_id'].str.split('_').str[0]
            
            # now we have eu_seq_forward/reverse_stacked2 which is the stacked version of the binned data
            # and eu_seq_forward/reverse_reg_expr which is the regular expression format data of summed nascent seq counts
            self.expression_df = eu_seq_reg_expr
            self.stacked_nascenet_expr = eu_seq_stacked2
            self.metadata_df = self.expression_df[['age']]
            self.meta_cols = self.metadata_df.columns.to_list()
        elif self.dataset == 'lu_2014':
            self.metadata = metadata_df
            self.expression_df = expression_df
            self.species = species
            # map affy probe to ensembl
            #self.map_affy_probe_to_ensembl()
        elif self.dataset == 'williams_2022':
            self.metadata = metadata_df
            self.expression_df = expression_df
            self.species = species
            # convert appropriate meta cols to float
            numeric_meta_cols = list(set(self.metadata_df) - set(['Strain3','Diet3','Sex', 'SacDate3']))
            self.metadata_df[numeric_meta_cols] = self.metadata_df[numeric_meta_cols].astype(float)
            # remove EarTagCurrent from expression bc will get merged in w metadata
            self.expression_df.drop(columns = 'EarTagCurrent', inplace=True)
            
            # map from transcript to gene
            dataset = Dataset(name='mmusculus_gene_ensembl', host='http://www.ensembl.org')
            self.transcript_to_gene = dataset.query(attributes=['ensembl_transcript_id', 'ensembl_gene_id'])
            self.transcript_to_gene.set_index('Transcript stable ID', inplace = True)
            self.transcipts_w_mapping = list(
                set(self.expression_df.columns).intersection(self.transcript_to_gene.index)
                )
            self.cols_in_transcipts_w_mapping = [
                x for x in self.expression_df.columns if x in self.transcipts_w_mapping
                ]
            mapped_expression_df = self.expression_df.copy(deep = True)
            mapped_expression_df = mapped_expression_df[self.cols_in_transcipts_w_mapping]
            mapped_expression_df.columns = self.transcript_to_gene.loc[
                self.cols_in_transcipts_w_mapping
                ]['Gene stable ID'].values
            summed_expression_df = mapped_expression_df.groupby(level=0, axis=1).sum()
            self.expression_df = summed_expression_df
            # merge in metadat 
            self.expression_df = self.expression_df.merge(
                self.metadata_df, left_index=True, right_index=True
                )
            self.meta_cols = self.metadata_df.columns.to_list()
        elif self.dataset == 'lu_2022':
            self.metadata = metadata_df
            self.expression_df = expression_df
            self.species = species
            # make expression columns uppercase
            self.expression_df.columns = self.expression_df.columns.str.upper()
            # make metadata columns uppercase
            self.metadata_df.columns = self.metadata_df.columns.str.upper()
            self.meta_cols = self.metadata_df.columns.to_list()
        elif self.dataset == 'synapse_rna_seq_harmonization':
            self.metadata = metadata_df
            self.expression_df = expression_df
            self.species = species
            # drop first 4 expr columns
            self.expression_df = self.expression_df.iloc[:, 4:]
            # average duplicate columns
            self.expression_df = self.expression_df.groupby(
                self.expression_df.columns, axis=1
                ).mean()
            # TODO: might be losing longitudinal info or something here
            # drop rows with duplicate index values
            self.metadata_df = self.metadata_df[~self.metadata_df.index.duplicated(keep='first')]
            # merge in metadata
            self.expression_df = self.expression_df.merge(
                self.metadata_df, left_index=True, right_index=True, how = 'left'
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            
            # fix columns
            # replace '90+' with 90 and convert to int
            self.expression_df['age_first_ad_dx'] = self.expression_df['age_first_ad_dx'].replace('90+', 90).astype(float)
            self.expression_df['age_death'] = self.expression_df['age_death'].replace('90+', 90).astype(float)
        elif self.dataset == 'liu_2021':
            self.metadata = metadata_df
            self.expression_df = expression_df
            self.species = species
            self.expression_df = self.expression_df.merge(
                self.metadata_df, left_index=True, right_index=True, how = 'left'
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            # keep only columns with value_type == 'FPKM'
            self.expression_df = self.expression_df.loc[self.expression_df['value_type'] == 'fpkm']
            self.metadata_df = self.metadata_df[self.metadata_df['value_type'] == 'fpkm']
        elif self.dataset == 'crosby_2022':
            self.metadata = metadata_df
            self.expression_df = expression_df
            self.species = species
            self.expression_df = self.expression_df.merge(
                self.metadata_df, left_index=True, right_index=True, how = 'left'
                )
            self.meta_cols = self.metadata_df.columns.to_list()
        elif self.dataset == 'petljak_2019':
            self.metadata = metadata_df
            self.expression_df = expression_df
            self.species = species
            self.expression_df = self.expression_df.merge(
                self.metadata_df, left_index=True, right_index=True, how = 'left'
                )
            self.meta_cols = self.metadata_df.columns.to_list()
            # convert all non-meta columns to float
            self.expression_df[self.expression_df.columns.difference(self.meta_cols)] = self.expression_df[self.expression_df.columns.difference(self.meta_cols)].astype(float)
        else:
            raise NotImplementedError(f"Dataset {self.dataset} not implemented")

    def map_illm_probe_to_ensembl(self):
        """Map illumina probe names to ensembl gene names"""
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
        
    def map_affy_probe_to_ensembl(self, ensembl_gene_mapper):
        """Map affy probe names to ensembl gene names"""
        if self.species == 'mouse':
            dataset = Dataset(name='mmusculus_gene_ensembl', host='http://www.ensembl.org')
            ensembl_gene_mapper = dataset.query(
                attributes=['ensembl_gene_id', 'external_gene_name','affy_mouse430_2']
                )
            self.expression_df.columns = self.expression_df.columns.map(
                    ensembl_gene_mapper.set_index('AFFY Mouse430 2 probe')['Gene stable ID'].to_dict()
                    )
            # drop columns with nan
            self.expression_df = self.expression_df.loc[:, ~self.expression_df.columns.isna()]
            # combine columns with duplicate names by summing
            self.expression_df = self.expression_df.groupby(
                    self.expression_df.columns, axis=1
                    ).sum()
        elif self.species == 'human':
            """dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
            ensembl_gene_mapper = dataset.query(
                attributes=['ensembl_gene_id', 'external_gene_name','affy_hg_u133_plus_2']
                )"""
            self.expression_df.columns = self.expression_df.columns.map(
                    ensembl_gene_mapper.set_index('AFFY HG U133 Plus 2 probe')['Gene stable ID'].to_dict()
                    )
            # drop columns with nan
            self.expression_df = self.expression_df.loc[:, ~self.expression_df.columns.isna()]
            # combine columns with duplicate names by summing
            self.expression_df = self.expression_df.groupby(
                    self.expression_df.columns, axis=1
                    ).sum()
                
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
    
    def convert_genes(self):
        """Convert DREAM genes to the species of the expression data"""
        dream_regulated_genes_names = self.dream_regulated_genes.index
        # convert them if need be
        if self.species != "human":
            need_to_convert_genes = True
        else:
            need_to_convert_genes = False
        # if we need to convert genes
        if need_to_convert_genes:
            if self.species == "mouse":
                self.gene_converter = pd.read_csv(
                    "/cellar/users/zkoch/dream/utilities/human_mouse_ensembl_genes.txt.gz",
                    sep="\t", index_col=0, 
                )
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
            elif self.species == 'mouse_transcripts':
                self.gene_converter = pd.read_csv(
                    "/cellar/users/zkoch/dream/utilities/human_mouse_ensembl_genes.txt.gz",
                    sep="\t", index_col=0, 
                )
                # get the mouse ensembl ids for the dream genes (e.g. ENSMUSG)
                dream_regulated_genes_names_converted = self.gene_converter.loc[
                    dream_regulated_genes_names, "Mouse gene stable ID"
                ]
                print("Converting genes to transcripts")
                dataset = Dataset(name='mmusculus_gene_ensembl', host='http://www.ensembl.org')
                
                transcript_to_gene = dataset.query(attributes=['ensembl_transcript_id', 'ensembl_gene_id'])
                # convert these gene ids to transcript ids
                dream_regulated_transcripts = transcript_to_gene[
                    transcript_to_gene['Gene stable ID'].isin(dream_regulated_genes_names_converted)
                    ]
                self.dream_regulated_genes_w_expression = list(
                    set(dream_regulated_transcripts['Transcript stable ID']).intersection(
                        set(self.expression_df.columns)
                    )
                )
                print("Converted DREAM genes to mouse transcripts")
            elif self.species == "rat":
                self.gene_converter = pd.read_csv(
                    "/cellar/users/zkoch/dream/utilities/human_rat_ensembl_genes.txt",
                    sep="\t", index_col=0, 
                )
                # convert to rat genes
                dream_regulated_genes_names_converted = self.gene_converter.loc[
                    dream_regulated_genes_names, "Rat gene stable ID"
                ]
                dream_regulated_genes_names_converted.dropna(inplace=True)
                # get expression of converted genes
                self.dream_regulated_genes_w_expression = list(
                    set(dream_regulated_genes_names_converted).intersection(
                        set(self.expression_df.columns)
                    )
                )
                print("Converted DREAM genes to rat genes")
            elif self.species == 'symbol':
                # use gene symbols instead of ensembl
                self.dream_regulated_genes.set_index('gene_name', inplace=True, drop = False)
                dream_regulated_genes_names = self.dream_regulated_genes.index
                self.dream_regulated_genes_w_expression = list(
                    set(dream_regulated_genes_names).intersection(
                        set(self.expression_df.columns)
                    )
                )
            elif self.species == 'worm':
                # convert from WB gene to human gene
                self.gene_converter = pd.read_csv(
                    "/cellar/users/zkoch/dream/utilities/human_worm_WB_ensembl_genes.tsv",
                    sep="\t", index_col=None, 
                )
                # set Human gene stable ID to index
                self.gene_converter.dropna(subset=['Human gene stable ID'], inplace=True)
                self.gene_converter.set_index('Human gene stable ID', inplace=True)
                dream_regulated_genes_names_in_converter = list(
                    set(dream_regulated_genes_names).intersection(
                        set(self.gene_converter.index)
                    )
                )
                # select the rows which match ensembl dream gene names and get the WB gene names (called 'Gene stable ID' in the gene_converter df)
                dream_regulated_genes_names_converted = self.gene_converter.loc[
                    dream_regulated_genes_names_in_converter, "Gene stable ID"
                ]
                dream_regulated_genes_names_converted.dropna(inplace=True)
                # get expression of converted genes
                self.dream_regulated_genes_w_expression = list(
                    set(dream_regulated_genes_names_converted).intersection(
                        set(self.expression_df.columns)
                    )
                )
                print("Converted DREAM genes to worm genes")
        else:
            self.dream_regulated_genes_w_expression = list(
                set(dream_regulated_genes_names).intersection(
                    set(self.expression_df.columns)
                )
            )
            print("Did not need to convert DREAM genes")
        print(f"Found {len(self.dream_regulated_genes_w_expression)} DREAM genes with expression")
        self.already_converted_genes = True
    
    def get_dream_gene_expression(
        self,
        use_new: bool = False
        ) -> pd.DataFrame:
        """Get expression of DREAM genes
        ### Parameters:
        use_new : bool
            Wether to use the new dream file
        ### Returns:
        None
        """
        print("using new dream file: ", use_new)
        self.dream_regulated_genes = utils.read_dream_files(use_new)
        print(f"Read in {len(self.dream_regulated_genes)} DREAM genes")
        # convert genes if necessary
        self.convert_genes()
        # select the dream regulated genes from the expression df
        dream_expression = self.expression_df[
            self.dream_regulated_genes_w_expression + self.meta_cols
            ].copy(deep = True)
        self.dream_expression = dream_expression
        self.dream_expression['mean_dream_reg_expr'] = self.dream_expression[
            self.dream_regulated_genes_w_expression
            ].mean(axis=1)
        # convert to floats
        self.dream_expression['mean_dream_reg_expr'] = self.dream_expression['mean_dream_reg_expr'].astype(float)
        
        self.scale_dream_by_seq_depth(col_name = 'mean_dream_reg_expr')

    def scale_dream_by_seq_depth(self, col_name: str, eq: str = None) -> None:
        """Scale the DREAM activity by sequencing depth
        ### Parameters:
        col_name : str
            Name of the column to scale
        eq : str
            Equation to use for scaling
        ### Returns:
        None
        """
        if eq is None:
            eq = f'{col_name} ~ total_seq_depth * n_genes_expressed'
        mut_ols = smf.ols(
            formula=eq,
            data=self.dream_expression
            ).fit()
        self.dream_expression[f'{col_name}_resid'] = mut_ols.resid
        # invert residuals making highest values low and lowest values high
        self.dream_expression[f'{col_name}_resid'] = -1 * self.dream_expression[f'{col_name}_resid']
        # then add the min value to make all values positive
        self.dream_expression[f'{col_name}_resid'] = self.dream_expression[f'{col_name}_resid'] - min(self.dream_expression[f'{col_name}_resid'])
        print(f"scaled {col_name} by sequence depth and created {col_name}_resid")
        
    def dream_enrichment_ssgsea(
        self,
        use_new: bool = False
        ) -> None:
        """Run ssgsea on the DREAM genes
        ### Parameters:
        use_new : bool
            Wether to use the new dream file
        ### Returns:
        None
        """
        print("using new dream file: ", use_new)
        if not hasattr(self, 'dream_regulated_genes'):
            self.dream_regulated_genes = utils.read_dream_files(use_new)
        # convert genes if necessary
        if self.already_converted_genes == False:
            self.convert_genes()
        # run ssgsea
        dream_gene_set_dict = {'dream_reg_genes' : self.dream_regulated_genes_w_expression}
        # select non-meta columns
        expr_df = self.expression_df.loc[:, ~self.expression_df.columns.isin(self.meta_cols)]
        ssgsea = gseapy.ssgsea(
            data=expr_df.T, gene_sets=dream_gene_set_dict,
            outdir=None, no_plot=True, max_size = 1000, verbose = True
            )
        self.ssgsea = ssgsea
        results_df = ssgsea.res2d
        #return results_df
        results_df.set_index('Name', inplace=True)
        results_df.drop(columns = ['Term'], inplace=True)
        results_df.rename(
            columns = {'NES':'DREAM_normalized_enrichment_score', 'ES': 'DREAM_enrichment_score'},
            inplace=True
            )
        # convert cols to float
        results_df = results_df.astype(float)
        # merge to dream expression on index
        self.dream_expression = self.dream_expression.merge(
            results_df, left_index=True, right_index=True, how='left'
            )
        self.dream_expression['DREAM_normalized_enrichment_score'] = self.dream_expression['DREAM_normalized_enrichment_score'].astype(float)
        self.dream_expression['DREAM_enrichment_score'] = self.dream_expression['DREAM_enrichment_score'].astype(float)
        self.dream_expression['total_seq_depth'] = self.dream_expression['total_seq_depth'].astype(float)
        
        self.scale_dream_by_seq_depth(col_name = 'DREAM_normalized_enrichment_score')
        self.scale_dream_by_seq_depth(col_name = 'DREAM_enrichment_score')
        
    def test_differential_dream_expression(
        self,
        class_col: str,
        treatment_classes: list,
        control_class: str,
        secondary_grouping_col: str = ''
        ) -> None:
        """
        Run GSEA comparing treatment_classes to control_class. Optionally, group by secondary_grouping like tissue or sex
        ### Parameters:
        class_col : str
            Column name in metadata_df to get classes from
        treatment_classes : list
            List of classes to compare to control_class
        control_class : str
            Class to compare treatment_classes to
        secondary_grouping_col : str
            Column name in metadata_df to group by
        ### Returns:
            None
        """
        dream_gene_set_dict = {'dream_reg_genes':self.dream_regulated_genes_w_expression}
        gs_results_dfs = []
        if secondary_grouping_col == '':
            for treatment_class in treatment_classes:
                # select this class's rows
                expr_df = self.expression_df.query(
                    f"{class_col} == '{treatment_class}' or {class_col} == '{control_class}'"
                    )
                # get class col
                classes = expr_df[class_col].to_list()
                # select non-meta columns
                expr_df = expr_df.loc[:, ~expr_df.columns.isin(self.meta_cols)]
                expr_df = expr_df.T
                # run
                gs_res = gseapy.GSEA(
                    data=expr_df, gene_sets=dream_gene_set_dict,
                    classes=classes, permutation_type = 'gene_set', method='s2n'
                    )
                gs_res.pheno_pos = treatment_class
                gs_res.pheno_neg = control_class
                gs_res.run()
                gs_results_df = gs_res.res2d
                gs_results_df['treatment_class'] = treatment_class
                gs_results_df['control_class'] = control_class
                gs_results_df['secondary_grouping'] = ''
                gs_results_dfs.append(gs_results_df)
                print("finished", treatment_class, "vs. ", control_class)
            # check if class has attribute all_gs_results_df
            if not hasattr(self, 'all_gs_results_df'):
                self.all_gs_results_df = pd.concat(gs_results_dfs)
                # make treatment_class and control_class the first columns
                self.all_gs_results_df = self.all_gs_results_df[
                    ['treatment_class', 'control_class', 'secondary_grouping'] + self.all_gs_results_df.columns.to_list()[:-3]
                ]
            else:
                all_gs_results_df = pd.concat(gs_results_dfs)
                all_gs_results_df = all_gs_results_df[
                    ['treatment_class', 'control_class', 'secondary_grouping'] + all_gs_results_df.columns.to_list()[:-3]
                ]
                # add to all_gs_results_df
                self.all_gs_results_df = pd.concat(
                    [self.all_gs_results_df, all_gs_results_df]
                    )
        else:
            for secondary_grouping_val in self.expression_df[secondary_grouping_col].unique():
                for treatment_class in treatment_classes:
                    # select this class's rows
                    expr_df = self.expression_df.query(
                        f"{secondary_grouping_col} == '{secondary_grouping_val}' and ({class_col} == '{treatment_class}' or {class_col} == '{control_class}')"
                        )
                    # get class col
                    classes = expr_df[class_col].to_list()
                    if len(set(classes)) != 2:
                        print("missing 2 classes")
                        print(secondary_grouping_val)
                        print(treatment_class)
                        print(control_class)
                        print(expr_df)
                        print(classes)
                        continue
                    # select non-meta columns
                    expr_df = expr_df.loc[:, ~expr_df.columns.isin(self.meta_cols)]
                    expr_df = expr_df.T
                    
                    gs_res = gseapy.GSEA(
                        data=expr_df, gene_sets=dream_gene_set_dict,
                        classes=classes, permutation_type = 'gene_set', method='s2n'
                        )
                    gs_res.pheno_pos = treatment_class
                    gs_res.pheno_neg = control_class
                    gs_res.run()
                    
                    gs_results_df = gs_res.res2d
                    gs_results_df['treatment_class'] = treatment_class
                    gs_results_df['control_class'] = control_class
                    gs_results_df['secondary_grouping'] = secondary_grouping_val
                    gs_results_dfs.append(gs_results_df)
                    print("finished", treatment_class, "vs. ", control_class, "in", secondary_grouping_val)
            # check if class has attribute all_gs_results_df
            if not hasattr(self, 'all_gs_results_df'):
                self.all_gs_results_df = pd.concat(gs_results_dfs)
                # make treatment_class and control_class the first columns
                self.all_gs_results_df = self.all_gs_results_df[
                    ['treatment_class', 'control_class', 'secondary_grouping'] + self.all_gs_results_df.columns.to_list()[:-3]
                ]
            else:
                all_gs_results_df = pd.concat(gs_results_dfs)
                all_gs_results_df = all_gs_results_df[
                    ['treatment_class', 'control_class', 'secondary_grouping'] + all_gs_results_df.columns.to_list()[:-3]
                ]
                # add to all_gs_results_df
                self.all_gs_results_df = pd.concat(
                    [self.all_gs_results_df, all_gs_results_df]
                    )
    
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
        self.expression_df['tissue'] = self.expression_df['tissue_or_organ_of_origin'].map(simple_tissue_map)
        # drop uknown, none, or nan
        self.expression_df = self.expression_df.query(
            "tissue != 'Unknown' and tissue != 'None'"
            )
        self.expression_df.dropna(subset=['tissue'], inplace=True)
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
        self.expression_df['cancer_type'] = self.expression_df['primary_diagnosis'].map(simple_cancer_map)
        # drop uknown, none, or nan
        self.expression_df = self.expression_df.query(
            "cancer_type != 'Unknown' and cancer_type != 'None'"
            )
        self.expression_df.dropna(subset=['cancer_type'], inplace=True)
        self.meta_cols.extend(['tissue', 'cancer_type'])
    
    def calc_total_seq_depth(self) -> None:
        """Calculate the total sequence depth of each sample
        ### Returns:
        None
        """
        all_columns = self.expression_df.columns
        # remove meta cols
        non_meta_cols = [x for x in all_columns if x not in self.meta_cols]
        # get sum of all expression columns
        self.expression_df['total_seq_depth'] = self.expression_df.loc[:, non_meta_cols].sum(axis=1)
        # claculate the number of genes with expression not equal to 0
        # check if any values in self.expression_df[non_meta_cols] are less than 0
        if (self.expression_df[non_meta_cols] < 0).any().any():
            # TODO: Keep in mind if results change a bunch in might be bc of this
            self.expression_df['n_genes_expressed'] = self.expression_df[non_meta_cols].astype(bool).sum(axis=1)
        # if none are then we can do gt(0)
        else:
            self.expression_df['n_genes_expressed'] = self.expression_df[non_meta_cols].gt(0).sum(axis=1)
        # check if  ['total_seq_depth', 'n_genes_expressed'] are in meta_cols
        if 'total_seq_depth' not in self.meta_cols and 'n_genes_expressed' not in self.meta_cols:
            self.meta_cols = self.meta_cols + ['total_seq_depth', 'n_genes_expressed']
        
        
    def scale_by_total_seq_depth(
        self,
        cpm: bool = False
        ) -> None:
        """Scale the expression df by the total sequence depth of each sample
        ### Parameters:
        cpm : bool
            Wether to convert to cpm by multiplying by 1e6
        ### Returns:
        None
        """
        # get total sequence depth and number of genes expressed
        self.calc_total_seq_depth()
        
        # remove meta cols
        non_meta_cols = [x for x in self.expression_df.columns if x not in self.meta_cols]
        scaled_expression = self.expression_df[non_meta_cols].div(
            self.expression_df['total_seq_depth'], axis=0
            )
        if cpm:
            # convert to cpm
            scaled_expression = scaled_expression * 1e6
        # set the scaled expression columns to the scaled values
        self.expression_df = pd.concat(
            [scaled_expression, self.expression_df[self.meta_cols]], axis=1
            )
            
    def log_scale_expr(self) -> None:
        """Log scale the expression df
        ### Returns:
        None
        """
        all_columns = self.expression_df.columns
        # remove meta cols
        non_meta_cols = [x for x in all_columns if x not in self.meta_cols]
        # log scale the numerical columns
        scaled_expression_cols = np.log2(
            self.expression_df[non_meta_cols] + 1
            )
        # min-max scale to be in range 0 to 10
        scaled_expression_cols = scaled_expression_cols / scaled_expression_cols.max()
        scaled_expression_cols = scaled_expression_cols * 10
        # reconstruct the df
        self.expression_df = pd.concat(
            [scaled_expression_cols, self.expression_df[self.meta_cols]], axis=1
        )
        
    def batch_correct(
        self, 
        batch_col : str
        ) -> None:
        """Batch correct the expression df"""
        # get the expression columns
        all_columns = self.expression_df.columns
        # remove meta cols
        non_meta_cols = [x for x in all_columns if x not in self.meta_cols]
        # get the expression columns
        expression = self.expression_df[non_meta_cols]
        # batch correct
        expression = sc.pp.combat(
            expression.T, batch=self.expression_df[batch_col]
            ).T
        
    def _read_gencode_v46(self) -> pd.DataFrame:
        """Read in gencode v46 grch38 annotation and process it"""
        gencode_v46 = pd.read_csv(
            "/cellar/users/zkoch/dream/utilities/gencode.v46.basic.annotation.gtf.gz",
            sep='\t', skiprows=5, header=None
            )
        gencode_v46.columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
        # split attributes column
        gencode_v46['gene_id'] = gencode_v46['attributes'].str.split('gene_id "').str[1].str.split('"').str[0]
        gencode_v46['gene_name'] = gencode_v46['attributes'].str.split('gene_name "').str[1].str.split('"').str[0]
        gencode_v46['gene_name'] = gencode_v46['gene_name'].str.upper()
        
        return gencode_v46
    
    def _calculate_shannon_entropy(self, row, zero_handling_method = 'drop'):
        """Calculate Shannon entropy for a row of gene expression values.
        
        Parameters:
        -----------
        row : pandas.Series
            Row of gene expression values
            
        Returns:
        --------
        float
            Shannon entropy value
        """
        # Convert to probabilities by normalizing
        total = row.sum()
        if total == 0:
            return 0
        probabilities = row / total
        
        EPSILON = 1e-10
        
        if zero_handling_method == 'drop':
            probabilities = probabilities[probabilities > 0]
            entropy = -np.sum(probabilities * np.log2(probabilities))
        elif zero_handling_method == 'pseudocount':
            # Laplace smoothing: Add small constant
            smoothed = probabilities + EPSILON
            # Renormalize to ensure sum = 1
            probabilities = smoothed / smoothed.sum()
            entropy = -np.sum(probabilities * np.log2(probabilities))
        elif zero_handling_method == 'clip':
            # Clip values to minimum threshold
            probabilities = np.clip(probabilities, EPSILON, 1.0)
            # Renormalize
            probabilities = probabilities / probabilities.sum()

            entropy = -np.sum(probabilities * np.log2(probabilities))
        elif zero_handling_method == 'limit':
            # Use the fact that lim(x->0) x*log(x) = 0
            # Manually set 0*log(0) = 0
            entropy = -np.sum(np.where(probabilities > 0,
                                probabilities * np.log2(probabilities),
                                0))
        return entropy

    def calculate_shannon_entropy(self, zero_handling_method = 'drop'):
        """Calculate Shannon entropy for the expression df"""
        entropy = self.expression_df.iloc[:, :-len(self.meta_cols)].apply(
            self._calculate_shannon_entropy, zero_handling_method=zero_handling_method, axis=1
        )
        return entropy

    def _calculate_cross_entropy(self, row1, row2, zero_handling_method='drop'):
        """Calculate cross entropy between two expression vectors.
        
        Parameters:
        -----------
        row1 : pandas.Series
            First row of gene expression values (true distribution)
        row2 : pandas.Series
            Second row of gene expression values (predicted distribution)
        zero_handling_method : str
            Method for handling zero values. Options are:
            - 'drop': Remove zero values
            - 'pseudocount': Add small constant
            - 'clip': Clip to minimum threshold
            - 'limit': Use mathematical limit
            
        Returns:
        --------
        float
            Cross entropy value between the two distributions
        """
        # Convert to probabilities by normalizing
        total1 = row1.sum()
        total2 = row2.sum()
        if total1 == 0 or total2 == 0:
            return 0
        p = row1 / total1  # true distribution
        q = row2 / total2  # predicted distribution
        
        EPSILON = 1e-10
        
        if zero_handling_method == 'drop':
            # Only keep positions where both are non-zero
            mask = (p > 0) & (q > 0)
            p = p[mask]
            q = q[mask]
            cross_entropy = -np.sum(p * np.log2(q))
        elif zero_handling_method == 'pseudocount':
            # Add small constant and renormalize
            p = p + EPSILON
            q = q + EPSILON
            p = p / p.sum()
            q = q / q.sum()
            cross_entropy = -np.sum(p * np.log2(q))
        elif zero_handling_method == 'clip':
            # Clip values to minimum threshold
            p = np.clip(p, EPSILON, 1.0)
            q = np.clip(q, EPSILON, 1.0)
            # Renormalize
            p = p / p.sum()
            q = q / q.sum()
            cross_entropy = -np.sum(p * np.log2(q))
        elif zero_handling_method == 'limit':
            # Use limit definition where needed
            cross_entropy = -np.sum(np.where((p > 0) & (q > 0),
                                           p * np.log2(q),
                                           0))
        return cross_entropy
    
    def calculate_cross_entropy(self, reference_sample: str, zero_handling_method: str = 'pseudocount') -> pd.Series:
        """Calculate cross entropy between a reference sample and all other samples.
        
        Parameters
        ----------
        reference_sample : str
            Sample ID to use as reference distribution (must be in expression_df index)
        zero_handling_method : str, optional
            Method for handling zero values, by default 'pseudocount'
            See cross_entropy() method for options
            
        Returns
        -------
        pd.Series
            Series with cross entropy values between reference sample and all other samples
        """
        if reference_sample not in self.expression_df.index:
            raise ValueError(f"Reference sample {reference_sample} not found in expression data")
            
        
        # remove meta cols and the reference sample row
        entropy_calc_df = self.expression_df.iloc[:, :-len(self.meta_cols)].copy(deep = True)
        # Get reference distribution
        ref_dist = entropy_calc_df.loc[reference_sample]
        # drop the reference sample row
        """entropy_calc_df.drop(
            index = reference_sample,
            inplace = True
            )"""
        
        # Calculate cross entropy for each sample compared to reference
        cross_entropies = entropy_calc_df.apply(
            lambda row: self._calculate_cross_entropy(ref_dist, row, zero_handling_method),
            axis=1
        )
        
        return cross_entropies