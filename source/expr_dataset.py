import pandas as pd
import scanpy as sc 
from pybiomart import Dataset
import statsmodels.formula.api as smf
import gseapy
import numpy as np


class ExpressionDataset:
    """ Class to represent an expression dataset """
    def __init__(
        self, 
        expression_df: pd.DataFrame, 
        species: str,
        metadata_df: pd.DataFrame,
        dataset: str
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
        ### Returns:
        None
        """
        self.expression_df = expression_df
        self.species = species
        self.metadata_df = metadata_df
        self.dataset = dataset
        # remove the trailing .N or .NN from genes if present
        self.expression_df.columns = self.expression_df.columns.str.replace(
            r'\.\d+$', '', regex=True
            )
        # do again incase multiple .N or .NN
        self.expression_df.columns = self.expression_df.columns.str.split(".").str[0]
        self.dream_regulated_genes_w_expression = None
        
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
            self.simplify_tissue_and_cancer_names()
        elif self.dataset == 'mSalt':
            self.meta_cols = self.metadata_df.columns.to_list()
            self.expression_df = self.expression_df.merge(
                self.metadata_df,
                left_index=True, right_on = 'sample_name'
                )
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
            # rename 
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
        
    def map_affy_probe_to_ensembl(self):
        """Map affy probe names to ensembl gene names"""
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
        if self.species == "human" or self.species == "symbol":
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
            elif self.species == 'symbol':
                # use gene symbols instead of ensembl
                self.dream_regulated_genes.set_index('gene_name', inplace=True)
                dream_regulated_genes_names = self.dream_regulated_genes.index
                self.dream_regulated_genes_w_expression = list(
                    set(dream_regulated_genes_names).intersection(
                        set(self.expression_df.columns)
                        )
                    )
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
        self.dream_expression['mean_dream_reg_expr'] = self.dream_expression[
            self.dream_regulated_genes_w_expression
            ].mean(axis=1)
        self.scale_dream_by_seq_depth(col_name = 'mean_dream_reg_expr')

    def scale_dream_by_seq_depth(self, col_name: str):
        mut_ols = smf.ols(
            formula=f'{col_name} ~ total_seq_depth * n_genes_expressed',
            data=self.dream_expression
            ).fit()
        self.dream_expression[f'{col_name}_resid'] = mut_ols.resid
        # invert residuals making highest values low and lowest values high
        self.dream_expression[f'{col_name}_resid'] = -1 * self.dream_expression[f'{col_name}_resid']
        # then add the min value to make all values positive
        self.dream_expression[f'{col_name}_resid'] = self.dream_expression[f'{col_name}_resid'] - min(self.dream_expression[f'{col_name}_resid'])
        print(f"scaled {col_name} by sequence depth and created {col_name}_resid")
        
    def dream_enrichment_ssgsea(self) -> None:
        """Run ssgsea on the DREAM genes"""
        if self.dream_regulated_genes_w_expression is None:
            # run get_dream_gene_expression
            self.get_dream_gene_expression()
        # run ssgsea
        dream_gene_set_dict = {'dream_reg_genes':self.dream_regulated_genes_w_expression}
        # select non-meta columns
        expr_df = self.expression_df.loc[:, ~self.expression_df.columns.isin(self.meta_cols)]
        ssgsea = gseapy.ssgsea(
            data=expr_df.T, gene_sets=dream_gene_set_dict,
            outdir=None, no_plot=True
            )
        results_df = ssgsea.res2d
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
            results_df, left_index=True, right_index=True
            )
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
                    # run
                    try:
                        gs_res = gseapy.GSEA(
                            data=expr_df, gene_sets=dream_gene_set_dict,
                            classes=classes, permutation_type = 'gene_set', method='s2n'
                            )
                        gs_res.pheno_pos = treatment_class
                        gs_res.pheno_neg = control_class
                        gs_res.run()
                    except:
                        print(secondary_grouping_val)
                        print(treatment_class)
                        print(control_class)
                        print(expr_df)
                        print(classes)
                        raise ValueError
                    
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
    
    def scale_by_total_seq_depth(self) -> None:
        """Scale the expression df by the total sequence depth of each sample
        ### Returns:
        None
        """
        all_columns = self.expression_df.columns
        # remove meta cols
        non_meta_cols = [x for x in all_columns if x not in self.meta_cols]
        # get sum of all expression columns
        self.expression_df['total_seq_depth'] = self.expression_df.loc[:, non_meta_cols].sum(axis=1)
        # claculate the number of genes with expression not equal to 0
        # TODO: Keep in mind if results change a bunch in might be bc of this
        self.expression_df['n_genes_expressed'] = self.expression_df[non_meta_cols].astype(bool).sum(axis=1)
        self.meta_cols = self.meta_cols + ['total_seq_depth', 'n_genes_expressed']
        
        # scale the numerical columns by the total sequence depth
        scaled_expression = self.expression_df[non_meta_cols].div(
            self.expression_df['total_seq_depth'], axis=0
            )
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