import pandas as pd
# my library code
from read_data import DatasetLoader
from typing import Union, List


class MetaExpressionDataset:
    """ Class to represent multiple expression datasets """
    def __init__(
        self,
        dataset_names : List[str]
        ) -> None:
        self.dataset_names = dataset_names
        self.datasets = {}
        # load each dataset
        for dataset_name in self.dataset_names:
            loader = DatasetLoader(dataset_name)
            expression_dset = loader.load_dataset()           
            self.datasets[dataset_name] = expression_dset
        
    def scale_by_total_seq_depth(self):
        for dataset_name in self.dataset_names:
            self.datasets[dataset_name].scale_by_total_seq_depth()
            
    def log_scale_expr(self):
        for dataset_name in self.dataset_names:
            self.datasets[dataset_name].log_scale_expr()
            
    def get_dream_gene_expression(self):
        for dataset_name in self.dataset_names:
            self.datasets[dataset_name].get_dream_gene_expression()
            
    def dream_enrichment_ssgsea(self):
        for dataset_name in self.dataset_names:
            self.datasets[dataset_name].dream_enrichment_ssgsea()
            
    def run_GSEA(self):
        all_gs_results_dfs = []
        for dataset_name in self.dataset_names:
            if dataset_name == 'tyshkovskiy':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query(
                        "condition != 'Control'"
                        )['condition'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'age'
                    )
                # since does not have tissue as secondary grouping col, ad tissue
                tissues = self.datasets[dataset_name].metadata_df['tissue'].unique()
                assert len(tissues) == 1
                self.datasets[dataset_name].all_gs_results_df['tissue'] = tissues[0]
            elif dataset_name == 'boutant_nestle':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'tissue'
                    )
                # since does have tissue as secondary grouping col, duplicate renaming it to tissue
                self.datasets[dataset_name].all_gs_results_df['tissue'] = self.datasets[
                    dataset_name
                    ].all_gs_results_df['secondary_grouping']
            elif dataset_name == 'martin_montalvo':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'tissue'
                    )
                # since does have tissue as secondary grouping col, duplicate renaming it to tissue
                self.datasets[dataset_name].all_gs_results_df['tissue'] = self.datasets[
                    dataset_name
                    ].all_gs_results_df['secondary_grouping']
            elif dataset_name == 'zhou_2012':
                self.datasets[dataset_name].expression_df['is_high_fat'] = self.datasets[dataset_name].expression_df['condition'].str.contains('HF')
                self.datasets[dataset_name].expression_df['is_high_fat'] = self.datasets[dataset_name].expression_df['is_high_fat'].map({True:'HF', False:'LF'})
                self.datasets[dataset_name].metadata_df['is_high_fat'] = self.datasets[dataset_name].metadata_df['condition'].str.contains('HF')
                self.datasets[dataset_name].metadata_df['is_high_fat'] = self.datasets[dataset_name].metadata_df['is_high_fat'].map({True:'HF', False:'LF'})
                self.datasets[dataset_name].expression_df['condition2'] = self.datasets[dataset_name].expression_df['condition'].map({
                    'LF':'Control', 'LF+Exercise':'Exercise','LF+CR':'CR', 'HF+Exercise':'Exercise', 'HF':'Control', 'HF+CR':'CR'
                })
                self.datasets[dataset_name].metadata_df['condition2'] = self.datasets[dataset_name].metadata_df['condition'].map({
                    'LF':'Control', 'LF+Exercise':'Exercise','LF+CR':'CR', 'HF+Exercise':'Exercise', 'HF':'Control', 'HF+CR':'CR'
                })
                self.datasets[dataset_name].meta_cols.append('condition2')
                self.datasets[dataset_name].meta_cols.append('is_high_fat')

                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition2',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition2 != 'Control'")['condition2'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'is_high_fat'
                    )
                # since does not have tissue as secondary grouping col, ad tissue
                tissues = self.datasets[dataset_name].metadata_df['tissue'].unique()
                assert len(tissues) == 1
                self.datasets[dataset_name].all_gs_results_df['tissue'] = tissues[0]
            elif dataset_name == 'fok_chronic_2014' or dataset_name == 'fok_short_term_2014':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'sex'
                    )
                # since does not have tissue as secondary grouping col, ad tissue
                tissues = self.datasets[dataset_name].metadata_df['tissue'].unique()
                assert len(tissues) == 1
                self.datasets[dataset_name].all_gs_results_df['tissue'] = tissues[0]
            elif dataset_name == 'fok_cr_2014':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control'
                )
                # since does not have tissue as secondary grouping col, add tissue
                tissues = self.datasets[dataset_name].metadata_df['tissue'].unique()
                assert len(tissues) == 1
                self.datasets[dataset_name].all_gs_results_df['tissue'] = tissues[0]
            elif dataset_name == 'yu_2012':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'tissue'
                )
                # since does have tissue as secondary grouping col, duplicate renaming it to tissue
                self.datasets[dataset_name].all_gs_results_df['tissue'] = self.datasets[
                    dataset_name
                    ].all_gs_results_df['secondary_grouping']
            elif dataset_name == 'mercken_2014':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'tissue'
                )
                # since does have tissue as secondary grouping col, duplicate renaming it to tissue
                self.datasets[dataset_name].all_gs_results_df['tissue'] = self.datasets[
                    dataset_name
                    ].all_gs_results_df['secondary_grouping']
            elif dataset_name == 'zhang_2023':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'tissue'
                )
                # since does have tissue as secondary grouping col, duplicate renaming it to tissue
                self.datasets[dataset_name].all_gs_results_df['tissue'] = self.datasets[
                    dataset_name
                    ].all_gs_results_df['secondary_grouping']
            elif dataset_name == 'neff_2013':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'tissue'
                )
                # since does have tissue as secondary grouping col, duplicate renaming it to tissue
                self.datasets[dataset_name].all_gs_results_df['tissue'] = self.datasets[
                    dataset_name
                    ].all_gs_results_df['secondary_grouping']
            elif dataset_name == 'eisenberg_2016':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control'
                )
                # since does not have tissue as secondary grouping col, add tissue
                tissues = self.datasets[dataset_name].metadata_df['tissue'].unique()
                assert len(tissues) == 1
                self.datasets[dataset_name].all_gs_results_df['tissue'] = tissues[0]
            elif dataset_name == 'aon_2020':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control', 
                )
                # since does not have tissue as secondary grouping col, add tissue
                tissues = self.datasets[dataset_name].metadata_df['tissue'].unique()
                assert len(tissues) == 1
                self.datasets[dataset_name].all_gs_results_df['tissue'] = tissues[0]
            elif dataset_name == 'barger_2008':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'tissue'
                )
                 # since does have tissue as secondary grouping col, duplicate renaming it to tissue
                self.datasets[dataset_name].all_gs_results_df['tissue'] = self.datasets[
                    dataset_name
                    ].all_gs_results_df['secondary_grouping']
            elif dataset_name == 'pearson_2008':
                self.datasets[dataset_name].test_differential_dream_expression(
                    class_col = 'condition',
                    treatment_classes = self.datasets[dataset_name].metadata_df.query("condition != 'Control'")['condition'].unique().tolist(),
                    control_class = 'Control', secondary_grouping_col = 'tissue'
                )
                 # since does have tissue as secondary grouping col, duplicate renaming it to tissue
                self.datasets[dataset_name].all_gs_results_df['tissue'] = self.datasets[
                    dataset_name
                    ].all_gs_results_df['secondary_grouping']
            else:
                raise NotImplementedError(f"run_GSEA not implemented for {dataset_name}")
            self.datasets[dataset_name].all_gs_results_df['dataset'] = dataset_name
            all_gs_results_dfs.append(self.datasets[dataset_name].all_gs_results_df)
        self.all_gs_results_df = pd.concat(all_gs_results_dfs)
        self.all_gs_results_df.reset_index(inplace=True, drop = True)
    