import pandas as pd
import utils
import gseapy
import statsmodels.formula.api as smf
from pylr2 import regress2
import numpy as np

class MethylationDataset:
    """ Class to represent a methylation dataset """
    def __init__(
        self, 
        methyl_df: pd.DataFrame, 
        metadata_df: pd.DataFrame,
        probe_map: pd.DataFrame,
        dataset: str,
        species: str
        ) -> None:
        """Constructor for MethylationDataset
        ### Parameters:
        methyl_df : pd.DataFrame
            Methylation dataframe, samples x probes
        metadata_df : pd.DataFrame
            Metadata dataframe
        probe_map : pd.DataFrame
            Mapping from probe names to gene names. index is probe name, first col is gene name
        dataset : str
            Dataset name
        species : str
            Species corresponding to gene names in probe_map
        ### Returns:
        None
        """
        self.methyl_df = methyl_df
        self.metadata_df = metadata_df
        self.probe_map = probe_map
        self.dataset = dataset
        self.species = species

        if self.dataset == 'mammalian_methylation_consort':
            self.metadata_df.set_index('SID', inplace=True)
            # remove first character of metadata_df index
            self.metadata_df.index = self.metadata_df.index.str[1:]
            # merge
            self.methyl_df = self.methyl_df.merge(
                self.metadata_df, left_index=True, right_index=True
            )
            self.meta_cols = self.metadata_df.columns.to_list()
            # read in species stats 
            self.species_stats = pd.read_excel(
                '/cellar/users/zkoch/dream/data/mammalian_methylation_consort/li_2024/sciadv.adm7273_tables_s1_to_s11.xlsx',
                sheet_name='S1. anAge Updated', skiprows=1
                )
            self.species_stats.rename(columns = {
                "Species Latin Name": "latin_bane", 
                "Common Name": "common_name",
                "Maximum Lifespan":"ML",
                "Average Age of Sexual Maturity (Years)": "ASM",
                "Body Mass (g)": "AW",
                }, inplace = True)
            
    
    def convert_genes(self):
        """Convert the DREAM genes to the species of the genes in the probe map"""
        dream_regulated_genes_names = self.dream_regulated_genes.index
        # convert them if need be
        if self.species != "human":
            need_to_convert_genes = True
        else:
            need_to_convert_genes = False
            
        # get the genes that have an associated measured methylation probe
        measured_probes = self.probe_map.index.intersection(self.methyl_df.columns)
        measured_genes = self.probe_map.loc[measured_probes].iloc[:, 0].dropna().drop_duplicates().values
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
                print("Converted DREAM genes to mouse genes")
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
                print("Converted DREAM genes to rat genes")
            elif self.species == 'symbol':
                # use gene symbols instead of ensembl
                self.dream_regulated_genes.set_index('gene_name', inplace=True, drop = False)
                dream_regulated_genes_names_converted = self.dream_regulated_genes.index
                print("Converted DREAM genes to gene symbols")
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
                print("Converted DREAM genes to worm genes")
            #get dream genes that have measured methylation 
            self.dream_regulated_genes_w_methyl = list(
                set(dream_regulated_genes_names_converted).intersection(
                    set(measured_genes)
                )
            )
        else:
            self.dream_regulated_genes_w_methyl = list(
                set(dream_regulated_genes_names).intersection(
                    set(self.expression_df.columns)
                )
            )
            print("Did not need to convert DREAM genes")
        print(f"Found {len(self.dream_regulated_genes_w_methyl)} DREAM genes with measured methylation")

    def get_dream_activity(
        self, 
        use_new : bool = False
        ) -> pd.DataFrame:
        """Calculate the DREAM complex activity 
        ### Parameters:
        use_new : bool
            Wether to use the new dream file
        ### Returns:
        None
        """
        print("Using new dream file: ", use_new)
        self.dream_regulated_genes = utils.read_dream_files(use_new)
        print(f"Read in {len(self.dream_regulated_genes)} DREAM genes")
        # convert genes if necessary
        self.convert_genes()
        # from these genes, get the dream_regulated_probes
        self.dream_regulated_probes = self.probe_map.loc[
            self.probe_map.iloc[:, 0].isin(self.dream_regulated_genes_w_methyl)
            ].index.to_list()
        # select dream regulated cpgs and metadata
        self.dream_methylation = self.methyl_df[
            self.dream_regulated_probes + self.meta_cols
            ].copy(deep=True)
        # naive measure
        self.dream_methylation['mean_dream_mf'] = self.dream_methylation[
            self.dream_regulated_probes
            ].mean(axis=1)
        
        self.scale_dream_activity_by_global_mf('mean_dream_mf')
        
    def get_dream_activity_ssgsea(
        self,
        use_new: bool = False
        ) -> None:
        """Calculate the DREAM activity using ssGSEA"""
        print("Using new dream file: ", use_new)
        self.dream_regulated_genes = utils.read_dream_files(use_new)
        print(f"Read in {len(self.dream_regulated_genes)} DREAM genes")
        # convert genes if necessary
        self.convert_genes()
        # from these genes, get the dream_regulated_probes
        self.dream_regulated_probes = self.probe_map.loc[
            self.probe_map.iloc[:, 0].isin(self.dream_regulated_genes_w_methyl)
            ].index.to_list()
        
        # run ssgsea
        dream_gene_set_dict = {'dream_reg_genes':self.dream_regulated_probes}
        # select non-meta columns
        methyl_df = self.methyl_df.loc[:, ~self.methyl_df.columns.isin(self.meta_cols)]
        ssgsea = gseapy.ssgsea(
            data=methyl_df.T, gene_sets=dream_gene_set_dict,
            outdir=None, no_plot=True, max_size = 1000
            )
        results_df = ssgsea.res2d
        results_df.set_index('Name', inplace=True)
        results_df.drop(columns = ['Term'], inplace=True)
        results_df.rename(
            columns = {'NES':'DREAM_normalized_enrichment_score', 'ES': 'DREAM_enrichment_score'},
            inplace=True
            )
        results_df = results_df.astype(float)
        self.dream_methylation = self.dream_methylation.merge(
            results_df, left_index=True, right_index=True
            )
        self.scale_dream_activity_by_global_mf('DREAM_normalized_enrichment_score')
        self.scale_dream_activity_by_global_mf('DREAM_enrichment_score')
        
        
    def calculate_global_stats(self) -> None:
        """Calculate global stats for the methylation"""
        
        # calculate the mean methylation fraction
        non_meta_cols = self.methyl_df.columns.difference(self.meta_cols)
        mean_global_mf = self.methyl_df[non_meta_cols].mean(axis=1)
        self.methyl_df['mean_mf'] = mean_global_mf
        self.meta_cols.append('mean_mf')
        
    def scale_dream_activity_by_global_mf(
        self, 
        col_name: str,
        eq: str = None
        ) -> None:
        """Scale the DREAM activity by the mean whole genome methylation fraction or otherwise
        ### Parameters:
        col_name : str
            Name of the column to scale
        eq : str
            Equation to use for scaling
        ### Returns:
        None
        """
        if eq is None:
            eq = f'{col_name} ~ mean_mf'
        mut_ols = smf.ols(
            formula=eq,
            data=self.dream_methylation
            ).fit()
        self.dream_methylation[f'{col_name}_resid'] = mut_ols.resid
        # no need to invert residuals, because high methylation means high DREAM activity...hopefully
        
        print(f"scaled {col_name} by sequence depth and created {col_name}_resid")  