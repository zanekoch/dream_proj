Downloaded from https://cellxgene.cziscience.com/collections/0b9d8a04-bb9d-44da-aa27-705bb65b54eb by clicking download next to the first results (All - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse) and then selecting the h5ad option which resulted in this link https://datasets.cellxgene.cziscience.com/2d6aefde-d7ab-4cd0-88d7-f4fba6f94504.h5ad

Raw data appears to be downloadable from AWS https://s3.console.aws.amazon.com/s3/buckets/czb-tabula-muris-senis/?region=us-west-2&bucketType=general&tab=objects

Following the link from the github to "processed data ready for use with scanpy" found mytself at this website (https://figshare.com/projects/Tabula_Muris_Senis/64982) then clicked procesesed files at the bottom and ran wget on the link to download all 5 files `wget https://figshare.com/ndownloader/articles/8273102/versions/3`. Then unzipped the files `unzip 3`

supplementary_table_8.xlsx is supplementary_table_8 from A single-cell transcriptomic atlas characterizes ageing tissues in the mouse 2020 (https://www.nature.com/articles/s41586-020-2496-1#Sec25)

adata_with_ercc_genecode_counts_for_gatk_with_metadata.h5ad: used in the creation of mutations in the notebook (https://github.com/czbiohub-sf/tabula-muris-senis/blob/master/1_tabula_muris_senis/13_figure_3/tabula-muris-senis-mutation-analysis.ipynb). Downloaded from aws bucket wget https://czb-tabula-muris-senis.s3.us-west-2.amazonaws.com/Data-to-reproduce-figures/mutation-analysis-objs/adata_with_ercc_genecode_counts_for_gatk_with_metadata.h5ad

adata_genecode_counts_for_gatk_with_metadata.h5ad: Downloaded from aws bucket wget
adata_with_ercc_genecode_counts_for_gatk_with_metadata.h5ad: Downloaded from aws bucket wget

# raw read count values
/cellar/users/zkoch/dream/data/tabula_muris_senis/tabula-muris-senis-facs-official-raw-obj.h5ad
# processed by TMS
/cellar/users/zkoch/dream/data/tabula_muris_senis/tabula-muris-senis-facs-processed-official-annotations.h5ad
# mutations from tms
/cellar/users/zkoch/dream/data/tabula_muris_senis/supplementary_table_8.xlsx