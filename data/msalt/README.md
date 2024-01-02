Data from (https://www.sciencedirect.com/science/article/pii/S0092867423004762?via%3Dihub#da0010)
# mice treated (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227358)
curl https://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227358/suppl/GSE227358%5FRNAseq%5Fraw%5Fcounts.csv.gz --output GSE227358%5FRNAseq%5Fraw%5Fcounts.csv.gz
gunzip GSE227358%5FRNAseq%5Fraw%5Fcounts.csv.gz
mv GSE227358%5FRNAseq%5Fraw%5Fcounts.csv treated_mice.csv
# get metadata
curl https://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227358/soft/GSE227358_family.soft.gz --output treated_mice_meta.soft.gz
gunzip treated_mice_meta.soft.gz


# longevity signature across species (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227359)
curl https://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227359/suppl/GSE227359%5FRNAseq%5Fraw%5Fcounts.csv.gz --output GSE227359%5FRNAseq%5Fraw%5Fcounts2.csv.gz
gunzip GSE227359%5FRNAseq%5Fraw%5Fcounts2.csv.gz
mv GSE227359%5FRNAseq%5Fraw%5Fcounts2.csv across_species.csv
# get metadata
curl https://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227359/soft/GSE227359_family.soft.gz --output across_species_meta.soft.gz
gunzip across_species_meta.soft.gz