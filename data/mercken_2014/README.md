### GSE49000_family.soft
Metadata and datafile (need to parse expression from here bc no unified df) downloaded from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49000) by clicking the .soft in download family
### GSE49000_expression.csv
Created by parsing expression from GSE49000_family.soft by removing the sequencing stuff up top and then running `grep '^ILMN' GSE49000_family.soft > GSE49000_expression.csv`