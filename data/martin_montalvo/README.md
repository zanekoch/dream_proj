### GSE40936_family.soft
Metadata and datafile (need to parse expression from here bc no unified df) downloaded from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40936) by clicking the .soft in download family
### GSE40936_expression.csv
Created by parsing expression from GSE40936_family.soft by removing the sequencing stuff up top and then running `grep '^ILMN' GSE40936_family.soft > GSE40936_expression.csv`