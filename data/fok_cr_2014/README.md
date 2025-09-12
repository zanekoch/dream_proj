### GSE40977_family.soft
Metadata and datafile (need to parse expression from here bc no unified df) downloaded from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40977) by clicking the .soft in download family
### GSE40977_expression.csv
Created by parsing expression from GSE40977_family.soft by removing the sequencing stuff up top and then running `grep '^ILMN' GSE40977_family_noSeqStuff.soft > GSE40977_expression.csv`