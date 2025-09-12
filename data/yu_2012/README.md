### GSE39313_family.soft
Metadata and datafile (need to parse expression from here bc no unified df) downloaded from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39313) by clicking the .soft in download family
### GSE39313_expression.csv
Created by parsing expression from GSE39313_family.soft by removing the sequencing stuff up top and then running `grep '^ILMN' GSE39313_family.soft > GSE39313_expression.csv`