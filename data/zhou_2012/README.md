### GSE36838_family.soft.gz
Metadata and datafile (need to parse expression from here bc no unified df) downloaded from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) by clicking the .soft in download family
### GSE36836_family.soft.gz 
Same but from (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi)

### GSE36838_expression.csv
Created by parsing expression from GSE36838_family.soft by removing the sequencing stuff up top and then running `zgrep -P '^\w+_at\b' GSE36838_noGeneTable.soft.gz >> GSE36838_expression.csv`
### GSE36836_expression.csv
Ditto


