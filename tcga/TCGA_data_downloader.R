library(TCGAbiolinks)
library(DT)


query_met.hg38 <- GDCquery(
  project= "TCGA-LGG", 
  data.category = "DNA Methylation", 
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450", 
  barcode = c("TCGA-HT-8111-01A-11D-2399-05","TCGA-HT-A5R5-01A-11D-A28N-05")
)
GDCdownload(query_met.hg38)
data.hg38 <- GDCprepare(query_met.hg38)
