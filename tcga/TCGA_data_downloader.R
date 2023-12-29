library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)

output_dir <- "/Users/zanekoch/Documents/ucsd/dream_proj/tcga/GDCdata"

datasets <- c("TCGA-BRCA", "TCGA-COAD")
data_categories <- c("Transcriptome Profiling", "DNA Methlation", "Simple Nucleotide Variation")
data_types  <- c("Gene Expression Quantification", "Methylation Beta Value", "Masked Somatic Mutation")
normal_sample_types <- c("Blood Derived Normal", "Solid Tissue Normal", "Buccal Cell Normal", "EBV Immortalized Normal", "Bone Marrow Normal")

# if methylation then platform = "Illumina Human Methylation 450" or "Illumina Methylation Epic"


# Loop through datasets, data categories, and data types
for (dataset in datasets) {
  for (category in data_categories) {
    # select dataypes
    if (category == "Transcriptome Profiling") {
      dtype <- data_types[1]
    } else if (category == "DNA Methlation") {
      dtype <- data_types[2]
      #platform <- c("Illumina Human Methylation 450", "Illumina Methylation Epic")
    } else { # simple nucleotide variation
      dtype <- data_types[3]
    }
    # Download data
    query <- GDCquery(project = dataset, 
                      data.category = category, 
                      data.type = dtype,
                      sample.type = normal_sample_types,
                      )
    
    GDCdownload(query, directory = output_dir)
    
    # Load the data
    data <- GDCprepare(query, directory = output_dir)
    
    # Display some information about the downloaded data
    cat("Dataset:", dataset, "\n")
    cat("Data Category:", category, "\n")
    cat("Data Type:", dtype, "\n")
    cat("Number of samples:", ncol(colData(data)), "\n")
    cat("\n")
  }
}





query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  experimental.strategy = "RNA-Seq", 
                  sample.type = "Solid Tissue Normal")
query$results

# choose the samples to download
query_met.hg38 <- GDCquery(
  project= "CPTAC-3", 
  data.category = "Transcriptome Profiling", 
  #data.type = "Methylation Beta Value",
  #platform = "Illumina Human Methylation 450", 
  #barcode = c("TCGA-HT-8111-01A-11D-2399-05","TCGA-HT-A5R5-01A-11D-A28N-05"),
  #sample.type = "Solid Tissue Normal"
  #sample.type = c("NB", "NT", "NBC", "NEBV", "NBM")
  sample.type = c("Blood Derived Normal", "Solid Tissue Normal", "Buccal Cell Normal", "EBV Immortalized Normal", "Bone Marrow Normal")
)
query_met.hg38$results

# download them to GDCdata dir
GDCdownload(query_met.hg38, directory = "/Users/zanekoch/Documents/ucsd/dream_proj/tcga/GDCdata")
# can  set summarizedExperiment to FALSE to get a df
data.hg38 <- GDCprepare(
  query_met.hg38, summarizedExperiment = FALSE,
  directory = "/Users/zanekoch/Documents/ucsd/dream_proj/tcga/GDCdata"
  )

colData(data.hg38)
# data matrix
assay(data.hg38)
# not sure what this is 
rowRanges(data.hg38)
