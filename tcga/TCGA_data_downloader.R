library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)

output_dir <- "/Users/zanekoch/Documents/ucsd/dream_proj/tcga/GDCdata"

datasets <- c("TCGA-COAD", "TCGA-BRCA")
data_categories <- c("DNA Methylation","Transcriptome Profiling",  "Simple Nucleotide Variation")
normal_sample_types <- c("Blood Derived Normal", "Solid Tissue Normal", "Buccal Cell Normal", "EBV Immortalized Normal", "Bone Marrow Normal")

# if methylation then platform = "Illumina Human Methylation 450" or "Illumina Methylation Epic"


# Loop through datasets, data categories, and data types
for (dataset in datasets) {
  for (category in data_categories) {
    # select dataypes
    if (category == "Transcriptome Profiling") {
      dtype <- "Gene Expression Quantification"
      # query
      query <- GDCquery(project = dataset, 
                        data.category = category, 
                        data.type = dtype,
                        sample.type = normal_sample_types,
                        )
    } else if (category == "DNA Methylation") {
      dtype <- "Methylation Beta Value"
      # TCGA is all 450k/27k but we are only going to use 450k
      if (startsWith(dataset, "TCGA")) {
        platform <- "Illumina Human Methylation 450"
      # while CPTAC and someothers are Epic
      } else {
        platform <- "Illumina Methylation Epic"
      }
      # query with platform
      query <- GDCquery(project = dataset, 
                        data.category = category, 
                        data.type = dtype,
                        sample.type = normal_sample_types,
                        platform = platform
                        )
    } else { # simple nucleotide variation
      dtype <- "Masked Somatic Mutation"
      # query
      query <- GDCquery(project = dataset, 
                        data.category = category, 
                        data.type = dtype,
                        sample.type = normal_sample_types,
                        )
    }
    
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


# CPTAC-3 mutations
query <- GDCquery(
  project = "CPTAC-3",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation", 
  #workflow.type = "STAR - Counts",
  #barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01")
)
View(getResults(query))

table(getResults(query)$sample_type)

# CPTAC-3 expr
query1 <- GDCquery(
  project = "CPTAC-3",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  #workflow.type = "STAR - Counts",
  #barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01")
)
table(getResults(query1)$sample_type)

# CPTAC-3 methyl
query2 <- GDCquery(
  project = "CPTAC-3",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value", 
  #workflow.type = "STAR - Counts",
  #barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01")
)
table(getResults(query2)$sample_type)


