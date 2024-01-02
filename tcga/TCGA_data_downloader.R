library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
# set download timeout limit to 100 minutes
options(timeout = 6000)

#datasets <- c("TCGA-COAD", "TCGA-BRCA")
#data_categories <- c("DNA Methylation","Transcriptome Profiling",  "Simple Nucleotide Variation")
#normal_sample_types <- c("Blood Derived Normal", "Solid Tissue Normal", "Buccal Cell Normal", "EBV Immortalized Normal", "Bone Marrow Normal")



# write a function that takes in a project and returns the samples with both mutation and expression data
# and then download the data for those samples
#get_mut_expr_samples <- function(project, output_dir="/cellar/users/zkoch/dream/tcga/GDCdata") {
  # get mutation data
output_dir<-"/cellar/users/zkoch/dream/tcga/GDCdata"
project <- "CPTAC-3"
query_mut <- GDCquery(
  project = project,
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation", 
  access = "open"
)
# get expression data
query_expr <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  access = "open",
)
# print how many samples of each type
print(paste("Mutation samples:", length(getResults(query_mut, cols = "cases"))))
print(paste("Expression samples:", length(getResults(query_expr, cols = "cases"))))
shared_samples <- intersect(
    substr(getResults(query_mut, cols = "cases"), 1, 12),
    substr(getResults(query_expr, cols = "cases"), 1, 12)
)
print(paste("Samples with both mutation and expression data:", length(shared_samples)))

GDCdownload(query_mut, directory = output_dir)
GDCdownload(query_expr, directory = output_dir)
print("Downloaded data")
# Load the data as a SummarizedExperiment
data_mut <- GDCprepare(query_mut, directory = output_dir, summarizedExperiment = TRUE)
data_expr <- GDCprepare(query_expr, directory = output_dir, summarizedExperiment = TRUE)
# return the data
return(list(data_mut, data_expr))
  #}

# run the function
data <- get_mut_expr_samples("CPTAC-3")
