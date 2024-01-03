### Usage ###
# Active r_go conda environment 
# Start R by typing R in the terminal
# run through this script line by line with command + enter

library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
# set download timeout limit to 100 minutes
options(timeout = 6000)

# get mutation data
output_dir<-"/cellar/users/zkoch/dream/data/tcga/GDCdata"
project <- "CPTAC-3"
query_mut <- GDCquery(
  project = project,
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation", 
  access = "open"
  # need to add workflow type???
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
print(paste("Expression samples:", length(table(getResults(query_expr, cols = "cases")))))
shared_samples <- intersect(
    substr(getResults(query_mut, cols = "cases"), 1, 12),
    substr(getResults(query_expr, cols = "cases"), 1, 12)
)
print(paste("Samples with both mutation and expression data:", length(shared_samples)))
# save file to sample mappings
file_mapping_mut = getResults(query_mut)
file_mapping_expr = getResults(query_expr)
#write.csv(file_mapping_mut, file = "/cellar/users/zkoch/dream/data/tcga/GDCdata/CPTAC-3/harmonized/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/file_mapping_mut.csv")
#write.csv(file_mapping_expr, file = "/cellar/users/zkoch/dream/data/tcga/GDCdata/CPTAC-3/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/file_mapping_expr.csv")

GDCdownload(query_mut, directory = output_dir)
GDCdownload(query_expr, directory = output_dir)
print("Downloaded data")
# Load the data as a SummarizedExperiment
data_mut <- GDCprepare(query_mut, directory = output_dir, summarizedExperiment = FALSE)
data_expr <- GDCprepare(query_expr, directory = output_dir, summarizedExperiment = FALSE)
# return the data
return(list(data_mut, data_expr))
  #}

# run the function
data <- get_mut_expr_samples("CPTAC-3")

tab <-  getSampleFilesSummary(project = "CPTAC-3")
