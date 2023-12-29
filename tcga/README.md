Downloaded TCGABIolinks and associated packages to local RStudio environemnt
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("TCGAbiolinks")
    install.packages('DT')
    BiocManager::install("sesameData")
    BiocManager::install("sesame")
    install.packages("arrow")
    
Then used TCGABiolinks documentation ("https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html")
- This page (https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html) has code for downloading each data type
  - the top of the page describes how to use colData, assay, and rowRanges to use the data that is downloaded (which comes in a summarizedData object)
- This page (https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/query.html) has the different value options for the arguments of query and other fxns


Went to tcga website (with filter: https://portal.gdc.cancer.gov/repository?facetTab=cases&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.samples.tissue_type%22%2C%22value%22%3A%5B%22normal%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.access%22%2C%22value%22%3A%5B%22open%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_category%22%2C%22value%22%3A%5B%22dna%20methylation%22%2C%22simple%20nucleotide%20variation%22%2C%22transcriptome%20profiling%22%5D%7D%7D%5D%7D&searchTableTab=cases)
Filtered for: 
    - tissue: normal
    - access: open
    - data type: gene expression, DNAm, or simple nucleotide variation
    - file type: everything except masked intensities (idat)
Downloaded metadata by clicking "Clinical button"
    - Uploaded by sftp
Downloaded manifest file
