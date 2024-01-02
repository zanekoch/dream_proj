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

TCGA_data_downloader.R: use this to download files

metadata_manual: Directory of manually downloaded TCGA metadata files, becuase the Biolinks prepare fxn does not seem to work well

TCGA_data_reader.py: NOT IMPLEMENTED YET. Use this to read in the data from the TCGA data files (many directories)