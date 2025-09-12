From the github (https://github.com/Pothof-Lab/Transcriptional-Stress/blob/master/Figure4/Main-fig-4_Chang.xlsx) of Gyenis et al. 2023 paper (https://www.nature.com/articles/s41588-022-01279-6#Abs1)
# Eu_20bins.xls
Sheets:
"forward": read counts from nascent RNA-seq in 20 bins for each gene. 6 sets which must correspond to 3 young and 3 old mice, but not sure how. (probably old1-3 then young1-3, because this is how it is in the "Sheet1"
"reverse": same but for reverse strand
"summary_average": the mean value of expression in each bin for each mouse and total diff between old and young

# Total_pol2_20bin.xls
Sheets:
"Totoal_p2_raw": raw read counts of RNAPII chip-seq across 20 gene bodies for each gene. 6 sets which must correspond to 3 young and 3 old mice, but not sure how. (probably old1-3 then young1-3)
    - has gene names in "1700010B13Rik" format
"pol2_log2FC": log2 fold change in each bin for each gene, presumably between old and young mice. 
    - ! HAS gene id to length mapping !
    - only has young vs. old because samples were pooled 
"eu-seq" : seems like log2 fold changes for eu-seq