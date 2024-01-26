# import scnorm
library(SCnorm)
library(zellkonverter)
library(dplyr)

tms_data_dir <- "/cellar/users/zkoch/dream/data/tabula_muris_senis"
facs_raw_fn <- "tabula-muris-senis-facs-official-raw-obj.h5ad"
# read in h5ad file with zellkonverter
facs_raw <- readH5AD(file.path(tms_data_dir, facs_raw_fn))
# choose only first 100 cells
facs_raw_small <- facs_raw[,1:5000]

# get conditions from tissue, converting to inter
#fake_conditions <- facs_raw$tissue
countDeptEst <- plotCountDepth(
    Data = facs_raw_small, Conditions = as.integer(facs_raw_small$tissue), 
    )
# convert countDeptEst to dataframe
countDeptEst_df <- as.data.frame(countDeptEst)
# describe X1.Slope within each X1.Group in countDeptEst_df
print(countDeptEst_df %>% group_by(X1.Group) %>% summarise(mean = mean(X1.Slope), sd = sd(X1.Slope)))

# do normalization
datanorm <- SCnorm(
    Data = facs_raw_small, 
    Conditions = as.integer(facs_raw_small$tissue),
    #PrintProgressPlots = TRUE,
    FilterCellNum = 10,
    reportSF = TRUE,
    #NCores = 1
    )
# write out normalized data
print("writing out normalized data")
writeH5AD(datanorm, file.path(tms_data_dir, "tabula-muris-senis-facs-official-raw-obj.SCnorm.h5ad"))


# recalc normalization
after_norm_countDeptEst <- plotCountDepth(
    Data = results(datanorm, "NormalizedData"), Conditions =fake_conditions, 
    )
after_norm_countDeptEst <- as.data.frame(after_norm_countDeptEst)
# describe X1.Slope within each X1.Group in countDeptEst_df
print(after_norm_countDeptEst %>% group_by(X1.Group) %>% summarise(mean = mean(X1.Slope), sd = sd(X1.Slope)))
