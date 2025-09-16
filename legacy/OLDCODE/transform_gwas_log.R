library(dplyr)
library(stringr)
# transform the results to log scale 
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/Association_analysis/"
args <- commandArgs(trailingOnly = TRUE)
ds <- args[1]
gwas_rslt <- read.table(paste0(DIR, ds, "/",ds,".assoc.logistic"), header=T)
gwas_rslt$BETA <- log(gwas_rslt$OR)
write.table(gwas_rslt, paste0(DIR, ds, "/",ds,".logOR"), sep = "\t", quote = FALSE,  row.names = FALSE)