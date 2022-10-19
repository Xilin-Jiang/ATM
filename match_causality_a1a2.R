# match A1/A2 for the logistica regression -- this is for LDSC rg analysis
require(dplyr)
library(stringr)
args <- commandArgs(trailingOnly = TRUE)                                  
ds_id <- args[1]
topic_id <- args[2]
setwd("/well/mcvean/xilin/Multimorbidity-biobank/")
DIR <- "/well/mcvean/xilin/Multimorbidity-biobank/causality_analysis/"
# using the function of package pdist
simple_assoc <- read.table(paste0(DIR,  ds_id,"_topic", topic_id,"/", ds_id,"_topic", topic_id, ".assoc"), header = T)
logistic_assoc <- read.table(paste0(DIR, ds_id,"_topic", topic_id,"/", ds_id,"_topic", topic_id,".assoc.logistic"), header = T)

logistic_a1_a2 <- logistic_assoc %>%
  left_join(select(simple_assoc, SNP, A1, A2), by = c("SNP", "A1"))
# save a file with A2
write.table(logistic_a1_a2, paste(DIR, ds_id,"_topic", topic_id,"/", ds_id,"_topic", topic_id, ".logistic.a2",sep=""), sep="\t", row.names = FALSE ,quote = F)




