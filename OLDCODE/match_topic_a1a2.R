# match A1/A2 for the logistica regression -- this is for LDSC rg analysis
require(dplyr)
library(stringr)
args <- commandArgs(trailingOnly = TRUE)                                  
topic_id <- args[1]
setwd("/well/mcvean/xilin/Multimorbidity-biobank/")
DIR <- "/well/mcvean/xilin/Multimorbidity-biobank/Association_analysis/"
# using the function of package pdist
ds_id <- "250.2"
simple_assoc <- read.table(paste0(DIR, ds_id, "/", ds_id, ".assoc"), header = T)
topic_assoc <- read.table(paste0(DIR, "association_results/topic", topic_id,"_loading_K10_rep10.assoc.linear"), header = T)
topic_a1_a2 <- topic_assoc %>%
  left_join(select(simple_assoc, SNP, A1, A2), by = c("SNP", "A1"))
# save a file with A2
write.table(topic_a1_a2, paste0(DIR, "association_results/topic", topic_id,"_loading.linear.a2"), sep="\t", row.names = FALSE ,quote = F)




