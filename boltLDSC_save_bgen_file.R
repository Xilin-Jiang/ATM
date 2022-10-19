library(dplyr)
args <- commandArgs(trailingOnly = TRUE)                                  
ds_id <- as.numeric(args[1])
topic_id <- args[2]

bgen_results <- read.table(paste0(ds_id,"_topic", topic_id,"/",ds_id,"_topic", topic_id, ".bolt.bgen.stat"), header = T, fill=TRUE, na.strings = c("", "NA") ) %>% 
  stats::na.omit() 

lmm_type <- names(bgen_results)[dim(bgen_results)[2]]
info_thre <- 0.9

if(lmm_type == "P_BOLT_LMM_INF" ){
  bgen_results <- bgen_results %>%
    rename(P = P_BOLT_LMM_INF, A1 = ALLELE1, A2 = ALLELE0) %>%
    filter(INFO > info_thre) 
}else{
  bgen_results <- bgen_results %>%
    rename(P = P_BOLT_LMM, A1 = ALLELE1, A2 = ALLELE0) %>%
    filter(INFO > info_thre)
}

bgen_results %>%
  mutate(INFO = as.numeric(INFO), BETA = as.numeric(BETA), SE = as.numeric(SE), P = as.numeric(P)) %>%
  stats::na.omit() %>%
  write.table(paste0(ds_id,"_topic", topic_id,"/",ds_id,"_topic", topic_id, ".filter_INFO.bolt.bgen.stat"), row.names = F, col.names = T, sep = "\t", quote = F)