# extract topic loading from the target results
# extract age distribution
setwd("/well/mcvean/xilin/Multimorbidity-biobank/")
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)                                  
ds_id <- args[1]
LOC <- paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/GxTopic/", ds_id, "/")

# save hits from two model
StratifyDisease_hit <- c()
GxTopic_hit <- c()
for(topic_id in 1:10){
  gxt_hits <- read.table(paste0(LOC, ds_id, "GxTopic_topic", topic_id,".assoc.logistic"), header = T) %>%
    filter(TEST == "ADDxCOV1")
  write.table(gxt_hits, paste0(LOC,ds_id, "GxTopic_topic", topic_id,".gxtopic"), sep="\t", col.names = T, row.names = FALSE, quote = F)
  GxTopic_hit[topic_id] <- gxt_hits %>% 
    filter(P < 5*10^(-8)) %>%
    tally() %>%
    pull
    
  # save stratified method
  strf_hits <- read.table(paste0(LOC, ds_id, "GxT_DiseaseStratify_topic", topic_id,".qassoc.gxe"), header = T) 
  StratifyDisease_hit[topic_id] <- strf_hits %>% 
    filter(P_GXE < 5*10^(-8)) %>%
    tally() %>%
    pull
}
data.frame(topic = 1:10, StratifyDisease_hit = StratifyDisease_hit, 
           GxTopic_hit = GxTopic_hit) %>% mutate(disease = ds_id) %>%
  write.table(paste0(LOC,ds_id, "_hits_number.txt"), sep="\t", col.names = T, row.names = FALSE, quote = F)
