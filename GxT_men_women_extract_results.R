# extract topic loading from the target results
# extract age distribution
setwd("/well/mcvean/xilin/Multimorbidity-biobank/")
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)                                  
ds_id <- args[1]
LOC <- paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/GxTopic/", ds_id, "/")

# save hits from two model
GxTopic_men_hit <- c()
GxTopic_women_hit <- c()
for(topic_id in 1:10){
  gxt_hits <- read.table(paste0(LOC, ds_id, "GxTopic_men_topic", topic_id,".assoc.logistic"), header = T) %>%
    filter(TEST == "ADDxCOV1")
  write.table(gxt_hits, paste0(LOC,ds_id, "GxTopic_men_topic", topic_id,".gxtopic"), sep="\t", col.names = T, row.names = FALSE, quote = F)
  GxTopic_men_hit[topic_id] <- gxt_hits %>% 
    filter(P < 5*10^(-8)) %>%
    tally() %>%
    pull
    
  # save stratified method
  gxt_hits <- read.table(paste0(LOC, ds_id, "GxTopic_women_topic", topic_id,".assoc.logistic"), header = T) %>%
    filter(TEST == "ADDxCOV1")
  write.table(gxt_hits, paste0(LOC,ds_id, "GxTopic_women_topic", topic_id,".gxtopic"), sep="\t", col.names = T, row.names = FALSE, quote = F)
  GxTopic_women_hit[topic_id] <- gxt_hits %>% 
    filter(P < 5*10^(-8)) %>%
    tally() %>%
    pull
}
data.frame(topic = 1:10, GxTopic_men_hit = GxTopic_men_hit, 
           GxTopic_women_hit = GxTopic_women_hit) %>% mutate(disease = ds_id) %>%
  write.table(paste0(LOC,ds_id, "_men_women_hits_number.txt"), sep="\t", col.names = T, row.names = FALSE, quote = F)
