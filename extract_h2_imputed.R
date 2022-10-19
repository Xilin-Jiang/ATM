library(dplyr)
library(stringr)
##########################################
# extract LDSC results
##########################################
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/LDSC_SEG/"
common_disease_within_topics <- read.table("/users/mcvean/xilin/xilin/Multimorbidity-biobank/BOLT_LMM_subtype_list.txt")
# need to add all disease incidence case for comparison
common_disease_within_topics <- common_disease_within_topics 

# everything we want to save
pheh2g <- matrix(NA, nrow = dim(common_disease_within_topics)[1], ncol = 1)
seh2g <- matrix(NA, nrow = dim(common_disease_within_topics)[1], ncol = 1)
for(i in 1:dim(common_disease_within_topics)[1]){
  ds_id <-  common_disease_within_topics$V1[i]
  topic_id <- common_disease_within_topics$V2[i]
  print(ds_id)
  try({
    # extract heritability for first subtype
    h2g.1 <- readLines(paste(DIR, ds_id,"_topic", topic_id,"/", ds_id, "_topic", topic_id,"_imputed.h2g.log",sep=""))
    pheh2g[i,1] <- as.numeric(str_split(h2g.1[length(h2g.1) - 6], "\\s+")[[1]][5])
    se_h2g <- str_split(h2g.1[length(h2g.1) - 6], "\\s+")[[1]][6]
    seh2g[i,1] <- as.numeric(substr(se_h2g, 2,nchar(se_h2g)-1))
  })
}
common_disease_within_topics$h2g <- pheh2g
common_disease_within_topics$seh2g <- seh2g
common_disease_within_topics <- common_disease_within_topics %>% 
  rename(disease = V1, topic = V2, age = V3, N = V4) %>%
  mutate(z_score = pheh2g/seh2g)
write.csv(common_disease_within_topics ,paste(DIR, "h2g_imputed.csv",sep=""), row.names = FALSE )

# elso extract the h2g for the topics
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/causality_analysis/"
# everything we want to save
num_topic <- 10
pheh2g <- matrix(NA, nrow = num_topic, ncol = 1)
seh2g <- matrix(NA, nrow = num_topic, ncol = 1)
for(topic_id in 1:num_topic){
  print(topic_id)
  # extract heritability for first subtype
  h2g.1 <- readLines(paste(DIR, "topic", topic_id,"/topic", topic_id, ".h2g.log",sep=""))
  pheh2g[topic_id,1] <- as.numeric(str_split(h2g.1[length(h2g.1) - 6], "\\s+")[[1]][5])
  se_h2g <- str_split(h2g.1[length(h2g.1) - 6], "\\s+")[[1]][6]
  seh2g[topic_id,1] <- as.numeric(substr(se_h2g, 2,nchar(se_h2g)-1))
}
topic_h2g <- data.frame(topic_id = 1:num_topic)
topic_h2g$h2g <- pheh2g
topic_h2g$seh2g <- seh2g
topic_h2g <- topic_h2g %>% 
  mutate(z_score = pheh2g/seh2g)
write.csv(topic_h2g ,paste(DIR, "topic_h2g.csv",sep=""), row.names = FALSE )



