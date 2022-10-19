library(dplyr)
library(stringr)
##########################################
# extract LDSC results
##########################################
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/causality_analysis/"
common_disease_within_topics <- read.table("/users/mcvean/xilin/xilin/Multimorbidity-biobank/all_disease_topic_list.txt")
hetero_ds <- read.table("/users/mcvean/xilin/xilin/Multimorbidity-biobank/top_hererogeneous_disease.txt")
load("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
# need to add all disease incidence case for comparison
common_disease_within_topics <- common_disease_within_topics %>%
  mutate(V2 = as.character(V2))
for(ds in hetero_ds$V1){
  j <- match(ds, para$list_above500occu$diag_icd10)
  common_disease_within_topics <- common_disease_within_topics %>%
    add_row(V1 = ds, V2 = "all", V3 =  
              mean(para$ds_list[[j]]$age_diag))
}

# everything we want to save
pheh2g <- matrix(NA, nrow = dim(common_disease_within_topics)[1], ncol = 1)
seh2g <- matrix(NA, nrow = dim(common_disease_within_topics)[1], ncol = 1)
for(i in 1:dim(common_disease_within_topics)[1]){
  ds_id <-  common_disease_within_topics$V1[i]
  topic_id <- common_disease_within_topics$V2[i]
  j <- match(ds_id, para$list_above500occu$diag_icd10)
  print(ds_id)
  try({
    # extract heritability for first subtype
    h2g.1 <- readLines(paste(DIR, ds_id,"_topic", topic_id,"/", ds_id, "_topic", topic_id,".h2g.log",sep=""))
    pheh2g[i,1] <- as.numeric(str_split(h2g.1 [length(h2g.1) - 6], "\\s+")[[1]][5])
    se_h2g <- str_split(h2g.1[length(h2g.1) - 6], "\\s+")[[1]][6]
    seh2g[i,1] <- as.numeric(substr(se_h2g, 2,nchar(se_h2g)-1))
  })
}
common_disease_within_topics$h2g <- pheh2g
common_disease_within_topics$seh2g <- seh2g
common_disease_within_topics <- common_disease_within_topics %>% 
  rename(disease = V1, topic = V2, age = V3) %>%
  mutate(z_score = pheh2g/seh2g)
write.csv(common_disease_within_topics ,paste(DIR, "h2g_causality.csv",sep=""), row.names = FALSE )
