library(dplyr)
library(stringr)
##########################################
# extract LDSC results
##########################################
args <- commandArgs(trailingOnly = TRUE)  
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/Association_analysis/"
ds_id <- args[1]
rep_id <- as.numeric(args[2])
set.seed(19940110 + rep_id)
if(rep_id == 1){
  permutation_Fst <- data.frame(weighted_Fst = as.numeric())
  write.csv(permutation_Fst, paste(DIR, ds_id,"/",ds_id,"_control_Fst.csv",sep=""), row.names = F)
}else{
  permutation_Fst <- read.csv(paste(DIR, ds_id,"/",ds_id,"_control_Fst.csv",sep=""))
  fst <- readLines(paste(DIR, ds_id,"/",ds_id,"_control_Fst.log",sep=""))
  permutation_Fst <- permutation_Fst %>%
    add_row(weighted_Fst = as.numeric(str_split(fst[length(fst) - 2], "\\s+")[[1]][4]))
  write.csv(permutation_Fst, paste(DIR, ds_id,"/",ds_id,"_control_Fst.csv",sep=""), row.names = F)
}
# first compute the range of topics weights
load("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
subtypes <- read.table(paste(DIR, ds_id,"/",ds_id,"subtypes.txt",sep=""), header = F)

# subtypes <- read.table("Association_analysis/153.2subtypes.txt", header = F)

topic_list <- unique(subtypes$V3)
patient_loadings <- sweep((para$alpha_z - 1), 1, rowSums(para$alpha_z -1), FUN="/") 
df_controls <- data.frame(eid = para$eid, loadings = patient_loadings[,topic_list]) %>%
  anti_join(subtypes, by=c("eid" = "V1"))
df_cases <- data.frame(eid = para$eid, loadings = patient_loadings[,topic_list]) %>%
  semi_join(subtypes, by=c("eid" = "V1")) %>%
  left_join(select(subtypes, V1, V3), by=c("eid" = "V1"))
control_Fst <- list() 
for(loc_id in 1:length(topic_list)){
  topic_id <- topic_list[loc_id]
  range_df <- df_cases %>%
    filter(V3 == topic_id) %>%
    pull(!!paste0("loadings.",loc_id))
  topic_quantiles <- quantile(range_df,c(0,0.25,0.5,0.75,1))
  size <- as.integer(length(range_df)/4)
  controls_topic <- list()
  controls_topic[[1]] <- df_controls %>% 
    filter(!!as.name(paste0("loadings.",loc_id))  > topic_quantiles[1] , !!as.name(paste0("loadings.",loc_id))  < topic_quantiles[2]) %>%
    sample_n(size)
  controls_topic[[2]]  <- df_controls %>% 
    filter(!!as.name(paste0("loadings.",loc_id))  > topic_quantiles[2] , !!as.name(paste0("loadings.",loc_id))  < topic_quantiles[3]) %>%
    sample_n(size)
  controls_topic[[3]]  <- df_controls %>% 
    filter(!!as.name(paste0("loadings.",loc_id))  > topic_quantiles[3] , !!as.name(paste0("loadings.",loc_id))  < topic_quantiles[4]) %>%
    sample_n(size)
  controls_topic[[4]]  <- df_controls %>% 
    filter(!!as.name(paste0("loadings.",loc_id))  > topic_quantiles[4] , !!as.name(paste0("loadings.",loc_id))  < topic_quantiles[5]) %>%
    sample_n(size)
  control_Fst[[loc_id]] <- bind_rows(controls_topic) %>%
    select(eid) %>%
    mutate(FIID = eid, topic = topic_id)
  # make sure there is no overlapping
  df_controls <- df_controls %>%
    anti_join(control_Fst[[loc_id]], by=c("eid"))
}
control_Fst <- bind_rows(control_Fst)
write.table(control_Fst, paste(DIR, ds_id,"/",ds_id,"_control_subtypes.txt",sep=""), sep="\t", col.names = F, row.names = FALSE ,quote = F)
