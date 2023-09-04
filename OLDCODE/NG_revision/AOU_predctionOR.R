# perform within population prediction
library(dplyr)
library(ATM)
########################################
########################################
args <- commandArgs(trailingOnly = TRUE) # args[1] defines the number of topics; args[2] is the repitation id
set.seed(19940110)

topic_number <- as.numeric(args[1]) 
cv_id  <- as.numeric(args[2]) 
# bootstrap samples for UKBB

print(paste0("running prediction OR for topic number: ", topic_number))

AOU_comorbidity_data <- read.csv("AOU_rec_data.csv")

cv_num <- 5
print(paste0("cross validation block: ", cv_id))

# df_ATM <- AOU_comorbidity_data # please use the data! XXXXX
# testing_ids <-  df_ATM %>% 
#   group_by(eid) %>% 
#   slice(1) %>% 
#   ungroup()
# fold_size <- floor(dim(testing_ids)[1]/cv_num) 
# if(cv_id == 5){
#   testing_ids <- testing_ids%>% 
#     slice(( (cv_id-1) * fold_size + 1): dim(testing_ids)[1])
# }else{
#   testing_ids <- testing_ids%>%  
#     slice(( (cv_id-1) * fold_size + 1): (cv_id * fold_size))
# }
# 
# training_data <- df_ATM %>% 
#   anti_join(testing_ids, by = "eid")
# testing_data <- df_ATM %>% semi_join(testing_ids, by = "eid")
# ATM_results <- wrapper_ATM(rec_data=training_data, topic_num = topic_number, CVB_num = 5)
# print(paste0("training finished"))
# save(ATM_results, file = paste0("../Results_NG_revision/", "AOU_cross_validation_training", topic_number, "cv", cv_id, ".Rdata") )

###############################################################
# the code has been modified by directly loading pre-trained ATM
# remove this block and uncomment to return to the original training
load(paste0("../Results_NG_revision/", "AOU_cross_validation_training", topic_number, "cv", cv_id, ".Rdata") )
testing_data <- AOU_comorbidity_data %>% 
  filter(!eid %in% ATM_results$patient_list)
###############################################################
###############################################################
print(paste0("Starting prediction OR"))
AOU_prediction_OR <- prediction_OR(testing_data = testing_data, ds_list = ATM_results$ds_list, topic_loadings = ATM_results$topic_loadings)
save(AOU_prediction_OR, file = paste0("../Results_NG_revision/", "AOU_predictionOR_topic_num", topic_number, "cv", cv_id, ".Rdata") )
