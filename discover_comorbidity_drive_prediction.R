source("topic_functions.R")
args <- commandArgs(trailingOnly = TRUE) # args[1] defines the number of topics; args[2] defines degree of freedom for the curves; args[3] is the repitation id

rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
# need to include survival age for PheRS
survive_age <- read.csv("survive_age.csv")

topic_num <- as.numeric(args[1]) 
degree_free_num <- as.numeric(args[2]) # degrees of freedom
rep_ID <- args[3]

# save the parameter
load(file = paste0("Results/","training_2rec_age_K",topic_num,"_P",degree_free_num,"_rep",rep_ID, ".RData"))

############################################
# preparing testing data: extract the training
############################################
all_eid <- rec_data %>%
  group_by(eid) %>%
  summarise() 
training_eid <- data.frame(eid = para$eid)
testing_eid <- all_eid %>% 
  anti_join(training_eid, by = "eid")

testing_data <- testing_eid %>%  
  left_join(rec_data, by = "eid")


# make predictions for single disease
PheRS_rslt <- prediction_PheRS_by_disease(testing_data, ds_list, para) 
# save the PheRS for each diseases
save(PheRS_rslt, file = paste0("Results/","single_diseae_prediction",para$K,"_P",para$P,"_rep",para$rep_ID, ".RData"))

max_predict <- 11 # predict the records until which step?
prediction_OR_rslt <- prediction_OR_onebyone(testing_data, ds_list, para, max_predict) 

# save the parameter
save(prediction_OR_rslt, file = paste0("Results/","OR_each_disease_K",para$K,"_P",para$P,"_rep",para$rep_ID, ".RData"))

# perform another task which saves a much larger set of data: how does the model make inference
save_decisions_data <- save_prediction_logics(testing_data, ds_list, para, max_predict)

# first step: save prediction accuracy data
collect_prediction_ranks <- save_decisions_data[[3]]
save(collect_prediction_ranks, file = paste0("Results/","rank_of_target_disease",para$K,"_P",para$P,"_rep",para$rep_ID, ".RData"))

# second step: save data connected to the decision procedure
collect_zn <- save_decisions_data[[1]]
collect_estimate_set <- save_decisions_data[[2]]
testing_data <- save_decisions_data[[4]]
# we also need the beta info to infer the posterior of the target disease
trajectory_data <- list(collect_zn, collect_estimate_set, testing_data, para$pi_beta_basis)
save(trajectory_data, file = paste0("Results/","prediction_trajectory_data_K",para$K,"_P",para$P,"_rep",para$rep_ID, ".RData"))




