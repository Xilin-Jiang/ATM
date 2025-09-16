source("topic_functions.R")
args <- commandArgs(trailingOnly = TRUE) # args[1] defines the number of topics; args[2] defines degree of freedom for the curves; args[3] is the repitation id

rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")

ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
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

max_predict <- 11 # predict the records until which step?
prediction_onebyone_rslt <- prediction_onebyone(testing_data, ds_list, para, max_predict) 

# save the parameter
save(prediction_onebyone_rslt, file = paste0("Results/","prediction_onebyone_age_K",para$K,"_P",para$P,"_rep",para$rep_ID, ".RData"))

