source("topic_functions.R")
args <- commandArgs(trailingOnly = TRUE) # args[1] defines the number of topics; args[2] is the repitation id

rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")

all_eid <- rec_data %>%
  group_by(eid) %>%
  summarise() 
training_eid <- all_eid %>% 
  sample_frac(size  = .8)
testing_eid <- all_eid %>% 
  anti_join(training_eid, by = "eid")

training_data <- training_eid %>%  
  left_join(rec_data, by = "eid")

testing_data <- testing_eid %>%  
  left_join(rec_data, by = "eid")

ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
topic_num <- as.numeric(args[1]) 
para <- topic_init_baseline(training_data, ds_list, topic_num)
para$rep_ID <- args[2]

############################
# start optimization
############################
# set the number of update 
para$max_itr <- 2000
para$alpha <- rep(1, para$K) 
para$lb <- data.frame("Iteration" = 0,"Lower_bound" = CVB_lb(para))
para$itr_check_lb <- 5
para$itr_save <- para$max_itr # don't need to save intermediate data
para$tol <- 10^(-6)
for(itr in 1:para$max_itr){
  print(paste0("Interation: ",itr))
  para <- CVB0_E_zn(para) # we choose CVB0 as papers shown it could converge quicker
  para <- update_beta_basic_lda(para)
  # para <- update_alpha(para) # we use non-informative alpha
  para$lb[nrow(para$lb) + 1,] <- c(itr, CVB_lb(para))
  if(itr %% para$itr_check_lb  ==0){
    curr_lb <- pull(filter(para$lb, Iteration == itr), Lower_bound) 
    prev_lb <- pull(filter(para$lb, Iteration == (itr -para$itr_check_lb )), Lower_bound) 
    print(paste0("Current Lower bound ", curr_lb, " at iteration: ",itr))
    try({
      if(is.finite((curr_lb - prev_lb)) & (curr_lb - prev_lb)/abs(prev_lb) < para$tol ){
        print(paste0("Optimization converged at step ", itr))
        break
      }
    })
  }
  if(itr %% para$itr_save ==0){
    save(para, file = paste0("Results/","training_baseline_K",para$K,"_rep",para$rep_ID,"_intermediate",itr, ".RData"))
  }
}

# save the parameter
save(para, file = paste0("Results/","training_baseline_K",para$K,"_rep",para$rep_ID, ".RData"))

# save a smaller dataset
order_by_post_z <- order(colMeans(para$alpha_z), decreasing = T)
ordered_beta <- para$beta[,order_by_post_z]
# save the ordered version of model_output
model_output <- list(ordered_beta, para$lb, para$D, para$M, para$K)
print(object.size(model_output), unit = "MB" , standard = "SI")
save(model_output, file = paste0("Results/","training_model_output_baseline_K",para$K,"_rep",para$rep_ID, ".RData"))

##########################
# perform prediction
##########################
max_predict <- 11 # predict the records until which step?
prediction_onebyone_rslt <- prediction_onebyone(testing_data, ds_list, para, max_predict) 
# save the parameter
save(prediction_onebyone_rslt, file = paste0("Results/","prediction_onebyone_baseline_K",para$K,"_rep",para$rep_ID, ".RData"))

# perform another task which saves a much larger set of data: how does the model make inference
save_decisions_data <- save_prediction_logics(testing_data, ds_list, para, max_predict)

# first step: save prediction accuracy data
collect_prediction_ranks <- save_decisions_data[[3]]
save(collect_prediction_ranks, file = paste0("Results/","baselineLDA_rank_of_target_disease",para$K,"_rep",para$rep_ID, ".RData"))





