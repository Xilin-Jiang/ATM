############################################
# the goal of this file is to:
# (1) create training and testing sets
# (2) estimate multiple CVB estimate of topics and pick the best one
###########################################################
source("topic_functions.R")
args <- commandArgs(trailingOnly = TRUE) # args[1] defines the number of topics; args[2] defines degree of freedom for the curves; args[3] is the repitation id
set.seed( 19940110+as.numeric(args[3]) )

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
degree_free_num <- as.numeric(args[2]) # degrees of freedom


############################
# start optimization
############################
# set the number of inference for getting the best topic estimation
CVB_num <- 5 
topics <- list()
ELBOs <- list()
for(cvb_rep in 1:CVB_num){
  print(paste0("CVB inference number: ", cvb_rep))
  para <- topic_init_age(training_data, ds_list, topic_num, degree_free_num)
  para$rep_ID <- args[3]
  # set the number of update 
  para$max_itr <- 2000
  para$alpha <- rep(1, para$K) 
  para$lb <- data.frame("Iteration" = 0,"Lower_bound" = CVB_lb(para))
  para$itr_beta <- 1
  para$itr_check_lb <- 5
  para$itr_save <- para$max_itr + 1 # don't need to save intermediate data
  para$tol <- 10^(-6)
  for(itr in 1:para$max_itr){
    print(paste0("Interation: ",itr))
    for(itr_inside in 1:para$itr_beta){ # in practice we perform quick steps a few times before move on.
      para <- CVB0_E_zn(para) # we choose CVB0 as papers shown it could converge quicker
      # para <- comp_E_lntheta(para)
    }
    para <- fast_update_age_depend_lda(para)
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
      save(para, file = paste0("Results/","training_2rec_age_K",para$K,"_P",para$P,"_rep",para$rep_ID,"_intermediate",itr, ".RData"))
    }
  }
  topics[[cvb_rep]] <- para$pi_beta_basis
  ELBOs[[cvb_rep]]  <- para$lb
  save(para, file = paste0("Results/","example_training_set_estimation.RData"))
}
multi_runs <- list(topics, ELBOs)
save(multi_runs, file = paste0("Results/","training_multiruns_ageLDA_RunNumber", CVB_num, "K",para$K,"_P",para$P,"_rep",para$rep_ID, ".RData"))
# find the best reps in the data
lb_rep <- data_frame(reps = as.integer(), lower_bound = as.numeric())
for(cvb_rep in 1:CVB_num){
  cvrg_lb <-  ELBOs[[cvb_rep]] %>% 
    filter(!is.na(Lower_bound)) %>%
    slice_tail %>% 
    pull(2)
  lb_rep <- lb_rep %>% 
    add_row(reps = cvb_rep, lower_bound = cvrg_lb)
}
best_id <- order(lb_rep$lower_bound, decreasing = T)[1]

# save a smaller dataset
model_output <- list(topics[[best_id]], ELBOs[[best_id]], para$D, para$M, para$K, para$P, testing_eid, lb_rep)
# need to save the testing_eid for analysis
save(model_output, file = paste0("Results/","training_best_output_AgeLDA_RunNumber", CVB_num, "K",para$K,"_P",para$P,"_rep",para$rep_ID, ".RData"))



