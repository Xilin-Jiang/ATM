source("topic_functions.R")
args <- commandArgs(trailingOnly = TRUE) # args[1] defines the number of topics; args[2] defines degree of freedom for the curves; args[3] is the repitation id

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
para <- topic_init_age(training_data, ds_list, topic_num, degree_free_num)
para$rep_ID <- args[3]

############################
# start optimization
############################
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

# save the parameter
save(para, file = paste0("Results/","training_2rec_age_K",para$K,"_P",para$P,"_rep",para$rep_ID, ".RData"))

# save a smaller dataset
order_by_post_z <- order(colMeans(para$alpha_z), decreasing = T)
ordered_pi_beta <- para$pi_beta_basis[,,order_by_post_z]
# save the ordered version of model_output
model_output <- list(ordered_pi_beta, para$lb, para$D, para$M, para$K, para$P)
print(object.size(model_output), unit = "MB" , standard = "SI")
save(model_output, file = paste0("Results/","training_model_output_2rec_age_K",para$K,"_P",para$P,"_rep",para$rep_ID, ".RData"))

############################################
# compute the likelihood for the testing set
############################################
# first order the incidences by age
testing_data <- testing_data %>%
  group_by(eid) %>%
  arrange(age_diag, .by_group = TRUE)

# choose the first half of individuals 
predicting_test_set <- testing_data %>%
  group_by(eid) %>%
  slice_tail(prop = .5) %>%
  ungroup

estimating_test_set <- testing_data %>%
  anti_join(predicting_test_set, by = c("eid", "diag_icd10")) %>%
  ungroup

para_testing <- topic_init_age(estimating_test_set, ds_list, topic_num, degree_free_num)
# assigning beta to the testing set
para_testing$beta_w_full <- apply(para$pi_beta_basis, 3, function(x) 
  x[as.matrix(select(para_testing$unlist_Ds_id, age_diag, Ds_id))]) 
para_testing$beta_w <- lapply(para_testing$patient_lst, function(x) para_testing$beta_w_full[x,,drop=F])

# updating Z_n
para_testing$max_itr <- 10
para_testing$alpha <- rep(1, para_testing$K) 
para_testing$lb <- data.frame("Iteration" = 0,"Lower_bound" = CVB_lb(para_testing))
para_testing$tol <- 10^(-6)
for(itr in 1:para_testing$max_itr){
  print(paste0("Interation: ",itr))
  para_testing <- CVB0_E_zn(para_testing) # we choose CVB0 as papers shown it could converge quicker
  para_testing$lb[nrow(para_testing$lb) + 1,] <- c(itr, CVB_lb(para_testing))
  curr_lb <- pull(filter(para_testing$lb, Iteration == itr), Lower_bound) 
  prev_lb <- pull(filter(para_testing$lb, Iteration == (itr - 1 )), Lower_bound) 
  print(paste0("Current Lower bound ", curr_lb, " at iteration: ",itr))
  try({
    if(is.finite((curr_lb - prev_lb)) & abs(curr_lb - prev_lb)/abs(prev_lb) < para_testing$tol ){
      print(paste0("Optimization converged at step ", itr))
      break
    }
  })
}

#########################################################################
# code below is used to compute prediction log likelihood for individual
# with different number of incidences
#########################################################################
predictions <- list()
testing_eid_included <- testing_data %>% # set incidence number threshold for inclusion in testing set
  group_by(eid) %>%
  summarise(n()) %>%
  filter(`n()` >= 2) %>%
  select(eid)
predict_subset_included <- testing_eid_included %>%
  left_join(predicting_test_set, by = "eid")

predictions[[1]] <- prediction_ageLDA(para_testing$w,predict_subset_included, para_testing$alpha_z, para_testing$eid, para$pi_beta_basis, para$list_above500occu$diag_icd10)[[1]]
# estimate the performance for individuals with different number of records
for(rec_num in 2:10){
  testing_eid_included <- testing_data %>% 
    group_by(eid) %>%
    summarise(n()) %>%
    filter(`n()` == rec_num) %>%
    select(eid)
  predict_subset_included <- testing_eid_included %>%
    left_join(predicting_test_set, by = "eid")
  
  predictions[[rec_num]] <- prediction_ageLDA(para_testing$w, predict_subset_included, para_testing$alpha_z, para_testing$eid, para$pi_beta_basis, para$list_above500occu$diag_icd10)[[1]]
}

# save the parameter
save(predictions, file = paste0("Results/","prediction_age_K",para$K,"_P",para$P,"_rep",para$rep_ID, ".RData"))

