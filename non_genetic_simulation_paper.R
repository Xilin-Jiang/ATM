setwd("/users/mcvean/xilin/xilin/Multimorbidity-biobank/")
source("simulation_functions.R")
library(reshape2)
library(topicmodels)
args <- commandArgs(trailingOnly = TRUE)
v_id <- as.integer(args[1])
s_id <- as.integer(args[2])
rep_id <- as.integer(args[3])

set.seed(19940110 + rep_id)
if(v_id == 1){
  sample_sz <- c(500, 2000, 5000,  10000, 20000)[s_id]
}else{
  sample_sz <- 20000 # 20000
}
topic_number <- 5 # 2
disease_number <- 50 # 30
degree_freedom <- 3

if(v_id == 2){
  inf_topic <- (2:15)[s_id] # number of topic number assumed in inference
}else{
  inf_topic <- 5 # 20000
}

if(v_id == 3){
  ds_per_idv <- c(2,4,6,8,10)[s_id] # number of diseaes per-individual
}else{
  ds_per_idv <- 6.1 # number of cases per individual
}


para_sim <- simulate_age_topic_data(sample_sz, topic_number, 
                                            disease_number, degree_freedom, 
                                            ds_per_idv = ds_per_idv)
# compute the disease distribution
true_topic <- sapply(1:topic_number, function(i) 
  colMeans(para_sim$true_beta[30:81,,i]))
# revove noise disease
true_topic[which(true_topic < 0.01)] <- 0

#############
# this is the truth for computing later parameters
true_topic_per_ds <- apply(true_topic, 1, function(x) ifelse(sum(x) == 0, 0, which.max(x)) )


rec_data <- para_sim$unlist_Ds_id %>% 
  rename(diag_icd10 = Ds_id) %>%
  select(eid, diag_icd10, age_diag) %>%
  group_by(eid, diag_icd10) %>% 
  arrange(age_diag, .by_group = T) %>% # keep only the first records of repeated diagnosis
  slice(1) %>%
  ungroup()

ds_list <- rec_data %>%
  group_by(diag_icd10) %>%
  summarise(occ = n())

##################################
# perform ageLDA inference
##################################
para_inf_age <- topic_init_age(rec_data, ds_list, inf_topic, degree_freedom)

para_inf_age$max_itr <- 1000
para_inf_age$alpha <- rep(1, para_inf_age$K)
para_inf_age$lb <- data.frame("Iteration" = 0,"Lower_bound" = comp_lda_lb(para_inf_age))
para_inf_age$itr_beta <- 1
para_inf_age$itr_check_lb <- 1
para_inf_age$tol <- 10^(-6)
cvb_tag <- T
cvb0_tag <- T
for(itr in 1:para_inf_age$max_itr){
  print(paste0("Interation: ",itr))
  for(itr_inside in 1:para_inf_age$itr_beta){ # in practice we perform quick steps a few times before move on.
    if(cvb_tag){
      if(cvb0_tag){
        para_inf_age <- CVB0_E_zn(para_inf_age)
      }else{
        para_inf_age <- CVB_E_zn(para_inf_age)
      }
    }else{
      para_inf_age <- comp_E_zn(para_inf_age)
      para_inf_age <- comp_E_lntheta(para_inf_age)
    }
  }
  para_inf_age <- fast_update_age_depend_lda(para_inf_age)
  # para_inf_age <- update_alpha(para_inf_age) # we use non-informative alpha
  if(cvb_tag){
    para_inf_age$lb[nrow(para_inf_age$lb) + 1,] <- c(itr, CVB_lb(para_inf_age))
  }else{
    para_inf_age$lb[nrow(para_inf_age$lb) + 1,] <- c(itr, comp_lda_lb(para_inf_age))
  }
  if(itr %% para_inf_age$itr_check_lb ==0){
    curr_lb <- pull(filter(para_inf_age$lb, Iteration == itr), Lower_bound) 
    prev_lb <- pull(filter(para_inf_age$lb, Iteration == (itr -para_inf_age$itr_check_lb )), Lower_bound) 
    print(paste0("Current Lower bound ", curr_lb, " at iteration: ",itr))
    try({
      if(is.finite((curr_lb - prev_lb)) & (curr_lb - prev_lb)/abs(prev_lb) < para_inf_age$tol ){
        print(paste0("Optimization converged at step ", itr))
        break
      }
    })
  }
}
loadings_per_ds <- sapply(1:para_inf_age$D, function(j) colMeans(para_inf_age$unlist_zn[para_inf_age$ds_list[[j]]$id,])) %>% t
ageLDA_topic_per_ds <- apply(loadings_per_ds, 1, which.max)

##################################
# perform basic LDA inference
##################################
para_inf_lda <- topic_init_baseline(rec_data, ds_list, inf_topic)
para_inf_lda$max_itr <- 500
para_inf_lda$alpha <- rep(1, para_inf_lda$K) 
para_inf_lda$lb <- data.frame("Iteration" = 0,"Lower_bound" = CVB_lb(para_inf_lda))
para_inf_lda$itr_check_lb <- 5
para_inf_lda$itr_save <- para_inf_lda$max_itr # don't need to save intermediate data
para_inf_lda$tol <- 10^(-6)
for(itr in 1:para_inf_lda$max_itr){
  print(paste0("Interation: ",itr))
  para_inf_lda <- CVB0_E_zn(para_inf_lda) # we choose CVB0 as papers shown it could converge quicker
  para_inf_lda <- update_beta_basic_lda(para_inf_lda)
  # para_inf_lda <- update_alpha(para_inf_lda) # we use non-informative alpha
  para_inf_lda$lb[nrow(para_inf_lda$lb) + 1,] <- c(itr, CVB_lb(para_inf_lda))
  if(itr %% para_inf_lda$itr_check_lb  ==0){
    curr_lb <- pull(filter(para_inf_lda$lb, Iteration == itr), Lower_bound) 
    prev_lb <- pull(filter(para_inf_lda$lb, Iteration == (itr -para_inf_lda$itr_check_lb )), Lower_bound) 
    print(paste0("Current Lower bound ", curr_lb, " at iteration: ",itr))
    try({
      if(is.finite((curr_lb - prev_lb)) & (curr_lb - prev_lb)/abs(prev_lb) < para_inf_lda$tol ){
        print(paste0("Optimization converged at step ", itr))
        break
      }
    })
  }
}
loadings_per_ds <- sapply(1:para_inf_lda$D, function(j) colMeans(para_inf_lda$unlist_zn[para_inf_lda$ds_list[[j]]$id,])) %>% t
LDA_topic_per_ds <- apply(loadings_per_ds, 1, which.max)


# k <- inf_topic
# control_LDA_VEM <-
#   list(estimate.alpha = F, alpha = 1, estimate.beta = TRUE,
#        verbose = 0, prefix = tempfile(), save = 0, keep = 0,
#        seed = as.integer(Sys.time()), nstart = 1, best = TRUE,
#        var = list(iter.max = 500, tol = 10^-6),
#        em = list(iter.max = 1000, tol = 10^-4),
#        initialize = "random")
# model_lda <- LDA(lda_x, k, method = "VEM", control = control_LDA_VEM)
# lda_Result <- posterior(model_lda)
# attributes(lda_Result)
# LDA_topic_per_ds <- apply(lda_Result$terms, 2, which.max)

##################################
# perform PCA inference
##################################
lda_x <- sapply(1:para_inf_lda$M, function(s) 1*(1:para_inf_lda$D %in% para_inf_lda$w[[s]]$Ds_id)) %>% t
pca_rslt <- prcomp(lda_x, scale = TRUE)
pca_per_ds <- pca_rslt$rotation[,1:inf_topic]
# get the substracted value
flip_pca <- apply(pca_per_ds, 2, function(x) ifelse(-min(x) > max(x), -1, 1))
pca_per_ds <- sapply(1:inf_topic, function(x) pca_per_ds[,x] * flip_pca[x])
PCA_topic_per_ds <- apply(pca_per_ds, 1, which.max)

####################################
# perform testing and training prediction
####################################
if(v_id == 2){
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
  
  # perform the training
  para_train_age <- topic_init_age(training_data, ds_list, inf_topic, degree_freedom)
  
  para_train_age$max_itr <- 1000
  para_train_age$alpha <- rep(1, para_train_age$K)
  para_train_age$lb <- data.frame("Iteration" = 0,"Lower_bound" = comp_lda_lb(para_train_age))
  para_train_age$itr_beta <- 1
  para_train_age$itr_check_lb <- 1
  para_train_age$tol <- 10^(-6)
  cvb_tag <- T
  cvb0_tag <- T
  for(itr in 1:para_train_age$max_itr){
    print(paste0("Interation: ",itr))
    for(itr_inside in 1:para_train_age$itr_beta){ # in practice we perform quick steps a few times before move on.
      if(cvb_tag){
        if(cvb0_tag){
          para_train_age <- CVB0_E_zn(para_train_age)
        }else{
          para_train_age <- CVB_E_zn(para_train_age)
        }
      }else{
        para_train_age <- comp_E_zn(para_train_age)
        para_train_age <- comp_E_lntheta(para_train_age)
      }
    }
    para_train_age <- fast_update_age_depend_lda(para_train_age)
    # para_train_age <- update_alpha(para_train_age) # we use non-informative alpha
    if(cvb_tag){
      para_train_age$lb[nrow(para_train_age$lb) + 1,] <- c(itr, CVB_lb(para_train_age))
    }else{
      para_train_age$lb[nrow(para_train_age$lb) + 1,] <- c(itr, comp_lda_lb(para_train_age))
    }
    if(itr %% para_train_age$itr_check_lb ==0){
      curr_lb <- pull(filter(para_train_age$lb, Iteration == itr), Lower_bound) 
      prev_lb <- pull(filter(para_train_age$lb, Iteration == (itr -para_train_age$itr_check_lb )), Lower_bound) 
      print(paste0("Current Lower bound ", curr_lb, " at iteration: ",itr))
      try({
        if(is.finite((curr_lb - prev_lb)) & (curr_lb - prev_lb)/abs(prev_lb) < para_train_age$tol ){
          print(paste0("Optimization converged at step ", itr))
          break
        }
      })
    }
  }
  
  
  max_predict <- 11 # predict the records until which step?
  prediction_onebyone_rslt <- prediction_onebyone(testing_data, ds_list, para_train_age, max_predict) 
  
}else{
  prediction_onebyone_rslt <- NA
}

nongenetic_tests <- list(true_topic_per_ds, ageLDA_topic_per_ds, LDA_topic_per_ds, PCA_topic_per_ds, para_inf_age, prediction_onebyone_rslt, para_sim, para_inf_lda)
save(nongenetic_tests, file = paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/simulations/non_genetic_simulation/","nongen_simulation_v",v_id,"s",s_id, "rep", rep_id, ".RData"))

