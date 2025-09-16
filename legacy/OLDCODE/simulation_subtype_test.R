setwd("/users/mcvean/xilin/xilin/Multimorbidity-biobank/")
source("simulation_functions.R")
args <- commandArgs(trailingOnly = TRUE)
rep_ID <- as.integer(args[1])
# add another layer to perform simulation over sample size and disease number
type_of_simulation <- as.integer(args[2])

set.seed(19940110+rep_ID)

degree_freedom <- 3

if(type_of_simulation == 1 | type_of_simulation == 5 | type_of_simulation == 6){ # change the mix ratio
  variable_list <- 0.1 * (0:9)

  sample_sz <- 10000 # 20000
  topic_number <- 2 # 2
  disease_number <- 20 # 30
  ds_per_idv <- 6.1 # number of disease per individual 
}else if(type_of_simulation == 2){ # change population size
  mix_ratio <- 0.4
  variable_list <- c(200, 500, 1000, 2000, 5000, 10000) # 20000
  topic_number <- 2 # 2
  disease_number <- 20 # 30
  ds_per_idv <- 6.1 # number of disease per individual 
}else if(type_of_simulation == 3){ # change topic number, each topic still have 10 disease
  mix_ratio <- 0.4
  sample_sz <- 10000 # 20000
  variable_list <- c(2,3,4,5,6) 
  ds_per_idv <- 6.1 # number of disease per individual 
}else if(type_of_simulation == 4){ # change the number of disease per individual
  mix_ratio <- 0.4
  sample_sz <- 10000 # 20000
  topic_number <- 2 # 2
  disease_number <- 20 # 30
  variable_list <- c(2,4,6,8,10) 
}

# infer topics on the simulated data set: varying the proportion of mixing between disease 1 v.s. disease 10 (same age distribution)
# & disease 1 & 11 (different age distribution)
accuracy_age <- matrix(NA, nrow = 3, ncol = length(variable_list))
recall_age <- matrix(NA, nrow = 3, ncol = length(variable_list))
precision_age <- matrix(NA, nrow = 3, ncol = length(variable_list))
auc_age <- matrix(NA, nrow = 3, ncol = length(variable_list))

accuracy_lda <- matrix(NA, nrow = 3, ncol = length(variable_list))
recall_lda <- matrix(NA, nrow = 3, ncol = length(variable_list))
precision_lda <- matrix(NA, nrow = 3, ncol = length(variable_list))
auc_lda <- matrix(NA, nrow = 3, ncol = length(variable_list))

for(rt in 1:length(variable_list)){
  if(type_of_simulation == 1 | type_of_simulation == 5 | type_of_simulation == 6){
    mix_ratio <- variable_list[rt]
  }else if(type_of_simulation == 2){
    sample_sz <- variable_list[rt] # 20000
  }else if(type_of_simulation == 3){
    topic_number <- variable_list[rt] # 2
    disease_number <- variable_list[rt] * 10 #10 diseaes for each topic
  }else if(type_of_simulation == 4){
    ds_per_idv <- variable_list[rt] 
  }
  
  
  # simulate a dataset for relabelling
  para_orginal_sim <- simulate_age_topic_data(sample_sz, topic_number, disease_number, degree_freedom, ds_per_idv = ds_per_idv)
  
  # need to make sure each diseaes only occur for once
  rec_data <- para_orginal_sim$unlist_Ds_id %>% 
    rename(diag_icd10 = Ds_id) %>%
    select(eid, diag_icd10, age_diag) %>%
    group_by(eid, diag_icd10) %>% 
    arrange(age_diag, .by_group = T) %>% # keep only the first records of repeated diagnosis
    slice(1) %>%
    ungroup() %>%
    mutate(rowid = 1:n())
  
  # add a 5th case to move the age difference to 10 years
  if(type_of_simulation == 5){
    rec_data <- rec_data %>%
      mutate(age_diag = ifelse(diag_icd10 == 11, age_diag - 10, age_diag))
  }else if(type_of_simulation == 6){
    rec_data <- rec_data %>%
      mutate(age_diag = ifelse(diag_icd10 == 11, age_diag - 15, age_diag))
  }
  
  new_rec_data <- list()
  new_ds_list <- list()
  # create first type of mixing 1&10 with similar age distribution
  ids_sample <- rec_data %>% 
    filter(diag_icd10 == 10) %>%
    pull(rowid )
  id_change <- sample(ids_sample, size = floor(length(ids_sample) *mix_ratio),replace = F)
  new_rec_data[[1]] <- rec_data %>%
    mutate(diag_icd10 = ifelse(rowid %in% id_change, 1, diag_icd10) )  %>%
    select(eid, diag_icd10, age_diag) # %>% # we need to keep the data order the same to know what are the true label for each disease
    # group_by(eid, diag_icd10) %>% 
    # arrange(age_diag, .by_group = T) %>% # keep only the first records of repeated diagnosis
    # slice(1) %>%
    # ungroup()
  new_ds_list[[1]] <- new_rec_data[[1]] %>%
    group_by(diag_icd10) %>%
    summarise(occ = n())
  
  # create second type of mixing 1&11 with different age distribution
  ids_sample <- rec_data %>% 
    filter(diag_icd10 == 11) %>%
    pull(rowid )
  id_change <- sample(ids_sample, size = floor(length(ids_sample) *mix_ratio),replace = F)
  new_rec_data[[2]] <- rec_data %>%
    mutate(diag_icd10 = ifelse(rowid %in% id_change, 1, diag_icd10) )  %>%
    select(eid, diag_icd10, age_diag) 
  new_ds_list[[2]] <- new_rec_data[[2]] %>%
    group_by(diag_icd10) %>%
    summarise(occ = n())
  
  # create third type of mixing 1&11; change label 1 -> 11 
  ids_sample <- rec_data %>% 
    filter(diag_icd10 == 1) %>%
    pull(rowid )
  id_change <- sample(ids_sample, size = floor(length(ids_sample) *mix_ratio),replace = F)
  new_rec_data[[3]] <- rec_data %>%
    mutate(diag_icd10 = ifelse(rowid %in% id_change, 11, diag_icd10) )  %>%
    select(eid, diag_icd10, age_diag) 
  new_ds_list[[3]] <- new_rec_data[[3]] %>%
    group_by(diag_icd10) %>%
    summarise(occ = n())
  
  for(type_mix in 1:3){
    para_mix <- topic_init_age(new_rec_data[[type_mix]], new_ds_list[[type_mix]], topic_number, degree_freedom)
    
    para_mix$max_itr <- 1000
    para_mix$alpha <- rep(1, para_mix$K)
    para_mix$lb <- data.frame("Iteration" = 0,"Lower_bound" = comp_lda_lb(para_mix))
    para_mix$itr_beta <- 1
    para_mix$itr_check_lb <- 1
    para_mix$tol <- 10^(-6)
    cvb_tag <- T
    cvb0_tag <- T
    for(itr in 1:para_mix$max_itr){
      print(paste0("Interation: ",itr))
      for(itr_inside in 1:para_mix$itr_beta){ # in practice we perform quick steps a few times before move on.
        if(cvb_tag){
          if(cvb0_tag){
            para_mix <- CVB0_E_zn(para_mix)
          }else{
            para_mix <- CVB_E_zn(para_mix)
          }
        }else{
          para_mix <- comp_E_zn(para_mix)
          para_mix <- comp_E_lntheta(para_mix)
        }
      }
      para_mix <- fast_update_age_depend_lda(para_mix)
      # para_mix <- update_alpha(para_mix) # we use non-informative alpha
      if(cvb_tag){
        para_mix$lb[nrow(para_mix$lb) + 1,] <- c(itr, CVB_lb(para_mix))
      }else{
        para_mix$lb[nrow(para_mix$lb) + 1,] <- c(itr, comp_lda_lb(para_mix))
      }
      if(itr %% para_mix$itr_check_lb ==0){
        curr_lb <- pull(filter(para_mix$lb, Iteration == itr), Lower_bound) 
        prev_lb <- pull(filter(para_mix$lb, Iteration == (itr -para_mix$itr_check_lb )), Lower_bound) 
        print(paste0("Current Lower bound ", curr_lb, " at iteration: ",itr))
        try({
          if(is.finite((curr_lb - prev_lb)) & (curr_lb - prev_lb)/abs(prev_lb) < para_mix$tol ){
            print(paste0("Optimization converged at step ", itr))
            break
          }
        })
      }
    }
    
    # inference using basic_LDA model
    para_mix_non_age <- topic_init_baseline(new_rec_data[[type_mix]], new_ds_list[[type_mix]], topic_number)
    para_mix_non_age$max_itr <- 1000
    para_mix_non_age$alpha <- rep(1, para_mix_non_age$K) 
    para_mix_non_age$lb <- data.frame("Iteration" = 0,"Lower_bound" = CVB_lb(para_mix_non_age))
    para_mix_non_age$itr_check_lb <- 5
    para_mix_non_age$itr_save <- para_mix_non_age$max_itr # don't need to save intermediate data
    para_mix_non_age$tol <- 10^(-6)
    for(itr in 1:para_mix_non_age$max_itr){
      print(paste0("Interation: ",itr))
      para_mix_non_age <- CVB0_E_zn(para_mix_non_age) # we choose CVB0 as papers shown it could converge quicker
      para_mix_non_age <- update_beta_basic_lda(para_mix_non_age)
      # para_mix_non_age <- update_alpha(para_mix_non_age) # we use non-informative alpha
      para_mix_non_age$lb[nrow(para_mix_non_age$lb) + 1,] <- c(itr, CVB_lb(para_mix_non_age))
      if(itr %% para_mix_non_age$itr_check_lb  ==0){
        curr_lb <- pull(filter(para_mix_non_age$lb, Iteration == itr), Lower_bound) 
        prev_lb <- pull(filter(para_mix_non_age$lb, Iteration == (itr -para_mix_non_age$itr_check_lb )), Lower_bound) 
        print(paste0("Current Lower bound ", curr_lb, " at iteration: ",itr))
        try({
          if(is.finite((curr_lb - prev_lb)) & (curr_lb - prev_lb)/abs(prev_lb) < para_mix_non_age$tol ){
            print(paste0("Optimization converged at step ", itr))
            break
          }
        })
      }
    }
    
    # compute the probability of assigning each incidence to the correct topic
    # three cases: 10 -> 1 or 11 -> 1 or 1 -> 11
    if(type_mix == 1){
      ds_late <- 10 # disease true label for the minor group
      ds_early <- 1 # disease true label for the major group
    }else if(type_mix == 2){
      ds_late <- 11
      ds_early <- 1
    }else if(type_mix == 3){
      ds_late <- 11
      ds_early <- 1
    }
    ##############
    # ageLDA model
    ##############
    ###############################################################
    # need to figure out how to compute precision in this case
    
    tp_ealry <- colMeans(para_mix$unlist_zn[which(para_mix$unlist_Ds_id$Ds_id == (ds_early + 1)),])  %>%
      which.max() # the topic label for the true disease
    tp_late <- colMeans(para_mix$unlist_zn[which(para_mix$unlist_Ds_id$Ds_id == (ds_late+ 1)),]) %>%
      which.max() 
    # the early and late group is identified in each condition
    if(sum(rec_data$diag_icd10 ==  ds_early) > sum(para_mix$unlist_Ds_id$Ds_id ==  ds_early)){
      early_subtype_list <- intersect(which(rec_data$diag_icd10 == ds_early),
                                      which(para_mix$unlist_Ds_id$Ds_id == ds_late))
      late_subtype_list <- intersect(which(rec_data$diag_icd10 == ds_late),
                                     which(para_mix$unlist_Ds_id$Ds_id == ds_late))
    }else{
      early_subtype_list <- intersect(which(rec_data$diag_icd10 == ds_early),
                                      which(para_mix$unlist_Ds_id$Ds_id == ds_early))
      late_subtype_list <- intersect(which(rec_data$diag_icd10 == ds_late),
                                     which(para_mix$unlist_Ds_id$Ds_id == ds_early))
    }
    prediction_early <- apply(para_mix$unlist_zn[early_subtype_list,], 1, which.max)
    prediction_late <- apply(para_mix$unlist_zn[late_subtype_list,], 1, which.max)

    accuracy_age[type_mix, rt] <- mean(c(prediction_early == tp_ealry, prediction_late == tp_late))
    if(type_mix == 3){
      recall_age[type_mix, rt]  <- mean(prediction_early == tp_ealry)
      precision_age[type_mix, rt]  <- sum(prediction_early == tp_ealry)/
        sum(c(prediction_early, prediction_late) == tp_ealry)
      # compute AUC for classification a minor subtypes
      try({
        pr_curve <- pr.curve(scores.class0 = para_mix$unlist_zn[early_subtype_list,tp_ealry] , scores.class1 = para_mix$unlist_zn[late_subtype_list,tp_ealry], curve = T)
        auc_age[type_mix, rt]  <- pr_curve$auc.integral
        # auc_age[type_mix, rt]  <- auc(c(rep(1,length(early_subtype_list)), rep(0,length(late_subtype_list))),  
        #                               c(para_mix$unlist_zn[early_subtype_list,tp_ealry], para_mix$unlist_zn[late_subtype_list,tp_ealry]) )
      })
    }else{
      recall_age[type_mix, rt]  <- mean(prediction_late == tp_late)
      precision_age[type_mix, rt]  <- sum(prediction_late == tp_late)/
        sum(c(prediction_early, prediction_late) == tp_late)
      # compute AUC for classification a minor subtypes
      try({
        pr_curve <- pr.curve(scores.class0 =para_mix$unlist_zn[late_subtype_list,tp_late],  scores.class1 = para_mix$unlist_zn[early_subtype_list,tp_late] , curve = T)
        auc_age[type_mix, rt]  <- pr_curve$auc.integral
        # auc_age[type_mix, rt]  <- auc(c(rep(1,length(early_subtype_list)), rep(0,length(late_subtype_list))),  
        #                               c(para_mix$unlist_zn[early_subtype_list,tp_ealry], para_mix$unlist_zn[late_subtype_list,tp_ealry]) )
      })
    }


    #################
    # basic LDA model
    #################
    tp_ealry <- colMeans(para_mix_non_age$unlist_zn[which(para_mix_non_age$unlist_Ds_id$Ds_id == (ds_early + 1)),])  %>%
      which.max() # the topic label for the true disease
    tp_late <- colMeans(para_mix_non_age$unlist_zn[which(para_mix_non_age$unlist_Ds_id$Ds_id == (ds_late+ 1)),]) %>%
      which.max() 
    # the early and late group is identified in each condition
    if(sum(rec_data$diag_icd10 ==  ds_early) > sum(para_mix_non_age$unlist_Ds_id$Ds_id ==  ds_early)){
      early_subtype_list <- intersect(which(rec_data$diag_icd10 == ds_early),
                                      which(para_mix_non_age$unlist_Ds_id$Ds_id == ds_late))
      late_subtype_list <- intersect(which(rec_data$diag_icd10 == ds_late),
                                     which(para_mix_non_age$unlist_Ds_id$Ds_id == ds_late))
    }else{
      early_subtype_list <- intersect(which(rec_data$diag_icd10 == ds_early),
                                      which(para_mix_non_age$unlist_Ds_id$Ds_id == ds_early))
      late_subtype_list <- intersect(which(rec_data$diag_icd10 == ds_late),
                                     which(para_mix_non_age$unlist_Ds_id$Ds_id == ds_early))
    }
    prediction_early <- apply(para_mix_non_age$unlist_zn[early_subtype_list,], 1, which.max)
    prediction_late <- apply(para_mix_non_age$unlist_zn[late_subtype_list,], 1, which.max)
    
    accuracy_lda[type_mix, rt] <- mean(c(prediction_early == tp_ealry, prediction_late == tp_late))
    if(type_mix == 3){
      recall_lda[type_mix, rt]  <- mean(prediction_early == tp_ealry)
      precision_lda[type_mix, rt]  <- sum(prediction_early == tp_ealry)/
        sum(c(prediction_early, prediction_late) == tp_ealry)
      # compute AUC for classification a minor subtypes
      try({
        pr_curve <- pr.curve(scores.class0 = para_mix_non_age$unlist_zn[early_subtype_list,tp_ealry] , scores.class1 = para_mix_non_age$unlist_zn[late_subtype_list,tp_ealry], curve = T)
        auc_lda[type_mix, rt]  <- pr_curve$auc.integral
        # auc_age[type_mix, rt]  <- auc(c(rep(1,length(early_subtype_list)), rep(0,length(late_subtype_list))),  
        #                               c(para_mix_non_age$unlist_zn[early_subtype_list,tp_ealry], para_mix_non_age$unlist_zn[late_subtype_list,tp_ealry]) )
      })
    }else{
      recall_lda[type_mix, rt]  <- mean(prediction_late == tp_late)
      precision_lda[type_mix, rt]  <- sum(prediction_late == tp_late)/
        sum(c(prediction_early, prediction_late) == tp_late)
      # compute AUC for classification a minor subtypes
      try({
        pr_curve <- pr.curve(scores.class0 =para_mix_non_age$unlist_zn[late_subtype_list,tp_late],  scores.class1 = para_mix_non_age$unlist_zn[early_subtype_list,tp_late] , curve = T)
        auc_lda[type_mix, rt]  <- pr_curve$auc.integral
        # auc_age[type_mix, rt]  <- auc(c(rep(1,length(early_subtype_list)), rep(0,length(late_subtype_list))),  
        #                               c(para_mix_non_age$unlist_zn[early_subtype_list,tp_ealry], para_mix_non_age$unlist_zn[late_subtype_list,tp_ealry]) )
      })
    }
    # this is the end for each loop over a specific mixing type (whether there are age difference or not)
  }
  # this is the end for each loop over a specific mixing ratio
}

subtype_tests <- list(accuracy_age, recall_age, precision_age, auc_age,
                      accuracy_lda, recall_lda, precision_lda, auc_lda)
save(subtype_tests, file = paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/simulations/non_genetic_simulation/subtype_tests_simulation_type", type_of_simulation, "_rep",rep_ID, ".RData"))

