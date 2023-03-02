##########################################
# Prediction functions: should be an abstraction for all three models
# for tree and baseline, we basically assume topics are constant over age
##########################################
# use the function below to estimate the probability that a disease appear in the top predicted disease list (1%, 2%, 5% 10%)
prediction_ageLDA <- function(estimate_set, predict_set, alpha_z, estimate_eid, beta, diag_icd10){
  predictions <- list()
  code2id <- function(x){
    return( match(x, diag_icd10))
  }
  # normalise alpha_z by row to get estimation of theta
  theta_mean <- sweep(alpha_z, 1, rowSums(alpha_z), FUN="/")

  # here I am rounding the disease time to year for computation efficiency
  predict_list_numeric <- predict_set %>%
    mutate(Ds_id = code2id(diag_icd10)) %>%
    select(-diag_icd10) %>%
    mutate(age_diag = round(age_diag)) %>%
    mutate(eid = match(eid, estimate_eid))

  theta_n <- theta_mean[predict_list_numeric$eid,]
  beta_n <- apply(beta, 3, function(x)
    x[as.matrix(select(predict_list_numeric, age_diag, Ds_id))])
  predictions$test_log_lik <- sum(log(rowSums(theta_n * beta_n)))

  # compute the average rank: for both the prediction and those used in estimation
  beta_rank_all <- sapply(1:dim(theta_n)[1], function(n) # compute the marginal probability for each disease
    rank(- tcrossprod(beta[predict_list_numeric$age_diag[n],,],theta_n[n,,drop = F]) )[c(predict_list_numeric$Ds_id[n],
                                                                                         estimate_set[[predict_list_numeric$eid[n]]]$Ds_id)] )

  # this is a complicated step: we want to remove the ICD-10 codes that have been used in estimation
  beta_rank <- beta_rank_all[1, ] - sapply(1:dim(beta_rank_all)[2], function(n)
    sum(beta_rank_all[1,n] > beta_rank_all[2:dim(beta_rank_all)[1],n]))

  # compute the specificity: first step remove the number of diseases used in estimation
  total_number <- length(diag_icd10) - mean(sapply(1:length(estimate_set), function(x) dim(estimate_set[[x]])[1]))
  predictions$mean_rank <- mean(beta_rank)/total_number
  predictions$mean_rank_top1 <- mean( beta_rank <= total_number/100)
  predictions$mean_rank_top2 <- mean( beta_rank <= total_number/50)
  predictions$mean_rank_top5 <- mean( beta_rank <= total_number/20)
  predictions$mean_rank_top10 <- mean( beta_rank <= total_number/10)

  return(list(predictions, beta_rank))
}
# using function below to perform onebyone prediction
prediction_onebyone <- function(testing_data, ds_list, para_training, max_predict){
  # first order the incidences by age
  testing_data <- testing_data %>%
    group_by(eid, diag_icd10) %>%
    arrange(age_diag, .by_group = T) %>% # keep only the first records of repeated diagnosis
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    group_by(eid) %>%
    arrange(age_diag, .by_group = TRUE)

  prediction_onebyone <- list()
  collect_prediction_ranks <- c()

  for(predict_num in 2:max_predict){
    testing_eid_included <- testing_data %>% # set incidence number threshold for inclusion in testing set
      group_by(eid) %>%
      summarise(n()) %>%
      filter(`n()` >= predict_num) %>%
      select(eid)

    testing_data_included <- testing_eid_included %>%
      left_join(testing_data, by = "eid")

    # for the last iteration, just predict all left
    if(predict_num == max_predict){
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        dplyr::slice(predict_num:n()) %>%
        dplyr::ungroup()
    }else{
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        dplyr::slice(predict_num) %>%
        dplyr::ungroup()
    }

    estimating_test_set <- testing_data_included %>%
      group_by(eid) %>%
      dplyr::slice(1:(predict_num-1)) %>%
      dplyr::ungroup()
    if(is.null(para_training$P)){ # using P to determine if age is included
      para_training$pi_beta_basis <- array( rep(para_training$beta, each = 81),dim =  c(81, para$D, para$K) )
      # for both treeLDA and baseline lda, we only need to initialise baseline case here
      para_testing <- topic_init_baseline(estimating_test_set, ds_list, para_training$K)
    }else{
      para_testing <- topic_init_age(estimating_test_set, ds_list, para_training$K, para_training$P)
    }
    # assigning beta to the testing set
    para_testing$beta_w_full <- apply(para_training$pi_beta_basis, 3, function(x)
      x[as.matrix(select(para_testing$unlist_Ds_id, age_diag, Ds_id))])
    para_testing$beta_w <- lapply(para_testing$patient_lst, function(x) para_testing$beta_w_full[x,,drop=F])

    # updating Z_n
    para_testing$max_itr <- 10
    para_testing$alpha <- para_training$alpha
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
    rslt_predict <- prediction_ageLDA(para_testing$w, predicting_test_set, para_testing$alpha_z, para_testing$eid, para_training$pi_beta_basis, para_training$list_above500occu$diag_icd10)
    prediction_onebyone[[predict_num]]  <- rslt_predict[[1]]
    collect_prediction_ranks <- c(collect_prediction_ranks, rslt_predict[[2]])
  }
  prediction_onebyone[[1]] <- list(mean(collect_prediction_ranks)/para_training$D,
                                   mean( collect_prediction_ranks <= para_training$D/100),
                                   mean( collect_prediction_ranks <= para_training$D/50),
                                   mean( collect_prediction_ranks <= para_training$D/20),
                                   mean( collect_prediction_ranks <= para_training$D/10) )
  return(prediction_onebyone)
}
# use the function to compute the OR for each disease
risk_each_disease <- function(estimate_set, predict_set, alpha_z, estimate_eid, beta, diag_icd10){
  code2id <- function(x){
    return( match(x, diag_icd10))
  }
  # normalise alpha_z by row to get estimation of theta
  theta_mean <- sweep(alpha_z, 1, rowSums(alpha_z), FUN="/")

  # here I am rounding the disease time to year for computation efficiency
  predict_list_numeric <- predict_set %>%
    mutate(Ds_id = code2id(diag_icd10)) %>%
    select(-diag_icd10) %>%
    mutate(age_diag = round(age_diag)) %>%
    mutate(eid = match(eid, estimate_eid))

  theta_n <- theta_mean[predict_list_numeric$eid,]
  beta_n <- apply(beta, 3, function(x)
    x[as.matrix(select(predict_list_numeric, age_diag, Ds_id))])

  # compute the marginal probability for each disease
  risk_all <- sapply(1:dim(theta_n)[1], function(n)
    tcrossprod(beta[predict_list_numeric$age_diag[n],,],theta_n[n,,drop = F]) )
  return(list(risk_all, predict_list_numeric$Ds_id, theta_n))
}
# using function below to compute per sd risk associated with odds-ratio
prediction_OR_onebyone <- function(testing_data, ds_list, para_training, max_predict){
  # first order the incidences by age
  testing_data <- testing_data %>%
    group_by(eid) %>%
    arrange(age_diag, .by_group = TRUE)

  collect_risk_set <- list() # save all risk profiles
  collect_ds_set <- c() # save all incidence index
  collect_loadings <- list()

  # use this for saving the final results
  OR_each_disease <- list()
  loadings_each_disease <- list()

  for(predict_num in 2:max_predict){
    testing_eid_included <- testing_data %>% # set incidence number threshold for inclusion in testing set
      group_by(eid) %>%
      summarise(n()) %>%
      filter(`n()` >= predict_num) %>%
      select(eid)

    testing_data_included <- testing_eid_included %>%
      left_join(testing_data, by = "eid")

    # for the last iteration, just predict all left
    if(predict_num == max_predict){
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        dplyr::slice(predict_num:n()) %>%
        dplyr::ungroup()
    }else{
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        dplyr::slice(predict_num) %>%
        dplyr::ungroup()
    }

    estimating_test_set <- testing_data_included %>%
      group_by(eid) %>%
      dplyr::slice(1:(predict_num-1)) %>%
      dplyr::ungroup()
    if(is.null(para_training$P)){ # using P to determine if age is included
      para_training$pi_beta_basis <- array( rep(para_training$beta, each = 81),dim =  c(81, para_training$D, para_training$K) )
      # for both treeLDA and baseline lda, we only need to initialise baseline case here
      para_testing <- topic_init_baseline(estimating_test_set, ds_list, para_training$K)
    }else{
      para_testing <- topic_init_age(estimating_test_set, ds_list, para_training$K, para_training$P)
    }
    # assigning beta to the testing set
    para_testing$beta_w_full <- apply(para_training$pi_beta_basis, 3, function(x)
      x[as.matrix(select(para_testing$unlist_Ds_id, age_diag, Ds_id))])
    para_testing$beta_w <- lapply(para_testing$patient_lst, function(x) para_testing$beta_w_full[x,,drop=F])

    # updating Z_n
    para_testing$max_itr <- 10
    para_testing$alpha <- para_training$alpha
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
    rslt_OR <- risk_each_disease(para_testing$w, predicting_test_set, para_testing$alpha_z, para_testing$eid, para_training$pi_beta_basis, para_training$list_above500occu$diag_icd10)
    collect_risk_set[[predict_num]]  <- rslt_OR[[1]]
    collect_ds_set[[predict_num]]  <- rslt_OR[[2]]
    collect_loadings[[predict_num]]  <- rslt_OR[[3]]

    # save results for each records
    risk <- t(rslt_OR[[1]])
    ds <- rslt_OR[[2]]
    tp_ld <- rslt_OR[[3]] # save loadings

    OR_rslt <- sapply(1:dim(risk)[2], function(j)
      mean(risk[which(ds == j),j])/mean(risk[-which(ds == j),j]))
    OR_rslt[is.na(OR_rslt)] <- 1
    OR_each_disease[[predict_num]] <- OR_rslt
    loadings_rslt <- sapply(1:dim(risk)[2], function(j)
      colMeans(tp_ld[which(ds == j), , drop = F]))
    loadings_rslt[is.na(loadings_rslt)] <- 1
    loadings_each_disease[[predict_num]] <- t(loadings_rslt)
  }
  risk_set <- do.call(cbind, collect_risk_set) %>% t
  ds_set <- unlist(collect_ds_set)
  loading_set  <- do.call(rbind, collect_loadings)

  OR_each_disease[[1]] <- sapply(1:dim(risk_set)[2], function(j)
    mean(risk_set[which(ds_set == j),j])/mean(risk_set[-which(ds_set == j),j]))

  loadings_each_disease[[1]] <- sapply(1:dim(risk_set)[2], function(j)
    colMeans(loading_set[which(ds_set == j), , drop = F])) %>% t

  # compute the logistic regression for each diseases: predictors are log risk ratio (sd = 1)
  logRR <- sapply(1:dim(risk_set)[2], function(x) log(risk_set[,x]/mean(risk_set[,x]))  )
  logRR <- sapply(1:dim(risk_set)[2], function(x) logRR[,x]/sd(logRR[,x]) )
  # use logistic regression to save the OR per sd logRR and p-value
  logstic_results <- sapply(1:dim(risk_set)[2], function(x)
    summary(glm((ds_set == x) ~ logRR[,x], family = binomial))$coefficients[2,] ) %>% t

  return(list(OR_each_disease, loadings_each_disease, para_training$pi_beta_basis, logstic_results))
}

# using function below save details on how the model makes decision
save_prediction_logics <- function(testing_data, ds_list, para_training, max_predict){
  # first order the incidences by age
  testing_data <- testing_data %>%
    group_by(eid) %>%
    arrange(age_diag, .by_group = TRUE)

  # collect zns
  collect_zn <- list()
  # collect history set
  collect_estimate_set <- list()
  collect_prediction_ranks <- c()

  for(predict_num in 2:max_predict){
    testing_eid_included <- testing_data %>% # set incidence number threshold for inclusion in testing set
      group_by(eid) %>%
      summarise(n()) %>%
      filter(`n()` >= predict_num) %>%
      select(eid)

    testing_data_included <- testing_eid_included %>%
      left_join(testing_data, by = "eid")

    # for the last iteration, just predict all left
    if(predict_num == max_predict){
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        dplyr::slice(predict_num:n()) %>%
        dplyr::ungroup()
    }else{
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        dplyr::slice(predict_num) %>%
        dplyr::ungroup()
    }

    estimating_test_set <- testing_data_included %>%
      group_by(eid) %>%
      dplyr::slice(1:(predict_num-1)) %>%
      dplyr::ungroup()
    if(is.null(para_training$P)){ # using P to determine if age is included
      para_training$pi_beta_basis <- array( rep(para$beta, each = 81),dim =  c(81, para$D, para$K) )
      # for both treeLDA and baseline lda, we only need to initialise baseline case here
      para_testing <- topic_init_baseline(estimating_test_set, ds_list, para_training$K)
    }else{
      para_testing <- topic_init_age(estimating_test_set, ds_list, para_training$K, para_training$P)
    }
    # assigning beta to the testing set
    para_testing$beta_w_full <- apply(para_training$pi_beta_basis, 3, function(x)
      x[as.matrix(select(para_testing$unlist_Ds_id, age_diag, Ds_id))])
    para_testing$beta_w <- lapply(para_testing$patient_lst, function(x) para_testing$beta_w_full[x,,drop=F])

    # updating Z_n
    para_testing$max_itr <- 10
    para_testing$alpha <- para_training$alpha
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

    collect_zn[[predict_num]]  <- para_testing$unlist_zn
    collect_estimate_set[[predict_num]]  <- estimating_test_set
    rslt_predict <- prediction_ageLDA(para_testing$w, predicting_test_set, para_testing$alpha_z, para_testing$eid, para_training$pi_beta_basis, para_training$list_above500occu$diag_icd10)
    collect_prediction_ranks <- c(collect_prediction_ranks, rslt_predict[[2]])
  }
  return(list(collect_zn, collect_estimate_set, collect_prediction_ranks, testing_data))
}

prediction_PheRS_by_disease <- function(testing_data, ds_list, para_training){
  # first order the incidences by age
  testing_data <- testing_data %>%
    group_by(eid) %>%
    arrange(age_diag, .by_group = TRUE)

  LogORestimates <- matrix(nrow = para_training$D, ncol = 4)
  AUC_methods <- matrix(nrow = para_training$D, ncol = 4)
  for(j in 1:para_training$D){
    ds_id <- para_training$list_above500occu$diag_icd10[j]
    print(ds_id)
    cases_eid <- testing_data %>% # select all the cases
      filter(diag_icd10 == ds_id) %>%
      rename(target_age = age_diag) %>%
      select(-diag_icd10)

    cases_data <- cases_eid %>% # select all disease that has a diagnosis that are earlier than target!
      left_join(testing_data, by = "eid") %>%
      filter(age_diag < target_age)

    # exclude cases that have no prior diagnosis
    # (i.e. target disease happen to be the first one; no way to predict!)
    targe_age_distribution <- cases_data %>%
      group_by(eid) %>%
      dplyr::slice(1)
    control_eid <- testing_data %>%
      anti_join(cases_eid, by = "eid") %>%
      group_by(eid) %>%
      dplyr::slice(1) %>%
      select(eid)
    control_eid$target_age <- sample(targe_age_distribution$target_age, dim(control_eid)[1], replace = T)

    control_data <- control_eid %>% # select all disease that has a diagnosis that are earlier than target!
      left_join(testing_data, by = "eid") %>%
      filter(age_diag < target_age)

    # save the case and control outcomes; save the age of the last event before target for computing area
    cases_eid <- cases_data %>%
      group_by(eid) %>%
      dplyr::slice_tail(n = 1) %>%
      mutate(outcome = 1) %>%
      select(eid, outcome, age_diag)
    control_eid <- control_data %>%
      group_by(eid) %>%
      dplyr::slice_tail(n = 1) %>%
      mutate(outcome = 0) %>%
      select(eid, outcome, age_diag)
    outcomes <- bind_rows(cases_eid, control_eid) %>%
      left_join(survive_age, by = "eid") %>%
      mutate(age_diag = round(age_diag), survive_year = round(survive_year))
    # we are underestimate the score for the cases, as some of them died of the disease
    outcomes %>% mutate(age_gap = survive_year - age_diag) %>% filter(outcome == 1) %>% dplyr::ungroup() %>% summarise(mean(age_gap))
    outcomes %>% mutate(age_gap = survive_year - age_diag) %>% filter(outcome == 0) %>% dplyr::ungroup() %>% summarise(mean(age_gap))
    # below is how we estimate the area under curve as well as the
    estimating_test_set <- bind_rows(cases_data, control_data) %>%
      select(- target_age) %>%
      dplyr::ungroup()

    if(is.null(para_training$P)){ # using P to determine if age is included
      para_training$pi_beta_basis <- array( rep(para_training$beta, each = 81),dim =  c(81, para_training$D, para_training$K) )
      # for both treeLDA and baseline lda, we only need to initialise baseline case here
      para_testing <- topic_init_baseline(estimating_test_set, ds_list, para_training$K)
    }else{
      para_testing <- topic_init_age(estimating_test_set, ds_list, para_training$K, para_training$P)
    }
    # assigning beta to the testing set
    para_testing$beta_w_full <- apply(para_training$pi_beta_basis, 3, function(x)
      x[as.matrix(select(para_testing$unlist_Ds_id, age_diag, Ds_id))])
    para_testing$beta_w <- lapply(para_testing$patient_lst, function(x) para_testing$beta_w_full[x,,drop=F])

    # updating Z_n
    para_testing$max_itr <- 10
    para_testing$alpha <- para_training$alpha
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
    alpha_z <- para_testing$alpha_z
    beta <- para_training$pi_beta_basis
    # normalise alpha_z by row to get estimation of theta
    theta_n <- sweep(alpha_z, 1, rowSums(alpha_z), FUN="/")
    beta_j <- beta[,j,]
    # there are different ways to compute risk score: we considered four: (1) mean till end (2) sum of area (3) risk at the prior event (4) three year average from prior event
    beta_n_mean <- sapply(1:dim(theta_n)[1], function(n) colMeans(beta_j[ outcomes$age_diag[n]:outcomes$survive_year[n], ,drop = F]) ) %>% t
    beta_n_sum <- sapply(1:dim(theta_n)[1], function(n) colSums(beta_j[ outcomes$age_diag[n]:outcomes$survive_year[n], ,drop = F]) ) %>% t
    beta_3year_mean <- sapply(1:dim(theta_n)[1], function(n) colSums(beta_j[ outcomes$age_diag[n]:pmin(81,outcomes$age_diag[n] + 3), ,drop = F]) ) %>% t
    beta_point <- sapply(1:dim(theta_n)[1], function(n) beta_j[ outcomes$age_diag[n], ,drop = F])  %>% t

    # compute the marginal probability for each disease
    risk_n_mean <- sapply(1:dim(theta_n)[1], function(n) tcrossprod(beta_n_mean[n,,drop = F],theta_n[n,,drop = F]) )
    risk_n_sum <- sapply(1:dim(theta_n)[1], function(n) tcrossprod(beta_n_sum[n,,drop = F],theta_n[n,,drop = F]) )
    risk_3year_mean <- sapply(1:dim(theta_n)[1], function(n) tcrossprod(beta_3year_mean[n,,drop = F],theta_n[n,,drop = F]) )
    risk_point <- sapply(1:dim(theta_n)[1], function(n) tcrossprod(beta_point[n,,drop = F],theta_n[n,,drop = F]) )
    # RS <- risk_disease/mean(risk_disease)
    # RS <- RS/sd(RS)
    outcomes$risk_n_mean <- risk_n_mean/sd(risk_n_mean)
    outcomes$risk_n_sum <- risk_n_sum/sd(risk_n_sum) # this should be the accurate case
    outcomes$risk_3year_mean <- risk_3year_mean/sd(risk_3year_mean)
    outcomes$risk_point <- risk_point/sd(risk_3year_mean)
    # using pROC package
    AUC_methods[j,1] <- pROC::auc(outcomes$outcome,outcomes$risk_n_sum )
    AUC_methods[j,2] <- pROC::auc(outcomes$outcome,outcomes$risk_n_mean )
    AUC_methods[j,3] <- pROC::auc(outcomes$outcome,outcomes$risk_3year_mean )
    AUC_methods[j,4] <- pROC::auc(outcomes$outcome,outcomes$risk_point )

    LogORestimates[j,] <- summary(glm(outcomes$outcome~ outcomes$risk_n_sum, family = binomial))$coefficients[2,]
  }
  return(list(LogORestimates, AUC_methods, outcomes$outcome, outcomes$risk_n_sum))
}

# now work on using LASSO for prediction
LASSO_predict <- function(rec_data, para){
  all_eid <- rec_data %>%
    group_by(eid) %>%
    summarise()
  training_eid <- data.frame(eid = para$eid)
  training_data <- training_eid %>%
    left_join(rec_data, by = "eid")

  testing_eid <- all_eid %>%
    anti_join(training_eid, by = "eid")

  testing_data <- testing_eid %>%
    left_join(rec_data, by = "eid") %>%
    group_by(eid) %>%
    arrange(age_diag, .by_group = TRUE)

  AUC_per_ds <- matrix(nrow = para$D, ncol = 1)
  coefficients <- list()
  for(j in 1:para$D){
    ds_id <- para$list_above500occu$diag_icd10[j]
    #####################
    # create training data
    #####################
    cases_training_eid <- training_data %>% # select all the cases
      filter(diag_icd10 == ds_id) %>%
      rename(target_age = age_diag) %>%
      select(-diag_icd10)
    training_cases_data <- cases_training_eid %>% # select all disease that has a diagnosis that are earlier than target!
      left_join(training_data, by = "eid") %>%
      filter(age_diag < target_age)

    training_targe_age_distribution <- training_cases_data %>%
      group_by(eid) %>%
      dplyr::slice(1)

    training_control_eid <- training_data %>%
      anti_join(cases_training_eid, by = "eid") %>%
      group_by(eid) %>%
      dplyr::slice(1) %>%
      select(eid)

    training_control_eid$target_age <- sample(training_targe_age_distribution$target_age, dim(training_control_eid)[1], replace = T)
    training_control_data <- training_control_eid %>% # select all disease that has a diagnosis that are earlier than target!
      left_join(training_data, by = "eid") %>%
      filter(age_diag < target_age)

    training_set <- bind_rows(training_cases_data, training_control_data)
    training_cases_eid <- training_cases_data %>%
      group_by(eid) %>%
      dplyr::slice_tail(n = 1) %>%
      mutate(outcome = 1) %>%
      select(eid, outcome, age_diag)
    training_control_eid <- training_control_data %>%
      group_by(eid) %>%
      dplyr::slice_tail(n = 1) %>%
      mutate(outcome = 0) %>%
      select(eid, outcome, age_diag)
    training_eid <- bind_rows(training_cases_eid, training_control_eid)

    #####################
    # create testing data
    #####################
    testing_cases_eid <- testing_data %>%
      filter(diag_icd10 == ds_id) %>%
      rename(target_age = age_diag) %>%
      select(-diag_icd10)

    # exclude cases that have no prior diagnosis
    # (i.e. target disease happen to be the first one; no way to predict!)
    testing_control_eid <- testing_data %>%
      anti_join(testing_cases_eid, by = "eid") %>%
      group_by(eid) %>%
      dplyr::slice(1) %>%
      select(eid)

    testing_cases_data <- testing_cases_eid %>% # select all disease that has a diagnosis that are earlier than target!
      left_join(testing_data, by = "eid") %>%
      filter(age_diag < target_age)

    targe_age_distribution <- testing_cases_data %>%
      group_by(eid) %>%
      dplyr::slice(1)

    testing_control_eid$target_age <- sample(targe_age_distribution$target_age, dim(testing_control_eid)[1], replace = T)

    testing_control_data <- testing_control_eid %>% # select all disease that has a diagnosis that are earlier than target!
      left_join(testing_data, by = "eid") %>%
      filter(age_diag < target_age)

    testing_set <- bind_rows(testing_cases_data, testing_control_data)
    # save the case and control outcomes; save the age of the last event before target for computing area
    testing_cases_eid <- testing_cases_data %>%
      group_by(eid) %>%
      dplyr::slice_tail(n = 1) %>%
      mutate(outcome = 1) %>%
      select(eid, outcome, age_diag)
    testing_control_eid <- testing_control_data %>%
      group_by(eid) %>%
      dplyr::slice_tail(n = 1) %>%
      mutate(outcome = 0) %>%
      select(eid, outcome, age_diag)
    testing_eid <- bind_rows(testing_cases_eid, testing_control_eid)

    #########################################
    # join and spread the data for glmnet
    #########################################

    # we want to join data together to make sure columns are aligned
    spread_data <- bind_rows(training_set, testing_set) %>%
      mutate(age_diag = 1) %>%
      select(-target_age) %>%
      dplyr::spread(key = diag_icd10, value = age_diag, fill = 0)
    # mutate(across(2:para$D, scale)) # normalise each column

    training_spread <- training_eid %>%
      left_join(spread_data, by = "eid")
    training_spread$age_diag <- scale(training_spread$age_diag)

    testing_spread <- testing_eid %>%
      left_join(spread_data, by = "eid")
    testing_spread$age_diag <- scale(testing_spread$age_diag)

    x_train <- as.matrix(training_spread[,3:dim(training_spread)[2]])
    y_train <- training_spread$outcome

    x_test <- as.matrix(testing_spread[,3:dim(testing_spread)[2]])
    y_test <- testing_spread$outcome

    num_itr <- floor(dim(x_train)[2] / 50)
    # create a data file to save controls
    scores <- matrix(nrow = dim(x_test)[1], ncol = num_itr + 1)
    coefs <- list()
    for(itr in 1:(num_itr+1)){
      print(paste0("Disease ", ds_id, " Itr ", itr))
      if(itr == (num_itr+1)){
        if(dim(x_train)[2] == 50 * num_itr) break
        x <- x_train[,(1 + 50*num_itr):ncol(x_train)]
        predict_x <- x_test[,(1 + 50*num_itr):ncol(x_train)]
      }else{
        x <- x_train[,(1+50*(itr-1)):(50*itr)]
        predict_x <- x_test[,(1+50*(itr-1)):(50*itr)]
      }

      cvfit_testing <- cv.glmnet(x, y_train, alpha = 1, family = "binomial")
      coefs[[itr]] <- as.matrix(coef(cvfit_testing, s = cvfit_testing$lambda.min)[2:(dim(x)[2]+1)])
      row.names(coefs[[itr]]) <- row.names(coef(cvfit_testing, s = cvfit_testing$lambda.min))[2:(dim(x)[2]+1)]
      scores[,itr] <- stats::predict(cvfit_testing, newx = predict_x, s = "lambda.min")
    }
    AUC_per_ds[j,1] <- pROC::auc(y_test, rowSums(scores))
    coefficients[[j]] <- do.call(rbind, coefs)
  }
  return(list(AUC_per_ds, coefficients))
}

# using function below to create an easy to use topic loadings to compute prediction accuracy in a test set
#' Title Compute prediction odds ratio for a testing data set using pre-training ATM topic loading. Note only diseases listed in the ds_list will be used.
#' The prediction odds ratio is the odds predicted by ATM versus a naive prediction using disease probability.
#'
#' @param testing_data A data set of the same format as ATM::HES_age_example; Note: for cross-validation, split the training and testing based on individuals (eid) instead of diagnosis to avoid using training data for testing.
#' @param ds_list The order of disease code that appears in the topic loadings. This is a required input as the testing data could miss some of the records.
#' @param topic_loadings A three dimension array of topic loading in the format of ATM::UKB_HES_10topics;
#' @param max_predict The logic of prediction is using 1,..N-1 records to predict the Nth diagnosis;
#'
#' @return The returned object has four components: OR_top1, OR_top2, OR_top5 is the prediction odds ratio using top 1%, top 2%, or top 5% of ATM predicted diseases as the target set;
#' the fourth component prediction_precision is as list, with first element saves the prediction probability for 1%, 2%, 5% and 10%; additional variables saves the percentile of target disease in the ATM predicted set; for example
#' 0.03 means the target disease ranked at 3% of the diseases ordered by ATM predicted probability.
#' @export
#'
#' @examples testing_data <- ATM::HES_age_example %>% slice(1:10000)
#' new_output <- prediction_OR(testing_data, ds_list = ATM::UKB_349_disease, topic_loadings =  ATM::UKB_HES_10topics, max_predict = 5)
prediction_OR <- function(testing_data, ds_list, topic_loadings, max_predict = 10){
  # first order the incidences by age
  testing_data <- testing_data %>%
    filter(diag_icd10 %in% ds_list$diag_icd10) %>% # only keep the diseases in ds_list
    group_by(eid, diag_icd10) %>%
    arrange(age_diag, .by_group = T) %>% # keep only the first records of repeated diagnosis
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    group_by(eid) %>%
    arrange(age_diag, .by_group = TRUE)

  prediction_onebyone <- list()
  collect_prediction_ranks <- c()

  for(predict_num in 2:max_predict){
    testing_eid_included <- testing_data %>% # set incidence number threshold for inclusion in testing set
      group_by(eid) %>%
      summarise(n()) %>%
      filter(`n()` >= predict_num) %>%
      select(eid)

    testing_data_included <- testing_eid_included %>%
      left_join(testing_data, by = "eid")

    # for the last iteration, just predict all left
    if(predict_num == max_predict){
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        dplyr::slice(predict_num:n()) %>%
        dplyr::ungroup()
    }else{
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        dplyr::slice(predict_num) %>%
        dplyr::ungroup()
    }

    estimating_test_set <- testing_data_included %>%
      group_by(eid) %>%
      dplyr::slice(1:(predict_num-1)) %>%
      dplyr::ungroup()

    # use topic_init_age to make sure order of eid and para_testing$w is the same
    if(length(dim(topic_loadings)) < 3){ # if topic loading is just a matrix
      topic_loadings <- array( rep(topic_loadings, each = 81),dim =  c(81, dim(topic_loadings)[1], dim(topic_loadings)[2]) )
      # for both treeLDA and baseline lda, we only need to initialise baseline case here
      para_testing <- topic_init_baseline(estimating_test_set, ds_list, para_training$K)
    }else{
      para_testing <- topic_init_age(estimating_test_set, ds_list, dim(ATM::UKB_HES_10topics)[3], degree_free_num = 5) # just use 5, it won't matter
    }
    # directly get the alpha_z
    test_set_topic_weights <- loading2weights(estimating_test_set, ds_list = ds_list, topic_loadings)
    alpha_z <- as.matrix(select(test_set_topic_weights$incidence_weight_sum, -eid)) + 1

    # to get prediction rank of the target disease, we use a method that exclude all diseases used in the estimation set
    # we need: estimate_set: the set of diagnosis in testing set used for inferring topic_weights;
    # predict_set: the set of diseases for predicting; estimate_eid: the order of individuals in the testing set (for retrieving the topic_weights)
    # alpha_z: topic weights; beta: topic loadings; diag_icd10: the disease list in the topic loading (in case some disease is not there)
    rslt_predict <- prediction_ageLDA(estimate_set = para_testing$w, predict_set = predicting_test_set,
                                      alpha_z = alpha_z, estimate_eid = para_testing$eid,
                                      beta = topic_loadings, diag_icd10 = ds_list$diag_icd10)
    prediction_onebyone[[predict_num]]  <- rslt_predict[[1]]
    collect_prediction_ranks <- c(collect_prediction_ranks, rslt_predict[[2]])
  }

  num_disease <- dim(ds_list)[1]
  prediction_onebyone[[1]] <- list(mean(collect_prediction_ranks)/num_disease,
                                   mean( collect_prediction_ranks <= num_disease/100),
                                   mean( collect_prediction_ranks <= num_disease/50),
                                   mean( collect_prediction_ranks <= num_disease/20),
                                   mean( collect_prediction_ranks <= num_disease/10) )

  # compute the odds just based on disease prevalence
  test_prevelance <- testing_data %>%
    group_by(diag_icd10) %>%
    summarise(occ = n())
  total_num <- sum(test_prevelance$occ)
  freq_top1 <- test_prevelance %>%
    arrange(desc(occ)) %>%
    slice(1:floor(num_disease/100)) %>%
    pull(occ) %>%
    sum
  freq_top2 <- test_prevelance %>%
    arrange(desc(occ)) %>%
    slice(1:floor(num_disease/50)) %>%
    pull(occ) %>%
    sum
  freq_top5 <- test_prevelance %>%
    arrange(desc(occ)) %>%
    slice(1:floor(num_disease/20)) %>%
    pull(occ) %>%
    sum
  Odds_freq_list <- c(freq_top1, freq_top2, freq_top5)/total_num
  Odds_freq_list <- Odds_freq_list/(1-Odds_freq_list)

  prdict_OR <- list()
  # compute the disease prevalence just using testing_data.
  prdict_OR$OR_top1 <- (prediction_onebyone[[1]][[2]]/(1 - prediction_onebyone[[1]][[2]]))/Odds_freq_list[1]
  prdict_OR$OR_top2 <- (prediction_onebyone[[1]][[3]]/(1 - prediction_onebyone[[1]][[3]]))/Odds_freq_list[2]
  prdict_OR$OR_top5 <- (prediction_onebyone[[1]][[4]]/(1 - prediction_onebyone[[1]][[4]]))/Odds_freq_list[3]
  prdict_OR$prediction_precision <- prediction_onebyone

  return(prdict_OR)
}



