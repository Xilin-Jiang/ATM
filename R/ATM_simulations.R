library(MendelianRandomization)
library(parallel)
library(dplyr)
library(tidyverse)
library(gtools)
library(maxLik)
library(Rcpp)
library(glmnet)
library(pROC)
library(gtools)


# source("R/ATM_inference.R")
####################################
# basic simulation without genetic components
####################################
simulate_basic_topic_data <- function(sample_sz, topic_number, disease_number, overlap = 2){
  para <- list()
  # set the data size
  para$M <- sample_sz
  para$Ns <- floor(rexp(para$M, rate = 1/(-1.5+6.1)) + 2)
  para$K <- topic_number # start with 3 component
  para$D <- disease_number # at the moment, only simulate 6 diseases

  # initiate the variables
  para$alpha <- rgamma(para$K, shape = 50, rate = 50)
  print(paste0("simulated alpha: ", para$alpha))

  off_set <- para$D *.2 - overlap
  # each column is a topic; D*K matrix
  para$beta <- matrix(c(rep(c(rep(0.1,floor(para$D *.2)),rep(0,floor(para$D*.8 - off_set))), para$K),
                        rep(1,(para$D - (floor(para$D *.2) + floor(para$D*.8 - off_set)) )*para$K) ), nrow = para$D, ncol = para$K)
  # normalise beta
  para$beta <- apply(para$beta, 2, function(x) x/sum(x))
  para$true_beta <- para$beta
  para$theta <- rdirichlet(para$M, para$alpha)
  para$zn <- list()
  para$w <-list()
  for(s in 1:para$M){
    para$zn[[s]] <- sample( 1:para$K, para$Ns[s], replace=TRUE, prob=para$theta[s,])
    para$w[[s]] <- data.frame(Ds_id = sapply(1:para$Ns[s], function(n) sample( 1:para$D, size= 1, replace=TRUE,
                                                                               prob=para$beta[,para$zn[[s]][n]]))) %>%
      mutate(eid = s)
  }
  # build dataframes for computation
  para$unlist_Ds_id <- bind_rows(para$w)
  para$ds_list <- para$unlist_Ds_id %>%
    mutate(id = row_number()) %>%
    group_by(Ds_id) %>%
    group_split()

  para$patient_lst <- para$unlist_Ds_id %>%
    mutate(id = row_number()) %>%
    select(eid, id) %>%
    group_by(eid) %>%
    group_split(keep = F) %>%
    lapply(pull)

  # this beta_w parameter save the beta for each word: it is a list of M elements and each contain a K*Ns matrix
  # initiate beta
  para$eta <- rgamma(para$D,shape = 100, rate = 100)
  # each column is a topic; D*K matrix
  para$beta <- t(rdirichlet(para$K, para$eta))
  para$beta_w <- list()
  para$beta_w <- lapply(para$w, function(w) para$beta[w$Ds_id,,drop=FALSE] )
  para$beta_w_full <- do.call(rbind, para$beta_w)
  # Matrix of M*K
  para$E_lntheta <- t(sapply(rgamma(para$M,shape = 100, rate = 100), function(x) x*(digamma(para$alpha) - digamma(sum(para$alpha))) ))

  # update E_zn: list of M; each element is matrix of Ns*K
  para$E_zn <- list()
  para <- comp_E_zn(para)
  return(para)
}
simulate_age_topic_data <- function(sample_sz, topic_number,
                                    disease_number, degree_freedom, overlap = 2,
                                    ds_per_idv = 6.1){
  para <- list()
  # set the data size
  para$M <- sample_sz
  para$Ns <- floor(rexp(para$M, rate = 1/(-1.5+ds_per_idv)) + 2)
  para$K <- topic_number
  para$D <- disease_number
  para$P <- degree_freedom

  # initiate the variables: simulated with balanced topics
  para$alpha <- rgamma(para$K, shape = 50, rate = 50)
  print(paste0("simulated alpha: ", para$alpha))
  # each column is a topic; D*K matrix
  num_ds_each_topic <- floor(para$D/para$K) - 1

  para$beta <- array(NA, dim =  c(81, para$D, para$K) )

  para$noise_ratio <- 0.1/para$K
  for(i in 1:para$K){
    units_topic <- c(rep(0, (i - 1) *num_ds_each_topic* (81 - 29)), rep(c(rep(1,30), rep(0,42)), num_ds_each_topic))
    # truncated the diseases that went into the random ones
    if(length(units_topic) > num_ds_each_topic*i * (81 - 29)) units_topic <- units_topic[1:(num_ds_each_topic*i * (81 - 29))]
    # each disease comes from only one topic, and the last para$K disease are noise (not assigned to any single topic)
    units_topic <- matrix(c(units_topic, rep(0, num_ds_each_topic*para$K * (81 - 29) - length(units_topic)) ,
                            para$noise_ratio * runif(para$K * (81-29))  ), nrow = 81 - 29, ncol = para$D)
    units_topic <- sapply(1:(81-29), function(t) units_topic[t,]/sum(units_topic[t,])) %>% t
    # add 0 to the rows between 1-29: there is no disease below 29 years old
    units_topic <- rbind(matrix(0, nrow = 29, ncol = para$D), units_topic)
    para$beta[,,i] <- units_topic
  }
  para$true_beta <- para$beta # save beta for plotting

  para$theta <- rdirichlet(para$M, para$alpha)
  para$zn <- list()
  para$w <-list()
  for(s in 1:para$M){
    para$zn[[s]] <- sample( 1:para$K, para$Ns[s], replace=TRUE, prob=para$theta[s,])
    age_diag <- sample(30:81,para$Ns[s], replace = T)
    para$w[[s]] <- data.frame(Ds_id = sapply(1:para$Ns[s], function(n) sample( 1:para$D, size= 1, replace=TRUE,
                                                                               prob=para$beta[age_diag[n],,para$zn[[s]][n]])),
                              age_diag = age_diag) %>%
      mutate(eid = s)
  }
  para$unlist_Ds_id <- bind_rows(para$w)
  para$ds_list <- para$unlist_Ds_id %>%
    select(-eid) %>%
    mutate(id = row_number()) %>%
    group_by(Ds_id) %>%
    group_split()
  para$patient_lst <- para$unlist_Ds_id %>%
    mutate(id = row_number()) %>%
    select(eid, id) %>%
    group_by(eid) %>%
    group_split(keep = F) %>%
    lapply(pull)

  # basis is for a age grid
  para$age_max <- 81
  para$age_min <- 30
  para$age_basis <- age_basis_spline(para$P, (1:ceiling(para$age_max))/ceiling(para$age_max), para$age_min/para$age_max, para$age_max/para$age_max)
  para$age_basis_discrte <- para$age_basis[min(para$unlist_Ds_id$age_diag):max(para$unlist_Ds_id$age_diag), ]
  para$basis_phi <- para$age_basis[para$unlist_Ds_id$age_diag,,drop = F] # combined age basis

  # each column is a topic; D*K matrix
  para$beta <-  array(rnorm(para$P * para$D * para$K,sd =0.1), dim = c(para$P, para$D, para$K))
  # initiate alpha, later we will set it to 1: here it is for random initialization of E_lntheta
  para$alpha <- rgamma(para$K, shape = 50, rate = 10)

  # this beta_w parameter save the beta for each word: it is a list of M elements and each contain a K*Ns matrix
  # first compute the whole age_beta_basis, for each topic it is T-by-D
  para$exp_age_beta_basis <- array( sapply(1:para$K,
                                           function(i) exp(para$age_basis %*% para$beta[,,i])),
                                    dim=c(dim(para$age_basis)[1], para$D, para$K))
  para$sum_exp_age_beta_basis <- apply(para$exp_age_beta_basis, c(1,3), sum)
  # this is the basis of softmax function
  para$pi_beta_basis <- apply(para$exp_age_beta_basis, c(1,3), function(x) x/sum(x)) %>%
    aperm(perm = c(2,1,3))
  para$beta_w_full <- apply(para$pi_beta_basis, 3, function(x)
    x[as.matrix(select(para$unlist_Ds_id, age_diag, Ds_id))])
  para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])

  # Matrix of M*K
  para$E_lntheta <- t(sapply(rgamma(para$M,shape = 100, rate = 100), function(x) x*(digamma(para$alpha) - digamma(sum(para$alpha))) ))

  # update E_zn: list of M; each element is matrix of Ns*K
  para <- comp_E_zn(para)
  return(para)
}
###

###########################################
# Part 1:
# simulation of topic: conditional on SNP, causal disease
###########################################
# put these functions to simulation_functions.R
generate_genetics <- function(maf, n){
  u <- runif(n)
  sapply(u, function(y) ifelse(y<maf^2, 2 , ifelse(y < (2*maf - maf^2), 1,0)))
}
# all genetic effect are on the first topic only
simulate_topics <- function(topic_number, num_snp=100, pop_sz=10000, disease2topic = 0, v2t = 20, snp2t = 0.04, snp2d = 0.15, liability_thre = 0.8){
  # set the data size
  para <- list()
  para$M <- pop_sz
  para$K <- topic_number

  allele_frequency <- runif(num_snp, max = 0.5) # assuming all MAF are risk
  genetics_population <- sapply(allele_frequency, function(x) generate_genetics(x, para$M))

  # simulate the causal disease
  lia_causal <- rnorm(para$M, mean = -3) +
    genetics_population[,1:20] %*% rep(snp2d, 20)
  causal_disease <- 1 * (quantile(lia_causal, liability_thre) < lia_causal)

  # simulate topic
  para$alpha <- rgamma(para$K, shape = 50, rate = 50)
  para$theta <- rdirichlet(para$M, para$alpha)
  # causal disease
  para$theta[,1] <- para$theta[,1] + disease2topic * causal_disease - mean(disease2topic * causal_disease)
  # genetic effects
  para$theta[,1] <- para$theta[,1] + genetics_population[,21:(20+v2t)] %*% rep(snp2t, v2t) -
    mean( genetics_population[,21:(20+v2t)] %*% rep(snp2t, v2t))

  # normalise para$theta
  para$theta[,1] <- (para$theta[,1] - min(para$theta[,1]))/
    (max(para$theta[,1]) - min(para$theta[,1]))
  para$theta <- sapply(1:para$M, function(j) para$theta[j,]/sum(para$theta[j,])) %>% t
  return(list(para, genetics_population, causal_disease))
}

###########################################
# Part 2:
# simulation diseases generated by topic
###########################################

# function below will out put rec_data and ds_list for inference
simulate_genetic_disease_from_topic <- function(para, genetics_population,causal_disease,
                                                disease_number, overlap = 2,
                                                ds_per_idv = 6.1, itr_effect = 0, topic2disease = 2, v2t = 20, liability_thre = 0.8){
  # set the data size
  para$Ns <- floor(rexp(para$M, rate = 1/(-1.5+ds_per_idv)) + 2)
  # note disease number is only for the non genetic diseases
  para$D <- disease_number

  # simulate disease interaction: using 41-60 SNP
  lia_interaction <- rnorm(para$M, mean = 0) +
                             para$theta[,1] * (genetics_population[,41:60] %*% rep(itr_effect, 20))
  interact_disease <- 1 * (quantile(lia_interaction, liability_thre) <  lia_interaction)
  # simulate a disease with both topic and genetic effect (pleiotropy)
  lia1 <- rnorm(para$M, mean = -3) + para$theta[,1] * topic2disease +
          genetics_population[,21:(20+v2t)] %*% rep(0.15, v2t)
  pleiotropy_disease <-  1 * (quantile(lia1, liability_thre) < lia1)

  lia2 <- rnorm(para$M, mean = -3) + para$theta[,1] * topic2disease +
                  genetics_population[,21:(20+v2t)] %*% rep(0.15, v2t) +
                  genetics_population[,21:(20+v2t)]^2 %*% rep(0.15, v2t)
  pleiotropy_g2_disease <-  1 * (quantile(lia2, liability_thre)  < lia2)

  lia3 <- rnorm(para$M, mean = -3) + (para$theta[,1] + para$theta[,1]^2) * topic2disease +
                  genetics_population[,21:(20+v2t)] %*% rep(0.15, v2t)
  pleiotropy_t2_disease <-  1 * (quantile(lia3, liability_thre) < lia3)
  # each column is a topic; D*K matrix
  num_ds_each_topic <- floor(para$D/para$K) - 1
  para$beta <- array(NA, dim =  c(81, para$D, para$K) )
  para$noise_ratio <- 0.1/para$K

  for(i in 1:para$K){
    units_topic <- c(rep(0, (i - 1) *num_ds_each_topic* (81 - 29)), rep(c(rep(1,30), rep(0,42)), num_ds_each_topic))
    # truncated the diseases that went into the random ones
    if(length(units_topic) > num_ds_each_topic*i * (81 - 29)) units_topic <- units_topic[1:(num_ds_each_topic*i * (81 - 29))]
    # each disease comes from only one topic, and the last para$K disease are noise (not assigned to any single topic)
    units_topic <- matrix(c(units_topic, rep(0, num_ds_each_topic*para$K * (81 - 29) - length(units_topic)) ,
                            para$noise_ratio * runif(para$K * (81-29))  ), nrow = 81 - 29, ncol = para$D)
    units_topic <- sapply(1:(81-29), function(t) units_topic[t,]/sum(units_topic[t,])) %>% t
    # add 0 to the rows between 1-29: there is no disease below 29 years old
    units_topic <- rbind(matrix(0, nrow = 29, ncol = para$D), units_topic)
    para$beta[,,i] <- units_topic
  }
  para$true_beta <- para$beta # save beta for plotting

  para$zn <- list()
  para$w <-list()
  for(s in 1:para$M){
    para$zn[[s]] <- sample( 1:para$K, para$Ns[s], replace=TRUE, prob=para$theta[s,])
    age_diag <- sample(30:81,para$Ns[s], replace = T)
    para$w[[s]] <- data.frame(Ds_id = sapply(1:para$Ns[s], function(n) sample( 1:para$D, size= 1, replace=TRUE,
                                                                               prob=para$beta[age_diag[n],,para$zn[[s]][n]])),
                              age_diag = age_diag) %>%
      mutate(eid = s)
  }
  para$unlist_Ds_id <- bind_rows(para$w)
  # add causal_disease (para$D + 1) and interact_disease (para$D + 2)
  causal_ds <- data.frame(eid = 1:para$M, Ds_id = causal_disease,
                          age_diag = sample(50:70, size = para$M, replace = T)) %>%
    filter(Ds_id == 1) %>%
    mutate(Ds_id = para$D + 1)
  interact_ds <-  data.frame(eid = 1:para$M, Ds_id = interact_disease,
                             age_diag = sample(50:70, size = para$M, replace = T)) %>%
    filter(Ds_id == 1) %>%
    mutate(Ds_id = para$D + 2)
  pleiotropy_ds <-  data.frame(eid = 1:para$M, Ds_id = pleiotropy_disease,
                             age_diag = sample(50:70, size = para$M, replace = T)) %>%
    filter(Ds_id == 1) %>%
    mutate(Ds_id = para$D + 3)
  pleiotropy_g2_ds <-  data.frame(eid = 1:para$M, Ds_id = pleiotropy_g2_disease,
                               age_diag = sample(50:70, size = para$M, replace = T)) %>%
    filter(Ds_id == 1) %>%
    mutate(Ds_id = para$D + 4)
  pleiotropy_t2_ds <-  data.frame(eid = 1:para$M, Ds_id = pleiotropy_t2_disease,
                                  age_diag = sample(50:70, size = para$M, replace = T)) %>%
    filter(Ds_id == 1) %>%
    mutate(Ds_id = para$D + 5)
  para$unlist_Ds_id <- bind_rows(para$unlist_Ds_id, causal_ds, interact_ds, pleiotropy_ds, pleiotropy_g2_ds, pleiotropy_t2_ds)

  # create the rec and ds_list  rec_data <- para_orginal_sim$unlist_Ds_id %>%
  rec_data <- para$unlist_Ds_id %>%
    rename(diag_icd10 = Ds_id) %>%
    select(eid, diag_icd10, age_diag) %>%
    group_by(eid, diag_icd10) %>%
    arrange(age_diag, .by_group = T) %>% # keep only the first records of repeated diagnosis
    slice(1) %>%
    ungroup()

  ds_list <- rec_data %>%
    group_by(diag_icd10) %>%
    summarise(occ = n())

  return(list(rec_data, ds_list, interact_disease, pleiotropy_disease))
}

inf_genetic_simulation <- function(rec_data, ds_list, topic_number){
  para <- topic_init_age(rec_data, ds_list, topic_number, 3)

  para$max_itr <- 500
  para$alpha <- rep(1, para$K)
  para$lb <- data.frame("Iteration" = 0,"Lower_bound" = comp_lda_lb(para))
  para$itr_beta <- 1
  para$itr_check_lb <- 1
  para$tol <- 10^(-6)
  cvb_tag <- T
  cvb0_tag <- T
  for(itr in 1:para$max_itr){
    print(paste0("Interation: ",itr))
    for(itr_inside in 1:para$itr_beta){ # in practice we perform quick steps a few times before move on.
      if(cvb_tag){
        if(cvb0_tag){
          para <- CVB0_E_zn(para)
        }else{
          para <- CVB_E_zn(para)
        }
      }else{
        para <- comp_E_zn(para)
        para <- comp_E_lntheta(para)
      }
    }
    para <- fast_update_age_depend_lda(para)
    # para <- update_alpha(para) # we use non-informative alpha
    if(cvb_tag){
      para$lb[nrow(para$lb) + 1,] <- c(itr, CVB_lb(para))
    }else{
      para$lb[nrow(para$lb) + 1,] <- c(itr, comp_lda_lb(para))
    }
    if(itr %% para$itr_check_lb ==0){
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
  }
  return(para)
}
##################################################
# test for SNP x topic interaction
##################################################
test_gxT <- function(patient_loadings, interact_disease, genetics_population){
  # perform int
  p_gxT <- rep(NA, dim(genetics_population)[2])
  for(g_id in 1:dim(genetics_population)[2]){
    gxT <- glm(interact_disease ~ patient_loadings * genetics_population[,g_id],
               family = binomial)
    p_gxT[g_id] <- summary(gxT)$coefficients[4,4]
  }
  return(p_gxT)
}
test_gxT_T2 <- function(patient_loadings, interact_disease, genetics_population){
  # perform int
  p_gxT <- rep(NA, dim(genetics_population)[2])
  for(g_id in 1:dim(genetics_population)[2]){
    non_linear_topic <- patient_loadings^2
    gxT <- glm(interact_disease ~ patient_loadings + non_linear_topic + genetics_population[,g_id] + patient_loadings : genetics_population[,g_id],
               family = binomial)
    p_gxT[g_id] <- summary(gxT)$coefficients[5,4]
  }
  return(p_gxT)
}

test_gxD <- function(patient_loadings, interact_disease, genetics_population){
  loading_case <- patient_loadings[which(interact_disease == 1)]
  loading_control <- patient_loadings[which(interact_disease == 0)]
  genetics_case <- genetics_population[which(interact_disease == 1), , drop = F]
  genetics_control <- genetics_population[which(interact_disease == 0), ,drop = F]
  z_gxT <- rep(NA, dim(genetics_population)[2])
  for(g_id in 1:dim(genetics_population)[2]){
    case_lm <- lm(loading_case ~ genetics_case[,g_id])
    control_lm <- lm(loading_control ~ genetics_control[,g_id])
    z_gxT[g_id] <- (summary(control_lm)$coefficients[2,1] - summary(case_lm)$coefficients[2,1])/
      sqrt(summary(control_lm)$coefficients[2,2]^2 + summary(case_lm)$coefficients[2,2]^2)
  }
  p_value <- 2*pnorm(q=abs(z_gxT), lower.tail=FALSE)
  return(p_value)
}

test_gxD_simple <- function(patient_loadings, interact_disease, genetics_population){
  p_gxD <- rep(NA, dim(genetics_population)[2])
  for(g_id in 1:dim(genetics_population)[2]){
    gxD <- lm(patient_loadings ~ interact_disease * genetics_population[,g_id])
    p_gxD[g_id] <- summary(gxD)$coefficients[4,4]
    }
  return(p_gxD)
}
##################################################
# MR analysis
##################################################
eggerMR <- function(loadings, disease, genetics, method = "ivw") {
  d_pvalue <- rep(NA, dim(genetics)[2])
  d_effects <- rep(NA, dim(genetics)[2])
  d_se <- rep(NA, dim(genetics)[2])
  t_pvalue <- rep(NA, dim(genetics)[2])
  t_effects <- rep(NA, dim(genetics)[2])
  t_se <- rep(NA, dim(genetics)[2])
  for(g_id in 1:dim(genetics)[2]){
    g2T <- lm(loadings ~ genetics[,g_id])
    t_pvalue[g_id] <- summary(g2T)$coefficients[2,4]
    t_effects[g_id] <- summary(g2T)$coefficients[2,1]
    t_se[g_id] <- summary(g2T)$coefficients[2,2]

    g2D <- glm(disease ~ genetics[,g_id],
               family = binomial)
    d_pvalue[g_id] <- summary(g2D)$coefficients[2,4]
    d_effects[g_id] <- summary(g2D)$coefficients[2,1]
    d_se[g_id] <- summary(g2D)$coefficients[2,2]
  }

  IV_topic <- which(p.adjust(t_pvalue,method = "fdr") < 0.05)
  IV_disease <- which(p.adjust(d_pvalue,method = "fdr") < 0.05)

  if(length(IV_topic) > 0){
    t2d_MRdata <- mr_input(bx = t_effects[IV_topic], bxse = t_se[IV_topic],
                           by = d_effects[IV_topic], byse = d_se[IV_topic])
    if(method == "egger"){
      topic2diseaseMR <- mr_egger(t2d_MRdata)
    }else if(method == "ivw"){
      topic2diseaseMR <- mr_ivw(t2d_MRdata)
    }

  }else{
    topic2diseaseMR <- NULL
  }

  if(length(IV_disease) > 0){
    d2t_MRdata <- mr_input(bx = d_effects[IV_disease], bxse = d_se[IV_disease],
                           by = t_effects[IV_disease], byse = t_se[IV_disease])
    if(method == "egger"){
      disease2topicMR <- mr_egger(d2t_MRdata)
    }else if(method == "ivw"){
      disease2topicMR <- mr_ivw(d2t_MRdata)
    }

  }else{
    disease2topicMR <- NULL
  }
  MR_results <- list(topic2diseaseMR, disease2topicMR)
  return(MR_results)
}

