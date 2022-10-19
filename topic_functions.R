library(parallel)
library(dplyr)
library(tidyverse)
library(gtools)
library(maxLik)
library(Rcpp)
library(glmnet)
library(pROC)
library(PRROC)
library(gtools)
########################################
# initialize parameters for the model
###########################################
topic_init_baseline <- function(rec_data, ds_list, topic_num){
  first_incidence_age <- rec_data %>%
    arrange(eid) 
  
  # plot the number distribution of indiviudal diseases
  df_number_records <- first_incidence_age %>%
    group_by(eid) %>%
    summarise(n())
  para <- list()
  
  para$eid <- df_number_records$eid
  para$list_above500occu <- ds_list 
  para$D <- dim(para$list_above500occu)[1] # disease number 
  para$M <- length(para$eid) # subject number 
  para$K <- topic_num # start with 10 component 
  para$Ns <- df_number_records$`n()`
  
  code2id <- function(x){
    return( match(x, para$list_above500occu$diag_icd10))
  }
  
  # here I am rounding the disease time to year for computation efficiency
  para$unlist_Ds_id <- first_incidence_age %>%
    mutate(Ds_id = code2id(diag_icd10)) %>%
    select(-diag_icd10) %>%
    mutate(age_diag = round(age_diag)) 
  
  # the patient_list provide the column index for efficiently breaking down matrix into list of matrices
  para$patient_lst <- para$unlist_Ds_id %>%
    mutate(id = row_number()) %>%
    select(eid, id) %>%
    group_by(eid) %>%
    group_split(keep = F) %>%
    lapply(pull)
  
  para$w <- para$unlist_Ds_id %>%
    group_by(eid) %>%
    group_split(keep = F)
  
  # this list is splitted by disease
  para$ds_list <- para$unlist_Ds_id %>%
    select(-eid) %>%
    mutate(id = row_number()) %>%
    group_by(Ds_id) %>%
    group_split()
  
  print(object.size(para$w), unit = "MB" , standard = "SI")
  
  # initiate beta
  para$eta <- rgamma(para$D,shape = 100, rate = 100)
  # each column is a topic; D*K matrix
  para$beta <- t(rdirichlet(para$K, para$eta))
  
  # initiate alpha
  para$alpha <- rgamma(para$K, shape = 50, rate = 10)
  
  # this beta_w parameter save the beta for each word: it is a list of M elements and each contain a K*Ns matrix
  para$beta_w_full <- para$beta[para$unlist_Ds_id$Ds_id,,drop=FALSE] 
  para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])
  
  # Matrix of M*K
  para$E_lntheta <- t(sapply(rgamma(para$M,shape = 100, rate = 100), function(x) x*(digamma(para$alpha) - digamma(sum(para$alpha))) ))
  
  # update E_zn: list of M; each element is matrix of Ns*K
  para <- comp_E_zn(para)
  
  # eventually reset alpha to be non-informative
  para$alpha <- rep(1, para$K)
  
  return(para)
}
topic_init_age <- function(rec_data, ds_list, topic_num, degree_free_num) {
  # arrange the data stack by individual is very important as we need to make sure the matrix could be rejoined into a single matrix
  first_incidence_age <- rec_data %>%
    arrange(eid) 
  
  # plot the number distribution of indiviudal diseases
  df_number_records <- first_incidence_age %>%
    group_by(eid) %>%
    summarise(n())

  para <- list()
  
  para$eid <- df_number_records$eid
  # add death to it
  para$list_above500occu <- ds_list 
  para$D <- dim(para$list_above500occu)[1] # disease number 
  para$M <- length(para$eid) # subject number 
  para$K <- topic_num # start with 10 component 
  para$P <- degree_free_num # degrees of freedom
  # also need to compute record number per individual for computing lower bound of cvb
  para$Ns <- df_number_records$`n()`
  
  code2id <- function(x){
    return( match(x, para$list_above500occu$diag_icd10))
  }
  
  # here I am rounding the disease time to year for computation efficiency
  para$unlist_Ds_id <- first_incidence_age %>%
    mutate(Ds_id = code2id(diag_icd10)) %>%
    select(-diag_icd10) %>%
    mutate(age_diag = round(age_diag)) 
  # the patient_list provide the column index for efficiently breaking down matrix into list of matrices
  para$patient_lst <- para$unlist_Ds_id %>%
    mutate(id = row_number()) %>%
    select(eid, id) %>%
    group_by(eid) %>%
    group_split(keep = F) %>%
    lapply(pull)
  
  para$w <- para$unlist_Ds_id %>%
    group_by(eid) %>%
    group_split(keep = F)
  
  # this list is column ID for each disease
  para$ds_list <- para$unlist_Ds_id %>%
    select(-eid) %>%
    mutate(id = row_number()) %>%
    group_by(Ds_id) %>%
    group_split()
  
  print(object.size(para$w), unit = "MB" , standard = "SI")
  
  # create an age matrix for each disease, it is a list same length as para$w, each of Ns-by-F matrix
  para$age_max <- first_incidence_age %>% 
    summarise(max(age_diag)) %>%
    pull
  para$age_min <- first_incidence_age %>% 
    summarise(min(age_diag)) %>%
    pull
  
  # I scale the magnitude of age to 1 (divided by para$age_max) to avoid numeric infinity in exponentials
  # basis is for a age grid
  para$age_basis <- age_basis_spline(para$P, (1:ceiling(para$age_max))/ceiling(para$age_max), para$age_min/para$age_max, para$age_max/para$age_max)
  # below is the set of age basis that are used for discrete computation
  para$age_basis_discrte <- para$age_basis[min(para$unlist_Ds_id$age_diag):max(para$unlist_Ds_id$age_diag), ]
  para$basis_phi <- para$age_basis[para$unlist_Ds_id$age_diag,,drop = F]
  
  # each column is a topic; D*K matrix
  para$beta <- array(rnorm(para$P * para$D * para$K,sd =0.1), dim = c(para$P, para$D, para$K))
  
  # initiate alpha, later we will set it to 1: here it is for random initialization of E_lntheta
  para$alpha <- rgamma(para$K, shape = 50, rate = 10)
  
   # first compute the whole age_beta_basis, for each topic it is T-by-D
  para$exp_age_beta_basis <- array( sapply(1:para$K, 
                                           function(i) exp(para$age_basis %*% para$beta[,,i])),
                                    dim=c(dim(para$age_basis)[1], para$D, para$K))
  para$sum_exp_age_beta_basis <- apply(para$exp_age_beta_basis, c(1,3), sum)
  # this is the basis of softmax function
  para$pi_beta_basis <- apply(para$exp_age_beta_basis, c(1,3), function(x) x/sum(x)) %>% 
    aperm(perm = c(2,1,3))
  # this beta_w parameter save the beta for each word: it is a list of M elements and each contain a K*Ns matri
  para$beta_w_full <- apply(para$pi_beta_basis, 3, function(x) 
    x[as.matrix(select(para$unlist_Ds_id, age_diag, Ds_id))]) 
  para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])
  
  # Matrix of M*K
  para$E_lntheta <- t(sapply(rgamma(para$M,shape = 100, rate = 100), function(x) x*(digamma(para$alpha) - digamma(sum(para$alpha))) ))
  
  # update E_zn: list of M; each element is matrix of Ns*K
  para <- comp_E_zn(para)
  
  # eventually reset alpha to be non-informative
  para$alpha <- rep(1, para$K)
  return(para)
}
topic_init_tree <- function(rec_data, ds_list, topic_num, tree_leveles){
  # arrange the data stack by individual is very important as we need to make sure the matrix could be rejoined into a single matrix
  first_incidence_age <- rec_data %>%
    arrange(eid) 
  
  # plot the number distribution of indiviudal diseases
  df_number_records <- first_incidence_age %>%
    group_by(eid) %>%
    summarise(n())
  
  para <- list()
  para$eid <- df_number_records$eid
  
  # add death to it
  para$list_above500occu <- ds_list
  # para$list_above500occu <- read.csv(paste0("listAbove500.csv"))
  para$D <- dim(para$list_above500occu)[1] # disease number 
  para$M <- length(para$eid) # subject number 
  para$K <- topic_num # start with 10 component 
  para$L <- tree_leveles # for ICD10 it is 4 layers 
  # also need to compute record number per individual for cvb
  para$Ns <- df_number_records$`n()`
  
  code2id <- function(x){
    return( match(x, para$list_above500occu$diag_icd10))
  }
  
  # here I am rounding the disease time to year for computation efficiency
  para$unlist_Ds_id <- first_incidence_age %>%
    mutate(Ds_id = code2id(diag_icd10)) %>%
    select(-diag_icd10) %>%
    mutate(age_diag = round(age_diag)) 
  
  # the patient_list provide the column index for efficiently breaking down matrix into list of matrices
  para$patient_lst <- para$unlist_Ds_id %>%
    mutate(id = row_number()) %>%
    select(eid, id) %>%
    group_by(eid) %>%
    group_split(keep = F) %>%
    lapply(pull)
  
  para$w <- para$unlist_Ds_id %>%
    group_by(eid) %>%
    group_split(keep = F)
  
  # this list is splited by disease
  para$ds_list <- para$unlist_Ds_id %>%
    select(-eid) %>%
    mutate(id = row_number()) %>%
    group_by(Ds_id) %>%
    group_split()
  
  print(object.size(para$w), unit = "MB" , standard = "SI")
  
  # initiate beta
  para$eta <- rgamma(para$D,shape = 100, rate = 100)
  
  # each column is a topic; D*K matrix
  para$beta <- t(rdirichlet(para$K, para$eta))
  
  # initiate alpha
  para$alpha <- rgamma(para$K, shape = 50, rate = 10)
  
  # this beta_w parameter save the beta for each word: it is a list of M elements and each contain a K*Ns matrix
  para$beta_w_full <- para$beta[para$unlist_Ds_id$Ds_id,,drop=FALSE] 
  para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])
  
  # Matrix of M*K
  para$E_lntheta <- t(sapply(rgamma(para$M,shape = 100, rate = 100), function(x) x*(digamma(para$alpha) - digamma(sum(para$alpha))) ))
  
  # update E_zn: list of M; each element is matrix of Ns*K
  para <- comp_E_zn(para)
  
  # eventually reset alpha to be non-informative
  para$alpha <- rep(1, para$K)
  
  # prepare the treeLDA
  para <- create_structure(para)
  return(para)
} 

#########################################
# functions for nominate comorbidities
#########################################
# function that making use of inferred topic to estimate individual weights 
topics2weights <- function(data, ds_list, degree_freedom, topics){
  para <- topic_init_age(data, ds_list, dim(topics)[length(dim(topics))], degree_freedom)
  # update beta_w: list of Ns-by-K 
  para$beta_w_full <- apply(topics, 3, function(x) 
    x[as.matrix(select(para$unlist_Ds_id, age_diag, Ds_id))]) 
  para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])
  
  # update z_n until convergence
  para$max_itr <- 10
  para$lb <- data.frame("Iteration" = 0,"Lower_bound" = CVB_lb(para))
  para$tol <- 10^(-6)
  for(itr in 1:para$max_itr){
    print(paste0("Interation: ",itr))
    para <- CVB0_E_zn(para) # we choose CVB0 as papers shown it could converge quicker
    para$lb[nrow(para$lb) + 1,] <- c(itr, CVB_lb(para))
    curr_lb <- pull(filter(para$lb, Iteration == itr), Lower_bound) 
    prev_lb <- pull(filter(para$lb, Iteration == (itr - 1 )), Lower_bound) 
    print(paste0("Current Lower bound ", curr_lb, " at iteration: ",itr))
    try({
      if(is.finite((curr_lb - prev_lb)) & abs(curr_lb - prev_lb)/abs(prev_lb) < para$tol ){
        print(paste0("Optimization converged at step ", itr))
        break
      }
    })
  }
  return(para)
}

###########################################
# implementing baseline LDA: mean-field VB
###########################################
comp_lda_lb <- function(para){
  # compute the lower bound for the whole dataset
  
  # terms for cross-entropy
  term1 <- para$M*(lgamma(sum(para$alpha)) - sum(lgamma(para$alpha)) ) +
    sum( para$E_lntheta %*% (para$alpha - 1) )
  term2 <- sapply(1:para$M,
                  function(s) sum( para$E_zn[[s]] %*% para$E_lntheta[s,] )) %>% sum
  # term3 <- sum(sapply(1:para$M, function(s) sum(para$E_zn[[s]] * log(para$beta_w[[s]]) ) ) )
  # term3 could be fully vectorized using para$unlist_zn and para$beta_w_full
  # term3 <- sum(para$unlist_zn * log(para$beta_w_full) ) 
  term3 <- sum( log(para$beta_w_full^para$unlist_zn) ) # use this method to avoid numeric issue
  
  # terms for entropy
  term4 <- sum(sapply(1:para$M,
                      function(s) lgamma(sum(para$alpha_z[s,])) - sum(lgamma(para$alpha_z[s,])))) +
    sum((para$alpha_z - 1)*para$E_lntheta)
  # term5 <- sum(sapply(1:para$M, function(s) sum( para$E_zn[[s]]*log(para$E_zn[[s]]) ) ) )
  # term5 could be fully vectorized using para$unlist_zn
  # term5 <- sum( para$unlist_zn*log(para$unlist_zn) ) 
  term5 <- sum( log(para$unlist_zn^para$unlist_zn) ) # similarly, avoid numeric issue
  
  return((term1 + term2 + term3 - term4 - term5))
}

comp_E_zn <- function(para){
  # compute E-step for zn
  # two layer of apply: outer layer is for documents s=1,..,M; saplly layer is for words in documents n=1,2,...Ns
  para$E_zn <-sapply(1:para$M, 
                     function(s)
                       (para$beta_w[[s]] %*% diag(exp(para$E_lntheta[s,]))  )/
                       (para$beta_w[[s]]) %*% t(exp(para$E_lntheta[rep(s,para$K),])),
                     simplify = FALSE)
  # update the variables to facilitate computation
  para$alpha_z <- sapply(1:para$M, function(s) para$alpha + colSums(para$E_zn[[s]])) %>% t
  para$unlist_zn <- do.call(rbind, para$E_zn)
  return(para)
}

comp_E_lntheta <- function(para){
  # compute E-step for ln(theta)
  para$E_lntheta <- sapply(1:para$M,
                           function(s) digamma(para$alpha_z[s,]) - 
                             digamma(sum(para$alpha_z[s,]))) %>% t 
  return(para)
}

update_beta_basic_lda <- function(para){
  # compute M-step for beta: basic case, direct maximize the upper bound
  para$beta <- sapply(1:para$D, function(j) colSums(para$unlist_zn[para$ds_list[[j]]$id,]) ) %>% t 
  # normalize beta
  para$beta <- sapply(1:para$K, function(i) para$beta[,i]/sum(para$beta[,i]))
  # this beta_w parameter save the beta for each word: it is a list of M elements and each contain a K*Ns matrix
  # para$beta_w <- lapply(para$w, function(w) para$beta[w$Ds_id,,drop=FALSE] )
  para$beta_w_full <- para$beta[para$unlist_Ds_id$Ds_id,,drop=FALSE]
  para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])
  return(para)
}

update_alpha <- function(para){
  # compute M-step for alpha; optimize the dirichlet with Newton-Raphson method
  para$lb_alpha <- function(alpha){
    para$M*(lgamma(sum(alpha)) - sum(lgamma(para$alpha)) ) +
      sum(colSums(para$E_lntheta) * (para$alpha - 1) )
  }
  
  para$grad_alpha <- function(alpha){
    para$M*(digamma(sum(alpha)) - digamma(para$alpha))  +
      colSums( para$E_lntheta )
  }
  
  para$hess_alpha <- function(alpha){
    para$M*( trigamma(sum(alpha)) - diag(trigamma(alpha)) )
  }
  # para$optim_alpha <- maxBFGS(para$lb_alpha, grad = para$grad_alpha, start = para$alpha, 
  #                             constraints=list(ineqA = diag(nrow = para$K), ineqB = matrix(0,nrow = para$K)) )
  para$optim_alpha <- maxNR(para$lb_alpha, grad = para$grad_alpha, hess = para$hess_alpha, start = para$alpha) # rep(10^(-5), para$K))
  
  para$alpha <- para$optim_alpha$estimate
  
  return(para)
}

##########################################################################
# implimenting collapse variational bayes
##########################################################################
# compute lower bound on the full likelihood
CVB_lb <- function(para){
  # compute the lower bound for the whole dataset
  
  # terms for cross-entropy
  term1 <- para$M*(lgamma(sum(para$alpha)) - sum(lgamma(para$alpha)) ) -
    sum( lgamma(sum(para$alpha) + para$Ns) )
  
  # we used approximation to obtain the expectation
  alpha_n_0 <- sapply(1:para$M, 
                      function(s) colSums(para$E_zn[[s]]) + para$alpha) %>% t
  # 2. compute the variance of z_n (excluding n)
  var_zn <- sapply(1:para$M, 
                   function(s) para$E_zn[[s]] * (1 - para$E_zn[[s]]),simplify = F)
  var_neg_zn <- sapply(1:para$M, 
                       function(s) colSums(var_zn[[s]])) %>% t
  term2 <- sum( (alpha_n_0 - .5)*log(alpha_n_0) - alpha_n_0 + 1/(12*alpha_n_0) + 
    var_neg_zn * (1/(2*alpha_n_0) + 1/(4*alpha_n_0^2) + 1/(12*alpha_n_0^3)) + 
    .5 *log(2*pi) )
   
  # term3 could be fully vectorized using para$unlist_zn and para$beta_w_full
  # term3 <- sum(para$unlist_zn * log(para$beta_w_full) ) 
  beta_w_full <- (para$beta_w_full + 10^(-256)) # to deal with numeric infinity 
  term3 <- sum( log(beta_w_full^para$unlist_zn) ) # use this method to avoid numeric issue
  
  # terms for entropy
  term4 <- sum( log(para$unlist_zn^para$unlist_zn) ) # similarly, avoid numeric issue
  
  return((term1 + term2 + term3 - term4))
}
# compute the zn update for collpased likelihood
CVB_E_zn <- function(para){
  # compute CVB E-step for zn
  # 1. compute the sum of z_n (excluding n)
  alpha_n_0 <- sapply(1:para$M, 
                function(s) t(colSums(para$E_zn[[s]]) + para$alpha - t(para$E_zn[[s]])),
                simplify = FALSE)
  # 2. compute the variance of z_n (excluding n)
  var_zn <- sapply(1:para$M, 
                   function(s) para$E_zn[[s]] * (1 - para$E_zn[[s]]),
                   simplify = FALSE)
  var_neg_zn <- sapply(1:para$M, 
                       function(s) t(colSums(var_zn[[s]]) - t(var_zn[[s]])),
                       simplify = FALSE)
  
  ##################
  # multiply everything together to get the distribution
  ###################
  para$E_zn <- sapply(1:para$M, 
                     function(s)
                       alpha_n_0[[s]] * para$beta_w[[s]] * exp(- var_neg_zn[[s]]/(2*alpha_n_0[[s]]^2)),
                     simplify = FALSE)
  # normalise the E_zn; add a small value to avoid NAs
  para$E_zn <- sapply(1:para$M, 
                      function(s) sweep(para$E_zn[[s]] + 10^(-256), 1, rowSums(para$E_zn[[s]]), FUN="/"),
                      simplify = FALSE)

  # update the variables to facilitate computation
  para$alpha_z <- sapply(1:para$M, function(s) para$alpha + colSums(para$E_zn[[s]])) %>% t
  para$unlist_zn <- do.call(rbind, para$E_zn)
  return(para)
}
# ingnore the second order information
# compute the zn update for collpased likelihood
CVB0_E_zn <- function(para){
  # compute CVB E-step for zn
  # compute the sum of z_n (excluding n)
  alpha_n_0 <- sapply(1:para$M, 
                           function(s) t(colSums(para$E_zn[[s]]) + para$alpha - t(para$E_zn[[s]])),
                           simplify = FALSE)
  para$E_zn <- sapply(1:para$M, 
                      function(s)
                        alpha_n_0[[s]] * para$beta_w[[s]],
                      simplify = FALSE)
  # normalise the E_zn; add a small value to avoid NAs
  para$E_zn <- sapply(1:para$M, 
                      function(s) sweep(para$E_zn[[s]] + 10^(-256), 1, rowSums(para$E_zn[[s]]), FUN="/"),
                      simplify = FALSE)
  ######################################
  # testing the new function!!!
  ######################################
  # update the variables to facilitate computation
  para$alpha_z <- sapply(1:para$M, function(s) para$alpha + colSums(para$E_zn[[s]])) %>% t
  para$unlist_zn <- do.call(rbind, para$E_zn)
  return(para)
}


###################################################################################################
# implementing age dependent topics
# for point estimate of beta, we don't need to change any other updates including the lower bound
####################################################################################################
age_basis_spline <- function(P, age_lst, age_min, age_max){
  # P is basis degree of freedom, age_lst is the ages that will be transform to the basis
  # age_min, age_max are the lowest and largest age in the data set which is used to compute
  X_base <- cbind(age_lst^0, age_lst^1, age_lst^2, age_lst^3)
  if(P <= 4){
    return(X_base[,1:P,drop = F])
  }else{
    X <- sapply(age_min+(1:(P-4))*(age_max-age_min)/(P-3), function(x) pmax((age_lst-x)^3,0))
    return(cbind(X_base, X))
  }
}

# this is the lower bound for term3 after introducing zeta and pi
lb_beta_zeta <- function(para){
  term <- sapply(1:para$K, 
                 function(i)sapply(1:para$D, function(j) 
                   t(para$basis_phi[para$ds_list[[j]]$id, ] %*% para$beta[, j, i] -
                                   (rowSums(exp(para$basis_phi[para$ds_list[[j]]$id, ] %*% para$beta[, , i]) ))/
                                   para$zeta_full[para$ds_list[[j]]$id, i] -
                                   log(para$zeta_full[para$ds_list[[j]]$id, i]) + 1) %*% 
                     para$unlist_zn[para$ds_list[[j]]$id, i] ))  %>% sum
  return(term)
}  

# compute the newton method for roots
fun_age_beta <- function(beta_ij,i_fn,j_fn, para){
  para$phi_z[[j_fn]][,i_fn] %*% beta_ij - 
    z_zeta_exp_phi_beta(para$basis_phi, beta_ij, para$z_zeta[, i_fn]) -
    1/200 * crossprod(beta_ij) # add a prior/regularisation, just using a standard normal to control the scale
}
grad_age_beta <- function(beta_ij,i_gd,j_gd, para){
  para$phi_z[[j_gd]][,i_gd] -
    phi_z_zeta_exp_phi_beta(x = para$basis_phi,beta = beta_ij, phi_z_zeta = para$phi_z_zeta[[i_gd]]) - 
    0.01 * beta_ij # add a prior/regularisation
}

update_age_depend_lda <- function(para){
  # first compute a full matrix then split: this is complicated but is much more efficient than directly lapply of para$w
  # compute M-step for beta: basic case, direct maximize the upper bound

  # # If I update zeta between each j, the optimization will be slower....
  # for(j in 1:para$D){
  #   # first compute the whole age_beta_basis, for each topic it is T-by-D
  #   para$exp_age_beta_basis <- array( sapply(1:para$K, 
  #                                            function(i) exp(para$age_basis %*% para$beta[,,i])),
  #                                     dim=c(dim(para$age_basis)[1], para$D, para$K))
  #   para$sum_exp_age_beta_basis <- apply(para$exp_age_beta_basis, c(1,3), sum)
  #   # compute the variational parameter zeta nsi
  #   para$zeta_full <- para$sum_exp_age_beta_basis[para$unlist_Ds_id$age_diag,,drop=F]
  #   para$zeta <- lapply(para$patient_lst, function(x) para$zeta_full[x,,drop=F]) # we actually only need zeta as a full list
  #   para$z_zeta <- para$unlist_zn/para$zeta_full # this term helps speed up
  #   
  #   para$beta[,j,] <- sapply(1:para$K, function(i) 
  #     optim(par = para$beta[,j,i],
  #           fn = function(x) fun_age_beta(x,i,j, para), 
  #           gr = function(x) grad_age_beta(x,i,j, para), method ="BFGS",
  #           control = list(fnscale = -1) )$par )
  # }
  para$zeta_full <- para$sum_exp_age_beta_basis[para$unlist_Ds_id$age_diag,,drop=F]
  # para$zeta <- lapply(para$patient_lst, function(x) para$zeta_full[x,,drop=F]) # we actually only need zeta as a full list
  para$z_zeta <- para$unlist_zn/para$zeta_full # this term helps speed up
  para$phi_z_zeta <- lapply(1:para$K, function(i) para$basis_phi*para$z_zeta[,i]) # help to speed up
  para$phi_z <- lapply(1:para$D, function(j) 
    crossprod(para$basis_phi[para$ds_list[[j]]$id, ,drop = F], para$unlist_zn[para$ds_list[[j]]$id, ]) )  # speed up
  
  # para$beta <- array(sapply(1:para$K, function(i) sapply(1:para$D, function(j) 
  #   optim(par = para$beta[,j,i],
  #         fn = function(x) fun_age_beta(x,i,j, para), 
  #         gr = function(x) grad_age_beta(x,i,j, para), method ="BFGS",
  #         control = list(fnscale = -1) )$par )), 
  #   dim = c(para$P, para$D, para$K))
  #############################################################
  # # using following piece to avoid memory blowing up
  #############################################################
  for(j in 1:para$D){
    # print(paste0("update beta: ", j))
    for(i in 1:para$K){
      fail_init <- (sum(!is.finite(c(fun_age_beta(para$beta[,j,i],i,j, para), grad_age_beta(para$beta[,j,i],i,j, para)) ) ) != 0)
      if(fail_init){ # handle the initialization case when gradient are numerically infeasible
        print("Warning: Handle infeasible initialization, should only be occuring at first few interations.")
        para$beta <-  array(rnorm(para$P * para$D * para$K,sd =0.1), dim = c(para$P, para$D, para$K))
        para$E_lntheta <- t(sapply(rgamma(para$M,shape = 100, rate = 100), function(x) x*(digamma(para$alpha) - digamma(sum(para$alpha))) ))
        break
      }
      para$beta[,j,i] <- optim(par = para$beta[,j,i],
              fn = function(x) fun_age_beta(x,i,j, para),
              gr = function(x) grad_age_beta(x,i,j, para), method ="BFGS",
              control = list(fnscale = -1) )$par
    }
    if(fail_init){break} # break from the iner circle as well 
  }
  
  # compute the beta_w: exponential divided by sum
  # first compute the whole age_beta_basis, for each topic it is T-by-D
  para$exp_age_beta_basis <- array( sapply(1:para$K, 
                                           function(i) exp(para$age_basis %*% para$beta[,,i])),
                                    dim=c(dim(para$age_basis)[1], para$D, para$K))
  para$sum_exp_age_beta_basis <- apply(para$exp_age_beta_basis, c(1,3), sum)
  # this is the basis of softmax function
  para$pi_beta_basis <- apply(para$exp_age_beta_basis, c(1,3), function(x) x/sum(x)) %>% 
    aperm(perm = c(2,1,3))
  
  # update beta_w: list of Ns-by-K 
  para$beta_w_full <- apply(para$pi_beta_basis, 3, function(x) 
    x[as.matrix(select(para$unlist_Ds_id, age_diag, Ds_id))]) 
  para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])

  return(para)
}

# discretise age for a fast computation
fun_age_beta_fast <- function(beta_ij,phi_z,z_zeta,age_bases){
  phi_z %*% beta_ij - 
    z_zeta %*% exp(age_bases %*% beta_ij) -
    1/200 * crossprod(beta_ij) # add a prior/regularisation, just using a standard normal to control the scale
}
grad_age_beta_fast <- function(beta_ij,phi_z,z_zeta,age_bases){
  phi_z -
    crossprod(age_bases, z_zeta * exp(age_bases %*% beta_ij) )  - 
    0.01 * beta_ij # add a prior/regularisation
}

fast_update_age_depend_lda <- function(para){
  # this is a much more 
  para$zeta_full <- para$sum_exp_age_beta_basis[para$unlist_Ds_id$age_diag,,drop=F]
  # para$zeta <- lapply(para$patient_lst, function(x) para$zeta_full[x,,drop=F]) # we actually only need zeta as a full list
  para$z_zeta <- para$unlist_zn/para$zeta_full # this term helps speed up
  para$z_zeta_sum_by_age <- sapply(min(para$unlist_Ds_id$age_diag):max(para$unlist_Ds_id$age_diag), 
                                   function(y) colSums(para$z_zeta[which(para$unlist_Ds_id$age_diag == y),,drop=F]) ) %>% t
  para$phi_z <- lapply(1:para$D, function(j) 
    crossprod(para$basis_phi[para$ds_list[[j]]$id, ,drop = F], para$unlist_zn[para$ds_list[[j]]$id, ]) )  # speed up
  
  #############################################################
  # # using following piece to avoid memory blowing up
  #############################################################
  for(j in 1:para$D){
    # print(paste0("update beta: ", j))
    for(i in 1:para$K){
      phi_z <- para$phi_z[[j]][,i]
      z_zeta <- para$z_zeta_sum_by_age[,i]
      age_bases <- para$age_basis_discrte
      fail_init <- (sum(!is.finite(c(fun_age_beta_fast(para$beta[,j,i],phi_z,z_zeta,age_bases), grad_age_beta_fast(para$beta[,j,i],phi_z,z_zeta,age_bases)) ) ) != 0)
      if(fail_init){ # handle the initialization case when gradient are numerically infeasible
        print("Warning: Handle infeasible initialization, should only be occuring at first few interations.")
        para$beta <-  array(rnorm(para$P * para$D * para$K,sd =0.1), dim = c(para$P, para$D, para$K))
        para$E_lntheta <- t(sapply(rgamma(para$M,shape = 100, rate = 100), function(x) x*(digamma(para$alpha) - digamma(sum(para$alpha))) ))
        break
      }
      para$beta[,j,i] <- optim(par = para$beta[,j,i],
                               fn = function(x) fun_age_beta_fast(x,phi_z,z_zeta,age_bases),
                               gr = function(x) grad_age_beta_fast(x,phi_z,z_zeta,age_bases), method ="BFGS",
                               control = list(fnscale = -1) )$par
    }
    if(fail_init){break} # break from the iner circle as well 
  }
  
  # compute the beta_w: exponential divided by sum
  # first compute the whole age_beta_basis, for each topic it is T-by-D
  para$exp_age_beta_basis <- array( sapply(1:para$K, 
                                           function(i) exp(para$age_basis %*% para$beta[,,i])),
                                    dim=c(dim(para$age_basis)[1], para$D, para$K))
  para$sum_exp_age_beta_basis <- apply(para$exp_age_beta_basis, c(1,3), sum)
  # this is the basis of softmax function
  para$pi_beta_basis <- apply(para$exp_age_beta_basis, c(1,3), function(x) x/sum(x)) %>% 
    aperm(perm = c(2,1,3))
  
  # update beta_w: list of Ns-by-K 
  para$beta_w_full <- apply(para$pi_beta_basis, 3, function(x) 
    x[as.matrix(select(para$unlist_Ds_id, age_diag, Ds_id))]) 
  para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])
  
  return(para)
}
##########################################
# Prediction functions: should be an abstraction for all three models
# for tree and baseline, we basically assume topics are constant over age
##########################################
# use the function below to estimate with half of the records and predict the other half
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
        slice(predict_num:n()) %>%
        ungroup
    }else{
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        slice(predict_num) %>%
        ungroup
    }
    
    estimating_test_set <- testing_data_included %>%
      group_by(eid) %>%
      slice(1:(predict_num-1)) %>%
      ungroup
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
# using function below to perform onebyone prediction
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
        slice(predict_num:n()) %>%
        ungroup
    }else{
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        slice(predict_num) %>%
        ungroup
    }
    
    estimating_test_set <- testing_data_included %>%
      group_by(eid) %>%
      slice(1:(predict_num-1)) %>%
      ungroup
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
        slice(predict_num:n()) %>%
        ungroup
    }else{
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        slice(predict_num) %>%
        ungroup
    }
    
    estimating_test_set <- testing_data_included %>%
      group_by(eid) %>%
      slice(1:(predict_num-1)) %>%
      ungroup
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
      slice(1) 
    control_eid <- testing_data %>% 
      anti_join(cases_eid, by = "eid") %>%
      group_by(eid) %>%
      slice(1) %>%
      select(eid)
    control_eid$target_age <- sample(targe_age_distribution$target_age, dim(control_eid)[1], replace = T)
    
    control_data <- control_eid %>% # select all disease that has a diagnosis that are earlier than target!
      left_join(testing_data, by = "eid") %>% 
      filter(age_diag < target_age)
    
    # save the case and control outcomes; save the age of the last event before target for computing area
    cases_eid <- cases_data %>%
      group_by(eid) %>%
      slice_tail(n = 1) %>%
      mutate(outcome = 1) %>%
      select(eid, outcome, age_diag)
    control_eid <- control_data %>%
      group_by(eid) %>%
      slice_tail(n = 1) %>%
      mutate(outcome = 0) %>%
      select(eid, outcome, age_diag)
    outcomes <- bind_rows(cases_eid, control_eid) %>%
      left_join(survive_age, by = "eid") %>%
      mutate(age_diag = round(age_diag), survive_year = round(survive_year))
    # we are underestimate the score for the cases, as some of them died of the disease
    outcomes %>% mutate(age_gap = survive_year - age_diag) %>% filter(outcome == 1) %>% ungroup %>% summarise(mean(age_gap))
    outcomes %>% mutate(age_gap = survive_year - age_diag) %>% filter(outcome == 0) %>% ungroup %>% summarise(mean(age_gap))
    # below is how we estimate the area under curve as well as the 
    estimating_test_set <- bind_rows(cases_data, control_data) %>%
      select(- target_age) %>%
      ungroup()
    
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
    AUC_methods[j,1] <- auc(outcomes$outcome,outcomes$risk_n_sum ) 
    AUC_methods[j,2] <- auc(outcomes$outcome,outcomes$risk_n_mean ) 
    AUC_methods[j,3] <- auc(outcomes$outcome,outcomes$risk_3year_mean ) 
    AUC_methods[j,4] <- auc(outcomes$outcome,outcomes$risk_point ) 
    
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
      slice(1)
    
    training_control_eid <- training_data %>% 
      anti_join(cases_training_eid, by = "eid") %>%
      group_by(eid) %>%
      slice(1) %>%
      select(eid)
    
    training_control_eid$target_age <- sample(training_targe_age_distribution$target_age, dim(training_control_eid)[1], replace = T)
    training_control_data <- training_control_eid %>% # select all disease that has a diagnosis that are earlier than target!
      left_join(training_data, by = "eid") %>% 
      filter(age_diag < target_age)
    
    training_set <- bind_rows(training_cases_data, training_control_data)
    training_cases_eid <- training_cases_data %>%
      group_by(eid) %>%
      slice_tail(n = 1) %>%
      mutate(outcome = 1) %>%
      select(eid, outcome, age_diag)
    training_control_eid <- training_control_data %>%
      group_by(eid) %>%
      slice_tail(n = 1) %>%
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
      slice(1) %>%
      select(eid)
    
    testing_cases_data <- testing_cases_eid %>% # select all disease that has a diagnosis that are earlier than target!
      left_join(testing_data, by = "eid") %>% 
      filter(age_diag < target_age)
    
    targe_age_distribution <- testing_cases_data %>%
      group_by(eid) %>%
      slice(1) 
    
    testing_control_eid$target_age <- sample(targe_age_distribution$target_age, dim(testing_control_eid)[1], replace = T)
    
    testing_control_data <- testing_control_eid %>% # select all disease that has a diagnosis that are earlier than target!
      left_join(testing_data, by = "eid") %>% 
      filter(age_diag < target_age)
    
    testing_set <- bind_rows(testing_cases_data, testing_control_data)
    # save the case and control outcomes; save the age of the last event before target for computing area
    testing_cases_eid <- testing_cases_data %>%
      group_by(eid) %>%
      slice_tail(n = 1) %>%
      mutate(outcome = 1) %>%
      select(eid, outcome, age_diag)
    testing_control_eid <- testing_control_data %>%
      group_by(eid) %>%
      slice_tail(n = 1) %>%
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
      spread(key = diag_icd10, value = age_diag, fill = 0) 
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
      scores[,itr] <- predict(cvfit_testing, newx = predict_x, s = "lambda.min")
    }
    AUC_per_ds[j,1] <- auc(y_test, rowSums(scores)) 
    coefficients[[j]] <- do.call(rbind, coefs)
  }
  return(list(AUC_per_ds, coefficients))
}



##############################################
# testing for subtypes functions: Mixture of dirichlet distribution
##############################################
comp_MixDir_pi <- function(parMixDir) {
  # update pi
  parMixDir$pi <- colSums(parMixDir$z_n)/sum(parMixDir$z_n)
  return(parMixDir)
}

comp_MixDir_zn <- function(parMixDir) {
  # update z_n
  z_n <- sapply(1:parMixDir$K, function(i) parMixDir$pi[i]*ddirichlet(parMixDir$X, parMixDir$alpha[i,]))
  parMixDir$z_n <- z_n/rowSums(z_n)
  return(parMixDir)
}

logL_full_MixDir <- function(alphal, parMixDir, l) {
  # term1 is the multivariate beta distribution multiplied by sum of gammas
  term1 <- (sum(lgamma(alphal)) - lgamma(sum(alphal))) * parMixDir$GammaSumN[l]
  term2 <- (alphal - 1) %*% parMixDir$GammaLogZn[l, ]
  # print(c(term2 , term1, lgamma(sum(alphal)), alphal) ) 
  return(term2 - term1)
}
grad_logL_mixDir <- function(alphal, parMixDir, l) {
  # term1 is the multivariate beta distribution multiplied by sum of gammas
  term1 <- (digamma(sum(alphal)) - digamma(alphal)) * parMixDir$GammaSumN[l]
  term2 <- parMixDir$GammaLogZn[l, ]
  # print(c(term2 , term1, lgamma(sum(alphal)), alphal) ) 
  return(term1 + term2)
}
########################
# maybe we will need gradient information
########################

comp_MixDir_alpha <- function(parMixDir){
  # update each alpha in turn, using full likelihood
  parMixDir$GammaSumN <- colSums(parMixDir$z_n)
  parMixDir$GammaLogZn <- crossprod(parMixDir$z_n, log(parMixDir$X))
  for(l in 1:parMixDir$K){
    parMixDir$alpha[l,] <- optim(par = parMixDir$alpha[l,],
                                 fn = function(x) logL_full_MixDir(x, parMixDir, l),
                                 gr = function(x) grad_logL_mixDir(x, parMixDir, l),
                                 control = list(fnscale = -1), 
                                 lower = rep(1, parMixDir$M),  upper = rep(10, parMixDir$M),
                                 method="L-BFGS-B")$par
  }
  return(parMixDir)
}

# compute likelihood for convergence
logL_MixDir <-  function(parMixDir) {
  terms <- 0
  for(l in 1:parMixDir$K){
    terms <- terms + parMixDir$pi[l] * ddirichlet(parMixDir$X,parMixDir$alpha[l,])
  }
  return(sum(log(terms)))
}

fit_MixDir <- function(loadings, K){
  parMixDir <- list()
  parMixDir$K <- K
  parMixDir$N <- dim(loadings)[1]
  parMixDir$M <- dim(loadings)[2]
  parMixDir$pi <- runif(parMixDir$K)
  parMixDir$pi <- parMixDir$pi/sum(parMixDir$pi)
  parMixDir$X <- loadings
  parMixDir$alpha <- matrix(runif(parMixDir$K * parMixDir$M), nrow = parMixDir$K, ncol = parMixDir$M)
  parMixDir$max_itr <- 200
  parMixDir$tol <- 10^(-6)
  parMixDir$lb <- data.frame("Iteration" = 0,"Lower_bound" = logL_MixDir(parMixDir))
  for(itr in 1:parMixDir$max_itr){
    # print(paste0("Interation: ",itr))
    parMixDir <- comp_MixDir_zn(parMixDir)
    # print(paste0("after zn", logL_MixDir(parMixDir)))
    parMixDir <- comp_MixDir_pi(parMixDir)
    # print(paste0("after pi", logL_MixDir(parMixDir)))
    parMixDir <- comp_MixDir_alpha(parMixDir)
    # print(paste0("after alpha", logL_MixDir(parMixDir)))
    # para <- update_alpha(para) # we use non-informative alpha
    parMixDir$lb[nrow(parMixDir$lb) + 1,] <- c(itr, logL_MixDir(parMixDir))
    if(itr %% 5 ==0){
      curr_lb <- pull(filter(parMixDir$lb, Iteration == itr), Lower_bound) 
      prev_lb <- pull(filter(parMixDir$lb, Iteration == (itr -para$itr_check_lb )), Lower_bound) 
      # print(paste0("Current Lower bound ", curr_lb, " at iteration: ",itr))
      try({
        if(is.finite((curr_lb - prev_lb)) & (curr_lb - prev_lb)/abs(prev_lb) < para$tol ){
          print(paste0("Optimization converged at step ", itr))
          break
        }
      })
    }
  }
  return(parMixDir)
}

#########################################################################################
# topics with Tree structures
# for point estimate of beta, we don't need to change any other updates including the lower bound
####################################################################################################
# The key is created a similar structure of the beta_w; each word is mapped to multiple layer.  
create_structure <- function(para){
  # there are three metrics to compute: 1. NodeDS index matrix of (L)-by-D: elements refereing parent for disease j at layer l
  # 2. a disease list: treeDS[[l]][[cl]] the records associated with each node, tree equivalent of para$ds_list (using the list of diseases code from 1:para$D that are under that node)
  # 3. the beta list: treeBeta[[l]] is the vector (matrix to include i) of betas at this layer 
  # 4. para$Cl[[l]] save the number of nodes at each layer
  # using Z for filling the empty
  filled_list <- as.character(para$list_above500occu$diag_icd10)
  # make sure death is separate from the tree
  if("Death" %in% filled_list){
    filled_list[which(filled_list == "Death")] <- "ZZZZ"}
  # healthy is also seperate from other branches
  if("Healthy" %in% filled_list){
    filled_list[which(filled_list == "Healthy")] <- "YYYY"}
  
  para$NodeDS <- matrix(nrow = para$L, ncol = para$D)
  para$treeDS <- list()
  para$treeBeta <- list()
  para$Cl <- list()
  for(l in 1:para$L){
    if(l == para$L){
      nodes_l <- unlist(lapply(filled_list, function(x) (substring(x, 1,l+4))))
    }else{
      nodes_l <- unlist(lapply(filled_list, function(x) (substring(x, 1,l)))) # longest ICD10 is 7 character
    }
    node_set <- unique(nodes_l)
    para$Cl[[l]] <- length(node_set) # save the number of children at each layer
    para$NodeDS[l,] <- match(nodes_l, node_set)
    para$treeDS[[l]] <- list()
    for(nd in 1:length(node_set)){
      para$treeDS[[l]][[nd]] <- bind_rows(para$ds_list[which(node_set[nd] == nodes_l)])
    }
  }
  # initialize treeBeta
  for(l in 1:para$L){
    para$treeBeta[[l]] <- matrix(nrow = para$Cl[[l]], ncol = para$K)
  }
  return(para)
}

update_beta_treeLDA <- function(para){
    # compute the beta
    for(l in 1:para$L){
      # compute M-step for beta: basic case, direct maximize the upper bound
      para$treeBeta[[l]] <- sapply(1:para$Cl[[l]], function(c) colSums(para$unlist_zn[para$treeDS[[l]][[c]]$id,]) ) %>% t 
      # the normaliztion is the key: we need to normalize with respect to common parent
      if(l == 1){
        para$treeBeta[[l]] <- sapply(1:para$K, function(i) para$treeBeta[[l]][,i]/sum(para$treeBeta[[l]][,i]))
      }else{
        for(c in 1:para$Cl[[l-1]]){
          parent_idx <- which(para$NodeDS[l-1,] == c) # the indices at the parent layer
          childeren_idx <-  unique(para$NodeDS[l, parent_idx])# the child under this c 
          para$treeBeta[[l]][childeren_idx, ] <- sapply(1:para$K, function(i) para$treeBeta[[l]][childeren_idx,i]/
                                                                 sum(para$treeBeta[[l]][childeren_idx,i]))
        }
      }
    }
    
    # compute the beta at leaves: using para$NodeDS going to go down the tree
    for(j in 1:para$D){
      beta_along_tree <- sapply(1:para$L, function(l) para$treeBeta[[l]][para$NodeDS[l,j],,drop = F])
      para$beta[j,] <- sapply(1:para$K, function(i) prod(beta_along_tree[i,]))
    }
  
    # this beta_w parameter save the beta for each word: it is a list of M elements and each contain a K*Ns matrix
    # para$beta_w <- lapply(para$w, function(w) para$beta[w$Ds_id,,drop=FALSE] )
    para$beta_w_full <- para$beta[para$unlist_Ds_id$Ds_id,,drop=FALSE]
    para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])
    return(para)
}

