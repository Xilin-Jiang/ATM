# save all the inference related functions here
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

##########################################################
# functions for inferring topic weights from known topic loadings
##########################################################
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
# ignore the second order information in the approximation of equation 7 in supplementary note
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
#########################
# topic loading inference
#########################
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

### the slower methods, not useful for large dataset but still intuitive for analytical purpose
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

wrapper_ATM <- function(rec_data, topic_num = 10, degree_free_num = 3, CVB_num = 5){
  dir.create("Results")
  ds_list <- rec_data %>%
    group_by(diag_icd10) %>%
    summarise(occ = n())
  topics <- list()
  ELBOs <- list()
  for(cvb_rep in 1:CVB_num){
    print(paste0("Male CVB inference number: ", cvb_rep))
    para <- topic_init_age(rec_data, ds_list, topic_num, degree_free_num)
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
    }
    topics[[cvb_rep]] <- para$pi_beta_basis
    ELBOs[[cvb_rep]]  <- para$lb
  }
  multi_runs <- list(topics, ELBOs)
  save(multi_runs, file = paste0("Results/","multirun", CVB_num, "K",para$K,"_P",para$P, ".RData"))
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
  model_output <- list(topics[[best_id]], ELBOs[[best_id]], para$D, para$M, para$K, para$P, lb_rep)
  save(model_output, file = paste0("Results/","best_output_AgeLDA_RunNumber", CVB_num, "K",para$K,"_P",para$P,"_rep",para$rep_ID, ".RData"))
}


