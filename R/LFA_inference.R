# save all the inference related functions here
topic_init_lfa_cvb <- function(rec_data, ds_list, topic_num){
  # arrange the data stack by individual is very important as we need to make sure the matrix could be rejoined into a single matrix
  first_incidence_age <- rec_data %>%
    arrange(eid) %>%
    select(eid, diag_icd10) %>%
    group_by(eid, diag_icd10) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

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

  para$unlist_Ds_id <- first_incidence_age %>%
    mutate(Ds_id = code2id(diag_icd10)) %>%
    select(-diag_icd10)

  # create a all-disease list
  para$unlist_Ds_id <- para$unlist_Ds_id %>%
    tidyr::expand(eid, Ds_id) %>%
    left_join(mutate(para$unlist_Ds_id, disease = 1), by = c("eid", "Ds_id")) %>%
    mutate(disease = tidyr::replace_na(disease, 0))

  # the patient_list provide the column index for efficiently breaking down matrix into list of matrices
  para$patient_lst <- para$unlist_Ds_id %>%
    mutate(id = dplyr::row_number()) %>%
    select(eid, id) %>%
    group_by(eid) %>%
    dplyr::group_split(.keep = F) %>%
    lapply(pull)

  # this list is splitted by disease
  para$ds_list <- para$unlist_Ds_id %>%
    select(-eid) %>%
    mutate(id = dplyr::row_number()) %>%
    group_by(Ds_id) %>%
    dplyr::group_split()

  # initiate beta
  para$eta <- rgamma(para$D,shape = 100, rate = 100)
  # each column is a topic; D*K matrix
  para$beta <- t(gtools::rdirichlet(para$K, para$eta))

  # compute beta using Bernoulli distribution: it is a list of M elements and each contain a K*D matrix
  para$beta_w_full <- para$beta[para$unlist_Ds_id$Ds_id,,drop=FALSE] * para$unlist_Ds_id$disease + (1-para$beta[para$unlist_Ds_id$Ds_id,,drop=FALSE]) * (1-para$unlist_Ds_id$disease)
  para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])


  # initiate alpha
  para$alpha <- rgamma(para$K, shape = 50, rate = 10)
  # Matrix of M*K
  para$E_lntheta <- t(sapply(rgamma(para$M,shape = 100, rate = 100), function(x) x*(digamma(para$alpha) - digamma(sum(para$alpha))) ))

  # update E_zn: list of M; each element is matrix of D*K
  para <- comp_E_zn(para)

  # eventually reset alpha to be non-informative
  para$alpha <- rep(1, para$K)

  print("LFA RAM occupation: ")
  print(object.size(para), unit = "GB" , standard = "SI")

  return(para)
}


# compute lower bound on the full likelihood
CVB_lb_lfa <- function(para){
  # compute the lower bound for the whole dataset

  # terms for cross-entropy
  term1 <- para$M*(lgamma(sum(para$alpha)) - sum(lgamma(para$alpha)) ) -
    sum( lgamma(sum(para$alpha) + para$D) )

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


update_beta_lfa <- function(para){
  # create a unlist_zn that are non-zero at only the diseases position to facilitate computation
  unlist_zn_disease <- para$unlist_zn * para$unlist_Ds_id$disease
  # compute M-step for beta: basic case, direct maximize the upper bound
  beta_nemerator <- sapply(1:para$D, function(j) colSums(unlist_zn_disease[para$ds_list[[j]]$id,]) ) %>% t
  para$beta <- sapply(1:para$D, function(j) colSums(para$unlist_zn[para$ds_list[[j]]$id,]) ) %>% t
  # normalize beta
  para$beta <- beta_nemerator/para$beta
  # this beta_w parameter save the beta for each word: it is a list of M elements and each contain a K*Ns matrix
  para$beta_w_full <- para$beta[para$unlist_Ds_id$Ds_id,,drop=FALSE] * para$unlist_Ds_id$disease + (1-para$beta[para$unlist_Ds_id$Ds_id,,drop=FALSE]) * (1-para$unlist_Ds_id$disease)

  para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])
  return(para)
}

# put a beta prior (shared for each bernouli probability)
update_beta_lfa_with_prior <- function(para){
  # create a unlist_zn that are non-zero at only the diseases position to facilitate computation
  unlist_zn_disease <- para$unlist_zn * para$unlist_Ds_id$disease
  # compute M-step for beta: basic case, direct maximize the upper bound
  beta_nemerator <- sapply(1:para$D, function(j) colSums(unlist_zn_disease[para$ds_list[[j]]$id,]) ) %>% t + (para$beta_prior_a - 1)
  para$beta <- sapply(1:para$D, function(j) colSums(para$unlist_zn[para$ds_list[[j]]$id,]) ) %>% t + (para$beta_prior_a + para$beta_prior_b - 2)
  # normalize beta
  para$beta <- beta_nemerator/para$beta
  # this beta_w parameter save the beta for each word: it is a list of M elements and each contain a K*Ns matrix
  para$beta_w_full <- para$beta[para$unlist_Ds_id$Ds_id,,drop=FALSE] * para$unlist_Ds_id$disease + (1-para$beta[para$unlist_Ds_id$Ds_id,,drop=FALSE]) * (1-para$unlist_Ds_id$disease)

  para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])
  return(para)
}


#' Run LFA on diagnosis data.
#'
#' Run LFA on diagnosis data to infer topic loadings and topic weights. Note one run of LFA on 100K individuals would take ~30min (defualt is 5 runs and pick the best fit);
#' if the data set is small and the goal is to infer patient-level topic weights (i.e. assign comorbidity profiles to individuals based on the disedases),
#' please use loading2weights.
#'
#' @param rec_data A diagnosis data frame with three columns; format data as HES_age_example; first column is individual ids, second column is the disease code;
#' third column is the age at diagnosis. Note for each individual, we only keep the first onset of each diseases. Therefore, if there are multiple incidences of the same disease
#' within each individual, the rest will be ignored.
#' @param topic_num Number of topics to infer.
#' @param CVB_num Number of runs with random initialization. The final output will be the run with highest ELBO value.
#' @param save_data A flag which determine whether full model data will be saved. If TRUE, a Results/ folder will be created and full model data will be saved. Default is set to be FALSE.
#' @param beta_prior_flag A flag if true, will use a beta prior on the topic loading. Default is set to be FALSE.
#' @param topic_weight_prior prior of individual topic weights, default is set to be a vector of one (non-informative)
#'
#' @return Return a list object with topic_loadings (of the best run), topic_weights (of the best run), ELBO_convergence (ELBO until convergence),
#' patient_list (list of eid which correspond to rows of topic_weights), ds_list (gives the ordering of diseases in the topic_loadings object), disease_number (number of total diseases), patient_number(total number of patients), topic_number (total number of topic),
#' ,multiple_run_ELBO_compare (ELBO of each runs).
#' @export
#'
#' @examples   HES_age_small_sample <- HES_age_example[1:100,]
#' inference_results <- wrapper_LFA(HES_age_small_sample, topic_num = 3, CVB_num = 1)
wrapper_LFA <- function(rec_data, topic_num, CVB_num = 5, save_data = F, beta_prior_flag = F, topic_weight_prior=NULL){
  ds_list <- rec_data %>%
    group_by(diag_icd10) %>%
    summarise(occ = n())
  topics <- list()
  ELBOs <- list()
  topic_weights <- list()
  for(cvb_rep in 1:CVB_num){
    print(paste0("CVB inference number: ", cvb_rep))
    para <- topic_init_lfa_cvb(rec_data, ds_list, topic_num)

    # set the beta priors
    para$beta_prior_a <- 1 # set it to be 1 to avoid numeric error
    para$beta_prior_b <- ( (para$D * para$M)/sum(para$unlist_Ds_id$disease) - 1) * para$beta_prior_a

    # set the number of update
    para$max_itr <- 2000
    if(is.null(topic_weight_prior)){
      para$alpha <- rep(1, para$K)
    }else{
      if(length(topic_weight_prior) != para$K){
        print("pre-specified prior of topic weight must be a postive vector of length K (number of topics)")
        para$alpha <- rep(1, para$K)
      }else{
        para$alpha <- topic_weight_prior
      }
    }
    para$lb <- data.frame("Iteration" = 0,"Lower_bound" = CVB_lb(para))
    para$itr_beta <- 1
    para$itr_check_lb <- 5
    para$itr_save <- para$max_itr + 1 # don't need to save intermediate data
    para$tol <- 10^(-7)
    for(itr in 1:para$max_itr){
      print(paste0("Interation: ",itr))
      for(itr_inside in 1:para$itr_beta){ # in practice we perform quick steps a few times before move on.
        para <- CVB_E_zn(para) # we choose CVB with second order approximation
        # para <- comp_E_lntheta(para)
      }
      if(beta_prior_flag){
        para <- update_beta_lfa_with_prior(para) # use prior for beta updates
      }else{
        para <- update_beta_lfa(para)
      }
      # para <- update_alpha(para) # we use non-informative alpha
      para$lb[nrow(para$lb) + 1,] <- c(itr, CVB_lb_lfa(para))
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
    topics[[cvb_rep]] <- para$beta
    ELBOs[[cvb_rep]]  <- para$lb
    topic_weights[[cvb_rep]] <- sweep((para$alpha_z - 1), 1, rowSums(para$alpha_z -1), FUN="/")

    if(save_data){
      print(paste0("Saving full model for CVB rep", cvb_rep))
      dir.create("Results")
      save(para, file = paste0("Results/","fullmodel_LFA_cvbrep_", cvb_rep, "K",para$K,"_P",para$P, ".RData"))
    }

  }
  multi_runs <- list(topics, ELBOs)
  # find the best reps in the data
  lb_rep <- tibble::tibble(reps = as.integer(), lower_bound = as.numeric())
  for(cvb_rep in 1:CVB_num){
    cvrg_lb <-  ELBOs[[cvb_rep]] %>%
      filter(!is.na(Lower_bound)) %>%
      dplyr::slice_tail(n = 1) %>%
      pull(2)
    lb_rep <- lb_rep %>%
      dplyr::add_row(reps = cvb_rep, lower_bound = cvrg_lb)
  }
  best_id <- order(lb_rep$lower_bound, decreasing = T)[1]

  # save a smaller dataset
  model_output <- list(topics[[best_id]], ELBOs[[best_id]], para$D, para$M, para$K, para$P, lb_rep)
  if(save_data){
    save(multi_runs, file = paste0("Results/","multirun", CVB_num, "K",para$K,"_P",para$P, ".RData"))
    save(model_output, file = paste0("Results/","best_output_LFA_RunNumber", CVB_num, "K",para$K,"_P",para$P, ".RData"))
  }
  output <- list()
  output$topic_loadings <- topics[[best_id]]
  output$topic_weights <- topic_weights[[best_id]]
  output$ELBO_convergence <- ELBOs[[best_id]]
  output$ds_list <- para$list_above500occu
  output$patient_list <- para$eid
  output$disease_number <-  para$D
  output$patient_number <- para$M
  output$topic_number <- para$K
  output$topic_configuration <-para$P
  output$multiple_run_ELBO_compare <-lb_rep
  return(output)
}

# plot topic loadings for LFA
#' Title plot topic loadings for LFA.
#'
#' @param disease_names the list of disease names, ordered as the topic.
#' @param beta disease topics, which should be a matrix of K-by-disease.
#' @param plot_title the title of the figure.
#'
#' @return a ggplot object of the topic loading.
#' @export
#'
#' @examples disease_list <- UKB_349_disease$diag_icd10[1:50]
#' topics <- matrix(rnorm(10*length(UKB_349_disease)), nrow = length(UKB_349_disease), ncol = 10)
#' plot_lfa_topics(disease_names = disease_list,
#'         beta = topics,
#'         plot_title = "Example noisy topics presentation")
plot_lfa_topics <- function(disease_names, beta,  plot_title = ""){
  longData <- reshape2::melt(beta) %>%
    mutate(Var1 =disease_names[Var1])
  longData <- longData %>%
    mutate(Var1 = factor(Var1, levels = disease_names))
  if(length(disease_names) <= 50){
    plt <- ggplot() +
      ggplot2::geom_tile(data = longData, aes(x = Var2, y = Var1, fill=red, alpha = value,width = 0.9)) +
      ggplot2::scale_alpha_continuous(range = c(0,1)) +
      ggplot2::labs(x="", y="", title="") +
      ggplot2::scale_x_discrete(expand=c(0,0)) +
      ggplot2::scale_y_discrete(expand=c(0,0)) +
      ggplot2::theme(axis.line=ggplot2::element_blank(),axis.text.x=ggplot2::element_blank(),
                     axis.ticks=ggplot2::element_blank(),
                     axis.title.x=ggplot2::element_blank(),
                     axis.title.y=ggplot2::element_blank(),legend.position="none",
                     panel.background=ggplot2::element_blank(),panel.border=ggplot2::element_blank(),panel.grid.major=ggplot2::element_blank(),
                     panel.grid.minor=ggplot2::element_blank(),plot.background=ggplot2::element_blank())
  }else{
    print("Too many diseases (>50), hide the disease name.")
    plt <- ggplot() +
      ggplot2::geom_tile(data = longData, aes(x = Var2, y = Var1, fill=red, alpha = value,width = 0.9)) +
      ggplot2::scale_alpha_continuous(range = c(0,1)) +
      ggplot2::labs(x="", y="", title="") +
      ggplot2::scale_x_discrete(expand=c(0,0)) +
      ggplot2::scale_y_discrete(expand=c(0,0)) +
      ggplot2::theme(axis.line=ggplot2::element_blank(),axis.text.x=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),
                     axis.title.x=ggplot2::element_blank(),
                     axis.title.y=ggplot2::element_blank(),legend.position="none",
                     panel.background=ggplot2::element_blank(),panel.border=ggplot2::element_blank(),panel.grid.major=ggplot2::element_blank(),
                     panel.grid.minor=ggplot2::element_blank(),plot.background=ggplot2::element_blank())
  }
  return(plt)
}




