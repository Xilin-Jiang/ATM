test_that("test lfa inference", {
  para <- topic_init_lfa_cvb(rec_data = HES_age_example, ds_list=UKB_349_disease, topic_num=10)
  # update the zn; each step should increase the lower bound
  lb1 <- CVB_lb_lfa(para)
  para <- CVB0_E_zn(para)
  lb2 <- CVB_lb_lfa(para)
  expect_gt(lb2, lb1)
  # update the beta
  para <- update_beta_lfa(para)
  lb3 <- CVB_lb_lfa(para)
  expect_gt(lb3, lb2)
})

test_that("lfa accuracy", {
  # simulate some data using Yidong's code
  # Parameters:
  K <- 4    # number of topics
  S <- 15     # total number of diseases (internal and terminal)
  D <- 5000   # number of individuals in the training data
  alpha <- rep(0.1,4)    # Dirichlet prior for topic weights


  # manually set the topics
  topics <- c( c(0.9, 0.9, 0.9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(0, 0, 0, 0.9, 0.9, 0.9, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(0, 0, 0, 0, 0, 0, 0.9, 0.9, 0.9, 0, 0, 0, 0, 0, 0),
               c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9) )

  topics <- matrix(topics,nrow=K,byrow=TRUE)/3

  theta <- gtools::rdirichlet(D, alpha)

  zn_sim <- matrix(NA, nrow =  D, ncol = S)
  w_sim <- matrix(NA, nrow =  D, ncol = S)
  for(s in 1:D){
    zn_sim[s, ] <- sample( 1:K, S, replace=TRUE, prob=theta[s,])
    w_sim[s, ] <- sapply(1:S, function(n) sample( c(0,1), size= 1, replace=TRUE,
                                                           prob=c(1-topics[zn_sim[s,n],n],topics[zn_sim[s,n],n]) ))

  }

  data <- as.matrix(w_sim)
  colnames(data) <- paste("D",1:ncol(data),sep="")

  # change the format:
  rownames(data) <- 1:nrow(data)
  data <- cbind(rownames(data),data)
  colnames(data)[1] <- "patient"

  data <- data.frame(data)

  data_long <- data %>%
    pivot_longer(names_to ="disease", values_to ="count", cols = !patient) %>%
    filter(count == 1) %>%
    rename(eid = patient, diag_icd10 = disease) %>%
    select(eid, diag_icd10)
  data_long_eid <- data_long %>%
    group_by(eid) %>%
    dplyr::tally() %>%
    filter(n>1) %>%
    dplyr::pull(eid)
  data_long <- data_long %>% filter(eid %in% data_long_eid)
  # Run lfa-atm: using below to verify the model
  # simu_result <- wrapper_LFA( data_long, topic_num = 4, CVB_num=5)
  # order_disease <- sapply(paste("D",1:S,sep=""), function(x )which(x == simu_result$ds_list$diag_icd10))
  # topic_ordered <-  simu_result$topic_loadings[order_disease, ]
  # # plot inferred topics:
  # plot_lfa_topics(paste0("D",1:S,sep=""), beta = topic_ordered,  plot_title = "LFA topics")
})
