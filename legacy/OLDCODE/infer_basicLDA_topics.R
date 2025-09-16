source("topic_functions.R")
args <- commandArgs(trailingOnly = TRUE) # args[1] defines the number of topics; args[2] defines degree of freedom for the curves; args[3] is the repitation id

rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
topic_num <- as.numeric(args[1]) 
para <- topic_init_baseline(rec_data, ds_list, topic_num)
para$rep_ID <- args[3]

############################
# start optimization
############################
para$max_itr <- 2000
para$alpha <- rep(1, para$K) 
para$lb <- data.frame("Iteration" = 0,"Lower_bound" = CVB_lb(para))
para$itr_check_lb <- 5
para$itr_save <- para$max_itr # don't need to save intermediate data
para$tol <- 10^(-6)
for(itr in 1:para$max_itr){
  print(paste0("Interation: ",itr))
  para <- CVB0_E_zn(para) # we choose CVB0 as papers shown it could converge quicker
  para <- update_beta_basic_lda(para)
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

# save the parameter
save(para, file = paste0("Results/","Run_BaselinLDA_Phecode_K",para$K,"_rep",para$rep_ID, ".RData"))

# save a smaller dataset
order_by_post_z <- order(colMeans(para$alpha_z), decreasing = T)
ordered_pi_beta <- para$pi_beta_basis[,,order_by_post_z]
# save the ordered version of model_output
model_output <- list(ordered_pi_beta, para$lb, para$D, para$M, para$K, para$P)
print(object.size(model_output), unit = "MB" , standard = "SI")
  save(model_output, file = paste0("Results/","BaselinLDA_model_output_PheCode_K",para$K,"_P",para$P,"_rep",para$rep_ID, ".RData"))
