## file order_topic_output.R
# order the topics by posterior 
args <- commandArgs(trailingOnly = TRUE) # args[1] defines the number of topics; args[2] defines degree of freedom for the curves; args[3] is the repitation id
K <- as.numeric(args[1]) 
P <- as.numeric(args[2]) # degrees of freedom
rep_id <- args[3]

load(paste0("Results/Run_A2N_age_dependent_K",K,"_P",P,"_rep",rep_id, ".RData"))
order_by_post_z <- order(colMeans(para$alpha_z), decreasing = T)
ordered_pi_beta <- para$pi_beta_basis[,,order_by_post_z]
model_output <- list(ordered_pi_beta, para$lb, para$D, para$M, para$K, para$P)
save(model_output, file = paste0("Results/","model_output_A2N_age_dependent_K",K,"_P",P,"_rep",rep_id, ".RData"))
