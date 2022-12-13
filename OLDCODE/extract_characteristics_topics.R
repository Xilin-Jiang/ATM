# in this file, I include code needed for extracting topic information from the results on the cluster. 
# there will be multiple tasks including: 1. extracting age & sex (and other individual information) distribution 
# of the topics; 2. extract topic loadings for GWAS. 3. for the prediction data: compute the single disease
# risk factors that could be evaluated. 

# extract age distribution
source("topic_functions.R")
args <- commandArgs(trailingOnly = TRUE) # args[1] defines the number of topics; args[2] defines degree of freedom for the curves; args[3] is the repitation id

topic_num <- as.numeric(args[1]) 
degree_free_num <- as.numeric(args[2]) # degrees of freedom
rep_ID <- args[3]

load(paste0("Results/","Run_2rec_PheCode_age_dependent_K",topic_num,"_P",degree_free_num,"_rep",rep_ID, ".RData"))
loadings <- sweep(para$alpha_z, 1, rowSums(para$alpha_z), FUN="/")
save_loadings <- data.frame(para$eid, para$eid, loadings) %>%
  rename(FID = para.eid, IID = para.eid.1)
LOC <- "Association_analysis/Phenotypes/"
write.table(save_loadings, paste0(LOC,"ageLDA_topic_loadings",para$K,"_P",para$P,"_rep",para$rep_ID, ".txt"), sep="\t", col.names = TRUE, row.names = FALSE, quote = F)

# second task: extract the loading specific age for each diseases
# first extract the loading distribution
assignment_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[para$ds_list[[j]]$id,]) ) %>% t
# then compute the average age for each loading
age_topic_assoc_per_ds <- matrix(NA, nrow = para$D, ncol = para$K)
for(j in 1:para$D){
  loadings <- para$unlist_zn[para$ds_list[[j]]$id,]
  sum_loadings <- colSums(loadings)
  loadings <- sweep(loadings, 2, sum_loadings, FUN = "/")
  age <- para$unlist_Ds_id[para$ds_list[[j]]$id,]$age_diag
  age_topic_assoc_per_ds[j,] <- sapply(1:para$K, function(i) age %*% loadings[,i])
}

topic_characteristic <- list(assignment_per_ds, age_topic_assoc_per_ds)
save(topic_characteristic, file = paste0("Results/","topic_characteristic_K",para$K,"_P",para$P,"_rep",para$rep_ID, ".RData"))

