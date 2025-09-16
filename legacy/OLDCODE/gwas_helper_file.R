# using this file to help extract clumped SNP for GCTA
library(dplyr)
library(stringr)
# DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/Association_analysis/association_results/"
# combine clump chr data
# K <- 10
# for(k in 1:K){
#   for(rep_id in 1:10){
#     print(paste0("topic number: ",k, " repeat ID: ", rep_id))
#     pt <- paste0("^clump_topic", k,"_ageLDA_K", K, "_rep",rep_id, "_chr.*clumped")
#     temp <- list.files(paste(DIR, sep=""), pattern=pt)
#     clumped <- list()
#     for(chr in 1:length(temp)){
#       try({
#         clumped[[chr]] <- read.table(paste0(DIR,temp[chr]), header = T)
#       })
#     }
#     clumped <- bind_rows(clumped)
#     write.table(select(clumped, SNP), paste0(DIR, "SNP_topic", k,"_ageLDA_K", K, "_rep",rep_id, ".txt"), sep="\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
#   }
# }

# transform the results to log scale 
# ds_list <- read.table("/users/mcvean/xilin/xilin/Multimorbidity-biobank/top50_hererogeneous_disease.txt")
# DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/Association_analysis/"
# args <- commandArgs(trailingOnly = TRUE) 
# ds <- args[1]
# gwas_rslt <- read.table(paste0(DIR, ds, "/",ds,".assoc"), header=T)
# gwas_rslt$BETA <- log(gwas_rslt$OR)
# write.table(gwas_rslt, paste0(DIR, ds, "/",ds,".logOR"), sep = "\t", quote = FALSE,  row.names = FALSE)

# # separate the cases to two groups 
# DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/Association_analysis/"
# ds50 <- read.table("/users/mcvean/xilin/xilin/Multimorbidity-biobank/top50_hererogeneous_disease.txt")
# load("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep2.RData")
# for(ds_id in ds50[[1]]){
#   j <- match(ds_id, para$list_above500occu$diag_icd10)
#   loadings <- para$unlist_zn[para$ds_list[[j]]$id,]
#   eid <- para$unlist_Ds_id[para$ds_list[[j]]$id,]$eid
#   kfit <- kmeans(loadings, centers=2)
#   typeds <- data.frame(kmCluster = kfit$cluster, eid)
#   
#   data_case_control <- read.table(paste(DIR, ds_id,"/",ds_id,"_pheno_over_age.txt",sep=""))
#   data_two_types <- data_case_control %>%
#     left_join(typeds, by = c("V1" = "eid")) %>%
#     mutate(V3 = if_else(kmCluster == 1, 1,-9),V4 = if_else(kmCluster == 2, 1,-9)) %>%
#     mutate(V3 = if_else(is.na(V3), 0,V3), V4 = if_else(is.na(V4), 0,V4)) %>%
#     select(V1, V2, V3, V4) %>%
#     rename(FID = V1, IID = V2)
#   write.table(data_two_types, paste(DIR, ds_id,"/",ds_id,"subtypes.txt",sep=""), sep="\t", col.names = TRUE, row.names = FALSE ,quote = F)
#   write.table(select(data_two_types, FID, IID), paste(DIR, ds_id,"/",ds_id,"keep.txt",sep=""), sep="\t", col.names = F, row.names = FALSE)
# }

##########################################
# extract BOLT-LMM/BOLT-REML results
##########################################
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/Association_analysis/"
ds50 <- read.table("/users/mcvean/xilin/xilin/Multimorbidity-biobank/top50_hererogeneous_disease.txt")
load("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep2.RData")

# everything we want to save
sz_pheno1 <- rep(NA, length(ds50[[1]])) # save also sample size
sz_pheno2 <- rep(NA, length(ds50[[1]]))
age_pheno1 <- rep(NA, length(ds50[[1]]))
age_pheno2 <- rep(NA, length(ds50[[1]]))
phe1h2g <- rep(NA, length(ds50[[1]]))
phe2h2g <- rep(NA, length(ds50[[1]]))
corr_lmm <- rep(NA, length(ds50[[1]]))
for(i in 1:length(ds50[[1]])){
  ds_id <-  ds50[[1]][i]
  j <- match(ds_id, para$list_above500occu$diag_icd10)
  print(ds_id)
  try({
    age_data <- para$unlist_Ds_id[para$ds_list[[j]]$id,]
    subtypes <- read.table(paste(DIR, ds_id,"/",ds_id,"subtypes.txt",sep=""), header = T)
    sz_pheno1[i] <- subtypes %>%
      filter(V3 == 1) %>%
      tally() %>% pull(1)
    sz_pheno2[i] <- subtypes %>%
      filter(V4 == 1) %>%
      tally() %>% pull(1)
      
    
    age_pheno1[i] <- subtypes %>%
      filter(V3 == 1) %>%
      left_join(age_data, by = c("FID" = "eid")) %>%
      pull(age_diag) %>%
      mean()
    age_pheno2[i]  <- subtypes %>%
      filter(V4 == 1) %>%
      left_join(age_data, by = c("FID" = "eid")) %>%
      pull(age_diag) %>%
      mean()
    
    # extract heritability
    h2g <- readLines(paste(DIR, ds_id,"/",ds_id,"pheno1.heritability",sep=""))
    phe1h2g[i] <- as.numeric(str_split(h2g[length(h2g) -2], " ")[[1]][5])
    
    h2g <- readLines(paste(DIR, ds_id,"/",ds_id,"pheno2.heritability",sep=""))
    phe2h2g[i] <- as.numeric(str_split(h2g[length(h2g) -2], " ")[[1]][5])
    
    # compute the correlation of effect sizes
    lmm_pheno1 <- read.table(paste(DIR, ds_id,"/",ds_id,"lmm.pheno1.stats",sep=""), header = T) %>%
      mutate(A11 = ALLELE1, A01 = ALLELE0,BETA1= BETA) %>%
      select(SNP, A11, A01, BETA1)
    lmm_pheno2 <- read.table(paste(DIR, ds_id,"/",ds_id,"lmm.pheno2.stats",sep=""), header = T) %>%
      mutate(A12 = ALLELE1, A02 = ALLELE0,BETA2 = BETA) %>%
      select(SNP, A12, A02, BETA2)
    joint_lmm <- lmm_pheno1 %>% 
      left_join(lmm_pheno2, by = "SNP") %>%
      mutate(BETA2 = if_else(A11 == A12, BETA2, -BETA2)) # turn the risk/protective effect
    corr_lmm[i] <- cor(joint_lmm$BETA1, joint_lmm$BETA2)
  })
}
df_lmm <- data.frame(ds50[[1]],sz_pheno1, sz_pheno2, age_pheno1,age_pheno2,phe1h2g, phe2h2g, corr_lmm)
write.table(df_lmm ,paste(DIR, "subtypes_boltlmm.txt",sep=""), sep="\t", col.names = TRUE, row.names = FALSE ,quote = F)


