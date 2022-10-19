library(dplyr)
library(stringr)
##########################################
# extract LDSC results
##########################################
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/Association_analysis/"
ds_target <- read.table("/users/mcvean/xilin/xilin/Multimorbidity-biobank/top_hererogeneous_disease.txt")
load("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")

# everything we want to save
Fst <- matrix(NA, nrow = length(ds_target[[1]]), ncol = 1)
p_fst <- matrix(NA, nrow = length(ds_target[[1]]), ncol = 1)
for(i in 1:length(ds_target[[1]])){
  ds_id <-  ds_target[[1]][i]
  j <- match(ds_id, para$list_above500occu$diag_icd10)
  print(ds_id)
  try({
    # also save the Fst
    fst <- readLines(paste(DIR, ds_id,"/",ds_id,"subtp_Fst.log",sep=""))
    target_fst <- as.numeric(str_split(fst[length(fst) - 2], "\\s+")[[1]][4])
    Fst[i] <- target_fst
    permutation_Fst <- read.csv(paste(DIR, ds_id,"/",ds_id,"_permutation_Fst.csv",sep=""))
    p_fst[i]  <- (sum(target_fst < permutation_Fst$weighted_Fst) + 1)/dim(permutation_Fst)[1]
  })
}
df_rg <- data.frame(disease = ds_target[[1]],Fst = Fst, p_fst = p_fst)
write.csv(df_rg ,paste(DIR, "subtypes_Fst.csv",sep=""), row.names = FALSE )

#####################################
# extract Fst for topic matched controls
######################################
Fst <- matrix(NA, nrow = length(ds_target[[1]]), ncol = 1)
p_fst <- matrix(NA, nrow = length(ds_target[[1]]), ncol = 1)
for(i in 1:length(ds_target[[1]])){
  ds_id <-  ds_target[[1]][i]
  j <- match(ds_id, para$list_above500occu$diag_icd10)
  print(ds_id)
  try({
    # also save the Fst
    fst <- readLines(paste(DIR, ds_id,"/",ds_id,"subtp_Fst.log",sep=""))
    target_fst <- as.numeric(str_split(fst[length(fst) - 2], "\\s+")[[1]][4])
    Fst[i] <- target_fst
    permutation_Fst <- read.csv(paste(DIR, ds_id,"/",ds_id,"_control_Fst.csv",sep=""))
    p_fst[i]  <- (sum(target_fst < permutation_Fst$weighted_Fst) + 1)/dim(permutation_Fst)[1]
  })
}
df_rg <- data.frame(disease = ds_target[[1]],Fst = Fst, p_fst = p_fst)
write.csv(df_rg ,paste(DIR, "subtypes_Fst_matched_topic.csv",sep=""), row.names = FALSE )


