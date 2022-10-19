library(dplyr)
library(stringr)
##########################################
# extract LDSC results
##########################################
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/Association_analysis/"
ds_target <- read.table("/users/mcvean/xilin/xilin/Multimorbidity-biobank/top_hererogeneous_disease.txt")
load("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")

# everything we want to save
age_pheno <- matrix(NA, nrow = length(ds_target[[1]]), ncol = 3)
pheh2g <- matrix(NA, nrow = length(ds_target[[1]]), ncol = 3)
seh2g <- matrix(NA, nrow = length(ds_target[[1]]), ncol = 3)
rg <- matrix(NA, nrow = length(ds_target[[1]]), ncol = 2)
serg <- matrix(NA, nrow = length(ds_target[[1]]), ncol = 2)
p_fst <- matrix(NA, nrow = length(ds_target[[1]]), ncol = 1)
for(i in 1:length(ds_target[[1]])){
  ds_id <-  ds_target[[1]][i]
  j <- match(ds_id, para$list_above500occu$diag_icd10)
  print(ds_id)
  try({
    age_data <- para$unlist_Ds_id[para$ds_list[[j]]$id,]
    subtypes <- read.table(paste(DIR, ds_id,"/",ds_id,"subtypes.txt",sep=""), header = F)
    
    age_pheno[i,1] <- subtypes %>%
      filter(V3 == 1) %>%
      left_join(age_data, by = c("V1" = "eid")) %>%
      pull(age_diag) %>%
      mean()
    age_pheno[i,2]  <- subtypes %>%
      filter(V3 == 2) %>%
      left_join(age_data, by = c("V1" = "eid")) %>%
      pull(age_diag) %>%
      mean()
    age_pheno[i,3]  <- subtypes %>%
      filter(V3 == 3) %>%
      left_join(age_data, by = c("V1" = "eid")) %>%
      pull(age_diag) %>%
      mean()
    
    # extract heritability for first subtype
    h2g.1 <- readLines(paste(DIR, ds_id,"/",ds_id,"_h2g_subtp_1.log",sep=""))
    pheh2g[i,1] <- as.numeric(str_split(h2g.1 [length(h2g.1) - 6], "\\s+")[[1]][5])
    se_h2g <- str_split(h2g.1[length(h2g.1) - 6], "\\s+")[[1]][6]
    seh2g[i,1] <- as.numeric(substr(se_h2g, 2,nchar(se_h2g)-1))
    # extract heritability for second subtype
    h2g.2 <- readLines(paste(DIR, ds_id,"/",ds_id,"_h2g_subtp_2.log",sep=""))
    pheh2g[i,2] <- as.numeric(str_split(h2g.2[length(h2g.2) - 6], "\\s+")[[1]][5])
    se_h2g <- str_split(h2g.2[length(h2g.2) - 6], "\\s+")[[1]][6]
    seh2g[i,2] <- as.numeric(substr(se_h2g, 2,nchar(se_h2g)-1))
    # extract genetic correlation (if there is any)
    h2g <- readLines(paste(DIR, ds_id,"/",ds_id,"subtp.log",sep=""))
    idx_table <- which(h2g == "Summary of Genetic Correlation Results" )
    rg[i,1]  <- as.numeric(str_split(h2g[idx_table + 2], "\\s+")[[1]][4])
    serg[i,1]  <- as.numeric(str_split(h2g[idx_table + 2], "\\s+")[[1]][5])
    if(sum(subtypes$V3 == 3)){
      rg[i,2] <- as.numeric(str_split(h2g[idx_table + 3], "\\s+")[[1]][4])
      serg[i,2] <- as.numeric(str_split(h2g[idx_table + 3], "\\s+")[[1]][5])
      # read h2g file for the third topic
      h2g.3 <- readLines(paste(DIR, ds_id,"/",ds_id,"_h2g_subtp_3.log",sep=""))
      pheh2g[i,3] <- as.numeric(str_split(h2g.3[length(h2g.3) - 6], "\\s+")[[1]][5])
      se_h2g <- str_split(h2g.3[length(h2g.3) - 6], "\\s+")[[1]][6]
      seh2g[i,3] <- as.numeric(substr(se_h2g, 2,nchar(se_h2g)-1))
    }
    
    # also save the Fst
    fst <- readLines(paste(DIR, ds_id,"/",ds_id,"subtp_Fst.log",sep=""))
    target_fst <- as.numeric(str_split(fst[length(fst) - 2], "\\s+")[[1]][4])
    permutation_Fst <- read.csv(paste(DIR, ds_id,"/",ds_id,"_permutation_Fst.csv",sep=""))
    p_fst[i]  <- (sum(target_fst < permutation_Fst$weighted_Fst) + 1)/dim(permutation_Fst)[1]
  })
}
df_rg <- data.frame(disease = ds_target[[1]],age = age_pheno, h2g = pheh2g, h2g_se = seh2g, rg = rg, serg = serg, p_fst = p_fst)
write.csv(df_rg ,paste(DIR, "subtypes_rg.csv",sep=""), row.names = FALSE )
