library(dplyr)
args <- commandArgs(trailingOnly = TRUE)                                  
ds_id <- as.numeric(args[1])

phecode_age <- read.csv(paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/DiseaseAbove1000occur_PheCode.csv"), header = T)
british_isle <- read.table("/users/mcvean/xilin/xilin/UK_biobank/keepBritish.txt", header = F)

BIA_cases <- phecode_age %>%
  filter(diag_icd10 == ds_id, eid  %in% british_isle$V1) %>%
  select(eid) %>%
  mutate(IID = eid) %>%
  rename(FID = eid) %>%
  mutate(phenotype = 1) 

BIA_controls <- british_isle %>%
  anti_join(BIA_cases, by = c("V1" = "IID")) %>%
  rename(IID = V1, FID = V2) %>%
  mutate(phenotype = 0)

# threshold: 10% of cases (as per BOLT-LMM mannual) or half of the population
# first case: sample half of the population
if(dim(BIA_cases)[1] > 0.05*dim(british_isle)[1] ){
  control_number <- floor(dim(british_isle)[1]/2) -  dim(BIA_cases)[1]
  BIA_controls <- BIA_controls %>%
    sample_n(control_number, replace = F)
}else{ # if there is too few cases just sample 9 control for each cases
  control_number <- 9 * dim(BIA_cases)[1]
  BIA_controls <- BIA_controls %>%
    sample_n(control_number, replace = F)
}

BIA_sample <- bind_rows(BIA_cases, BIA_controls)

# save all the data & cross validation set 
cases_idx <- sample(1:nrow(BIA_cases)) 
control_idx <- sample(1:nrow(BIA_controls)) 

cases_fold_size <- floor(length(cases_idx)/5)  
control_fold_size <- floor(length(control_idx)/5)  

for(csvld in 1:5){
  if(csvld == 5){
    BIA_case_testing <- BIA_cases %>%
      slice(cases_idx[( (csvld-1) * cases_fold_size + 1): length(cases_idx)])
    BIA_ctrl_testing <- BIA_controls %>%
      slice(control_idx[( (csvld-1) * control_fold_size + 1):  length(control_idx)])
  }else{
    BIA_case_testing <- BIA_cases %>%
      slice(cases_idx[( (csvld-1) * cases_fold_size + 1): (csvld * cases_fold_size)])
    BIA_ctrl_testing <- BIA_controls %>%
      slice(control_idx[( (csvld-1) * control_fold_size + 1): (csvld * control_fold_size)])
  }
  
  BIA_training <- BIA_sample %>%
    anti_join(BIA_case_testing, by = "IID") %>%
    anti_join(BIA_ctrl_testing, by = "IID")
  write.table(BIA_training, paste0(ds_id,"/",ds_id, "_csvld_", csvld, ".pheno"), row.names = F, col.names = T, sep = "\t", quote = F)
}
BIA_sample %>%
  select(-phenotype) %>%
  write.table(paste0(ds_id,"/",ds_id, "keep.txt"), row.names = F, col.names = F, sep = "\t", quote = F)



