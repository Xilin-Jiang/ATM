library(dplyr)
args <- commandArgs(trailingOnly = TRUE)                                  
ds_id <- as.numeric(args[1])
topic_id <- args[2]

dir.create(paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/LDSC_SEG/", ds_id,"_topic", topic_id))

phecode_age <- read.csv(paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/DiseaseAbove1000occur_PheCode.csv"), header = T)
british_isle <- read.table("/users/mcvean/xilin/xilin/UK_biobank/keepBritish.txt", header = F)
load("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")

BIA_cases <- phecode_age %>%
  filter(diag_icd10 == ds_id, eid  %in% british_isle$V1) %>%
  select(eid) %>%
  mutate(IID = eid) %>%
  rename(FID = eid) %>%
  mutate(phenotype = 1) 

# make sure the controls does not include cases from other subtypes
BIA_controls <- british_isle %>%
  anti_join(BIA_cases, by = c("V1" = "IID")) %>%
  rename(IID = V1, FID = V2) %>%
  mutate(phenotype = 0)

if(topic_id == "all"){
  BIA_cases_subtype <- BIA_cases
}else{
  topic_id <- as.numeric(topic_id)
  j <- match(ds_id, para$list_above500occu$diag_icd10)
  topic_specific_id <- which(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max) == topic_id)
  topic_specific_set <- para$unlist_Ds_id[para$ds_list[[j]]$id,][topic_specific_id,]
  BIA_cases_subtype <- BIA_cases %>%
    semi_join(topic_specific_set, by = c("FID" = "eid"))
}
 

# threshold: 10% of cases (as per BOLT-LMM mannual) or half of the population
# first case: sample half of the population
if(dim(BIA_cases_subtype)[1] > 0.05*dim(british_isle)[1] ){
  control_number <- floor(dim(british_isle)[1]/2) -  dim(BIA_cases_subtype)[1]
  BIA_controls <- BIA_controls %>%
    sample_n(control_number, replace = F)
}else{ # if there is too few cases just sample 9 control for each cases
  control_number <- 9 * dim(BIA_cases_subtype)[1]
  BIA_controls <- BIA_controls %>%
    sample_n(control_number, replace = F)
}

BIA_sample <- bind_rows(BIA_cases_subtype, BIA_controls)

BIA_sample %>%
  write.table(paste0(ds_id,"_topic", topic_id,"/",ds_id,"_topic", topic_id, "_keep.txt"), row.names = F, col.names = T, sep = "\t", quote = F)



