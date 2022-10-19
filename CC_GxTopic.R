library(ggplot2)
require(dplyr)
library(stringr)
library(pdist)
args <- commandArgs(trailingOnly = TRUE)                                  
ds_id <- args[1]
setwd("/well/mcvean/xilin/Multimorbidity-biobank/")
# using the function of package pdist
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")

# create subtype files
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/Association_analysis/"
subtype_diseases_list <- read.csv("/users/mcvean/xilin/xilin/Multimorbidity-biobank/subtype_disesea_list.csv")
load("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")

sex_PC <- read.table("/users/mcvean/xilin/xilin/UK_biobank/covariates_for_asso_over_time.txt") %>%
  select(c(1,3:13))

j <- match(ds_id, para$list_above500occu$diag_icd10)
subtype <- read.table(paste0(DIR, ds_id,"/", ds_id, "topic_id_subtp.txt"))
patient_loadings <- sweep((para$alpha_z - 1), 1, rowSums(para$alpha_z -1), FUN="/") 
if(dim(subtype)[1] > 1){
  try({
    for(ccgwas in 1:(dim(subtype)[1])){
      assoc.linear <- read.table(paste0(DIR, ds_id,"/", ds_id, "linear_ccgwas_", ccgwas, ".assoc.linear"), header = T)
      assoc.linear %>% 
        filter(P < 5 * 10^(-8)) %>%
        select(SNP) %>%
        write.table(paste(DIR, ds_id,"/",ds_id,"_hits_ccgwas_", ccgwas,".txt",sep=""), sep="\t", col.names = F, row.names = FALSE ,quote = F)
      data.frame(eid = para$eid, loadings = patient_loadings[,subtype$V3[ccgwas]]) %>%
        mutate(FID = eid, IID = eid) %>%
        select(FID, IID, loadings) %>%
        left_join(sex_PC, by = c("FID" = "V1")) %>%
        write.table(paste(DIR, ds_id,"/",ds_id,"_GxTopic_covariates_", ccgwas,".txt",sep=""), sep="\t", col.names = F, row.names = FALSE ,quote = F)
    }
  })
}


