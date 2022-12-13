library(ggplot2)
require(dplyr)
library(stringr)
library(pdist)
args <- commandArgs(trailingOnly = TRUE)                                  
ds_id <- args[1]
setwd("/users/mcvean/xilin/xilin/Multimorbidity-biobank/")
dir.create(paste0("Association_analysis/", ds_id))
# using the function of package pdist
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")

# create subtype files
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/Association_analysis/"
subtype_diseases_list <- read.csv("/users/mcvean/xilin/xilin/Multimorbidity-biobank/subtype_disesea_list.csv")
load("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
  
j <- match(ds_id, para$list_above500occu$diag_icd10)
sample_sz <- dim(para$unlist_zn[para$ds_list[[j]]$id,])[1]
cases_eid <- list()
for(topic_id in 1:para$K){
  topic_specific_id <- which(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max) == topic_id)
  topic_specific_set <- para$unlist_Ds_id[para$ds_list[[j]]$id,][topic_specific_id,]
  if(dim(topic_specific_set)[1] > 500 & dim(topic_specific_set)[1]/sample_sz > 0.05){
    cases_eid[[topic_id]] <- topic_specific_set %>%
      mutate(topic = topic_id) %>%
      mutate(Ds_id = para$list_above500occu[Ds_id, 1])
  }
}
cases_eid <- bind_rows(cases_eid)

cases_eid %>%
  mutate(FID = eid, IID = eid) %>%
  select(FID, IID, topic) %>%
  write.table(paste(DIR, ds_id,"/",ds_id,"subtypes.txt",sep=""), sep="\t", col.names = F, row.names = FALSE ,quote = F)

# discrete data type: save phenotype case control data for CC-gwas 
topic_list <- unique(cases_eid$topic)
phenotypes <- cases_eid %>%
  mutate(FID = eid, IID = eid) %>%
  select(FID, IID, topic) 
if(length(topic_list) > 1){
  for(topic_id in topic_list[2:length(topic_list)]){
    phenotypes <- phenotypes %>%
      mutate(!!paste0(topic_id) := if_else(topic == topic_id, 2, 1))
  }
}
phenotypes %>%
  select(-topic) %>%
  write.table(paste(DIR, ds_id,"/",ds_id,"cc_gwas.txt",sep=""), sep="\t", col.names = F, row.names = FALSE ,quote = F)
data.frame(ds_id, topic_list) %>%
  write.table(paste(DIR, ds_id,"/",ds_id,"topic_id_subtp.txt",sep=""), sep="\t", col.names = F, row.names = T ,quote = F)

# continuous data type: save patient loading (mode of dirichlet distribution) for  CC-gwas 
patient_loadings <- sweep((para$alpha_z - 1), 1, rowSums(para$alpha_z -1), FUN="/") 
topic_list <- unique(cases_eid$topic)

df_umap_patient <- data.frame(eid = para$eid, loadings = patient_loadings[,topic_list]) 

subtype_loadings <- cases_eid %>%
  mutate(FID = eid, IID = eid) %>%
  select(FID, IID) %>%
  left_join(df_umap_patient, by = c("FID" = "eid"))

subtype_loadings %>%
  write.table(paste(DIR, ds_id,"/",ds_id,"loadings_cc_gwas.txt",sep=""), sep="\t", col.names = F, row.names = FALSE ,quote = F)


