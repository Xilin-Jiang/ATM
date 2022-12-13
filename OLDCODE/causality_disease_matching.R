library(ggplot2)
require(dplyr)
library(stringr)
library(pdist)
args <- commandArgs(trailingOnly = TRUE)                                  
ds_id <- args[1]
topic_id <- args[2]

dir.create(paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/causality_analysis/", ds_id,"_topic", topic_id))
# using the function of package pdist
rec_data <- read.csv("/users/mcvean/xilin/xilin/Multimorbidity-biobank/rec2subjectAbove1000occur_include_death_PheCode.csv")

# create subtype files
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/causality_analysis/"
common_disease_within_topics <- read.table("/users/mcvean/xilin/xilin/Multimorbidity-biobank/all_disease_topic_list.txt", header = F)
load("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
  
j <- match(ds_id, para$list_above500occu$diag_icd10)
loadings <- para$unlist_zn[para$ds_list[[j]]$id,]

DsData <- rec_data %>%
  filter(diag_icd10 == ds_id)

birthyear <- read.csv("/users/mcvean/xilin/xilin/UK_biobank/Year_of_birth.csv") %>%
  select(eid, X31.0.0, X34.0.0, X23104.0.0)

data <- read.table("/users/mcvean/xilin/xilin/UK_biobank/covariates_for_asso_over_time.txt", sep="\t") %>%
  select(-V3) %>%
  rename(eid = V1)

data_pca <- data %>%
  left_join(birthyear, by=c("eid"="eid"))  %>%
  mutate(SEX = replace(X31.0.0, is.na(X31.0.0), sample(1:2, 1, replace = T)), 
         year_of_birth = replace(X34.0.0, is.na(X34.0.0), mean(X34.0.0, na.rm = T)),
         bmi = replace(X23104.0.0, is.na(X23104.0.0), mean(X23104.0.0, na.rm = T))) %>% 
  select(-X31.0.0 , - X34.0.0, - X23104.0.0) # death data also need to be removed.

pca_comp <- data_pca %>%
  select(-eid, -V2) %>%
  prcomp()

# get the control group set
data_pca <- data_pca %>% 
  select(eid) %>%
  cbind(pca_comp[["x"]])

CaseId <- DsData %>% 
  select(eid)
# this step makes sure only samples that are not of any subtypes will be considered as controls
ControlPCA <- data_pca %>%
  anti_join(CaseId, by = "eid") 
#################################
# filter only the disease of specific topic
#################################
topic_specific_id <- which(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max) == topic_id)
topic_specific_set <- para$unlist_Ds_id[para$ds_list[[j]]$id,][topic_specific_id,]
CasePCA <- data_pca %>%
  semi_join(topic_specific_set, by = "eid")

num_itr <- floor(dim(CasePCA)[1] / 100)
# create a data file to save controls
control <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("eid"))
for(itr in 1:(num_itr+1)){
  if(itr == (num_itr+1)){
    if(dim(CasePCA)[1] == 100 * num_itr) break
    dist_mat <- as.matrix(pdist(data.matrix(ControlPCA[,2:44]), 
                                data.matrix(CasePCA[(1 + 100*num_itr):nrow(CasePCA),2:44])))
  }else
    dist_mat <- as.matrix(pdist(data.matrix(ControlPCA[,2:44]), 
                                data.matrix(CasePCA[(1+100*(itr-1)):(100*itr),2:44])))
  
  dist_mat <- ControlPCA %>%
    select(eid) %>%
    cbind(dist_mat)
  for(ColId in 1:(dim(dist_mat)[2] - 1)){
    control_add <- dist_mat %>%
      select(1, ColId+1) %>%
      anti_join(control, by = "eid") %>%
      top_n(-4) %>%
      select(eid)
    control <- rbind(control, control_add)
  }
}
##########################
# need the diagnose age
# new file is created in gwas_25
case <- CasePCA %>% 
  select(eid) %>% 
  mutate(FID = eid, phenotype = 2)

control_for_gwas <- control %>%
  mutate(FID = eid, phenotype = 1)

DataGWAS <- case %>% rbind(control_for_gwas) 

write.table(control_for_gwas, paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/causality_analysis/", ds_id,"_topic", topic_id,"/", ds_id, "_topic", topic_id, "_control.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
write.table(case, paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/causality_analysis/", ds_id,"_topic", topic_id,"/", ds_id, "_topic", topic_id, "_case.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
write.table(DataGWAS, paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/causality_analysis/", ds_id,"_topic", topic_id,"/", ds_id, "_topic", topic_id, "_pheno_over_age.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)


