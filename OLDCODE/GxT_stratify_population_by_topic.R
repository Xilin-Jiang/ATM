library(ggplot2)
require(dplyr)
library(stringr)
library(pdist)
args <- commandArgs(trailingOnly = TRUE)                                  
ds_id <- args[1]

#########################
# testing data loading start
# ds_id <- "250.2" # "401.1"
# rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
# load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
# birthyear <- read.csv("Year_of_birth.csv") %>%
#   select(eid, X31.0.0, X34.0.0, X23104.0.0)
# 
# data <- read.table("../../genetics_longitudinal_data/longitudinal_data/covariates_for_asso_over_time.txt", sep="\t") %>%
#   select(-V3) %>%
#   rename(eid = V1)

# coding data loading end
#########################

# use 4 controls except for hypertension:
if(ds_id == "401.1"){
  control_num <- 1
}else{
  control_num <- 4
}

dir.create(paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/GxTopic/", ds_id))
# using the function of package pdist
rec_data <- read.csv("/users/mcvean/xilin/xilin/Multimorbidity-biobank/rec2subjectAbove1000occur_include_death_PheCode.csv")

# create subtype files
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/GxTopic/"
load("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
  
# j <- match(ds_id, para$list_above500occu$diag_icd10)
# loadings <- para$unlist_zn[para$ds_list[[j]]$id,]

DsData <- rec_data %>%
  filter(diag_icd10 == ds_id)

birthyear <- read.csv("/users/mcvean/xilin/xilin/UK_biobank/Year_of_birth.csv") %>%
  select(eid, X31.0.0, X34.0.0, X23104.0.0)

data <- read.table("/users/mcvean/xilin/xilin/UK_biobank/covariates_for_asso_over_time.txt", sep="\t") %>%
  select(-V3) %>%
  rename(eid = V1)

# compute the patient loadings 
patient_loadings <- sweep((para$alpha_z-1), 1, rowSums(para$alpha_z - 1), FUN="/")
patient_loadings <- bind_cols(eid = para$eid, loadings = patient_loadings)

data_pca <- data %>%
  select(1:12) %>%
  right_join(patient_loadings, by = "eid") %>%
  left_join(birthyear, by=c("eid"="eid"))  %>%
  mutate(SEX = replace(X31.0.0, is.na(X31.0.0), sample(1:2, 1, replace = T)), 
         year_of_birth = replace(X34.0.0, is.na(X34.0.0), mean(X34.0.0, na.rm = T)),
         bmi = replace(X23104.0.0, is.na(X23104.0.0), mean(X23104.0.0, na.rm = T))) %>% 
  mutate(across(3:12, ~ replace(.x, is.na(.x), mean(.x, na.rm = T))))  %>% 
  select(-X31.0.0 , - X34.0.0, - X23104.0.0) # death data also need to be removed.

pca_comp <- data_pca %>%
  select(-eid, -V2) %>%
  prcomp()

# get the control group set
data_pca <- data_pca %>% 
  select(eid) %>%
  cbind(pca_comp[["x"]])

num_covariates <- dim(data_pca)[2]

CaseId <- DsData %>% 
  select(eid)
# this step makes sure only samples that are not of any subtypes will be considered as controls
ControlPCA <- data_pca %>%
  anti_join(CaseId, by = "eid") 
# disease PCA
CasePCA <- data_pca %>%
  semi_join(CaseId, by = "eid")

num_itr <- floor(dim(CasePCA)[1] / 100)
# create a data file to save controls
casectrl_data <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("case", "control"))
for(itr in 1:(num_itr+1)){
  if(itr == (num_itr+1)){
    if(dim(CasePCA)[1] == 100 * num_itr) break
    dist_mat <- as.matrix(pdist(data.matrix(ControlPCA[,2:num_covariates]), 
                                data.matrix(CasePCA[(1 + 100*num_itr):nrow(CasePCA),2:num_covariates])))
    case_list <- CasePCA$eid[(1 + 100*num_itr):nrow(CasePCA)]
  }else{
    dist_mat <- as.matrix(pdist(data.matrix(ControlPCA[,2:num_covariates]), 
                                data.matrix(CasePCA[(1+100*(itr-1)):(100*itr),2:num_covariates])))
    case_list <- CasePCA$eid[(1+100*(itr-1)):(100*itr)]
  }

  dist_mat <- ControlPCA %>%
    select(eid) %>%
    cbind(dist_mat)
  
  for(ColId in 1:(dim(dist_mat)[2] - 1)){
    control_add <- dist_mat %>%
      select(1, ColId+1) %>%
      anti_join(casectrl_data, by = c("eid" = "control")) %>%
      top_n(-control_num) %>%
      select(eid) %>%
      mutate(case = case_list[ColId]) %>%
      rename(control = eid)
    casectrl_data <- rbind(casectrl_data, control_add)
  }
}
write.table(casectrl_data, paste0(DIR, ds_id,"/", ds_id, "_matching_topic_casectrl.txt"), sep="\t", col.names = T, row.names = FALSE)


