require(dplyr)
library(stringr)
setwd("/well/mcvean/xilin/Multimorbidity-biobank/")
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/GxTopic/"
# using the function of package pdist
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
birthyear <- read.csv("/users/mcvean/xilin/xilin/UK_biobank/Year_of_birth.csv") %>%
  select(eid, X31.0.0, X34.0.0, X23104.0.0)
rec_data %>%
  group_by(eid) %>% 
  slice(1) %>% 
  ungroup %>%
  select(eid) %>%
  mutate(IID = eid) %>%
  left_join(select(birthyear, eid, X23104.0.0), by = "eid") %>%
  write.table(paste(DIR, "BMI_pheno.txt",sep=""), sep="\t", col.names = F, row.names = FALSE ,quote = F)

# create topic files
load("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
subtype <- read.table(paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Association_analysis/278.1/278.1topic_id_subtp.txt"))
sex_PC <- read.table("/users/mcvean/xilin/xilin/UK_biobank/covariates_for_asso_over_time.txt") %>%
  select(c(1,3:13))
patient_loadings <- sweep((para$alpha_z - 1), 1, rowSums(para$alpha_z -1), FUN="/") 

for(topic_id in 1:para$K){
  data.frame(eid = para$eid, loadings = patient_loadings[,topic_id]) %>%
    mutate(FID = eid, IID = eid) %>%
    select(FID, IID, loadings) %>%
    left_join(sex_PC, by = c("FID" = "V1")) %>%
    write.table(paste(DIR, "topic_covariates/GxTopic_covariates_topic", topic_id,".txt",sep=""), sep="\t", col.names = F, row.names = FALSE ,quote = F)
  
  data.frame(eid = para$eid, loadings = patient_loadings[,topic_id], 
             loadingsqure = patient_loadings[,topic_id]^2) %>%
    mutate(FID = eid, IID = eid) %>%
    select(FID, IID, loadings, loadingsqure) %>%
    left_join(sex_PC, by = c("FID" = "V1")) %>%
    write.table(paste(DIR, "topic_covariates/GxTopic_nonlinear_topic", topic_id,".txt",sep=""), sep="\t", col.names = F, row.names = FALSE ,quote = F)
}



