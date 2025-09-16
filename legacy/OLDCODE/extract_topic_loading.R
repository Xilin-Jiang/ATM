# extract topic loading from the target results
# extract age distribution
setwd("/well/mcvean/xilin/Multimorbidity-biobank/")
library(dplyr)
topic_num <- 10
degree_free_num <- 5 # degrees of freedom
rep_ID <- 10

load(paste0("Results/","Run_2rec_PheCode_age_dependent_K",topic_num,"_P",degree_free_num,"_rep",rep_ID, ".RData"))
loadings <- sweep((para$alpha_z - 1), 1, rowSums(para$alpha_z -1), FUN="/")
save_loadings <- data.frame(para$eid, para$eid, loadings) %>%
  rename(FID = para.eid, IID = para.eid.1)
LOC <- "Association_analysis/Phenotypes/"
write.table(save_loadings, paste0(LOC,"ageLDA_topic_loadings",para$K,"_P",para$P,"_rep",para$rep_ID, ".txt"), sep="\t", col.names = FALSE, row.names = FALSE, quote = F)

