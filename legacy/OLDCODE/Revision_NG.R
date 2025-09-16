##########################################################
# preparing disease list for AllOfUs & CVD-COVID-UK
##########################################################

# step 1: get the list of ICD-10 codes (as opposed to phecodes)
ICD_list1 <- ATM::phecode_icd10cm %>% 
  select(ICD10, phecode) %>%
  filter(phecode %in% ATM::UKB_349_disease$diag_icd10)

ICD_list2 <- ATM::phecode_icd10 %>% 
  select(ICD10, PheCode) %>%
  filter(PheCode %in% ATM::UKB_349_disease$diag_icd10) %>%
  rename(phecode = PheCode)

ICD_list3 <- ATM::short_icd10cm %>% 
  select(ICD10, parent_phecode) %>%
  filter(!(ICD10 %in% ICD_list1$ICD10), !(ICD10 %in% ICD_list2$ICD10) ) %>%
  filter(parent_phecode %in% ATM::UKB_349_disease$diag_icd10)%>%
  rename(phecode = parent_phecode)

phe_phecode <- read.csv("info_phenotype_phecode.csv")
ICD_list <- bind_rows(list(ICD_list1, ICD_list2, ICD_list3)) %>%
  group_by(ICD10) %>%
  slice(1) %>%
  ungroup() %>%
  select(ICD10, phecode) %>%
  left_join(select(phe_phecode, phecode, phenotype), by = "phecode") 
 

write.csv(ICD_list, file = "icd10_list_baseline_UKBB.csv", row.names = F)
CVD_covid_list <- ICD_list %>%
  filter(stringr::str_detect(ICD10, "^[A-N]")) %>%
  filter(stringr::str_length(ICD10) < 5) %>%
  mutate(phenotype = str_replace_all(phenotype, ";", "")) %>%
  mutate(phenotype = str_replace_all(phenotype, ",", ""))
write.table(CVD_covid_list, file = "CVD_COVID_UK_icd10_list_348_phecodes.txt", row.names = F, quote = F, sep = "\t", col.names = F)


# load the AllOfUs data
topic_number <- 13
load(paste0("Revision_NG/ATM-SNOMED-349_topic_num-", topic_number,"_dof-5_CVB-5.Rdata"))
AllOfUs_results <- ATM
disease_list <- AllOfUs_results$ds_list %>%
  dplyr::left_join(mutate(ATM::disease_info_phecode_icd10, phecode = as.character(phecode)), by = c("diag_icd10"="phecode")) %>%
  dplyr::pull(phenotype)
topic_id <- 13
ATM::plot_age_topics(disease_names = disease_list, AllOfUs_results$topic_loadings[15:85,,topic_id], plot_title = paste0("topic ", topic_id),
                     top_ds = 10,start_age = 15)
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")

#################################
# revision Figure: AllOfUs results
#################################
# plot cross-population prediction odds ratio AllOfUs to UKBB
predict_OR_ukbb_aouTopicLoading <- list()
for(topic_number in 5:15){
  print(paste0("running prediction OR for topic number: ", topic_number))
  load(paste0("Revision_NG/ATM-SNOMED-349_topic_num-", topic_number,"_dof-5_CVB-5.Rdata"))
  AllOfUs_results <- ATM
  predict_OR_ukbb_aouTopicLoading[[topic_number]] <- ATM::prediction_OR(rec_data, ds_list = AllOfUs_results$ds_list, topic_loadings = AllOfUs_results$topic_loadings, max_predict = 10)
}
save(predict_OR_ukbb_aouTopicLoading, file = "Revision_NG/prediction_AllOfUs_to_UKBB.Rdata")
predict_OR_across_populations <- tibble(OR_top1= as.numeric(), OR_top2= as.numeric(), OR_top5= as.numeric(), topic_number = as.integer()) 
for(topic_number in 5:15){
  predict_OR_across_populations <- predict_OR_across_populations %>% 
    add_row(OR_top1= predict_OR_ukbb_aouTopicLoading[[topic_number]]$OR_top1, 
            OR_top2= predict_OR_ukbb_aouTopicLoading[[topic_number]]$OR_top2, 
            OR_top5= predict_OR_ukbb_aouTopicLoading[[topic_number]]$OR_top5, topic_number = topic_number)
}
predict_OR_across_populations <- predict_OR_across_populations %>%
  pivot_longer(cols = - topic_number, names_to = "Target_size", values_to = "Odds_Ratio") %>%
  mutate(topic_number = factor(topic_number, levels = 5:15))
plt_predic_OR <- ggplot(data=predict_OR_across_populations,aes(x=topic_number, y=Odds_Ratio, color = Target_size)) +
  geom_point(size = 3, alpha=0.6) +
  # geom_jitter(width=0.1, size = 0.2, alpha=0.4) + 
  labs(x = "Number of topics", y = "Predicting Odds Ratio") + 
  scale_color_manual(name = "Prediction set",values=c("OR_top1" = blue, "OR_top2" = green, "OR_top5" = red )) + 
  theme(panel.background=element_blank()) + 
  geom_hline(yintercept = 1, linetype='dashed', color = red)
ggsave(filename = "Revision_NG/Prediction_OR_AOU_on_UKB.png", plot = plt_predic_OR, width = 5, height = 5)

# plot jackknife prediction odds ratio
jkf_predict_OR_across_populations <- tibble(OR_top1= as.numeric(), OR_top2= as.numeric(), OR_top5= as.numeric(), topic_number = as.integer()) 
for(topic_number in 5:15){
  print(paste0("running prediction OR for topic number: ", topic_number))
  load(paste0("Results_NG_revision/prediction_AllOfUs_to_UKBB_topic_num", topic_number,".Rdata"))
  for(jkf_id in 1:length(jkf_predict_OR_ukbb_aouTopicLoading)){
    jkf_predict_OR_across_populations <- jkf_predict_OR_across_populations %>% 
      add_row(OR_top1= jkf_predict_OR_ukbb_aouTopicLoading[[jkf_id]]$OR_top1, 
              OR_top2= jkf_predict_OR_ukbb_aouTopicLoading[[jkf_id]]$OR_top2, 
              OR_top5= jkf_predict_OR_ukbb_aouTopicLoading[[jkf_id]]$OR_top5, topic_number = topic_number)
  }

}
jkf_predict_OR_across_populations <- jkf_predict_OR_across_populations %>%
  mutate(topic_number = as.factor(topic_number))
  
# plot p-value for the 13 topic optimal model
jkf_13aou_topic <- jkf_predict_OR_across_populations %>%
  filter(topic_number == 13) %>%
  summarise(se_top2 =sqrt( 9/10 * sum((OR_top2-mean(OR_top2))^2) ) , mean_top2 = mean(OR_top2))
1- pnorm( ( jkf_13aou_topic$mean_top2 - 1)/jkf_13aou_topic$se_top2 )

plt_OR <- ggplot(data=jkf_predict_OR_across_populations,aes(x=topic_number, y=OR_top2 )) +
  geom_boxplot(fill = red, alpha = 0.6, outlier.shape = NA, width=.5,position=position_dodge(width=0.85)) +
  geom_jitter(size = 0.5, alpha=0.6) + 
  # geom_jitter(width=0.1, size = 0.2, alpha=0.4) + 
  labs(x = "Number of topics", y = "Cross population Prediction OR") + 
  theme(panel.background=element_blank()) + 
  geom_hline(yintercept = 1, linetype='dashed', color = red)
ggsave(filename = "Revision_NG/AOU_cross_prediction_jkf.png", plot = plt_OR, width = 5, height = 5)


# comparing the ELBO and cross population prediction accuracy for 5-15 topic numbers. 
run_ELBO <- list()
for(topic_number in 5:15){
  load(paste0("Revision_NG/ATM-SNOMED-349_topic_num-", topic_number,"_dof-5_CVB-5.Rdata"))
  AllOfUs_results <- ATM
  run_ELBO[[topic_number]] <- AllOfUs_results$multiple_run_ELBO_compare %>%
    mutate(topic_number = topic_number)
}
run_ELBO <- bind_rows(run_ELBO) %>%
  mutate(topic_number = as.factor(topic_number))
plt_ELBO <- ggplot(data=run_ELBO,aes(x=topic_number, y=lower_bound )) +
  geom_boxplot(fill = red, alpha = 0.6, outlier.shape = NA, width=.8,position=position_dodge(width=0.85)) +
  geom_jitter(size = 0.5, alpha=0.6) + 
  # geom_jitter(width=0.1, size = 0.2, alpha=0.4) + 
  labs(x = "Number of topics", y = "ELBO") + 
  theme(panel.background=element_blank()) 
ggsave(filename = "Revision_NG/AOU_ELBO.png", plot = plt_ELBO, width = 5, height = 5)
#########################################################
# compute the cosine similarity for AOU topics versus UKB
#########################################################
run_ELBO %>% arrange(desc(lower_bound))
topic_number <- 13
load(paste0("Revision_NG/ATM-SNOMED-349_topic_num-", topic_number,"_dof-5_CVB-5.Rdata"))
AllOfUs_results <- ATM
AOU_loading_mean <- apply( AllOfUs_results$topic_loadings[30:80, , ], c(2,3), mean)
order_topic_diseases <- sapply(1:dim(AOU_loading_mean)[1] , function(i) which(AllOfUs_results$ds_list$diag_icd10[i] ==  ATM::UKB_349_disease$diag_icd10))
UKB_loading_mean <- ATM::UKB_HES_10topics[30:80,,]
UKB_loading_mean <- apply( UKB_loading_mean[, order_topic_diseases, ], c(2,3), mean)
correlation_UKB_AOU <- cor(UKB_loading_mean, AOU_loading_mean)

apply(correlation_UKB_AOU, 1, max)
order_aou <- c()
for(UKB_id in 1:10){
  new_element <- setdiff(order(correlation_UKB_AOU[UKB_id,], decreasing = T), order_aou)[1]
  order_aou <- c(order_aou, new_element)
}
order_aou <- c(order_aou, setdiff(1:topic_number, order_aou))
correlation_UKB_AOU <- correlation_UKB_AOU[, order_aou]

# plot the correlation
longData<-reshape2::melt(correlation_UKB_AOU) 
plt <- ggplot() + geom_tile(data = longData, aes(x = Var2, y = Var1,  fill = value, width = 0.9, height = 0.9)) + 
  scale_fill_gradient2(low=blue, mid="white", high=red,  midpoint=0)+
  labs(x="", y="", title="") +
  scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
# plt <- ggplot() + geom_tile(data = longData, aes(x = Var2, y = Var1,  alpha = value, width = 0.9, height = 0.9), fill=red) + 
#   scale_alpha_continuous(range = c(0,1)) +
#   labs(x="", y="", title="") +
#   scale_x_discrete(expand=c(0,0)) + 
#   scale_y_discrete(expand=c(0,0)) +
#   theme(axis.line=element_blank(),axis.text.x=element_blank(),
#         axis.text.y=element_blank(),axis.ticks=element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),legend.position="none",
#         panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),plot.background=element_blank())
ggsave("Revision_NG/topic_correlation.png", plt, width = 12, height = 10)

# make a legend 
longData<-reshape2::melt(matrix(seq(from = min(correlation_UKB_AOU), to = max(correlation_UKB_AOU), length.out = 100),ncol = 1)) 
plt <- ggplot() + 
  geom_tile(data = longData, aes(x = Var2, y = Var1,fill = value,width = 1)) + 
  scale_fill_gradient2(low=blue, mid="white", high=red,  midpoint=0)+
  labs(x="", y="", title="") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
# below for using alpha
# longData<-reshape2::melt(matrix(c(1:100/100),ncol = 1)) 
# plt <- ggplot() + 
#   geom_tile(data = longData, aes(x = Var2, y = Var1,alpha = value,width = 0.9), fill=red) + 
#   scale_alpha_continuous(range = c(0,1)) +
#   labs(x="", y="", title="") +
#   theme(axis.line=element_blank(),axis.text.x=element_blank(),
#         axis.text.y=element_blank(),axis.ticks=element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),legend.position="none",
#         panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),plot.background=element_blank())
ggsave(paste0("Revision_NG/correlation_color_legends.png"), plt, width = 1, height = 5)
# legends and text for the figure
sum(apply(correlation_UKB_AOU, 1, max) > 0.5) 
mean(apply(correlation_UKB_AOU, 1, max)) 
max(correlation_UKB_AOU)
min(correlation_UKB_AOU)

topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")

# save the topic correlation 
rownames(correlation_UKB_AOU) <- topic_name
colnames(correlation_UKB_AOU) <- sapply(1:13, function(x) paste0("AOU.topic.", x))
write.csv(correlation_UKB_AOU, file = "Revision_NG/correlation_aou_ukb_topic_loading.csv")

#########################################################
# compute the cosine similarity for AOU topics versus UKB
# at age slice 50, 60, and 70
#########################################################
for(age in c(50,60,70)){
  AOU_loading_50 <- AllOfUs_results$topic_loadings[age, , ]
  order_topic_diseases <- sapply(1:dim(AOU_loading_mean)[1] , function(i) which(AllOfUs_results$ds_list$diag_icd10[i] ==  ATM::UKB_349_disease$diag_icd10))
  UKB_loading_50 <- ATM::UKB_HES_10topics[age,,]
  UKB_loading_50 <- UKB_loading_50[order_topic_diseases, ]
  correlation_UKB_AOU_ageslice <- cor(UKB_loading_50, AOU_loading_50)
  print(apply(correlation_UKB_AOU_ageslice, 1, max))
  print(c(max(correlation_UKB_AOU_ageslice), min(correlation_UKB_AOU_ageslice)))
  correlation_UKB_AOU_ageslice <- correlation_UKB_AOU_ageslice[, order_aou]
  longData<-reshape2::melt(correlation_UKB_AOU_ageslice) 
  plt <- ggplot() + geom_tile(data = longData, aes(x = Var2, y = Var1,  fill = value, width = 0.9, height = 0.9)) + 
    scale_fill_gradient2(low=blue, mid="white", high=red,  midpoint=0)+
    labs(x="", y="", title="") +
    scale_x_discrete(expand=c(0,0)) + 
    scale_y_discrete(expand=c(0,0)) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  ggsave(paste0("Revision_NG/topic_correlation_age_slice", age, ".png"), plt, width = 12, height = 10)
}

# save topic loadings as a csv
topic_number <- 13
load(paste0("Revision_NG/ATM-SNOMED-349_topic_num-", topic_number,"_dof-5_CVB-5.Rdata"))
AllOfUs_results <- ATM
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD", "unmatched 1", "unmatched 2", "unmatched 3")
ds.system <- AllOfUs_results$ds_list %>%
  mutate(diag_icd10 = as.numeric(diag_icd10)) %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") ) %>%
  select(-exclude_range, -occ, -exclude_name) %>%
  rename(Phecode = diag_icd10, Description=phenotype)
topic_loadings <- list()
for(topic_id in 1:topic_number){
  topic_id <- order_aou[topic_id]
  trajs <- AllOfUs_results$topic_loadings[20:85,,topic_id] %>%
    t() %>%
    data.frame()# trajectories
  names(trajs) <- 20:85
  topic_loadings[[topic_id]] <- ds.system %>%
    mutate(Topic = topic_name[topic_id]) %>%
    select(Topic, Phecode, ICD10, Description) %>%
    bind_cols(trajs)
}
topic_loadings %>% bind_rows() %>%
  write.csv("Revision_NG/AOU_topic_loadings.csv", row.names = F)
####################################
# save a prevalence sup table: AOU vs. UKB
####################################
phe_phecode <- read.csv("info_phenotype_phecode.csv")
topic_number <- 13
load(paste0("Revision_NG/ATM-SNOMED-349_topic_num-", topic_number,"_dof-5_CVB-5.Rdata"))
AllOfUs_results <- ATM

AllOfUs_results$ds_list %>% 
  mutate(AOU_occ = occ/AllOfUs_results$patient_number, diag_icd10 = as.numeric(diag_icd10)) %>%
  select(-occ) %>%
  left_join(ATM::UKB_349_disease, by = "diag_icd10") %>%
  mutate(UKB_occ = occ/282957) %>%
  select(-occ) %>%
  arrange(diag_icd10) %>%
  rename(phecode = diag_icd10) %>%
  left_join(select(phe_phecode, ICD10, phecode, phenotype), by = "phecode") %>%
  select(phecode, ICD10, phenotype, AOU_occ, UKB_occ) %>%
  write.csv("Revision_NG//AOU_UKB_prevelence_comparison.csv", row.names = F)


##################################################
# plot all topics loadings comparison 
##################################################
phe_phecode <- read.csv("info_phenotype_phecode.csv")
topic_number <- 13
load(paste0("Revision_NG/ATM-SNOMED-349_topic_num-", topic_number,"_dof-5_CVB-5.Rdata"))
AllOfUs_results <- ATM
AOU_loading_mean <- apply( AllOfUs_results$topic_loadings[30:80, , ], c(2,3), mean)
order_topic_diseases <- sapply(1:dim(AOU_loading_mean)[1] , function(i) which(AllOfUs_results$ds_list$diag_icd10[i] ==  ATM::UKB_349_disease$diag_icd10))
UKB_loading_mean <- ATM::UKB_HES_10topics[30:80,,]
UKB_loading_mean <- apply( UKB_loading_mean[, order_topic_diseases, ], c(2,3), mean)
correlation_UKB_AOU <- cor(UKB_loading_mean, AOU_loading_mean)

best_aou_topic <- apply(correlation_UKB_AOU, 1, which.max)
apply(correlation_UKB_AOU, 1, max)
disease_list <- AllOfUs_results$ds_list %>%
  dplyr::left_join(mutate(ATM::disease_info_phecode_icd10, phecode = as.character(phecode)), by = c("diag_icd10"="phecode")) %>%
  dplyr::pull(phenotype)
for(ukb_id in 1:10){
  aou_id <- best_aou_topic[ukb_id]
  UKB_loadings <- ATM::UKB_HES_10topics[, order_topic_diseases, ]
  # plot UKB topic
  plt_UKB <- ATM::plot_age_topics(disease_names = disease_list, UKB_loadings[30:80,,ukb_id], plot_title = paste0("UKB topic ", topic_name[ukb_id]),
                                  top_ds = 7,start_age = 30)
  ggsave(paste0("Revision_NG/topic_loadings/UKB_topic_", topic_name[ukb_id], ".png"), plt_UKB, width = 8, height = 8)
  # plot AOU topic
  plt_AOU <- ATM::plot_age_topics(disease_names = disease_list, AllOfUs_results$topic_loadings[15:85,,aou_id], plot_title = paste0("AOU topic corrresponding to ", topic_name[ukb_id]),
                                  top_ds = 7,start_age = 15)
  ggsave(paste0("Revision_NG/topic_loadings/AOU_topic_", topic_name[ukb_id], ".png"), plt_AOU, width = 8, height = 8)
}

# per requirement of the NG editor, remove all the grids
source("plotting_functions.R")
pal_age <- colorRampPalette(c(blue, red, orange, purple, green))
pal_age_vector <- pal_age(length(disease_list))
for(ukb_id in 1:10){
  aou_id <- best_aou_topic[ukb_id]
  UKB_loadings <- ATM::UKB_HES_10topics[, order_topic_diseases, ]
  # plot AOU topic
  plt_AOU <- plot_age_topics(icd10 = disease_list, AllOfUs_results$topic_loadings[15:85,,aou_id], pal_age_vector, plot_title = paste0("AOU topic corrresponding to ", topic_name[ukb_id]), start_age = 15,top_ds = 7)
  
  # plt_AOU <- ATM::plot_age_topics(disease_names = disease_list, AllOfUs_results$topic_loadings[15:85,,aou_id], plot_title = paste0("AOU topic corrresponding to ", topic_name[ukb_id]),
  #                                 top_ds = 7,start_age = 15)
  ggsave(paste0("Revision_NG/topic_loadings/AOU_topic_", topic_name[ukb_id], ".png"), plt_AOU, width = 8, height = 8)
}


#####################################################
# code for Arun: running five fold cross validation
#####################################################
cv_id <- 1 # cross validation id, should be 1:5
topic_number <- 5 # should be from 5 to 16, and 20
df_ATM <- HES_age_example # please use the data! XXXXX
testing_ids <-  df_ATM %>% 
  group_by(eid) %>% 
  slice(1) %>% 
  ungroup()
fold_size <- floor(dim(testing_ids)[1]/5) 
if(cv_id == 5){
  testing_ids <- testing_ids%>% 
    slice(( (cv_id-1) * fold_size + 1): dim(testing_ids)[1])
}else{
  testing_ids <- testing_ids%>%  
    slice(( (cv_id-1) * fold_size + 1): (cv_id * fold_size))
}

training_data <- HES_age_example %>% 
  anti_join(testing_ids, by = "eid")
testing_data <- HES_age_example %>% semi_join(testing_ids, by = "eid")
ATM_results <- wrapper_ATM(rec_data=training_data, topic_num = topic_number, CVB_num = 5)
testing_prediction_OR <- prediction_OR(testing_data = testing_data, ds_list = ATM_results$ds_list, topic_loadings = ATM_results$topic_loadings)
########
# save both testing_prediction_OR, in case there is things wrong, perhaps could also save ATM_results  
########



#############################################
# plot prediction odds ratio of AOU
#############################################
df_predict_lik_P_K <- data_frame(df_K = as.integer(), 
                                 OR_top1 = as.numeric(),OR_top2 = as.numeric(),OR_top5 = as.numeric())
rep_number <- 5
for(K in c(5:20)){
    for(rep_id in 1:rep_number){
      try({
        load(paste0("Results_NG_revision/cross_validation_AOU_prediction_OR/AOU_predictionOR_topic_num", K, "cv", rep_id, ".Rdata"))
        df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
          add_row(df_K = K, 
                  OR_top1 = AOU_prediction_OR$OR_top1, 
                  OR_top2 = AOU_prediction_OR$OR_top2, 
                  OR_top5 = AOU_prediction_OR$OR_top5, 
          )
      })
  }
}
df_boxplot <- df_predict_lik_P_K %>%
  filter(df_K  < 21, df_K > 4) %>% # help visualise one or two Ks
  mutate( df_K = as.factor(df_K)) 
plt <- ggplot(data=df_boxplot,aes(x=df_K, y=OR_top2)) +
  geom_boxplot(fill = red, alpha = 0.6, outlier.shape = NA, width=.8,position=position_dodge(width=0.85)) +
  geom_jitter(size = 0.5, alpha=0.6) + 
  # geom_jitter(width=0.1, size = 0.2, alpha=0.4) + 
  labs(x = "Number of topics", y = "Predicting Odds Ratio") + 
  theme(panel.background=element_blank()) 
  
ggsave(paste0("Revision_NG/AOU_prediction_onebyone_OR_topic_number.png"), plt, width = 6, height = 6)


##############################################################################
# plotting the subtype distribution; similar to Fig 4A. (need to get Arun's help)
##############################################################################
# first get the best runs
topic_number <- 13
load(paste0("Results_NG_revision/Full_run_para_AOU/multirun5K", topic_number,"_P5.Rdata"))
best_id <- which.max(unlist(lapply(multi_runs[[2]], function(x) dplyr::pull(dplyr::slice_tail(x), Lower_bound))))
load(paste0("Results_NG_revision/Full_run_para_AOU/fullmodel_ATM_cvbrep_", best_id, "K", topic_number, "_P5.Rdata"))
AOU2UKB_rec_data <- para$unlist_Ds_id %>% 
  mutate(diag_icd10 = para$list_above500occu$diag_icd10[Ds_id]) %>%
  select(-Ds_id)
write.csv(AOU2UKB_rec_data, "Results_NG_revision/AOU_rec_data.csv", row.names = F)

#################################################################
# save the para for the best AOU run, which is 13 topics
#################################################################
source("topic_functions.R")
topic_number <- 13
load(paste0("Revision_NG/ATM-SNOMED-349_topic_num-", topic_number,"_dof-5_CVB-5.Rdata"))
AllOfUs_results <- ATM
ds_list <- AllOfUs_results$ds_list
topics <- AllOfUs_results$topic_loadings
data <- read.csv("Results_NG_revision/AOU_rec_data.csv") %>%
  filter(diag_icd10 %in% ds_list$diag_icd10)
para <- topic_init_age(data, ds_list, dim(topics)[length(dim(topics))], degree_free_num = 5) # for internal note: degree_free_num doesn't really matter in this case
# update beta_w: list of Ns-by-K
para$beta_w_full <- apply(topics, 3, function(x)
  x[as.matrix(select(para$unlist_Ds_id, age_diag, Ds_id))])
para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])

# update z_n until convergence
para$max_itr <- 100
para$lb <- data.frame("Iteration" = 0,"Lower_bound" = CVB_lb(para))
para$tol <- 10^(-7)
for(itr in 1:para$max_itr){
  print(paste0("Interation: ",itr))
  para <- CVB0_E_zn(para) # we choose CVB0 as papers shown it could converge quicker
  para$lb[nrow(para$lb) + 1,] <- c(itr, CVB_lb(para))
  curr_lb <- pull(filter(para$lb, Iteration == itr), Lower_bound)
  prev_lb <- pull(filter(para$lb, Iteration == (itr - 1 )), Lower_bound)
  print(paste0("Current Lower bound ", curr_lb, " at iteration: ",itr))
  try({
    if(is.finite((curr_lb - prev_lb)) & abs(curr_lb - prev_lb)/abs(prev_lb) < para$tol ){
      print(paste0("Optimization converged at step ", itr))
      break
    }
  })
}
para$pi_beta_basis <- topics
save(para, file = paste0("Results_NG_revision/full_para_ATM-SNOMED-349_topic_num-", topic_number,"_dof_5_best_model.Rdata"))
best_AOU_model_path <- "Results_NG_revision/full_para_ATM-SNOMED-349_topic_num-13_dof_5_best_model.Rdata"


# get the para for AOU
load(best_AOU_model_path)
para_aou <- para
AOU_loading_mean <- apply( para_aou$pi_beta_basis[30:80, , ], c(2,3), mean)
order_topic_diseases <- sapply(1:dim(AOU_loading_mean)[1] , function(i) which(para_aou$list_above500occu$diag_icd10[i] ==  ATM::UKB_349_disease$diag_icd10))
UKB_loading_mean <- ATM::UKB_HES_10topics[30:80,,]
UKB_loading_mean <- apply( UKB_loading_mean[, order_topic_diseases, ], c(2,3), mean)
correlation_UKB_AOU <- cor(UKB_loading_mean, AOU_loading_mean)

order_aou <- c()
for(UKB_id in 1:10){
  new_element <- setdiff(order(correlation_UKB_AOU[UKB_id,], decreasing = T), order_aou)[1]
  order_aou <- c(order_aou, new_element)
}
order_aou <- c(order_aou, setdiff(1:topic_number, order_aou))

topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")

### making the figure similar to Fig 5A
subtypes_list <- read.table( "subtype_disease_topic_list.txt", sep="\t", header = F)
aou_subtype_diseases_list <- subtypes_list %>%
  group_by(V1) %>%
  slice(1) %>%
  filter(V1 %in% para_aou$list_above500occu$diag_icd10)

library(colorBlindness)
displayAvailablePalette(color="white")

# for the subtypes, map the 13 AOU topics to 10 UKB topics
aou2ukb_mapping <- apply(correlation_UKB_AOU, 2, which.max)
# # small function mapps aou_idx when the order_aou has been computed
# map_aou2ukb_topic <- function(aou_idx, order_aou){
#   return(sapply(aou_idx, function(x) which(x == order_aou)))
# }

for(idx in 1:dim(aou_subtype_diseases_list)[1]){
  ds_id <- aou_subtype_diseases_list$V1[idx]
  j <- match(as.numeric(ds_id), para_aou$list_above500occu$diag_icd10)
  cases_eid <- list()
  for(topic_id in 1:para_aou$K){
    topic_specific_id <- which(apply(para_aou$unlist_zn[para_aou$ds_list[[j]]$id,], 1, which.max) == topic_id)
    topic_specific_set <- para_aou$unlist_Ds_id[para_aou$ds_list[[j]]$id,][topic_specific_id,]
    
    cases_eid[[topic_id]] <- topic_specific_set %>%
      mutate(topic = topic_id) %>%
      mutate(Ds_id = para_aou$list_above500occu[Ds_id, 1])
  }
  cases_eid <- bind_rows(cases_eid)
  ###############################################################
  ###############################################################
  ###############################################################
  # change the topic index that matches UKBB: this is not straightforward!
  # cases_eid <- cases_eid %>%
  #   mutate(topic = map_aou2ukb_topic(topic, order_aou = order_aou)) 
  cases_eid <- cases_eid %>%
    mutate(topic = aou2ukb_mapping[topic]) 
  
  colors_topic <- c(cbPalette, paletteMartin[c(8,11)])
  names(colors_topic) <- as.character(1:10)
  # plot histogram for each group
  plt <- ggplot(cases_eid,aes(x = age_diag, fill = factor(topic, levels = as.character(order_ds)), 
                              colour = factor(topic, levels = as.character(order_ds)))) +
    geom_histogram(alpha = 0.5, position = "stack", binwidth = 1) + 
    scale_fill_manual(values = colors_topic) + 
    scale_color_manual(values =colors_topic) +
    theme_bw(base_size = 15) + 
    theme(legend.position = "None") + 
    labs(x="Age (years)", y="incidence count", title=paste0("Disease: ", aou_subtype_diseases_list$V4[idx])) 
  ggsave(paste0("Revision_NG/Subtype_histgram/","aou_subtypes_age_distribution_",ds_id,".png"), plt, width = 6, height = 4)
}

##########################################################################################
# use comorbidity correlation to compare disease subtypes between AOU and UKB
##########################################################################################
# parameter from previous section
best_AOU_model_path <- "Results_NG_revision/full_para_ATM-SNOMED-349_topic_num-13_dof_5_best_model.Rdata"
load(best_AOU_model_path)
para_aou <- para
AOU_loading_mean <- apply( para_aou$pi_beta_basis[30:80, , ], c(2,3), mean)
order_topic_diseases <- sapply(1:dim(AOU_loading_mean)[1] , function(i) which(para_aou$list_above500occu$diag_icd10[i] ==  ATM::UKB_349_disease$diag_icd10))
UKB_loading_mean <- ATM::UKB_HES_10topics[30:80,,]
UKB_loading_mean <- apply( UKB_loading_mean[, order_topic_diseases, ], c(2,3), mean)
correlation_UKB_AOU <- cor(UKB_loading_mean, AOU_loading_mean)

subtypes_list <- read.table( "subtype_disease_topic_list.txt", sep="\t", header = F)
aou_subtype_diseases_list <- subtypes_list %>%
  group_by(V1) %>%
  slice(1) %>%
  filter(V1 %in% para_aou$list_above500occu$diag_icd10)

variance_ukb_aou <- pmax(correlation_UKB_AOU,0)^2
aou2ubk_topic_mapping <- sweep(variance_ukb_aou, 2, colSums(variance_ukb_aou), FUN="/")
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
para_ukb <- para
loadings_per_ds_ukb <- matrix(NA, nrow = dim(aou_subtype_diseases_list)[1], ncol = para_ukb$K)
loadings_per_ds_aou <- matrix(NA, nrow = dim(aou_subtype_diseases_list)[1], ncol = para_ukb$K)
for(idx in 1:dim(aou_subtype_diseases_list)[1]){
  ds_id <- aou_subtype_diseases_list$V1[idx]
  j <- match(as.numeric(ds_id), para_aou$list_above500occu$diag_icd10)
  loadings_per_ds_aou[idx, ] <- aou2ubk_topic_mapping %*% colMeans(para_aou$unlist_zn[para_aou$ds_list[[j]]$id,])
  
  j <- match(as.numeric(ds_id), para_ukb$list_above500occu$diag_icd10)
  loadings_per_ds_ukb[idx, ] <- colMeans(para_ukb$unlist_zn[para_ukb$ds_list[[j]]$id,])
  print(cor(loadings_per_ds_aou[idx, ], loadings_per_ds_ukb[idx, ]))
}
cor_matrix <- cor(t(loadings_per_ds_aou), t(loadings_per_ds_ukb) )
background_disease_cor <- data.frame(correlation = cor_matrix[col(cor_matrix) != row(cor_matrix)]) %>%
  mutate(type = "controls")
same_disease_css_ukb_aou <- data.frame(correlation = diag(cor_matrix)) %>%
  mutate(type = "disease_pair") # comorbidity sharing score between ukb and aou results

t.test(background_disease_cor$correlation, same_disease_css_ukb_aou$correlation, alternative = "two.sided", var.equal = FALSE)
mean_grey <- mean(background_disease_cor$correlation)
mean_red <- mean(same_disease_css_ukb_aou$correlation)
disease_pair_correlation <- bind_rows(background_disease_cor, same_disease_css_ukb_aou)
plt_cor <- ggplot(disease_pair_correlation) + 
  geom_density(aes(x = correlation, fill = type), alpha = 0.5) + 
  scale_fill_manual(values = c("disease_pair" = red, "controls" = grey)) + 
  theme_bw(base_size = 15) +
  theme(legend.position = "None") + 
  labs(x="Correlations between UKB and AOU topic assignments") +
  geom_vline(xintercept=mean_grey, linetype="dashed", color = grey, size = 1) + 
  geom_vline(xintercept=mean_red, linetype="dashed", color = red, size = 1) 
ggsave(paste0("Revision_NG/subtype_correlation_aou_ukb.png"), plt_cor, width = 6, height = 5)

# also plot the heatmap
colnames(cor_matrix) <- aou_subtype_diseases_list$V4
rownames(cor_matrix) <- aou_subtype_diseases_list$V4
longData<-reshape2::melt(cor_matrix) 
plt_subtpcor <- ggplot() + geom_tile(data = longData, aes(x = Var2, y = Var1,  fill = value, width = 0.9, height = 0.9)) + 
  geom_tile(data = filter(longData, Var2 == Var1), aes(x = Var2, y = Var1), color= "black",size = 1, alpha = 0) +  
  scale_fill_gradient2(low=blue, mid="white", high=red,  midpoint=0)+
  labs(x="", y="", title="") +
  scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.line=element_blank(),axis.text.x=element_text(angle = 90, vjust = 0.7, hjust=1, size=13),
        axis.text.y=element_text(size=13), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
ggsave(paste0("Revision_NG/heatmap_subtype_correlation_aou_ukb.png"), plt_subtpcor, width = 12, height = 12)

longData<-reshape2::melt(matrix(seq(from = min(cor_matrix), to = max(cor_matrix), length.out = 100),ncol = 1)) 
plt <- ggplot() + 
  geom_tile(data = longData, aes(x = Var2, y = Var1,fill = value,width = 1)) + 
  scale_fill_gradient2(low=blue, mid="white", high=red,  midpoint=0)+
  labs(x="", y="", title="") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
ggsave(paste0("Revision_NG/subtpcorrelation_legends.png"), plt, width = 1, height = 5)

# running the comorbidity similarity across 233 diseases
phe_phecode <- read.csv("info_phenotype_phecode.csv") 

AOU_ds_list <- para_aou$list_above500occu %>%
  mutate(V1 = as.numeric(diag_icd10)) %>%
  left_join(phe_phecode, by = c("V1" = "phecode")) %>%
  rename(V4 = phenotype) %>%
  arrange(V1)
  
loadings_per_ds_ukb <- matrix(NA, nrow = dim(AOU_ds_list)[1], ncol = para_ukb$K)
loadings_per_ds_aou <- matrix(NA, nrow = dim(AOU_ds_list)[1], ncol = para_ukb$K)
for(idx in 1:dim(AOU_ds_list)[1]){
  ds_id <- AOU_ds_list$V1[idx]
  j <- match(as.numeric(ds_id), para_aou$list_above500occu$diag_icd10)
  loadings_per_ds_aou[idx, ] <- aou2ubk_topic_mapping %*% colMeans(para_aou$unlist_zn[para_aou$ds_list[[j]]$id,])
  
  j <- match(as.numeric(ds_id), para_ukb$list_above500occu$diag_icd10)
  loadings_per_ds_ukb[idx, ] <- colMeans(para_ukb$unlist_zn[para_ukb$ds_list[[j]]$id,])
  print(cor(loadings_per_ds_aou[idx, ], loadings_per_ds_ukb[idx, ]))
}
cor_matrix <- cor(t(loadings_per_ds_aou), t(loadings_per_ds_ukb) )
background_disease_cor <- data.frame(correlation = cor_matrix[col(cor_matrix) != row(cor_matrix)]) %>%
  mutate(type = "controls")
same_disease_css_ukb_aou <- data.frame(correlation = diag(cor_matrix)) %>%
  mutate(type = "disease_pair") # comorbidity sharing score between ukb and aou results

t.test(background_disease_cor$correlation, same_disease_css_ukb_aou$correlation, alternative = "two.sided", var.equal = FALSE)
mean_grey <- mean(background_disease_cor$correlation)
mean_red <- mean(same_disease_css_ukb_aou$correlation)
disease_pair_correlation <- bind_rows(background_disease_cor, same_disease_css_ukb_aou)
plt_cor <- ggplot(disease_pair_correlation) + 
  geom_density(aes(x = correlation, fill = type), alpha = 0.5) + 
  scale_fill_manual(values = c("disease_pair" = red, "controls" = grey)) + 
  theme_bw(base_size = 15) +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = BLACK)
  ) + 
  labs(x="Correlations between UKB and AOU topic assignments") +
  geom_vline(xintercept=mean_grey, linetype="dashed", color = grey, size = 1) + 
  geom_vline(xintercept=mean_red, linetype="dashed", color = red, size = 1) 
ggsave(paste0("Revision_NG/233ds_subtype_correlation_aou_ukb.png"), plt_cor, width = 6, height = 5)

colnames(cor_matrix) <- AOU_ds_list$V4
rownames(cor_matrix) <- AOU_ds_list$V4
write.csv(cor_matrix, file = paste0("Revision_NG/233ds_subtype_correlation_aou_ukb.csv")) # row is AOU, colums are UKB

########################################################
# make table 2 subtype AOU vs UKB
########################################################
# use cor_matrix from above
subtype_cor <- reshape2::melt(cor_matrix) %>%
  filter(Var1 == Var2)

load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
subtypes_list <- read.table( "subtype_disease_topic_list.txt", sep="\t", header = F)
subtype_diseases_list <- subtypes_list %>%
  group_by(V1) %>%
  slice(1)
library(colorBlindness)
displayAvailablePalette(color="white")
subtype_info <- list()
for(idx in 1:dim(subtype_diseases_list)[1]){
  ds_id <- subtype_diseases_list$V1[idx]
  j <- match(as.numeric(ds_id), para$list_above500occu$diag_icd10)
  cases_eid <- list()
  for(topic_id in 1:para$K){
    topic_specific_id <- which(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max) == topic_id)
    topic_specific_set <- para$unlist_Ds_id[para$ds_list[[j]]$id,][topic_specific_id,]
    
    cases_eid[[topic_id]] <- topic_specific_set %>%
      mutate(topic = topic_id) %>%
      mutate(Ds_id = para$list_above500occu[Ds_id, 1])
  }
  cases_eid <- bind_rows(cases_eid)
  
  subtype_info[[idx]] <- cases_eid %>% 
    group_by(topic) %>%
    summarise(age = mean(age_diag), subtype_count = n()) %>%
    mutate(subtype_percentage = subtype_count/sum(subtype_count), disease = ds_id)
  
}
subtype_info <- bind_rows(subtype_info)
subtype_count <- subtype_info %>% 
  mutate(topic = factor(topic, levels = order_ds)) %>%
  left_join(select(ds.system, diag_icd10, phenotype), by=c("disease" = "diag_icd10")) %>%
  arrange(topic) %>%
  pivot_wider(id_cols = c("disease", "phenotype"), names_from = "topic", values_from = "subtype_count") %>%
  replace(is.na(.), 0) %>%
  arrange(disease)
topic_ordered_namte <- c("NRI", "CER","SRD", "CVD", "UGI", "LGI", "FGND",  "MGND", "MDS", "ARP")
names(subtype_count) <- c("Phecode", "Description",topic_ordered_namte)

subtype_count %>% 
  pivot_longer(cols = NRI:ARP, names_to = "topic", values_to = "counts") %>%
  filter(counts >= 500) %>%
  group_by(Phecode) %>% 
  mutate(Subtypes = paste(topic, collapse = ", ")) %>%
  slice(1) %>%
  select( -topic, - counts) %>%
  left_join(select(subtype_cor, - Var2), by = c("Description" = "Var1")) %>%
  rename(AOU_UKB_cor=value) %>% 
  write.csv("~/Desktop/comorbidity/Multi-morbidity_biobank/Revision_NG/Table2_AOU_UKB_Subtype_correlation.csv", row.names = F)


#############################################################################################
# compute additional validation on PRS association with topic subtypes in the controls.
#############################################################################################
load(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/bolt_topic_regression.Rdata" ))
case_pvalues <- subtype_regression[[2]]
case_coefficients <- subtype_regression[[3]]
case_se_coefficients <- subtype_regression[[4]]

# perform the analysis on the controls
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
patient_loadings <- sweep((para$alpha_z - 1), 1, rowSums(para$alpha_z -1), FUN="/") 

rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv") 
topic_individuals <- rec_data %>%
  group_by(eid) %>%
  slice(1)

ds_list <- read.table("BOLT_LMM_disease_list.txt")
# playing with the code
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("V1" = "phecode"))
coefficients <- matrix(nrow = length(ds_list$V1), ncol = 10)
se_coef <- matrix(nrow = length(ds_list$V1), ncol = 10)
p_values <- matrix(nrow = length(ds_list$V1), ncol = 10)
for(idx in 1:length(ds_list$V1)){
  ds_id <- as.character(ds_list$V1[idx])
  print(ds_id)
  
  disease_topic_data <- rec_data %>%
    filter(diag_icd10 == ds_id)
  
  # topic_id <- 7
  for(topic_id in 1:10){
    # df_umap_patient <- data.frame(eid = para$eid, loadings = patient_loadings[,topic_id])  %>%
    #   mutate(loadings = loadings/sd(loadings))
    # could directly use topic weights as unit
    df_umap_patient <- data.frame(eid = para$eid, loadings = patient_loadings[,topic_id])  
    
    prediction_set <- list()
    for(cv_id in 1:5){
      cv_profile <- read.table(paste0("bolt_prs_result/", ds_id, "_csvld_", cv_id, ".bolt.prs.profile"), header = T)
      training_set <- read.table(paste0("bolt_prs_result/", ds_id, "_csvld_", cv_id, ".pheno"), header = T) 
      prediction_set[[cv_id]] <- read.table(paste0("bolt_prs_result/", ds_id, "keep.txt"), header = F) %>%
        anti_join(training_set, by = c("V1" = "FID")) %>%
        semi_join(topic_individuals,by = c("V1" = "eid")) %>%
        mutate(phenotype = if_else(V1 %in% disease_topic_data$eid, 1, 0)) %>%
        left_join(select(cv_profile, FID, SCORE), by= c("V1" = "FID"))
    }
    prediction_set <- bind_rows(prediction_set) %>%
      mutate(SCORE = SCORE/sd(SCORE)) %>%
      filter(phenotype == 0)
    
    regression_df <- prediction_set %>% 
      left_join(df_umap_patient, by = c("V1" = "eid" )) 
    regression_md <- lm(SCORE ~ loadings, data = regression_df)
    coefficients[idx, topic_id] <- summary(regression_md)$coefficients[2,1]
    se_coef[idx, topic_id] <- summary(regression_md)$coefficients[2,2]
    p_values[idx, topic_id] <- summary(regression_md)$coefficients[2,4]
  }
  
}

pasted_rslt <- matrix(mapply(function(x,y) paste0(as.character(x), " (P = ", as.character(y), ")"), 
                             round(coefficients,digits = 4), round(p_values, digits = 4)), dim(ds_list)[1])
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
colnames(pasted_rslt) <- topic_name
df_subtype_regression <- data.frame(Phecode = ds_list$V1, Phenotype = ds_list$phenotype, pasted_rslt) 
write.csv(df_subtype_regression, paste0("Revision_NG/controls_bolt_topic_regression.csv" ))
controls_subtype_regression <- list(ds_list, p_values, coefficients, se_coef)
save(controls_subtype_regression, file = paste0("Revision_NG/controls_bolt_topic_regression.Rdata" ))

load(paste0("Revision_NG/controls_bolt_topic_regression.Rdata" ))

plot(x = -log10(p_values), y  = -log10(case_pvalues) )

plot(x = coefficients, y  = case_coefficients )
sum(p.adjust(p_values) < 0.05)
sum(p.adjust(case_pvalues) < 0.05)

ids_case_significant <- which(p.adjust(case_pvalues) < 0.05)
# plot all the coefficient estimation for those significant hits
data_compare_PRS_subtype_controls <- data.frame(case_coef =case_coefficients[ids_case_significant],
                                                cose_se = case_se_coefficients[ids_case_significant], 
                                                control_coef = coefficients[ids_case_significant],
                                                control_se = se_coef[ids_case_significant])
plt_compare_PRS_subtype_controls <- ggplot(data_compare_PRS_subtype_controls) +
  geom_pointrange(aes(x = case_coef, y = control_coef, 
                      ymin =control_coef -1.96*control_se,
                      ymax = control_coef +1.96*control_se) , color = red) + 
  geom_pointrange(aes(x = case_coef, y = control_coef, 
                      xmin =case_coef -1.96*cose_se,
                      xmax = case_coef +1.96*cose_se), color = red) + 
  geom_abline(slope = 1, color = grey, linetype = "dashed") +
  labs(x = "Excess PRS in cases", y = "Excess PRS in controls") + 
  lims(x = c(-0.45, 0.45), y = c(-0.45, 0.45)) + 
  geom_smooth(data = data_compare_PRS_subtype_controls, mapping=aes(x=case_coef, y=control_coef), color = red, method='lm') +
  theme_bw(base_size = 20) +
  theme(legend.position = "None",panel.background=element_blank()) 
md <- lm(control_coef~case_coef, data_compare_PRS_subtype_controls)
summary(md)
ggsave(paste0("Revision_NG/compare_PRS_subtype_controls.png" ), plt_compare_PRS_subtype_controls, width = 5, height = 5)


#############################################
# for Yidong: histogram of multiple counting
#############################################
all_data_ukb <- read.csv(file = "../UKBB_interim_files/hes_records_application130418.csv")
HES_data <- read.table(paste0("~/Desktop/comorbidity/UKBB_interim_files/cambridge_application/cambridge_local_rec_data.txt"), sep = "\t", header = T)
Yidong_codes <- read.table("UKB_data_yidong/top100_code.txt", header = T, sep = "\t")
HES_data_primary <- HES_data %>%
  filter(code_type == "ICD-10", cause_type == "primary")

HES_data_secondary <- HES_data %>%
  filter(code_type == "ICD-10", cause_type == "secondary")
# goal, plot disease per-person for primary and secondary disease
hist_primary <- HES_data_primary %>%
  mutate(short_code = substr(code, 1, 3) ) %>%
  filter(short_code %in% Yidong_codes$code) %>%
  group_by(eid, short_code) %>%
  summarise(counts=n()) %>%
  mutate(cause_type = "primary")

hist_secondary <- HES_data_secondary %>%
  mutate(short_code = substr(code, 1, 3) ) %>%
  filter(short_code %in% Yidong_codes$code) %>%
  group_by(eid, short_code) %>%
  summarise(counts=n()) %>%
  mutate(cause_type = "secondary")

hist_per_person_records <- bind_rows(hist_primary, hist_secondary)

plt_yidong_primary_secondary_compare <- ggplot(data = hist_per_person_records, aes(x = counts)) +
  geom_histogram(aes(y=..density.., fill = cause_type), alpha=0.5, position = "dodge", binwidth = 1) + 
  scale_fill_manual(name = "Cause",values = c("primary" = red, "secondary"=blue)) + 
  scale_x_continuous(breaks = seq(1, 10, 1), limits = c(0,10)) + 
  labs(title = "Number of recurrent diseae code") +
  theme_bw(base_size = 15)

ggsave(filename = "UKB_data_yidong/counts_primary_secondary.png", plt_yidong_primary_secondary_compare, width = 5, height = 5)

################################################
# replication of main figure 3; figure 4B
################################################
# need to fix the ordering!
library(reshape2)
best_AOU_model_path <- "Results_NG_revision/full_para_ATM-SNOMED-349_topic_num-13_dof_5_best_model.Rdata"
load(best_AOU_model_path)
loadings_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[para$ds_list[[j]]$id,])) %>% t
# do it separately for younger and older
younger60 <- lapply(para$ds_list, function(x) filter(x, age_diag <=60))
loadings_young_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[younger60[[j]]$id,])) %>% t
older60 <- lapply(para$ds_list, function(x) filter(x, age_diag >60))
loadings_old_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[older60[[j]]$id,,drop=FALSE])) %>% t
loadings_old_per_ds[is.na(loadings_old_per_ds)] <- 0

phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds.system <- para$list_above500occu %>%
  mutate(diag_icd10 = as.numeric(diag_icd10)) %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") ) 
# combine the last four systems to other:congenital anomalies, symptoms, injuries & poisonings, <NA>  
ds.system <- ds.system %>%
  mutate(exclude_name = if_else(exclude_name %in% c("congenital anomalies", "symptoms", "injuries & poisonings","", NA), "others", exclude_name))
# fix the ordering: Phecode should be ordered by numeric number!
ds_rev_aou_order <- sort(as.numeric(ds.system$diag_icd10), decreasing = T)
system_aou_order <- ds.system %>%
  mutate(diag_icd10 = as.numeric(diag_icd10)) %>%
  arrange(diag_icd10) %>%
  pull(exclude_name)

ds.system <- ds.system %>%
  mutate(exclude_name = factor(exclude_name, levels = unique(system_aou_order)))

systems_associated <- list()
disease_number <- dim(loadings_per_ds)[1]
topic_number <- dim(loadings_per_ds)[2]
for(i in 1:topic_number){
  ds.system$loadings <- loadings_per_ds[, i]
  systems_associated[[i]] <- ds.system %>% 
    group_by(exclude_name) %>%
    summarise(mean_loading = mean(loadings), sum_loading = sum(loadings)) %>%
    filter(mean_loading > .5 | sum_loading > 10) %>%
    pull(exclude_name)
}
order_ds <- order(sapply(systems_associated, function(x) x[1]))
systems_associated <- systems_associated[order_ds]
ds.system$loadings <- loadings_per_ds[, order_ds]

loadings_age_diff <- matrix(NA, nrow = disease_number, ncol = 2*topic_number)
for(i in 1:length(order_ds)){
  loadings_age_diff[,2*i -1] <- loadings_young_per_ds[,order_ds[i]]
  loadings_age_diff[,2*i] <- loadings_old_per_ds[,order_ds[i]]
}

longData<-melt(loadings_age_diff) %>%
  mutate(Var1 =ds.system$diag_icd10[Var1])

ds.groups <- ds.system %>% 
  select(diag_icd10, exclude_name) %>%
  mutate(Var1 = diag_icd10, group_range = (as.integer(exclude_name) %% 2)) %>%
  mutate(group_range = ifelse(is.na(group_range), 1, group_range)) %>%
  mutate(group_range = as.factor(group_range))
longData <- longData %>%
  left_join(ds.groups, by = c("Var1")) %>%
  mutate(Var1 = factor(Var1, levels = ds_rev_aou_order))
plt <- ggplot() + 
  geom_tile(data = longData, aes(x = Var2, y = Var1, fill=group_range, alpha = value,width = 0.9)) + 
  scale_alpha_continuous(range = c(0,1)) +
  scale_fill_manual(values = c("0" = blue, "1" = red))+
  labs(x="", y="", title="") +
  scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
# add the vertical separation line
df_segments <- data.frame(x = 0:topic_number * 2 + 0.5, xend =  0:topic_number * 2+ 0.5, y = rep(0,(topic_number + 1)), yend = rep(disease_number,(topic_number + 1))) 
plt <- plt + geom_segment(data=df_segments, aes(x,y,xend=xend, yend=yend), size=.5, alpha = 0.3, inherit.aes=F)
ggsave("Revision_NG/AOU_topic_disease_highlight.png",plt, width = 6, height = 10)

categories <- ds.system$exclude_name
longData<-melt(as.matrix(categories)) %>%
  mutate(Var1 = factor(ds.system$diag_icd10[Var1], levels = ds_rev_aou_order))
cols <- setNames(rep(c(red, blue), 9)[1:14], unique(system_aou_order))
plt_pallete <- ggplot() + 
  geom_tile(data = longData, aes(x = Var2, y = Var1, fill=value), alpha = 0.5) + 
  scale_fill_manual(values = cols) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
ggsave("Revision_NG/AOU_topic_pallete.png",plt_pallete, width = 1, height = 10)

# figure 4 B-C
best_AOU_model_path <- "Results_NG_revision/full_para_ATM-SNOMED-349_topic_num-13_dof_5_best_model.Rdata"
load(best_AOU_model_path)
loadings_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[para$ds_list[[j]]$id,])) %>% t
top_1_value <- sapply(1:dim(loadings_per_ds)[1], function(x) sort(loadings_per_ds[x,], decreasing = T)) %>% t
df_box <- melt(top_1_value) %>%
  rename(topic_order = Var2, weights = value, record_number = Var1) %>%
  mutate(topic_order = factor(topic_order, levels = 1:13))


plt_all <- ggplot(data=df_box) +
  geom_boxplot(mapping=aes(x=topic_order, y=weights, fill = red), color = red,   width=0.5, alpha = 0.6, outlier.shape = NA) +
  # geom_jitter(mapping=aes(x=topic_order, y=weights),width=0.2, size = 0.1, alpha=0.4) + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = "Disease topic loadings rank", y = "Disease topic loadings") + 
  theme_bw(base_size = 20) +
  theme(legend.position = "None",panel.background=element_blank()) + 
  geom_hline(aes(yintercept = 1), color = "red", linetype = "dashed")
ggsave("Revision_NG/AOU_fig3_disease_sparsity.png", plt_all, width = 5, height = 5)
print(paste0("proportion of max topic value > 0.95: ", mean(top_1_value[,1] > 0.95) ) )


patient_loadings <- sweep((para$alpha_z-1), 1, rowSums(para$alpha_z - 1), FUN="/")
top_1_value <- sapply(1:dim(patient_loadings)[1], function(x) sort(patient_loadings[x,], decreasing = T)) %>% t
df_box <- melt(top_1_value) %>%
  rename(topic_order = Var2, weights = value, record_number = Var1) %>%
  mutate(topic_order = factor(topic_order, levels = 1:13))

plt <- ggplot(data=df_box) +
  geom_boxplot(mapping=aes(x=topic_order, y=weights, fill = red), color = red,  width=0.5, alpha = 0.6, outlier.shape = NA) +
  # geom_jitter(mapping=aes(x=topic_order, y=weights),width=0.2, size = 0.1, alpha=0.4) + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = "Patient topic weight rank", y = "Patient topic weight") + 
  theme_bw(base_size = 20) +
  theme(legend.position = "None",panel.background=element_blank()) + 
  geom_hline(aes(yintercept = 1), color = "red", linetype = "dashed")
ggsave("Revision_NG/AOU_fig3_patient_sparsity.png", plt, width = 5, height = 5)



##################################################
# Supp Figuure 22* perform simulation using realistic MAF information
##################################################
source("simulation_functions.R")

# extract the MAF from the GWAS SNPs
SNP_files <- list.files("Results_NG_revision/GWAS_SNP_list", pattern="", full.names=TRUE)
SNP_list <- list()
for(f_id in 1:length(SNP_files)){
  try(SNP_list[[f_id]] <- read.table(SNP_files[[f_id]]))
}
SNP_list <- bind_rows(SNP_list)
SNP_list <- unique(SNP_list$V1)
MAF_ukb <- read.table("~/Desktop/PublicRepository/LIBRARY/Jiang_2022_ATM/GxTopic_simulation/genotypes_maf.frq", header = T) %>%
  filter(SNP %in% SNP_list)
hist(MAF_ukb$MAF)
write.csv(MAF_ukb, file = "~/Desktop/PROJECTS/health_history_prediction/EHR_prediction/NG_revision/maf_ukb_GWAS.csv")

# plot the genetic simulation figures
# 1. SNP -> T -> D
v_id <- 1
pw_thr <- 0.05 # threshold of power 
df_1 <- data.frame( causal_v2t = as.numeric(), gxT_mean= as.numeric(), gxT_se = as.numeric(),
                    gT2_mean= as.numeric(), gT2_se = as.numeric())
for(s_id in 1:10){
  print(s_id)
  load(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/Results_NG_revision/genetic_simulation/", 
              "realmaf_gen_simulation_v",v_id,"s", s_id,".RData" ))
  p_gxT <- NULL
  p_gT2 <- NULL
  p_gxD <- NULL
  for(rep_id in 1:length(genetic_simulation)){
    genetics_population <- genetic_simulation[[rep_id]][[1]]
    patient_loadings <- genetic_simulation[[rep_id]][[2]] 
    disease_matrix <- genetic_simulation[[rep_id]][[3]] 
    sim_para <- genetic_simulation[[rep_id]][[4]]
    # find which topic is genetically modulated
    tp_genetic <- colMeans(patient_loadings[which(disease_matrix[,1] == 1),]) %>% which.max
    loadings <- patient_loadings[,tp_genetic]
    for(ds_id in 1:((sim_para$disease_number-2)/2)){
      # for(ds_id in 1:1){
      p_gxT <- c(p_gxT, test_gxT(loadings, disease_matrix[,ds_id], genetics_population[,21:(20+sim_para$cont_v2t)]))
      p_gT2 <- c(p_gT2, test_gxT_T2(loadings, disease_matrix[,ds_id], genetics_population[,21:(20+sim_para$cont_v2t)]))
      # p_gxD <- c(p_gxD,test_gxD_simple(loadings, disease_matrix[,ds_id], genetics_population[,21:(20+sim_para$cont_v2t)]))
    }
  }
  gxT_mean <- mean(p_gxT < pw_thr,na.rm = T)
  gxT_se <- sqrt(gxT_mean*(1-gxT_mean)/length(p_gxT))
  gT2_mean <- mean(p_gT2 < pw_thr,na.rm = T)
  gT2_se <- sqrt(gT2_mean*(1-gT2_mean)/length(p_gT2))
  # gT2_mean <- mean(p_gxD < pw_thr,na.rm = T)
  # gT2_se <- sqrt(gT2_mean*(1-gxT_mean)/length(p_gxD))
  df_1 <- df_1 %>%
    add_row(causal_v2t = sim_para$cont_v2t,  gxT_mean= gxT_mean, gxT_se = gxT_se,
            gT2_mean= gT2_mean, gT2_se = gT2_se)
}
plt_g2t2d <- ggplot(df_1,aes(x = causal_v2t)) +
  geom_line(aes(y = gxT_mean), color = red, linetype = "dashed") + 
  geom_line(aes(y =  gT2_mean), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = gxT_mean, ymin=gxT_mean-1.96*gxT_se, ymax = gxT_mean + 1.96*gxT_se, color = "G x Topic model"), size = .3) + 
  geom_pointrange(aes(y = gT2_mean, ymin=gT2_mean-1.96*gT2_se, ymax = gT2_mean + 1.96*gT2_se, color = "G x Disease model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "G x Disease model" = green)) + 
  labs(x ="Number of causal SNP to topic", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
ggsave(paste0("Revision_NG/","realmaf_interaction_g2t2d.png"), plt_g2t2d, width = 4, height = 4)


# 2. causal disease SNP->D->T
v_id <- 2
pw_thr <- 0.05 # threshold of power 
df_2 <- data.frame( disease2topic = as.numeric(), gxT_mean= as.numeric(), gxT_se = as.numeric(),
                    gT2_mean= as.numeric(), gT2_se = as.numeric())
for(s_id in 1:10){
  print(s_id)
  load(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/Results_NG_revision/genetic_simulation/", 
              "realmaf_gen_simulation_v",v_id,"s", s_id,".RData" ))
  p_gxT <- NULL
  p_gT2 <- NULL
  for(rep_id in 1:length(genetic_simulation)){
    genetics_population <- genetic_simulation[[rep_id]][[1]]
    patient_loadings <- genetic_simulation[[rep_id]][[2]] 
    disease_matrix <- genetic_simulation[[rep_id]][[3]] 
    sim_para <- genetic_simulation[[rep_id]][[4]]
    # find which topic is genetically modulated
    tp_genetic <- colMeans(patient_loadings[which(disease_matrix[,1] == 1),]) %>% which.max
    loadings <- patient_loadings[,tp_genetic]
    ds_id <- sim_para$disease_number + 1
    p_gxT <- c(p_gxT, test_gxT(loadings, disease_matrix[,ds_id], genetics_population[,1:20]))
    p_gT2 <- c(p_gT2,test_gxT_T2(loadings, disease_matrix[,ds_id], genetics_population[,1:20]))
  }
  gxT_mean <- mean(p_gxT < pw_thr,na.rm = T)
  gxT_se <- sqrt(gxT_mean*(1-gxT_mean)/length(p_gxT))
  gT2_mean <- mean(p_gT2 < pw_thr,na.rm = T)
  gT2_se <- sqrt(gT2_mean*(1-gxT_mean)/length(p_gT2))
  df_2 <- df_2 %>%
    add_row(disease2topic = sim_para$disease2topic,  gxT_mean= gxT_mean, gxT_se = gxT_se,
            gT2_mean= gT2_mean, gT2_se = gT2_se)
}
plt_g2d2t <- ggplot(df_2,aes(x = disease2topic)) +
  geom_line(aes(y = gxT_mean), color = red, linetype = "dashed") + 
  geom_line(aes(y =  gT2_mean), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = gxT_mean, ymin=gxT_mean-1.96*gxT_se, ymax = gxT_mean + 1.96*gxT_se, color = "G x Topic model"), size = .3) + 
  geom_pointrange(aes(y = gT2_mean, ymin=gT2_mean-1.96*gT2_se, ymax = gT2_mean + 1.96*gT2_se, color = "G x Disease model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "G x Disease model" = green)) + 
  labs(x ="Causal disease effect on topic", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
ggsave(paste0("Revision_NG/","realmaf_interaction_g2d2t.png"), plt_g2d2t, width = 4, height = 4)

######################################################
# compute per-sd change in disease to topic

# 3. SNP*T->D
v_id <- 3
pw_thr <- 0.05 # threshold of power 
df_3 <- data.frame( itr_effect = as.numeric(), gxT_mean= as.numeric(), gxT_se = as.numeric(),
                    gT2_mean= as.numeric(), gT2_se = as.numeric())
for(s_id in 1:10){
  print(s_id)
  load(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/Results_NG_revision/genetic_simulation/", 
              "realmaf_gen_simulation_v",v_id,"s", s_id,".RData" ))
  p_gxT <- NULL
  p_gT2 <- NULL
  for(rep_id in 1:length(genetic_simulation)){
    genetics_population <- genetic_simulation[[rep_id]][[1]]
    patient_loadings <- genetic_simulation[[rep_id]][[2]] 
    disease_matrix <- genetic_simulation[[rep_id]][[3]] 
    sim_para <- genetic_simulation[[rep_id]][[4]]
    # find which topic is genetically modulated
    tp_genetic <- colMeans(patient_loadings[which(disease_matrix[,1] == 1),]) %>% which.max
    loadings <- patient_loadings[,tp_genetic]
    ds_id <- sim_para$disease_number + 2
    p_gxT <- c(p_gxT, test_gxT(loadings, disease_matrix[,ds_id], genetics_population[,41:60]))
    p_gT2 <- c(p_gT2,test_gxT_T2(loadings, disease_matrix[,ds_id], genetics_population[,41:60]))
    ############################
    # important observation
    # when simulating multiple interaction and testing for one -- the power is lost!
    ############################
    # summary(glm(disease_matrix[,ds_id] ~ genetics_population[,41:60]*loadings, family = binomial))
  }
  gxT_mean <- mean(p_gxT < pw_thr,na.rm = T)
  gxT_se <- sqrt(gxT_mean*(1-gxT_mean)/length(p_gxT))
  gT2_mean <- mean(p_gT2 < pw_thr,na.rm = T)
  gT2_se <- sqrt(gT2_mean*(1-gxT_mean)/length(p_gT2))
  df_3 <- df_3 %>%
    add_row(itr_effect = sim_para$itr_effect, gxT_mean= gxT_mean, gxT_se = gxT_se,
            gT2_mean= gT2_mean, gT2_se = gT2_se)
}
plt_gxt <- ggplot(df_3,aes(x = itr_effect)) +
  geom_line(aes(y = gxT_mean), color = red, linetype = "dashed") + 
  geom_line(aes(y =  gT2_mean), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = gxT_mean, ymin=gxT_mean-1.96*gxT_se, ymax = gxT_mean + 1.96*gxT_se, color = "G x Topic model"), size = .3) + 
  geom_pointrange(aes(y = gT2_mean, ymin=gT2_mean-1.96*gT2_se, ymax = gT2_mean + 1.96*gT2_se, color = "G x Disease model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "G x Disease model" = green)) + 
  labs(x ="Effect size of GxTopic interaction", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
ggsave(paste0("Revision_NG/","realmaf_interaction_gxt.png"), plt_gxt, width = 4, height = 4)

# 4. SNP->D; SNP->T->D
v_id <- 4
pw_thr <- 0.05 # threshold of power 
df_4 <- data.frame( topic2disease = as.numeric(), gxT_mean= as.numeric(), gxT_se = as.numeric(),
                    gT2_mean= as.numeric(), gT2_se = as.numeric())
for(s_id in 1:10){
  print(s_id)
  load(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/Results_NG_revision/genetic_simulation/", 
              "realmaf_gen_simulation_v",v_id,"s", s_id,".RData" ))
  p_gxT <- NULL
  p_gT2 <- NULL
  for(rep_id in 1:length(genetic_simulation)){
    genetics_population <- genetic_simulation[[rep_id]][[1]]
    patient_loadings <- genetic_simulation[[rep_id]][[2]] 
    disease_matrix <- genetic_simulation[[rep_id]][[3]] 
    sim_para <- genetic_simulation[[rep_id]][[4]]
    # find which topic is genetically modulated
    tp_genetic <- colMeans(patient_loadings[which(disease_matrix[,1] == 1),]) %>% which.max
    loadings <- patient_loadings[,tp_genetic]
    ds_id <- sim_para$disease_number + 3
    p_gxT <- c(p_gxT, test_gxT(loadings, disease_matrix[,ds_id], genetics_population[,21:40]))
    p_gT2 <- c(p_gT2,test_gxT_T2(loadings, disease_matrix[,ds_id], genetics_population[,21:40]))
  }
  gxT_mean <- mean(p_gxT < pw_thr,na.rm = T)
  gxT_se <- sqrt(gxT_mean*(1-gxT_mean)/length(p_gxT))
  gT2_mean <- mean(p_gT2 < pw_thr,na.rm = T)
  gT2_se <- sqrt(gT2_mean*(1-gxT_mean)/length(p_gT2))
  df_4 <- df_4 %>%
    add_row(topic2disease = sim_para$topic2disease, gxT_mean= gxT_mean, gxT_se = gxT_se,
            gT2_mean= gT2_mean, gT2_se = gT2_se)
}
plt_pleiotropy <- ggplot(df_4,aes(x = topic2disease)) +
  geom_line(aes(y = gxT_mean), color = red, linetype = "dashed") + 
  geom_line(aes(y =  gT2_mean), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = gxT_mean, ymin=gxT_mean-1.96*gxT_se, ymax = gxT_mean + 1.96*gxT_se, color = "G x Topic model"), size = .3) + 
  geom_pointrange(aes(y = gT2_mean, ymin=gT2_mean-1.96*gT2_se, ymax = gT2_mean + 1.96*gT2_se, color = "G x Disease model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "G x Disease model" = green)) + 
  labs(x ="Effect of topic on disease", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
ggsave(paste0("Revision_NG/","realmaf_pleitropy_gxt.png"), plt_pleiotropy, width = 4, height = 4)

# additional cases with non-linear genetic effect
# 4. SNP + SNP^2->D; SNP->T->D
v_id <- 4
pw_thr <- 0.05 # threshold of power 
df_5 <- data.frame( topic2disease = as.numeric(), gxT_mean= as.numeric(), gxT_se = as.numeric(),
                    gT2_mean= as.numeric(), gT2_se = as.numeric())
for(s_id in 1:10){
  print(s_id)
  load(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/Results_NG_revision/genetic_simulation/", 
              "realmaf_gen_simulation_v",v_id,"s", s_id,".RData" ))
  p_gxT <- NULL
  p_gT2 <- NULL
  for(rep_id in 1:length(genetic_simulation)){
    genetics_population <- genetic_simulation[[rep_id]][[1]]
    patient_loadings <- genetic_simulation[[rep_id]][[2]] 
    disease_matrix <- genetic_simulation[[rep_id]][[3]] 
    sim_para <- genetic_simulation[[rep_id]][[4]]
    # find which topic is genetically modulated
    tp_genetic <- colMeans(patient_loadings[which(disease_matrix[,1] == 1),]) %>% which.max
    loadings <- patient_loadings[,tp_genetic]
    ds_id <- sim_para$disease_number + 4
    p_gxT <- c(p_gxT, test_gxT(loadings, disease_matrix[,ds_id], genetics_population[,21:40]))
    p_gxD <- c(p_gxD,test_gxT_T2(loadings, disease_matrix[,ds_id], genetics_population[,21:40]))
  }
  gxT_mean <- mean(p_gxT < pw_thr,na.rm = T)
  gxT_se <- sqrt(gxT_mean*(1-gxT_mean)/length(p_gxT))
  gT2_mean <- mean(p_gxD < pw_thr,na.rm = T)
  gT2_se <- sqrt(gT2_mean*(1-gxT_mean)/length(p_gxD))
  df_5 <- df_5 %>%
    add_row(topic2disease = sim_para$topic2disease, gxT_mean= gxT_mean, gxT_se = gxT_se,
            gT2_mean= gT2_mean, gT2_se = gT2_se)
}
plt_nonlinrG <- ggplot(df_5,aes(x = topic2disease)) +
  geom_line(aes(y = gxT_mean), color = red, linetype = "dashed") + 
  geom_line(aes(y =  gT2_mean), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = gxT_mean, ymin=gxT_mean-1.96*gxT_se, ymax = gxT_mean + 1.96*gxT_se, color = "G x Topic model"), size = .3) + 
  geom_pointrange(aes(y = gT2_mean, ymin=gT2_mean-1.96*gT2_se, ymax = gT2_mean + 1.96*gT2_se, color = "G x Disease model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "G x Disease model" = green)) + 
  labs(x ="Effect of topic on disease", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
ggsave(paste0("Revision_NG/","realmaf_interaction_nonlinear.png"), plt_nonlinrG, width = 4, height = 4)


# D = T + T^2 non linear case
v_id <- 4
pw_thr <- 0.05 # threshold of power 
df_6 <- data.frame( topic2disease = as.numeric(), gxT_mean= as.numeric(), gxT_se = as.numeric(),
                    gT2_mean= as.numeric(), gT2_se = as.numeric())
for(s_id in 1:10){
  print(s_id)
  load(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/Results_NG_revision/genetic_simulation/", 
              "realmaf_gen_simulation_v",v_id,"s", s_id,".RData" ))
  p_gxT <- NULL
  p_gT2 <- NULL
  for(rep_id in 1:length(genetic_simulation)){
    genetics_population <- genetic_simulation[[rep_id]][[1]]
    patient_loadings <- genetic_simulation[[rep_id]][[2]] 
    disease_matrix <- genetic_simulation[[rep_id]][[3]] 
    sim_para <- genetic_simulation[[rep_id]][[4]]
    # find which topic is genetically modulated
    tp_genetic <- colMeans(patient_loadings[which(disease_matrix[,1] == 1),]) %>% which.max
    loadings <- patient_loadings[,tp_genetic]
    ds_id <- sim_para$disease_number + 5
    p_gxT <- c(p_gxT, test_gxT(loadings, disease_matrix[,ds_id], genetics_population[,21:40]))
    p_gT2 <- c(p_gT2,test_gxT_T2(loadings, disease_matrix[,ds_id], genetics_population[,21:40]))
  }
  gxT_mean <- mean(p_gxT < pw_thr,na.rm = T)
  gxT_se <- sqrt(gxT_mean*(1-gxT_mean)/length(p_gxT))
  gT2_mean <- mean(p_gT2 < pw_thr,na.rm = T)
  gT2_se <- sqrt(gT2_mean*(1-gxT_mean)/length(p_gT2))
  df_6 <- df_6 %>%
    add_row(topic2disease = sim_para$topic2disease, gxT_mean= gxT_mean, gxT_se = gxT_se,
            gT2_mean= gT2_mean, gT2_se = gT2_se)
}
plt_nonlinrG <- ggplot(df_6,aes(x = topic2disease)) +
  geom_line(aes(y = gxT_mean), color = red, linetype = "dashed") + 
  geom_line(aes(y =  gT2_mean), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = gxT_mean, ymin=gxT_mean-1.96*gxT_se, ymax = gxT_mean + 1.96*gxT_se, color = "G x Topic model"), size = .3) + 
  geom_pointrange(aes(y = gT2_mean, ymin=gT2_mean-1.96*gT2_se, ymax = gT2_mean + 1.96*gT2_se, color = "G x Disease model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "G x Disease model" = green)) + 
  labs(x ="Effect of topic on disease", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
ggsave(paste0("Revision_NG/","realmaf_interaction_t2_nonlinear.png"), plt_nonlinrG, width = 4, height = 4)









