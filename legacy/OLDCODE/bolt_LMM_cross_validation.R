# create cross validation samples for BOLT_LMM

# British isle ancestry
british_isle <- read.table("../UKBB_interim_files/keep_britishIsle.txt", header = F)

# this data is larger than "rec2subjectAbove1000occur_include_death_PheCode.csv" as it includes the individual who have only one records 
phecode_age <- read.csv(paste0("DiseaseAbove1000occur_PheCode.csv"), header = T)

# this is the group of people to finally generate the PGS
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv") 

phecode_age %>% filter(diag_icd10 == 250.2) %>% dim

rec_data %>% filter(diag_icd10 == 250.2) %>% dim

ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
subtypes_data <- read.table("subtype_disease_topic_list.txt", sep="\t" )
# decide what are the diseases to work on
h2g_disease_topics <- read.csv("causality_analysis/h2g_causality.csv")
h2g_disease_topics %>% 
  arrange(desc(z_score)) %>%
  filter(z_score > 4) %>%
  group_by(disease) %>%
  slice(1) %>%
  ungroup %>%
  left_join(ds_list, by = c("disease" = "diag_icd10")) %>%
  filter(disease %in% subtypes_data$V1) %>% # keep only diseases that have subtypes
  arrange(desc(z_score)) %>%
  slice(1:10) %>% # I could also include all the diseases (it is only 19)
  select(disease) %>%
  write.table("BOLT_LMM_disease_list.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)
############################################
# create cross validation file 
# file: bolt_create_cv_phenotypes.R 
############################################

args <- commandArgs(trailingOnly = TRUE)                                  
ds_id <- as.numeric(args[1])

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
  write.table(paste0(ds_id,"/",ds_id, "_csvld_", csvld, ".pheno"), row.names = F, col.names = T, sep = "\t")
}
  

########################################
# compute the cross validation PRSxTopic 
########################################
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
patient_loadings <- sweep((para$alpha_z - 1), 1, rowSums(para$alpha_z -1), FUN="/") 

rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv") 
topic_individuals <- rec_data %>%
  group_by(eid) %>%
  slice(1)


# ds_id <- "250.2"
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
    df_umap_patient <- data.frame(eid = para$eid, loadings = patient_loadings[,topic_id])  %>%
      mutate(loadings = loadings/sd(loadings))
    
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
      filter(phenotype == 1)
    
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
write.csv(df_subtype_regression, paste0("~/Desktop/comorbidity/paper_writing/Production_figures/bolt_topic_regression.csv" ))
subtype_regression <- list(ds_list, p_values, coefficients, se_coef)
save(subtype_regression, file = paste0("~/Desktop/comorbidity/paper_writing/Production_figures/bolt_topic_regression.Rdata" ))

# adjust for multiple testing
ds_list <- subtype_regression[[1]]
p_values <- subtype_regression[[2]]
p_values <- bind_cols(ds_list$phenotype, data.frame(p_values))
colnames(p_values) <-  c("phenotype","MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
p_values <- p_values %>%
  pivot_longer(cols = !phenotype, names_to = "topic", values_to = "p_value") %>%
  filter(!is.na(p_value))
coefficients <- subtype_regression[[3]]
coefficients <- bind_cols(ds_list$phenotype, data.frame(coefficients))
colnames(coefficients) <- c("phenotype","MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
coefficients <- coefficients %>%
  pivot_longer(cols = !phenotype, names_to = "topic", values_to = "coefficients") %>%
  filter(!is.na(coefficients))
regression_PRS <- p_values %>%
  left_join(coefficients, by = c("phenotype", "topic")) 
regression_PRS$FDR <- p.adjust(regression_PRS$p_value)
regression_PRS %>% select(phenotype, topic, coefficients, p_value, FDR) %>%
  filter(FDR < 0.1) %>%
  write.csv("~/Desktop/comorbidity/paper_writing/Production_figures/bolt_PRS_regression.csv", row.names = F)

# Fig 4 sub figure
# first step: add one analysis of the ARP topic for hypertension
load("~/Desktop/comorbidity/paper_writing/Production_figures/bolt_topic_regression.Rdata" )
p_values <- subtype_regression[[2]]
coefficients <- subtype_regression[[3]]
se_coef <-  subtype_regression[[4]]

library(reshape2)
idx_plot <- which(ds_list$V1 %in% c(401.1, 272.11,  495, 250.2))
p_plot <- p_values[idx_plot, ]
colnames(p_plot) <- topic_name
rownames(p_plot) <- c( "Essential hypertension", "Hypercholesterolemia", "Asthma", "Type 2 diabetes")
p_plot <- melt(p_plot[,c(3,5,6,7)]) %>%
  rename(P = value)
coef_plot <- coefficients[idx_plot, ]
colnames(coef_plot) <- topic_name
rownames(coef_plot) <- c( "Essential hypertension", "Hypercholesterolemia", "Asthma", "Type 2 diabetes")
coef_plot <- melt(coef_plot[,c(3,5,6,7)])%>%
  rename(coefficient = value)
se_plot <- se_coef[idx_plot, ]
colnames(se_plot) <- topic_name
rownames(se_plot) <- c( "Essential hypertension", "Hypercholesterolemia", "Asthma", "Type 2 diabetes")
se_plot <- melt(se_plot[,c(3,5,6,7)])%>%
  rename(se_coef = value)

df_plot <- coef_plot %>%
  left_join(p_plot, by = c("Var1", "Var2")) %>%
  left_join(se_plot, by = c("Var1", "Var2")) %>%
  mutate(P_levels = -log10(P) ) 
# use color that are different
paletteMartin[[4]]
Brown2Blue12Steps[[10]]
plt <- ggplot(df_plot) +
  geom_point(aes(x = Var2, y = Var1, fill = coefficient, size = P_levels), shape = 21) +
  scale_fill_gradient2(low = "#32E3FF", mid = "white", high = "#ff6db6") + 
  scale_size_area(name = expression(-log[10](P)), max_size = 20) + 
  labs(x = "", y = "") +
  theme_bw(base_size = 20)+
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/Subtype_liability_interaction.png", plt, width = 6, height = 4.5)
plt <- ggplot(df_plot) +
  geom_point(aes(x = Var2, y = Var1, fill = coefficient, size = P_levels), shape = 21) +
  scale_fill_gradient2(low = "#32E3FF", mid = "white", high = "#ff6db6") + 
  scale_size_area(name = expression(-log[10](P)), max_size = 40) + 
  labs(x = "", y = "") +
  theme_bw(base_size = 20)
legend_plt <- cowplot::get_legend(plt)
grid.newpage()
legend_plt <- grid.draw(legend_plt)

# plot bar plot
topic_name_pallete <- as.vector(colors_topic)
names(topic_name_pallete) <- topic_name
plt_bar <- ggplot(df_plot, mapping = aes(x = Var1, y = coefficient, fill = Var2)) +
  geom_col( position=position_dodge(0.8), width=.7, alpha = 0.8) + 
  geom_errorbar(aes(ymin=coefficient- se_coef, ymax=coefficient+ se_coef), width=.2, position=position_dodge(0.8)) +
  scale_fill_manual(values = topic_name_pallete) +
  coord_flip() + 
  theme_bw(base_size = 20)+
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) + labs(x = "", y = "") 
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/PRS_topic_bar.png", plt_bar, width = 6, height = 4.5)
plt <- ggplot(df_plot, mapping = aes(x = Var1, y = coefficient, fill = Var2)) +
  geom_col( position=position_dodge(0.8), width=.7, alpha = 0.5) + 
  geom_errorbar(aes(ymin=coefficient- se_coef, ymax=coefficient+ se_coef), width=.2, position=position_dodge(0.8)) +
  scale_fill_manual(values = topic_name_pallete) +
  coord_flip() + 
  theme_bw(base_size = 20)
legend_plt <- cowplot::get_legend(plt)
grid.newpage()
legend_plt <- grid.draw(legend_plt)

# show the significant stars 
fdr_bolt <- read.csv("~/Desktop/comorbidity/paper_writing/Production_figures/bolt_PRS_regression.csv")
fdr_bolt %>%
  filter(topic %in% c("CVD", "MGND", "CER", "FGND"), FDR < 0.01)

##############
# supplementary figures: bar plot for 10 diseases with 10 topics
##############
load("~/Desktop/comorbidity/paper_writing/Production_figures/bolt_topic_regression.Rdata" )
p_values <- subtype_regression[[2]]
coefficients <- subtype_regression[[3]]
se_coef <-  subtype_regression[[4]]

library(reshape2)
colnames(p_values) <- topic_name
p_plot <- data.frame(Phecode = ds_list$V1, Phenotype = ds_list$phenotype, p_values) 
p_plot <- p_plot %>% 
  pivot_longer(cols = MDS:SRD, names_to = "Topic", values_to = "P")

colnames(se_coef) <- topic_name
se_plot <- data.frame(Phecode = ds_list$V1, Phenotype = ds_list$phenotype, se_coef) 
se_plot <- se_plot %>% 
  pivot_longer(cols = MDS:SRD, names_to = "Topic", values_to = "se_coef")

colnames(coefficients) <- topic_name
coef_plot <- data.frame(Phecode = ds_list$V1, Phenotype = ds_list$phenotype, coefficients) 
coef_plot <- coef_plot %>% 
  pivot_longer(cols = MDS:SRD, names_to = "Topic", values_to = "coefficient")

topic_ordered_namte <- c("NRI", "CER","SRD", "CVD", "UGI", "LGI", "FGND",  "MGND", "MDS", "ARP")
df_plot <- coef_plot %>%
  left_join(p_plot, by = c("Phecode", "Phenotype", "Topic")) %>%
  left_join(se_plot, by = c("Phecode", "Phenotype","Topic")) %>%
  mutate(Topic = factor(Topic, levels = rev(topic_ordered_namte)), P_levels = -log10(P) ) 
df_plot$FDR <- p.adjust(df_plot$P)
df_plot <- df_plot %>%
  mutate(stars = if_else(FDR < 0.1, 1, 0))
plt_bar <- ggplot(df_plot, mapping = aes(x = Phenotype, y = coefficient, fill = Topic)) +
  geom_col( position=position_dodge(0.8), width=.7, alpha = 0.8) + 
  geom_errorbar(aes(ymin=coefficient- se_coef, ymax=coefficient+ se_coef), width=.2, position=position_dodge(0.8)) +
  scale_fill_manual(values = topic_name_pallete) +
  geom_point(aes(y = coefficient + sign(coefficient) * (se_coef + 0.005), alpha = stars), shape = 8, position = position_dodge(width=0.8))  +
  scale_alpha(range = c(0, 1)) + 
  coord_flip() + 
  theme_bw(base_size = 20)+
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) + labs(x = "", y = "") 
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/all_PRS_topic_bar.png", plt_bar, width = 8, height = 10)
plt <- ggplot(df_plot, mapping = aes(x = Phenotype, y = coefficient, fill = Topic)) +
  geom_col( position=position_dodge(0.8), width=.7, alpha = 0.8) + 
  geom_errorbar(aes(ymin=coefficient- se_coef, ymax=coefficient+ se_coef), width=.2, position=position_dodge(0.8)) +
  scale_fill_manual(values = topic_name_pallete) +
  coord_flip() + 
  theme_bw(base_size = 20)
legend_plt <- cowplot::get_legend(plt)
grid.newpage()
legend_plt <- grid.draw(legend_plt)

# plot PRS enrichment with respect to sample size.
library(colorBlindness)
displayAvailablePalette(color="white")
subtype_info <- list()
for(idx in 1:dim(ds_list)[1]){
  ds_id <- ds_list$V1[idx]
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
  filter(subtype_percentage > 0.1) %>%
  mutate(topic = topic_name[topic]) %>%
  left_join(select(df_plot, Phecode, Topic, coefficient, se_coef), by = c("disease" ="Phecode", "topic" = "Topic"))
  
plt <- ggplot(subtype_count) +
  geom_pointrange(aes(x = subtype_percentage, y = coefficient, ymin =coefficient -1.96*se_coef, ymax = coefficient + 1.96*se_coef), color = red) + 
  labs(x = "Subtype proportion", y = "PRS enrichment") + 
  scale_x_continuous(labels = scales::percent) + 
  theme_bw(base_size = 20) +
  theme(legend.position = "None",panel.background=element_blank()) 
ggsave("~/Desktop/comorbidity/paper_writing/supplementary_files/PRS_vs_sampe_size.png", plt, width = 6, height = 5)

plt <- ggplot(subtype_count) +
  geom_pointrange(aes(x = age, y = coefficient, ymin =coefficient -1.96*se_coef, ymax = coefficient + 1.96*se_coef), color = red) + 
  labs(x = "Age (Years)", y = "PRS enrichment") + 
  theme_bw(base_size = 20) +
  theme(legend.position = "None",panel.background=element_blank()) 
ggsave("~/Desktop/comorbidity/paper_writing/supplementary_files/PRS_vs_age.png", plt, width = 6, height = 5)

 ###########################################
# plot the percentile plot for the four diseases
###########################################
british_isle <- read.table("../UKBB_interim_files/keep_britishIsle.txt", header = F)
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
DIR <- "~/Desktop/comorbidity/association_results/"
ds_list <- read.table("BOLT_LMM_disease_list.txt")
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
# playing with the code
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("V1" = "phecode"))

colors_topic <- c(cbPalette, paletteMartin[c(8,11)])
names(colors_topic) <- as.character(1:10)

for(ds_id in c("401.1", "272.11", "250.2", "495")){
  idx <- which(ds_list$V1 == ds_id)
  
  # topic_list <- c(3, 5,6,7)
  topic_list <- c(6,7)
  # topic_list <- c(5, 7)
  # check the PRS data are indeed correct
  ds_eid <- rec_data %>% 
    filter(diag_icd10  == ds_id) %>%
    select(eid)

  # load the cross validation risk profile
  prediction_set <- list()
  for(cv_id in 1:5){
    cv_profile <- read.table(paste0("bolt_prs_result/", ds_id, "_csvld_", cv_id, ".bolt.prs.profile"), header = T)
    training_set <- read.table(paste0("bolt_prs_result/", ds_id, "_csvld_", cv_id, ".pheno"), header = T) 
    prediction_set[[cv_id]] <- read.table(paste0("bolt_prs_result/", ds_id, "keep.txt"), header = F) %>%
      anti_join(training_set, by = c("V1" = "FID")) %>%
      semi_join(topic_individuals,by = c("V1" = "eid")) %>%
      mutate(phenotype = if_else(V1 %in% ds_eid$eid, 1, 0)) %>%
      left_join(select(cv_profile, FID, SCORE), by= c("V1" = "FID"))
  }
  PRS_profile <- bind_rows(prediction_set) 
  
  # PRS_profile <- read.table(paste0(DIR, ds_id, ".profile"), header =T)
  
  # get the disease of different groups
  j <- match(as.numeric(ds_id), para$list_above500occu$diag_icd10)
  cases_eid <- list()
  for(topic_id in 1:para$K){
    topic_specific_id <- which(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max) == topic_id)
    topic_specific_set <- para$unlist_Ds_id[para$ds_list[[j]]$id,][topic_specific_id,]
    if(dim(topic_specific_set)[1] > 500 & 
       dim(topic_specific_set)[1]/dim(ds_eid)[1] > 0.05 & 
       (topic_id %in% topic_list)){
      cases_eid[[topic_id]] <- topic_specific_set %>%
        mutate(topic = topic_id) %>%
        mutate(Ds_id = para$list_above500occu[Ds_id, 1])
    }
  }
  cases_eid <- bind_rows(cases_eid) %>%
    semi_join(british_isle, by = c("eid" = "V1"))
  
  # plot percentile PRS with incidence rate
  tiles_num <- 50
  percentile <- PRS_profile %>% mutate(percentiles = ntile(SCORE,tiles_num))
  prevalence <- rep(NA, tiles_num)
  prevalence_df <- list()
  for(k in topic_list){
    eid_group <- cases_eid %>%
      filter(topic == k)
    prevalence_subgroup <- rep(NA, tiles_num)
    for(pct in 1:tiles_num){
      pct_eid <- percentile %>%
        filter(percentiles == pct) %>%
        pull(1)
      # compute log(Relative risk)
      prevalence_subgroup[pct] <- tiles_num * sum(pct_eid %in% eid_group$eid )/length(eid_group$eid)
    }
    prevalence_df[[k]] <- data.frame(percentile = 1:tiles_num/tiles_num, topic = k, prevalence_subgroup = prevalence_subgroup)
  }
  
  # compute background (all cases)
  prevalence_all <- rep(NA, tiles_num)
  for(pct in 1:tiles_num){
    pct_eid <- percentile %>%
      filter(percentiles == pct) %>%
      pull(1)
    # compute log(Relative risk)
    prevalence_all[pct] <- tiles_num * sum(pct_eid %in% cases_eid$eid )/length(cases_eid$eid)
  }
  prevalence_all <- data.frame(percentile = 1:tiles_num/tiles_num, topic = 0, prevalence_all = prevalence_all)
  
  prevalence_df <- bind_rows(prevalence_df)
  plt <- ggplot(prevalence_df) +
    geom_point(aes(x = percentile, y = log(prevalence_subgroup), color = as.character(topic)), alpha = 0.7) +
    geom_smooth( aes(x = percentile,  y = log(prevalence_subgroup), color = as.character(topic), fill = as.character(topic) ), method = "lm", alpha = 0.5) + 
    scale_color_manual(values =colors_topic) + 
    scale_fill_manual(values =colors_topic) + 
    scale_x_continuous(labels = scales::percent) + 
    theme_bw(base_size = 15) + 
    theme(legend.position = "None") + 
    labs(x="PRS percentile", y="Log relative risk", title=paste0("Disease: ", ds_list$phenotype[idx])) 
    
  ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/","subtypesExample_PRS_",ds_id,".png"), plt, width = 6, height = 4)
  
  plt_lines <- ggplot(prevalence_df) +
    geom_smooth(aes(x = percentile, y = prevalence_subgroup, color = as.character(topic)), alpha = 0.7) +
    scale_color_manual(values =colors_topic) + 
    scale_x_continuous(labels = scales::percent) + 
    theme_bw(base_size = 15) + 
    theme(legend.position = "None") + 
    labs(x="PRS percentile", y="Relative risk", title=paste0("Disease: ", ds_list$phenotype[idx])) 
  ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/","subtypesExample_PRS_smoothlines",ds_id,".png"), plt_lines, width = 6, height = 4)
  
}

 

