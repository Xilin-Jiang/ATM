########################
# perform LDSC-seg analysis
########################
h2g_disease_topics <- read.csv("causality_analysis/h2g_causality.csv")
filtered_heritable <- h2g_disease_topics %>% 
  filter(z_score > 4 | disease == 250.2) 
phe_phecode <- read.csv("info_phenotype_phecode.csv") 

# there are totally 11 diseases that have subtypes of h2g z-score > 4
# ds_of_interest <- filtered_heritable  %>% 
#   group_by(disease) %>% 
#   summarise(subnum = n()) %>% 
#   filter(subnum > 1) %>%
#   left_join(phe_phecode, by= c("disease" = "phecode"))

all_disease <- read.table("all_disease_topic_list.txt", header =  F) 
ds_of_interest <- all_disease %>%
  group_by(V1) %>% 
  summarise(subnum = n()) %>% 
  filter(subnum > 1) %>%
  rename(disease = V1) %>%
  left_join(phe_phecode, by= c("disease" = "phecode"))
# compute the sample size 
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
thre_pick <- 1000
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- para$list_above500occu %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") )
icd10 <- ds_list$phenotype
# using the larger results file: ordering are changed in "model_output"
betas <- para$pi_beta_basis
topic_assign <- c()
disease_idx <- c()
mean_age <- c()
number_case <- c()
for(topic_id in 1:K){
  # get the average age of incidences
  for(j in 1:para$D){
    topic_specific_id <- which(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max) == topic_id)
    if(length(topic_specific_id) > thre_pick){
      number_case <- c(number_case, length(topic_specific_id))
      mean_age <- c(mean_age, mean(para$ds_list[[j]][topic_specific_id, ]$age_diag))
      topic_assign <- c(topic_assign, topic_id)
      disease_idx <- c(disease_idx, para$list_above500occu$diag_icd10[j])
    }
  }
}
common_disease_within_topics <- data.frame(disease = disease_idx, topic = topic_assign, mean_age = mean_age, number_case = number_case) %>%
  mutate(topic = as.character(topic))
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")

h2g_disease_topics %>% 
  filter(disease %in% ds_of_interest$disease) %>%
  left_join(select(common_disease_within_topics, disease, topic, number_case), by = c("disease" = "disease", "topic" = "topic")) %>%
  left_join(ds_list, by = c("disease" = "diag_icd10")) %>%
  mutate(number_case = if_else(is.na(number_case), occ, number_case)) %>%
  write.table("ldsc_SEG_disease_topic_list.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)

##########################################
# extract tissue specific analysis results
##########################################
CTS_dir <- "~/Desktop/comorbidity/Multi-morbidity_biobank/CTS_results/"

disease_ldsc_seg <- read.table("ldsc_SEG_disease_topic_list.txt", header =  F)

names(disease_ldsc_seg) <- c("disease", "topic", "age", "h2g", "h2g_se", "h2g_zscore", "subtype_size", "disease_size")

# keep diseases where at least one subtype have h2g > 4 
focused_disease <- disease_ldsc_seg %>%
  #filter(topic != "all")  %>% 
  filter(h2g_zscore > 4) %>%
  pull(disease) %>%
  unique()

disease_ldsc_seg <- disease_ldsc_seg %>% 
  #filter(topic != "all")  %>% 
  filter(disease %in% focused_disease)


cts_ldsc <- list()
for(id in 1:dim(disease_ldsc_seg)[1]){
  ds_id <- disease_ldsc_seg$disease[id]
  topic_id <- disease_ldsc_seg$topic[id]
  
  cts_ldsc[[id]] <- read.table(paste0(CTS_dir, ds_id, "_topic", topic_id, ".cell_type_results.txt"), header = T) %>%
    mutate(disease = ds_id, topic = topic_id)
}

cts_ldsc <- bind_rows(cts_ldsc)
plot_cts <- cts_ldsc

# plot_cts <- cts_ldsc %>%
#   filter(disease == 300.1, topic == 1)

plt_qq <- data.frame(y = sort(-log10(plot_cts$Coefficient_P_value)), x = sort(-log10(runif(length(plot_cts$Coefficient_P_value)))))

ggplot(plt_qq) + 
  geom_point(aes(x = x, y = y), size = 0.5, alpha = 0.5) + 
  geom_abline(slope = 1, linetype = "dashed", color = red) +
  theme(legend.position = "none",panel.background=element_blank()) + 
  xlab(expression("Expected" * -log[10](P)) ) + ylab(expression("Observed" * -log[10](P))) +
  ggtitle(paste0(ds_name, " x all topic: "))


plot_cts$fdr <- p.adjust(plot_cts$Coefficient_P_value, method = "fdr")
plot_cts %>% arrange(fdr)

# plot similar figure as Finucane et al. 2018 NG
tissue_category <- read.csv("tissue_category.csv") %>%
  mutate(Tissue.category.for.display = if_else(Tissue.category.for.display == "Musculoskeletal/connective",
                                               "Musculoskeletal/Connective", Tissue.category.for.display)) %>%
  mutate(Tissue = str_replace_all(Tissue, " ", "_")) # remove all the space with lower slash
order_tissues_catogory <- unique(tissue_category$Tissue.category.for.display)
tissue_category <- tissue_category %>%
  mutate(Tissue.category.for.display = factor(Tissue.category.for.display, levels = order_tissues_catogory)) %>%
  arrange(Tissue.category.for.display)
order_tissues_type <- unique(tissue_category$Tissue)
tissue_category <- tissue_category %>%
  mutate(Tissue = factor(Tissue, levels = order_tissues_type)) %>%
  select(Tissue, Tissue.category.for.display)

plot_cts <- plot_cts %>%
  left_join(tissue_category, by = c("Name" = "Tissue")) %>%
  mutate(Name = factor(Name, levels = order_tissues_type))
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CAD", "UGI", "LGI", "SRD")
disease_ldsc_seg <- disease_ldsc_seg %>%
  left_join(phe_phecode, by = c("disease" = "phecode") )
for(id in 1:dim(disease_ldsc_seg)[1]){
  ds_id <- disease_ldsc_seg$disease[id]
  topic_id <- disease_ldsc_seg$topic[id]
  plot_df <- plot_cts %>%
    filter(disease == ds_id, topic == topic_id)
  if(topic_id == "all"){
    subtype_name <- "all"
  }else{
    subtype_name <- topic_name[as.numeric(topic_id)]
  }
  plt <- ggplot(plot_df) +
    geom_point(aes(x = Name, y = -log10(Coefficient_P_value), color = Tissue.category.for.display)) + 
    labs(x = "Tissue/Cell Type", y =  expression(-log[10](P)),
         title = paste0(disease_ldsc_seg$phenotype[id], " \n ",
                        subtype_name, "\n",
                        "sample size: ", disease_ldsc_seg$subtype_size[id])) + 
    theme_bw(base_size = 20)+
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.ticks.x = element_blank(), axis.text.x = element_blank()
    )
  ggsave(paste0("~/Desktop/comorbidity/figures/LDSC_SEG/", disease_ldsc_seg$phenotype[id], " x " ,  subtype_name, ".png" ), plt, width = 6, height = 4)
    
}

########################
# analysis procedure
########################
# get all the significants for each trait with all
disease_ldsc_seg <- read.table("ldsc_SEG_disease_topic_list.txt", header =  F)
names(disease_ldsc_seg) <- c("disease", "topic", "age", "h2g", "h2g_se", "h2g_zscore", "subtype_size", "disease_size")
all_ds <- disease_ldsc_seg 
cts_fdr_ldsc <- list()
for(id in 1:dim(all_ds)[1]){
  ds_id <- all_ds$disease[id]
  topic_id <- all_ds$topic[id]
  
  rslt_cts <- read.table(paste0(CTS_dir, ds_id, "_topic", topic_id, ".cell_type_results.txt"), header = T) %>%
    mutate(disease = ds_id, topic = topic_id)
  
  rslt_cts$fdr <- p.adjust(rslt_cts$Coefficient_P_value, method = "fdr")
  
  # cts_fdr_ldsc[[id]] <- rslt_cts %>%
  #   filter(fdr < 0.1)
  cts_fdr_ldsc[[id]] <- rslt_cts %>%
       filter(Coefficient_P_value < 0.01)
  
}
cts_fdr_ldsc <- bind_rows(cts_fdr_ldsc)

ds_list <- Coefficient_P_value %>%
ds_id <- 244.4 # 495

cts_specific <- cts_fdr_ldsc %>%
  filter(disease == ds_id, topic == "all") %>%
  pull(Name)

specific_ds <- cts_ldsc %>%
  filter(disease == ds_id) %>%
  filter(topic != "all")

topics_sub <- specific_ds$topic %>%
  unique()

subtp_1 <- cts_ldsc %>%
  filter(Name %in% cts_specific, topic == topics_sub[1], disease == ds_id)
subtp_2 <- cts_ldsc %>%
  filter(Name %in% cts_specific, topic == topics_sub[2], disease == ds_id)

joint_data <- subtp_1 %>%
  select(Name, Coefficient, Coefficient_std_error) %>%
  left_join(select(subtp_2, Name, Coefficient, Coefficient_std_error), by ="Name") %>%
  mutate(p_diff = (1-pnorm(abs(Coefficient.x  - Coefficient.y )/
                             sqrt(Coefficient_std_error.x^2 + Coefficient_std_error.y^2) ))*2)
  


#######################################
# BOLT-LMM analysis to improve power
#######################################
# use all subtypes (>1000, no more constraint on threshold of 0.5 for per-case topic weights)
subtypes <- common_disease_within_topics %>% 
  filter(disease %in% focused_disease)
diseases_nosubtype <- disease_ldsc_seg %>%
  filter(topic == "all") %>%
  rename(mean_age = age, number_case= subtype_size) %>%
  select(disease, topic, mean_age, number_case)
bolt_diseaes_subtype <- bind_rows(subtypes, diseases_nosubtype)
# use this for BOLT_LMM analysis
write.table(bolt_diseaes_subtype, "BOLT_LMM_subtype_list.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)

###################################
# Plot BOLT-LMM results
###################################
CTS_dir <- "~/Desktop/comorbidity/Multi-morbidity_biobank/CTS_results/"
disease_ldsc_seg <- read.table("BOLT_LMM_subtype_list.txt", header =  F) %>%
  arrange(V1)
names(disease_ldsc_seg) <- c("disease", "topic", "age",  "subtype_size")

cts_ldsc <- list()
for(id in 1:dim(disease_ldsc_seg)[1]){
  ds_id <- disease_ldsc_seg$disease[id]
  topic_id <- disease_ldsc_seg$topic[id]
  try({
    cts_ldsc[[id]] <- read.table(paste0(CTS_dir, ds_id, "_topic", topic_id, "_imputed.cell_type_results.txt"), header = T) %>%
      mutate(disease = ds_id, topic = topic_id)
  })

}

cts_ldsc <- bind_rows(cts_ldsc)
plot_cts <- cts_ldsc 

plot_cts <- cts_ldsc %>%
  filter(disease == 250.2, topic ==  "5")

plt_qq <- data.frame(y = sort(-log10(plot_cts$Coefficient_P_value)), x = sort(-log10(runif(length(plot_cts$Coefficient_P_value)))))

ggplot(plt_qq) + 
  geom_point(aes(x = x, y = y), size = 0.5, alpha = 0.5) + 
  geom_abline(slope = 1, linetype = "dashed", color = red) +
  theme(legend.position = "none",panel.background=element_blank()) + 
  xlab(expression("Expected" * -log[10](P)) ) + ylab(expression("Observed" * -log[10](P))) +
  ggtitle(paste0(ds_name, " x all topic: "))


plot_cts$fdr <- p.adjust(plot_cts$Coefficient_P_value, method = "fdr")
plot_cts %>%
  filter(disease == 250.2, topic != "all") %>%
  arrange(fdr)

# plot similar figure as Finucane et al. 2018 NG
tissue_category <- read.csv("tissue_category.csv") %>%
  mutate(Tissue.category.for.display = if_else(Tissue.category.for.display == "Musculoskeletal/connective",
                                               "Musculoskeletal/Connective", Tissue.category.for.display)) %>%
  mutate(Tissue = str_replace_all(Tissue, " ", "_")) # remove all the space with lower slash
order_tissues_catogory <- unique(tissue_category$Tissue.category.for.display)
tissue_category <- tissue_category %>%
  mutate(Tissue.category.for.display = factor(Tissue.category.for.display, levels = order_tissues_catogory)) %>%
  arrange(Tissue.category.for.display)
order_tissues_type <- unique(tissue_category$Tissue)
tissue_category <- tissue_category %>%
  mutate(Tissue = factor(Tissue, levels = order_tissues_type)) %>%
  select(Tissue, Tissue.category.for.display)

plot_cts <- plot_cts %>%
  left_join(tissue_category, by = c("Name" = "Tissue")) %>%
  mutate(Name = factor(Name, levels = order_tissues_type))
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CAD", "UGI", "LGI", "SRD")
disease_ldsc_seg <- disease_ldsc_seg %>%
  left_join(phe_phecode, by = c("disease" = "phecode") ) 
for(id in 1:dim(disease_ldsc_seg)[1]){
  ds_id <- disease_ldsc_seg$disease[id]
  topic_id <- disease_ldsc_seg$topic[id]
  plot_df <- plot_cts %>%
    filter(disease == ds_id, topic == topic_id)
  if(topic_id == "all"){
    subtype_name <- "all"
  }else{
    subtype_name <- topic_name[as.numeric(topic_id)]
  }
  plt <- ggplot(plot_df) +
    geom_point(aes(x = Name, y = -log10(Coefficient_P_value), color = Tissue.category.for.display)) + 
    labs(x = "Tissue/Cell Type", y =  expression(-log[10](P)),
         title = paste0(disease_ldsc_seg$phenotype[id], " \n ",
                        subtype_name, "\n",
                        "sample size: ", disease_ldsc_seg$subtype_size[id])) + 
    theme_bw(base_size = 20)+
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.ticks.x = element_blank(), axis.text.x = element_blank()
    )
  ggsave(paste0("~/Desktop/comorbidity/figures/LDSC_SEG/", disease_ldsc_seg$phenotype[id], " x " ,  subtype_name, "_BOLT_imputed.png" ), plt, width = 6, height = 4)
  
}

##############################################
# focusing on only those with > 5000 samples 
#############################################
CTS_dir <- "~/Desktop/comorbidity/Multi-morbidity_biobank/CTS_results/"
disease_ldsc_seg <- read.table("BOLT_LMM_subtype_list.txt", header =  F) %>%
  arrange(V1)
names(disease_ldsc_seg) <- c("disease", "topic", "age",  "subtype_size")

# all_ldsc_seg <- read.table("ldsc_SEG_disease_topic_list.txt", header =  F)
# names(all_ldsc_seg) <- c("disease", "topic", "age", "h2g", "h2g_se", "h2g_zscore", "subtype_size", "disease_size")
all_ldsc_seg <- read.csv("h2g_imputed.csv", header =  T)
names(all_ldsc_seg) <- c("disease", "topic", "age", "subtype_size", "h2g", "h2g_se", "h2g_zscore")

all_ldsc_seg %>% 
  filter(h2g_zscore > 5, topic != "all") %>% 
  arrange(disease) %>%
  left_join(select(phe_phecode, phecode, phenotype), by = c("disease" = "phecode") )  

filtered_results <- disease_ldsc_seg %>%
  left_join(select(all_ldsc_seg, disease, topic, h2g, h2g_se, h2g_zscore), by =c("disease", "topic")) %>%
  #filter(topic != "all")
  filter( h2g_zscore > 5, topic != "all")
cts_ldsc <- list()
for(id in 1:dim(filtered_results)[1]){
  ds_id <- filtered_results$disease[id]
  topic_id <- filtered_results$topic[id]
  try({
    cts_ldsc[[id]] <- read.table(paste0(CTS_dir, ds_id, "_topic", topic_id, "_imputed.cell_type_results.txt"), header = T) %>%
      mutate(disease = ds_id, topic = topic_id)
  })
  
}

cts_ldsc <- bind_rows(cts_ldsc)


#######################################
# plot the P-value of all cases v.s. subtypes
#######################################
source("topic_functions.R")
CTS_dir <- "~/Desktop/comorbidity/Multi-morbidity_biobank/CTS_results/"
all_ds_results <- disease_ldsc_seg %>%
  left_join(select(all_ldsc_seg, disease, topic, h2g, h2g_se, h2g_zscore), by =c("disease", "topic")) %>%
  filter(topic == "all")
 # filter(subtype_size > 5000, h2g_zscore > 3, topic == "all")
all_cts_ldsc <- list()
for(id in 1:dim(all_ds_results)[1]){
  ds_id <- all_ds_results$disease[id]
  topic_id <- all_ds_results$topic[id]
  try({
    all_cts_ldsc[[id]] <- read.table(paste0(CTS_dir, ds_id, "_topic", topic_id, "_imputed.cell_type_results.txt"), header = T) %>%
      mutate(disease = ds_id, topic = topic_id)
  })
  
}

all_cts_ldsc <- bind_rows(all_cts_ldsc) %>%
  rename(P_all_case = Coefficient_P_value)
matched_all_sub_cts_ldsc <- cts_ldsc %>%
  left_join(select(all_cts_ldsc, Name, disease, P_all_case), by = c("disease", "Name") ) %>%
  filter(Coefficient_P_value < 0.1 | P_all_case < 0.1)
  # filter(Coefficient_P_value < 0.1 | P_all_case < 0.1)

ggplot(matched_all_sub_cts_ldsc) + 
  geom_point(aes(x =-log10(P_all_case), y = -log10(Coefficient_P_value))) + 
  geom_abline(slope = 1) +
  lims(x= c(0,9), y = c(0, 9))

cts_ldsc$fdr <- p.adjust(cts_ldsc$Coefficient_P_value, method = "fdr")
cts_ldsc %>%
  arrange(fdr) %>%
  left_join(select(all_cts_ldsc, Name, disease, P_all_case), by = c("disease", "Name") ) %>%
  filter(fdr < 0.2, topic != "all") %>% # , disease != 244.4)
  rename(Tissue = Name) %>%
  left_join(select(phe_phecode, phecode, phenotype), by = c("disease" = "phecode") )  %>%
  mutate(topic = topic_name[as.numeric(topic)]) %>%
  select(-disease) %>%
  arrange(phenotype) %>%
  rename(Topic = topic, FDR = fdr, Disease = phenotype) %>%
  write.csv("~/Desktop/comorbidity/paper_writing/supplementary_files/LDSC_SEG_table.csv",row.names = F)
  
##########################################################
# use the jackknife samples to compute difference p-value
##########################################################
CTS_dir <- "~/Desktop/comorbidity/Multi-morbidity_biobank/CTS_results/CTS_results/"
disease_ldsc_seg <- read.table("BOLT_LMM_subtype_list.txt", header =  F) %>%
  arrange(V1)
names(disease_ldsc_seg) <- c("disease", "topic", "age",  "subtype_size")

cts_names <- read.table(paste0(CTS_dir, "153.2_topic4_imputed.cell_type_results.txt"), header = T) %>%
  pull(Name)


topic_id <- 7
disease_id <- 250.2

strong_genetic_subtype <- read.csv("Association_analysis/subtypes_Fst_matched_topic.csv") %>%
  filter(p_fst < 0.011) %>%
  pull(disease)

all_ldsc_seg <- read.csv("h2g_imputed.csv", header =  T)
names(all_ldsc_seg) <- c("disease", "topic", "age", "subtype_size", "h2g", "h2g_se", "h2g_zscore")

ds_high_h2g <- all_ldsc_seg %>% 
  filter(disease %in% strong_genetic_subtype, topic != "all", h2g_zscore > 5) %>%
  left_join(select(phe_phecode, phecode, phenotype), by = c("disease" = "phecode") )  
  

# ds_high_h2g <- all_ldsc_seg %>% 
#   filter(h2g_zscore > 5, topic != "all") %>% 
#   arrange(disease) %>%
#   left_join(select(phe_phecode, phecode, phenotype), by = c("disease" = "phecode") )  

LDSC_seg_subtp_diff <- list()
for(ds_tp_id in 1:dim(ds_high_h2g)[1]){
  topic_id <- ds_high_h2g$topic[ds_tp_id]
  disease_id <- ds_high_h2g$disease[ds_tp_id]
  print(disease_id)
  
  cts_mean <- c()
  cts_se <- c()
  cts_p <- c()
  for(num_cts in 1:length(cts_names)){
    cts_id <- cts_names[num_cts]
    subtp_cts <- read.table(paste0(CTS_dir, disease_id, "_topic", topic_id,"_imputed.", cts_id, ".part_delete")) %>%
      pull(V54)
    
    all_cts <- read.table(paste0(CTS_dir, disease_id, "_topicall_imputed.", cts_id, ".part_delete")) %>%
      pull(V54)
    Jacknife_samples <- (subtp_cts - all_cts)
    cts_mean[num_cts] <- mean(Jacknife_samples)
    cts_se[num_cts] <- sqrt(var(Jacknife_samples) * (length(Jacknife_samples) - 1)^2/ length(Jacknife_samples) )
    cts_p[num_cts] <- (1-pnorm(cts_mean[num_cts]/cts_se[num_cts]))
  }
  LDSC_seg_subtp_diff[[ds_tp_id]] <- data.frame(Name = cts_names, diff_mean = cts_mean, diff_se = cts_se, P = cts_p) %>%
    mutate(topic = topic_id, disease = disease_id)
  
}
LDSC_seg_subtp_diff <- bind_rows(LDSC_seg_subtp_diff)
LDSC_seg_subtp_diff$FDR <- p.adjust(LDSC_seg_subtp_diff$P , method = "fdr")
LDSC_seg_subtp_diff %>% 
  arrange(P)


plot_ds <- LDSC_seg_subtp_diff 
plt_qq <- data.frame(y = sort(-log10(plot_ds$P)), x = sort(-log10(runif(length(plot_ds$P)))))

ggplot(plt_qq) + 
  geom_point(aes(x = x, y = y), size = 0.5, alpha = 0.5) + 
  geom_abline(slope = 1, linetype = "dashed", color = red) +
  theme(legend.position = "none",panel.background=element_blank()) + 
  xlab(expression("Expected" * -log[10](P)) ) + ylab(expression("Observed" * -log[10](P))) +
  ggtitle(paste0("LDSC: Subtypes v.s. all disease"))





    
    
    
  

