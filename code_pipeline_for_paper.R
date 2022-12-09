# this file serve as the full pipeline for results in the paper
# (please put the exploration analysis in other files!)

source("topic_functions.R")
source("plotting_functions.R")

########################################################
# data preparation: Phecode and prevalence Threshold
########################################################

birthyear <- read.csv("~/Desktop/genetics_longitudinal_data/longitudinal_data/Year_of_birth.csv")

hes_diag <- read.table(file = '~/Desktop/genetics_longitudinal_data/longitudinal_data/hesin_diag10_130418.tsv',quote = '', sep = '\t', header = TRUE, fill=TRUE) %>%
  mutate(record_id = as.factor(record_id))
hes <- read.table(file = '~/Desktop/genetics_longitudinal_data/longitudinal_data/hesin_130418.tsv',quote = '', sep = '\t', header = TRUE, fill=TRUE)

# need to include all data from hes & hes_diag
new_data <- hes_diag %>%
  left_join(select(hes, record_id, epistart), by = "record_id") %>%
  select(eid, record_id, diag_icd10, epistart) %>%
  rbind(select(hes, eid, record_id, diag_icd10, epistart), by="record_id") %>% # get the union of the data
  mutate(epistart = as.Date(epistart)) %>%
  # we will lost a few patient as some of them don't have birthyear info
  merge(birthyear, by="eid") %>%
  mutate(birth_year=as.Date(paste(X34.0.0,X52.0.0,"01", sep="-"))) %>%
  mutate(age_diag = difftime(epistart, birth_year, units = "days")/365) %>%
  filter(!is.na(age_diag))

# save a survive age
survive_age <- birthyear %>% # "2016-03-01" is the date of the last records: use it as the end of the censoring
  mutate(eid = as.character(eid), X40000.0.0 = as.character(X40000.0.0)) %>%
  mutate(birth_year=as.Date(paste(X34.0.0, X52.0.0, "01", sep="-"))) %>%
  mutate(survive_date=as.Date(replace(X40000.0.0, X40000.0.0 == "", "2016-03-01"))) %>%
  mutate(survive_year = difftime(survive_date ,birth_year,units = "days")/365) %>%
  select(eid, survive_year) %>%
  filter(!is.na(survive_year))

write.csv(survive_age, paste("survive_age.csv", sep = ""), row.names = FALSE)


# mapping all ICD-10 to phecodes; need to first exclude the non-mapping (one icd10 to different phecode)
non_one2one_map <- read.csv("Phecode_map_v1_2_icd10cm_beta.csv") %>%
  group_by(phecode) %>%
  summarise(occ = n())

phecode_icd10cm <- read.csv("Phecode_map_v1_2_icd10cm_beta.csv") %>%
  left_join(non_one2one_map, by = "phecode") %>%
  group_by(icd10cm) %>%
  arrange( desc(occ), .by_group = T, ) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(ICD10 = sub("[.]","",icd10cm)) %>% # remove the dots
  select(ICD10, phecode,exclude_range, exclude_name)

short_icd10cm <- phecode_icd10cm %>%
  mutate(ICD10 = substring(ICD10, 1,4)) %>%
  left_join(non_one2one_map, by = "phecode") %>%
  group_by(ICD10) %>%
  arrange( desc(occ), .by_group = T, ) %>%
  slice(1) %>%
  ungroup() %>%
  rename(parent_phecode = phecode)

phecode_icd10 <- read.csv("phecode_icd10.csv") %>%
  left_join(non_one2one_map, by = c("PheCode" = "phecode")) %>%
  group_by(ICD10) %>%
  arrange( desc(occ), .by_group = T, ) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(ICD10 = sub("[.]","",ICD10)) %>% # remove the dots
  select(ICD10, PheCode, Excl..Phecodes, Excl..Phenotypes)

# save a phecode classifications for future plotting
read.csv("Phecode_map_v1_2_icd10cm_beta.csv") %>%
  group_by(phecode) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(ICD10 = sub("[.]","",icd10cm)) %>% # remove the dots
  select(ICD10, phecode, phecode_str, exclude_range, exclude_name) %>%
  rename(phenotype = phecode_str) %>%
  write.csv("info_phenotype_phecode.csv", row.names = F)

# map icd10 to phecodes
new_data <- new_data %>%
    select(eid, diag_icd10, age_diag) %>%
    filter(str_detect(diag_icd10, "^[A-N]")) %>%
    left_join(phecode_icd10cm, by = c("diag_icd10" = "ICD10")) %>%
    mutate(diag_icd10 = substring(diag_icd10, 1,4)) %>%
    left_join(phecode_icd10, by = c("diag_icd10" = "ICD10"))  %>%
    left_join(short_icd10cm, by = c("diag_icd10" = "ICD10"))

# mapping rate: 0.997053
# new_data %>% filter(is.na(phecode), is.na(PheCode), is.na(parent_phecode)) %>% dim

# finish the mapping
new_data <- new_data %>%
  mutate(phecode = if_else(is.na(phecode), parent_phecode, phecode)) %>%
  mutate(PheCode = if_else(is.na(PheCode), phecode, PheCode)) %>%
  filter(!is.na(PheCode)) %>%
  select(eid, PheCode, diag_icd10, age_diag)


# find all of the ICD-10 that is above 500 occurence (> 0.1%)
ds_occ_thre <- 1000 # 500
list_above500occu <- new_data %>%
  group_by(eid, PheCode) %>%
  slice(1) %>%
  group_by(PheCode) %>%
  summarise(occ = n()) %>%
  filter(occ > ds_occ_thre )

# only keep the first incidence of disease
ptm <- proc.time()
first_incidence_age <- new_data %>%
  select(eid, PheCode, age_diag) %>%
  filter(PheCode %in% list_above500occu$PheCode) %>%
  group_by(eid, PheCode) %>%
  # slice(1) %>% # just pick the first record and usually it is ranked with age
  filter(n() == 1 | age_diag == min(age_diag) ) %>% # this row is highly optimized, a lot faster the slice_min ### don't change
  slice(1) %>% # avoid the ties in min
  ungroup()
print(proc.time() - ptm) # 166.086s to run


# create another disease record where each individual has to have at least two records
individual_two_records <- first_incidence_age %>%
  group_by(eid) %>%
  filter(n() > 1)

# change all the names of PheCode to diag_icd10 to make it consistent
list_above500occu <- list_above500occu %>%
  rename(diag_icd10 = PheCode)
first_incidence_age <- first_incidence_age %>%
  rename(diag_icd10 = PheCode)
individual_two_records <- individual_two_records %>%
  rename(diag_icd10 = PheCode)

##################################
# save the data for future use
##################################
# A2N is the range of codes considered here
write.csv(first_incidence_age, paste0("DiseaseAbove",ds_occ_thre,"occur_PheCode.csv"), row.names = F)

write.csv(list_above500occu, paste0("listAbove",ds_occ_thre,"_PheCode.csv"), row.names = F)

# save the data with death event
death_age <- read.csv("death_age.csv") %>%  # include death as one of the event
  mutate(eid = as.character(eid)) %>%
  mutate(diag_icd10 = 0.01)

first_incidence_age %>%
  mutate(age_diag = as.double(age_diag)) %>%
  bind_rows(death_age) %>%
  arrange(eid) %>% # important step!
  write.csv(paste0("DiseaseAbove",ds_occ_thre,"occur_include_death_PheCode.csv"), row.names = F)

list_above500occu %>%
  add_row(diag_icd10 = 0.01, occ = dim(death_age)[1]) %>%
  write.csv(paste0("listAbove",ds_occ_thre,"include_deaths_PheCode.csv"), row.names = F)

# also save the data for individuals with at least one records
write.csv(individual_two_records, paste0("rec2subjectAbove",ds_occ_thre,"occur_PheCode.csv"), row.names = F)

# only include death events to those who have two disease records
eid_list <- individual_two_records %>%
  group_by(eid) %>%
  summarise()

death_age <- read.csv("death_age.csv") %>%  # include death as one of the event
  mutate(eid = as.character(eid)) %>%
  semi_join(eid_list, by = "eid") %>%
  mutate(diag_icd10 = 0.01)

individual_two_records %>%
  mutate(age_diag = as.double(age_diag)) %>%
  bind_rows(death_age) %>%
  arrange(eid) %>% # important step!
  write.csv(paste0("rec2subjectAbove",ds_occ_thre,"occur_include_death_PheCode.csv"), row.names = F)

#################################################
# Fig 1: schematic small figures
#################################################
K <- 10
df_P <- 5
DIR <- "~/Desktop/comorbidity/Results/"
# pt <- paste0("CVB0_model_output_A2N_age_dependent_K", K,"_P",df_P, "_rep*")
pt <- paste0("^rec2CVB0_model_output_PheCode_age_dependent_K", K,"_P",df_P, "_rep*")
temp <- list.files(paste(DIR, sep=""), pattern=pt)
lb_rep <- data_frame(df_P = as.integer(), lower_bound = as.numeric())
for(rep_id in 1:length(temp)){
  load(paste0(DIR,temp[rep_id]))
  cvrg_lb <-  model_output[[2]] %>%
    filter(!is.na(Lower_bound)) %>%
    slice_tail %>%
    pull(2)
  lb_rep <- lb_rep %>%
    add_row(df_P = df_P, lower_bound = cvrg_lb)
}
rep_id <- order(lb_rep$lower_bound, decreasing = T)[1]

load(paste0(DIR,temp[rep_id]))
### save the rec2CVB0_model_output_PheCode_age_dependent_K10_P5_rep10.RData as UKB_HES_10topics
UKB_HES_10topics <- model_output[[1]]
usethis::use_data(UKB_HES_10topics)

para$D <- model_output[[3]]
pal_age <- colorRampPalette(c(blue, red, orange, purple, green))
pal_age_vector <- pal_age(para$D)
thre_pick <- 10/para$D
phe_phecode <- read.csv("info_phenotype_phecode.csv")
ds_list <- para$list_above500occu %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") )
icd10 <- ds_list$phenotype
icd10 <- rep("", length(ds_list$phenotype))
# using the larger results file: ordering are changed in "model_output"
betas <- para$pi_beta_basis
dominant_ds_id <-  match(716.9, para$list_above500occu$diag_icd10) # hypertension 401.1 250.2 174.11 716.9
col_ds <- pal_age_vector[dominant_ds_id]
for(topic_id in 1:K){
  trajs <- betas[30:80,dominant_ds_id,topic_id, drop = F] # trajectories
  # plot_title <- paste0("Topic: ", topic_id)
  plot_title <- paste0("")
  print(paste0("Number of diseases selected: ",length(dominant_ds_id)))
  start_age <- 30
  df_topics <- data.frame(age=start_age:(start_age + dim(trajs)[1] - 1), inferred_topics =  trajs)
  plt <- ggplot(data = df_topics, aes(x = age)) +
    theme_bw(base_size = 40) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_line(aes(x = age, y = inferred_topics), color = col_ds, size = 2) +
    labs(x="Age", y="Topic loadings", title=plot_title) +
    ylim(c(0,0.22))
  ggsave(paste0("~/Desktop/comorbidity/figures/topics/rep", "Fig1_ds_", dominant_ds_id, "topic", topic_id,".png"),
         plt, width = 8, height = 8)
}



##################################
# compare collapsed VB v.s. mean-field VB
##################################
rep_number <- 10
df_P <- 5
df_lb_P_K <- data_frame(df_P = as.integer(),df_K = as.integer(), lower_bound = as.numeric(), inference = as.character())
for(K in c(5:20) ){
  lb_lst <- list()
    for(rep_id in 1:rep_number){
      try({
        load(paste0("~/Desktop/comorbidity/Results/rec2CVB0_model_output_A2N_age_dependent_K", K,"_P",df_P, "_rep",rep_id,".RData"))
        cvrg_lb <-  model_output[[2]] %>%
          filter(!is.infinite(Lower_bound)) %>%
          slice_tail %>%
          pull(2)
        df_lb_P_K  <- df_lb_P_K %>%
          add_row(df_P = df_P, df_K = K, lower_bound = cvrg_lb, inference = "CVB")
        # also load the mean-field results
        load(paste0("~/Desktop/comorbidity/Results/model_output_A2N_age_dependent_K", K,"_P",df_P, "_rep",rep_id,".RData"))
        cvrg_lb <-  model_output[[2]] %>%
          filter(!is.infinite(Lower_bound)) %>%
          slice_tail %>%
          pull(2)
        df_lb_P_K  <- df_lb_P_K %>%
          add_row(df_P = df_P, df_K = K, lower_bound = cvrg_lb, inference = "VB")
      })
    }
}

# plot a box plot for selecting best topic number and degree of freedom
df_boxplot <- df_lb_P_K %>%
  # filter(df_K  < 21, df_K > 4) %>% # help visualise one or two Ks
  mutate(df_P = as.character(df_P), df_K = as.factor(df_K)) %>%
  filter(lower_bound < 0)

df_boxplot %>% arrange(desc(lower_bound))
P_chosen <- 5
K_chosen <- 10

plt <- ggplot(data=df_boxplot,aes(x=df_K, y=lower_bound, fill = inference)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width=.8,position=position_dodge(width=0.85)) +
  geom_point(position=position_jitterdodge(jitter.width = .02, dodge.width=0.85), size = 0.3, alpha=0.6) +
  # geom_jitter(width=0.1, size = 0.2, alpha=0.4) +
  labs(x = "Number of topics", y = "Lower bound") +
  scale_fill_manual(values=cbPalette[2:7]) +
  theme(panel.background=element_blank())
ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/CVB_VB_lb_compare.png"), plt, width = 10, height = 6)


##################################
# selection of best model
##################################
# step 1: basicLDA lower bound analysis
rep_number <- 10
basicLDA_lb_K <- data_frame(df_K = as.integer(), lower_bound = as.numeric())
for(K in c(9:12) ){
  lb_lst <- list()
    for(rep_id in 1:rep_number){
      try({
        load(paste0("~/Desktop/comorbidity/Results/BaselinLDA_model_output_PheCode_K", K,"_P_rep",rep_id,".RData"))
        cvrg_lb <-  model_output[[2]] %>%
          filter(!is.infinite(Lower_bound)) %>%
          slice_tail %>%
          pull(2)
        basicLDA_lb_K  <- basicLDA_lb_K %>%
          add_row( df_K = K, lower_bound = cvrg_lb)
      })
  }
}

# plot a box plot for selecting best topic number and degree of freedom
basicLDA_lb_K <- basicLDA_lb_K %>%
  mutate( df_K = as.factor(df_K))
ggplot(data=basicLDA_lb_K,aes(x=df_K, y=lower_bound)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width=.5,position=position_dodge(width=0.85), fill = red) +
  geom_jitter(width=0.1, size = 1, alpha=0.4) +
  labs(x = "Number of topics", y = "Lower bound") +
  theme(panel.background=element_blank())

# step 1: get the best fitting results for ageLDA
rep_number <- 10
maxP <- 6
df_lb_P_K <- data_frame(df_P = as.integer(),df_K = as.integer(), lower_bound = as.numeric())
for(K in c(5:20, 25, 30, 35, 40, 45, 50) ){
  lb_lst <- list()
  for(df_P in 2:maxP){
    for(rep_id in 1:rep_number){
      try({
        load(paste0("~/Desktop/comorbidity/Results/rec2CVB0_model_output_PheCode_age_dependent_K", K,"_P",df_P, "_rep",rep_id,".RData"))
        cvrg_lb <-  model_output[[2]] %>%
          filter(!is.infinite(Lower_bound)) %>%
          slice_tail %>%
          pull(2)
        df_lb_P_K  <- df_lb_P_K %>%
          add_row(df_P = df_P, df_K = K, lower_bound = cvrg_lb)
      })
    }
  }
}

# plot a box plot for selecting best topic number and degree of freedom
df_boxplot <- df_lb_P_K %>%
  # filter(df_K  < 21, df_K > 4) %>% # help visualise one or two Ks
  mutate(df_P = as.character(df_P), df_K = as.factor(df_K))

df_boxplot %>% arrange(desc(lower_bound))
P_chosen <- 5
K_chosen <- 10

plt <- ggplot(data=df_boxplot,aes(x=df_K, y=lower_bound, fill = df_P)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width=.8,position=position_dodge(width=0.85)) +
  geom_point(position=position_jitterdodge(jitter.width = .02, dodge.width=0.85), size = 0.3, alpha=0.6) +
  # geom_jitter(width=0.1, size = 0.2, alpha=0.4) +
  labs(x = "Number of topics", y = "Lower bound") +
  scale_fill_manual(values=cbPalette[2:7]) +
  theme(panel.background=element_blank())
ggsave(paste0("~/Desktop/comorbidity/figures/CVB0_lb_change_with_topic_number_degree_freedom.png"), plt, width = 10, height = 4)
#######################################################
# prediction accuracy: using OR for top 1%, 2% or 5%
#######################################################
# Denominator of OR is the odds of choosing the disease by frequency
source("topic_functions.R")
rep_number <- 10
maxP <- 7
df_predict_lik_P_K <- data_frame(df_P = as.integer(),df_K = as.integer(),
                                 OR_top1 = as.numeric(),OR_top2 = as.numeric(),OR_top5 = as.numeric())
# compute risk when randomly choosing disease (pick topic frequency)
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
total_num <- sum(ds_list$occ)
freq_top1 <- ds_list %>%
  arrange(desc(occ)) %>%
  slice(1:floor(para$D/100)) %>%
  pull(occ) %>%
  sum
freq_top2 <- ds_list %>%
  arrange(desc(occ)) %>%
  slice(1:floor(para$D/50)) %>%
  pull(occ) %>%
  sum
freq_top5 <- ds_list %>%
  arrange(desc(occ)) %>%
  slice(1:floor(para$D/20)) %>%
  pull(occ) %>%
  sum
Odds_freq_list <- c(freq_top1, freq_top2, freq_top5)/total_num
Odds_freq_list <- Odds_freq_list/(1-Odds_freq_list)
for(K in 5:20){
  lb_lst <- list()
  for(df_P in 2:maxP){
    for(rep_id in 1:rep_number){
      try({
        load(paste0("../Results/prediction_onebyone_age_K",K,"_P",df_P,"_rep",rep_id, ".RData"))
        df_predict_lik_P_K  <- df_predict_lik_P_K %>%
          add_row(df_P = df_P, df_K = K,
                  OR_top1 = (prediction_onebyone_rslt[[1]][[2]]/(1 - prediction_onebyone_rslt[[1]][[2]]))/Odds_freq_list[1],
                  OR_top2 = (prediction_onebyone_rslt[[1]][[3]]/(1 - prediction_onebyone_rslt[[1]][[3]]))/Odds_freq_list[2],
                  OR_top5 = (prediction_onebyone_rslt[[1]][[4]]/(1 - prediction_onebyone_rslt[[1]][[4]]))/Odds_freq_list[3],
                  )
      })
    }
  }
}

# plot a box plot for selecting best topic number and degree of freedom
df_boxplot <- df_predict_lik_P_K %>%
  filter(df_K  < 21, df_K > 4) %>% # help visualise one or two Ks
  mutate(df_P = as.character(df_P), df_K = as.factor(df_K))
plt <- ggplot(data=df_boxplot,aes(x=df_K, y=OR_top1, fill = df_P)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width=.8,position=position_dodge(width=0.85)) +
  geom_point(position=position_jitterdodge(jitter.width = .02, dodge.width=0.85), size = 0.3, alpha=0.6) +
  # geom_jitter(width=0.1, size = 0.2, alpha=0.4) +
  labs(x = "Number of topics", y = "Predicting Odds Ratio") +
  scale_fill_manual(name = "d.f.",values=cbPalette[2:7]) +
  theme(panel.background=element_blank())
ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/prediction_onebyone_OR_topic_number_degree_freedom.png"), plt, width = 10, height = 6)

# prediction OR: campare ageLDA v.s. basic LDA
rep_number <- 10
df_P <- 5
df_predict_lik_P_K <- data_frame(df_P = as.integer(),df_K = as.integer(),
                                 OR_top1 = as.numeric(),OR_top2 = as.numeric(),OR_top5 = as.numeric(), model = as.character())
# compute risk when randomly choosing disease (pick topic frequency)
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
total_num <- sum(ds_list$occ)
freq_top1 <- ds_list %>%
  arrange(desc(occ)) %>%
  slice(1:floor(para$D/100)) %>%
  pull(occ) %>%
  sum
freq_top2 <- ds_list %>%
  arrange(desc(occ)) %>%
  slice(1:floor(para$D/50)) %>%
  pull(occ) %>%
  sum
freq_top5 <- ds_list %>%
  arrange(desc(occ)) %>%
  slice(1:floor(para$D/20)) %>%
  pull(occ) %>%
  sum
Odds_freq_list <- c(freq_top1, freq_top2, freq_top5)/total_num
Odds_freq_list <- Odds_freq_list/(1-Odds_freq_list)
for(K in 5:20){
  lb_lst <- list()
    for(rep_id in 1:rep_number){
      try({
        load(paste0("../Results/prediction_onebyone_age_K",K,"_P",df_P,"_rep",rep_id, ".RData"))
        df_predict_lik_P_K  <- df_predict_lik_P_K %>%
          add_row(df_P = df_P, df_K = K,
                  OR_top1 = (prediction_onebyone_rslt[[1]][[2]]/(1 - prediction_onebyone_rslt[[1]][[2]]))/Odds_freq_list[1],
                  OR_top2 = (prediction_onebyone_rslt[[1]][[3]]/(1 - prediction_onebyone_rslt[[1]][[3]]))/Odds_freq_list[2],
                  OR_top5 = (prediction_onebyone_rslt[[1]][[4]]/(1 - prediction_onebyone_rslt[[1]][[4]]))/Odds_freq_list[3],
                  model = "ATM"
          )
        # basic LDA
        load(paste0("../Results/prediction_onebyone_baseline_K",K,"_rep",rep_id, ".RData"))
        df_predict_lik_P_K  <- df_predict_lik_P_K %>%
          add_row(df_P = df_P, df_K = K,
                  OR_top1 = (prediction_onebyone_rslt[[1]][[2]]/(1 - prediction_onebyone_rslt[[1]][[2]]))/Odds_freq_list[1],
                  OR_top2 = (prediction_onebyone_rslt[[1]][[3]]/(1 - prediction_onebyone_rslt[[1]][[3]]))/Odds_freq_list[2],
                  OR_top5 = (prediction_onebyone_rslt[[1]][[4]]/(1 - prediction_onebyone_rslt[[1]][[4]]))/Odds_freq_list[3],
                  model = "LDA"
          )
      })
    }
}

# plot a box plot for selecting best topic number and degree of freedom
df_boxplot <- df_predict_lik_P_K %>%
  filter(df_K  < 21, df_K > 4) %>% # help visualise one or two Ks
  mutate(df_P = as.character(df_P), df_K = as.factor(df_K))
plt <- ggplot(data=df_boxplot,aes(x=df_K, y=OR_top1, fill = model)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width=.8,position=position_dodge(width=0.85)) +
  geom_point(position=position_jitterdodge(jitter.width = .02, dodge.width=0.85), size = 0.3, alpha=0.6) +
  # geom_jitter(width=0.1, size = 0.2, alpha=0.4) +
  labs(x = "Number of topics", y = "Predicting Odds Ratio") +
  scale_fill_manual(name = "Models",values=cbPalette[2:7]) +
  theme(panel.background=element_blank())
ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/OR_ageLDA_LDA_comparison.png"), plt, width = 10, height = 6)



############################################
# plotting record number per individual
############################################
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
first_incidence_age <- rec_data %>%
  arrange(eid)

# compute the age span of individual
survive_age <- read.csv(paste("survive_age.csv", sep = ""))
age_span <- rec_data %>%
  group_by(eid) %>%
  summarise(min_age = min(age_diag)) %>%
  left_join(survive_age, by = "eid") %>%
  mutate(age_span = survive_year - min_age)
mean(age_span$min_age)
mean(age_span$survive_year)
mean(age_span$age_span)

# plot the number distribution of indiviudal diseases
df_number_records <- first_incidence_age %>%
  group_by(eid) %>%
  summarise(n())
df_simu_pois <- data.frame(num_records = floor(2+rexp(282957, 1/(-1.5+6.1))))
ggplot(df_number_records) +
  geom_histogram(aes(x = `n()`, fill = "true"), alpha = 1, binwidth = 1) +
  geom_histogram(data = df_simu_pois, aes(x = num_records, fill = "simulated"), alpha = .5,, binwidth = 1) +
  lims(x = c(0,40)) +
  scale_fill_manual(values = c("true" = grey, "simulated" = red)) +
  theme(legend.position = c(.8,.8),panel.background=element_blank())

# print the disease distribution over age
onset_by_year <- first_icidence_age %>%
  mutate(age_diag = floor(age_diag)) %>%
  group_by(age_diag) %>%
  summarise(record_per_age = n())

ggplot(data = onset_by_year) +
  geom_line(aes(x = age_diag, y = record_per_age))

# sd of diagnosis for each disease
std_age_ds <-rec_data %>%
  group_by(diag_icd10) %>%
  summarise(std_age_ds = sd(age_diag))

##################################
# Step 2: visualise all the topics
##################################

##################################
# Figure 3: plot topic weights based on disease assignment to topic
# reorder the topic by Phecodes
##################################

load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
loadings_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[para$ds_list[[j]]$id,])) %>% t
# do it separately for younger and older
younger60 <- lapply(para$ds_list, function(x) filter(x, age_diag <=60))
loadings_young_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[younger60[[j]]$id,])) %>% t
older60 <- lapply(para$ds_list, function(x) filter(x, age_diag >60))
loadings_old_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[older60[[j]]$id,,drop=FALSE])) %>% t
loadings_old_per_ds[is.na(loadings_old_per_ds)] <- 0

phe_phecode <- read.csv("info_phenotype_phecode.csv")
ds.system <- para$list_above500occu %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") )
# combine the last four systems to other:congenital anomalies, symptoms, injuries & poisonings, <NA>
ds.system <- ds.system %>%
  mutate(exclude_name = if_else(exclude_name %in% c("congenital anomalies", "symptoms", "injuries & poisonings","", NA), "others", exclude_name))
ds.system <- ds.system %>%
  mutate(exclude_name = factor(exclude_name, levels = unique(ds.system$exclude_name) ))
####################################################
# Figure 3: Topic overview: topic weight distribution
####################################################
systems_associated <- list()
for(i in 1:para$K){
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

loadings_age_diff <- matrix(NA, nrow = para$D, ncol = 2*para$K)
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
  mutate(Var1 = factor(Var1, levels = rev(ds.system$diag_icd10)))
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
df_segments <- data.frame(x = 0:10 * 2 + 0.5, xend =  0:10 * 2+ 0.5, y = rep(0,11), yend = rep(para$D,11))
plt <- plt + geom_segment(data=df_segments, aes(x,y,xend=xend, yend=yend), size=.5, alpha = 0.3, inherit.aes=F)

# these code could be used to add rectangular to the heatmap
# df_rect <- data.frame(xmin = as.numeric(), xmax = as.numeric(), ymin = as.numeric(), ymax = as.numeric() )
# for(i in 1:para$K){
#   for(j in 1:length(systems_associated[[i]])){
#     type <- systems_associated[[i]][j]
#     if(!is.na(type) & length(which(rev(ds.system$exclude_name)== type)) > 10){
#       df_rect <- df_rect %>%
#         add_row(xmin = i - 0.48, xmax = i + 0.48,
#                 ymin = which(rev(ds.system$exclude_name) == type)[1],
#                 ymax = which(rev(ds.system$exclude_name) == type)[length(which(rev(ds.system$exclude_name)== type))])
#     }
#
#   }
# }
#
# plt <- plt + geom_rect(data = df_rect,
#                        aes(xmin = xmin,
#                            xmax = xmax,
#                            ymin = ymin,
#                            ymax = ymax),
#                        fill = "transparent", color = blue, size = 1)


ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/topic_disease_highlight.png",plt, width = 6, height = 10)

# making the color legends
longData<-melt(matrix(c(1:100/100, 1:100/100),ncol = 2))
df_subplot <- longData %>%
  mutate(group_range = as.factor(Var2 %% 2))

plt <- ggplot() +
  geom_tile(data = df_subplot, aes(x = Var2, y = Var1, fill=group_range, alpha = value,width = 0.9)) +
  scale_alpha_continuous(range = c(0,1)) +
  scale_fill_manual(values = c("0" = blue, "1" = red))+
  labs(x="", y="", title="") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/color_legends.png"), plt, width = 2, height = 5)

# save the data table
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
phe_phecode <- read.csv("info_phenotype_phecode.csv")
ds.system <- para$list_above500occu %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") ) %>%
  select(-exclude_range) %>%
  rename(Category=exclude_name)
ds.system$loadings <- loadings_age_diff

ds.system %>%
  write.csv("~/Desktop/comorbidity/paper_writing/supplementary_files/Figure3DataTable.csv", row.names = F)
# save colored x for making figures
categories <- ds.system$exclude_name
categories[346] <- categories[347] # fixing one error in Phecode coding
longData<-melt(as.matrix(categories)) %>%
  mutate(Var1 = factor(ds.system$diag_icd10[Var1], levels = rev(ds.system$diag_icd10)))
cols <- setNames(rep(c(red, blue), 9)[1:18], unique(ds.system$exclude_name))
plt_pallete <- ggplot() +
  geom_tile(data = longData, aes(x = Var2, y = Var1, fill=value), alpha = 0.5) +
  scale_fill_manual(values = cols) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/pallete.png",plt_pallete, width = 1, height = 10)

# for figure 3 panel B save specific part of the graph
longData<-melt(loadings_per_ds[, order_ds]) %>%
  mutate(Var1 =ds.system$diag_icd10[Var1])
ds.groups <- ds.system %>%
  select(diag_icd10, exclude_name) %>%
  mutate(Var1 = diag_icd10, group_range = (as.integer(exclude_name) %% 2)) %>%
  mutate(group_range = ifelse(is.na(group_range), 1, group_range)) %>%
  mutate(group_range = as.factor(group_range))
longData <- longData %>%
  left_join(ds.groups, by = c("Var1")) %>%
  mutate(Var1 = factor(Var1, levels = rev(ds.system$diag_icd10)))
# filter down to only two topics and a few systems
topic_subset <- c(7,8)
topics_of_interest <- order_ds[topic_subset]
systems_of_interest <- c()
for(i in topics_of_interest){
  ds.system$loadings <- loadings_per_ds[, i]
  system_subset <- ds.system %>%
    group_by(exclude_name) %>%
    summarise(mean_loading = mean(loadings), sum_loading = sum(loadings)) %>%
    filter(sum_loading > 5) %>%
    pull(exclude_name)
  systems_of_interest <- c(systems_of_interest, as.character(system_subset))
}
df_subplot <- longData %>%
  filter(exclude_name %in% systems_of_interest, Var2 %in% topic_subset)
df_subplot <- df_subplot %>%
  mutate(exclude_name = factor(exclude_name, levels = unique(as.character(df_subplot$exclude_name)))) %>%
  mutate(group_range = (as.integer(exclude_name) %% 2)) %>%
  mutate(group_range = as.factor(group_range), Var2 = as.factor(Var2))

plt <- ggplot() +
  geom_tile(data = df_subplot, aes(x = Var2, y = Var1, fill=group_range, alpha = value,width = 0.9)) +
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
ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/topic",topic_subset[1],topic_subset[2],"subset.png"), plt, width = 2, height = 5)

categories <- ds.system$exclude_name
categories[346] <- categories[347] # fixing one error in Phecode coding
longData<-melt(as.matrix(categories)) %>%
  mutate(Var1 = factor(ds.system$diag_icd10[Var1], levels = rev(ds.system$diag_icd10))) %>%
  filter(value %in% systems_of_interest)
cols <- setNames(rep(c(red, blue), 9)[1:length(unique(df_subplot$exclude_name))], unique(as.character(df_subplot$exclude_name)))
plt_pallete <- ggplot() +
  geom_tile(data = longData, aes(x = Var2, y = Var1, fill=value), alpha = 0.5) +
  scale_fill_manual(values = cols) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures//pallete",topic_subset[1],topic_subset[2],"subset.png") ,plt_pallete, width = 1, height = 5)

# save legneds:
plt_legend <- ggplot() +
  geom_tile(data = longData, aes(x = Var2, y = Var1, fill=group_range, alpha = value,width = 0.9))  + theme_bw()
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/legends_topic_distribution.png",plt_legend, width = 6, height = 10)

# save topic properties for acronym
info_topics <- data.frame(topic_id = order_ds, systems = rep(0, length(order_ds)))
for(i in 1:length(order_ds)){
  info_topics$systems[i] <- paste(as.character(systems_associated[[i]]), collapse = ", ")
}
# topic association analysis with BMI, sex, survive age
birthyear <- read.csv("~/Desktop/genetics_longitudinal_data/longitudinal_data/Year_of_birth.csv") %>%
  rename(sex = X31.0.0, BMI = X23104.0.0, birth_year = X34.0.0, ds_num = X52.0.0) %>%
  select(eid, sex, BMI,birth_year, ds_num)
df_birthyear <- data.frame(eid = para$eid) %>%
  left_join(birthyear,  by = "eid")
patient_loadings <- sweep((para$alpha_z - 1), 1, rowSums(para$alpha_z -1), FUN="/")
BMI_r2 <- rep(NA, para$K)
BMI_coef <- rep(NA, para$K)
sex_r2 <- rep(NA, para$K)
sex_coef <- rep(NA, para$K)
birthyear_r2 <- rep(NA, para$K)
birthyear_coef <- rep(NA, para$K)
dsnum_r2 <- rep(NA, para$K)
dsnum_coef <- rep(NA, para$K)
for(i in 1:para$K){
  fit.bmi <- lm(patient_loadings[, order_ds[i]] ~ df_birthyear$BMI)
  BMI_r2[i] <- summary(fit.bmi)$r.squared * 100
  BMI_coef[i] <- summary(fit.bmi)$coefficients[2,1]

  fit.sex <- lm(patient_loadings[, order_ds[i]] ~ df_birthyear$sex)
  sex_r2[i] <- summary(fit.sex)$r.squared * 100
  sex_coef[i] <- summary(fit.sex)$coefficients[2,1]

  fit.birthyear <- lm(patient_loadings[, order_ds[i]] ~ df_birthyear$birth_year)
  birthyear_r2[i] <- summary(fit.birthyear)$r.squared * 100
  birthyear_coef [i] <- summary(fit.birthyear)$coefficients[2,1]

  fit.ds_num <- lm(patient_loadings[, order_ds[i]] ~ df_birthyear$ds_num)
  dsnum_r2[i] <- summary(fit.ds_num)$r.squared * 100
  dsnum_coef[i] <- summary(fit.ds_num)$coefficients[2,1]
}
info_topics$BMI_r2 <- BMI_r2
info_topics$BMI_coef <- BMI_coef
info_topics$sex_r2 <- sex_r2
info_topics$sex_coef <- sex_coef
info_topics$birthyear_r2 <- birthyear_r2
info_topics$birthyear_coef <- birthyear_coef

# save the top three systems for each disease
systems_top <- list()
disease_top <- list()
num_disease <- c()
for(i in 1:para$K){
  ds.system$loadings <- loadings_per_ds[, i]
  systems_top[[i]] <- ds.system %>%
    group_by(exclude_name) %>%
    summarise(mean_loading = mean(loadings), sum_loading = sum(loadings)) %>%
    dplyr::arrange(desc(sum_loading)) %>%
    slice(1:3) %>%
    pull(exclude_name)
  # get the top three loadings
  # filter(mean_loading > .5 | sum_loading > 10) %>%

  # step 2, get the top 5 diseases
  disease_top[[i]] <- ds.system %>%
    filter(occ > 5000, loadings > 0.3) %>%
    dplyr::arrange(desc(loadings)) %>%
    slice(1:10) %>%
    pull(phenotype)
  # step 3: compute the number of diseases associated with each topic
  num_disease[i] <- ds.system %>%
    filter(loadings > 0.5) %>%
    pull(diag_icd10)  %>%
    length()
}

systems_top <- systems_top[order_ds]
disease_top <- disease_top[order_ds]
info_topics$top_systems <- rep(0, length(order_ds))
info_topics$disease_top <- rep(0, length(order_ds))
for(i in 1:length(order_ds)){
  info_topics$top_systems[i] <- paste(as.character(systems_top[[i]]), collapse = ", ")
  info_topics$disease_top[i] <- paste(as.character(disease_top[[i]]), collapse = "; ")
}

info_topics$num_disease <- num_disease[order_ds]

# average topic weight
info_topics$avg_topic_weight_per_disease <- colMeans(para$unlist_zn)[order_ds]
patient_loadings <- sweep((para$alpha_z-1), 1, rowSums(para$alpha_z - 1), FUN="/")
info_topics$avg_topic_weight_per_individual <- colMeans(patient_loadings)[order_ds]
# weighted average of age
topic_age <- c()
for(i in 1:para$K){
  topic_age[i] <- (para$unlist_zn[,i] %*% para$unlist_Ds_id$age_diag)/sum(para$unlist_zn[,i])
}
info_topics$topic_age <- topic_age[order_ds]

# downloading the heritability
h2g_topic <- read.csv("topic_h2g.csv")
pasted_rslt <- matrix(mapply(function(x,y) paste0(as.character(x), " (s.e.=", as.character(y), ")"),
                             round(h2g_topic$h2g,digits = 3), round(h2g_topic$seh2g, digits = 3)), para$K,ncol = 1)

info_topics$h2g <- pasted_rslt[order_ds]

write.csv(info_topics, file = "~/Desktop/comorbidity/paper_writing/Production_figures/Topic_information.csv", row.names = F)
############### save this into a csv file

##################################################
# only compute the weights within topics (alternative, not chosen)
##################################################
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
loadings_per_ds <- sapply(1:para$D, function(j)  colMeans(para$pi_beta_basis[31:80,j,]))%>% t
loadings_per_ds <- apply(loadings_per_ds, 2, function(x) x/max(x))
# do it separately for younger and older
loadings_young_per_ds <- sapply(1:para$D, function(j)  colMeans(para$pi_beta_basis[31:60,j,]))%>% t
loadings_young_per_ds <- apply(loadings_young_per_ds, 2, function(x) x/max(x))
loadings_old_per_ds <- sapply(1:para$D, function(j)  colMeans(para$pi_beta_basis[61:80,j,]))%>% t
loadings_old_per_ds <- apply(loadings_old_per_ds, 2, function(x) x/max(x))

phe_phecode <- read.csv("info_phenotype_phecode.csv")
ds.system <- para$list_above500occu %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") )
ds.system <- ds.system %>%
  mutate(exclude_name = factor(exclude_name, levels = unique(ds.system$exclude_name) ))

systems_associated <- list()
for(i in 1:para$K){
  ds.system$loadings <- loadings_per_ds[, i]
  systems_associated[[i]] <- ds.system %>%
    group_by(exclude_name) %>%
    summarise(mean_loading = mean(loadings), sum_loading = sum(loadings)) %>%
    filter(mean_loading > .1 | sum_loading > 5) %>%
    pull(exclude_name)
}
# order the topics by the Phecode structure
order_ds <- order(sapply(systems_associated, function(x) x[1]))
systems_associated <- systems_associated[order_ds]
ds.system$loadings <- loadings_per_ds[, order_ds]

loadings_age_diff <- matrix(NA, nrow = para$D, ncol = 2*para$K)
for(i in 1:length(order_ds)){
  loadings_age_diff[,2*i -1] <- loadings_young_per_ds[,order_ds[i]]
  loadings_age_diff[,2*i] <- loadings_old_per_ds[,order_ds[i]]
}

longData<-melt(loadings_per_ds[, order_ds]) %>%
  mutate(Var1 =ds.system$diag_icd10[Var1])

longData<-melt(loadings_age_diff) %>%
  mutate(Var1 =ds.system$diag_icd10[Var1])

ds.groups <- ds.system %>%
  select(diag_icd10, exclude_name) %>%
  mutate(Var1 = diag_icd10, group_range = (as.integer(exclude_name) %% 2)) %>%
  mutate(group_range = ifelse(is.na(group_range), 1, group_range)) %>%
  mutate(group_range = as.factor(group_range))
longData <- longData %>%
  left_join(ds.groups, by = c("Var1")) %>%
  mutate(Var1 = factor(Var1, levels = rev(ds.system$diag_icd10)))
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

ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/topic_compression.png",plt, width = 6, height = 10)

###############################################
# Figure 3: Topic overview: save topics (not normalised)
################################################
# first find the best rep
K <- 10
df_P <- 5
DIR <- "~/Desktop/comorbidity/Results/"
# pt <- paste0("CVB0_model_output_A2N_age_dependent_K", K,"_P",df_P, "_rep*")
pt <- paste0("^rec2CVB0_model_output_PheCode_age_dependent_K", K,"_P",df_P, "_rep*")
temp <- list.files(paste(DIR, sep=""), pattern=pt)
lb_rep <- data_frame(df_P = as.integer(), lower_bound = as.numeric())
for(rep_id in 1:length(temp)){
  load(paste0(DIR,temp[rep_id]))
  cvrg_lb <-  model_output[[2]] %>%
    filter(!is.na(Lower_bound)) %>%
    slice_tail %>%
    pull(2)
  lb_rep <- lb_rep %>%
    add_row(df_P = df_P, lower_bound = cvrg_lb)
}
rep_id <- order(lb_rep$lower_bound, decreasing = T)[1]

# load(paste0(DIR,temp[rep_id]))

load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
para$D <- model_output[[3]]
pal_age <- colorRampPalette(c(blue, red, orange, purple, green))
pal_age_vector <- pal_age(para$D)
thre_pick <- 10/para$D
phe_phecode <- read.csv("info_phenotype_phecode.csv")
ds_list <- para$list_above500occu %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") )
icd10 <- ds_list$phenotype
# icd10 <- rep("", length(ds_list$phenotype))
# using the larger results file: ordering are changed in "model_output"
# betas <- model_output[[1]]
betas <- para$pi_beta_basis
for(topic_id in 1:K){
  trajs <- betas[30:80,,topic_id] # trajectories
  # plot_title <- paste0("Topic: ", topic_id)
  plot_title <- paste0("")
  plt <- plot_age_topics(icd10, trajs, pal_age_vector, plot_title, start_age = 30,top_ds = 7)
  ggsave(paste0("~/Desktop/comorbidity/figures/topics/rep", 10, "K",K,"P",df_P,"age_topics",topic_id,".png"),
         plt, width = 8, height = 8)
}

# save the topic loadings for all topics
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
ds.system <- para$list_above500occu %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") ) %>%
  select(-exclude_range, -occ, -exclude_name) %>%
  rename(Phecode = diag_icd10, Description=phenotype)
topic_loadings <- list()
for(topic_id in 1:para$K){
  trajs <- betas[30:80,,topic_id] %>%
    t() %>%
    data.frame()# trajectories
  names(trajs) <- 30:80
  topic_loadings[[topic_id]] <- ds.system %>%
    mutate(Topic = topic_name[topic_id]) %>%
    select(Topic, Phecode, ICD10, Description) %>%
    bind_cols(trajs)
}
topic_loadings %>% bind_rows() %>%
  write.csv("~/Desktop/comorbidity/paper_writing/supplementary_files/Topic_loadings.csv", row.names = F)

###############################################
# Figure 3+: normalised topic weights
################################################
# first find the best rep
K <- 10
df_P <- 5
DIR <- "~/Desktop/comorbidity/Results/"
# pt <- paste0("CVB0_model_output_A2N_age_dependent_K", K,"_P",df_P, "_rep*")
pt <- paste0("^rec2CVB0_model_output_PheCode_age_dependent_K", K,"_P",df_P, "_rep*")
temp <- list.files(paste(DIR, sep=""), pattern=pt)
lb_rep <- data_frame(df_P = as.integer(), lower_bound = as.numeric())
for(rep_id in 1:length(temp)){
  load(paste0(DIR,temp[rep_id]))
  cvrg_lb <-  model_output[[2]] %>%
    filter(!is.na(Lower_bound)) %>%
    slice_tail %>%
    pull(2)
  lb_rep <- lb_rep %>%
    add_row(df_P = df_P, lower_bound = cvrg_lb)
}
rep_id <- order(lb_rep$lower_bound, decreasing = T)[1]

# load(paste0(DIR,temp[rep_id]))

load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
para$D <- model_output[[3]]
pal_age <- colorRampPalette(c(blue, red, orange, purple, green))
pal_age_vector <- pal_age(para$D)
phe_phecode <- read.csv("info_phenotype_phecode.csv")
ds_list <- para$list_above500occu %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") )
icd10 <- ds_list$phenotype
# icd10 <- rep("", length(ds_list$phenotype))
# using the larger results file: ordering are changed in "model_output"
# betas <- model_output[[1]]
betas <- para$pi_beta_basis

loadings_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[para$ds_list[[j]]$id,])) %>% t
disease_show <- 7
for(topic_id in 1:K){
  ds_list <- order(loadings_per_ds[,topic_id], decreasing = T)[1:disease_show]
  trajs <- betas[30:80, ,topic_id] # trajectories
  # plot_title <- paste0("Topic: ", topic_id)
  plot_title <- paste0("")
  plt <- plot_age_topics_specify_id(icd10, trajs, ds_list, pal_age_vector, plot_title, start_age = 30)
  ggsave(paste0("~/Desktop/comorbidity/figures/topics/rep", 10, "K",K,"P",df_P,"normalised_age_topics",topic_id,".png"),
         plt, width = 8, height = 8)
}


##########################################################
# supplementary figure: male/female topics
##########################################################
# compare topics loading between male/female
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
loadings_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[para$ds_list[[j]]$id,])) %>% t

# first order the disease topics
phe_phecode <- read.csv("info_phenotype_phecode.csv")
ds.system <- para$list_above500occu %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") )
ds.system <- ds.system %>%
  mutate(exclude_name = factor(exclude_name, levels = unique(ds.system$exclude_name) ))
systems_associated <- list()
for(i in 1:para$K){
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

# load male and female data
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
keep_woman <- read.table(file="keep_women.txt", header=FALSE, sep=" ") %>%
  mutate(eid = V1)
male_data <- rec_data %>%
  anti_join(keep_woman, by = "eid")
female_data <- rec_data %>%
  semi_join(keep_woman, by = "eid")
# extract individual weights
K <- 10
df_P <- 5
DIR <- "~/Desktop/comorbidity/Results/"

load(paste0(DIR,paste0("female_best_output_AgeLDA_RunNumber5K", K,"_P",df_P, "_rep1.RData")))
female_topic <- model_output[[1]]
load(paste0(DIR,paste0("male_best_output_AgeLDA_RunNumber5K", K,"_P",df_P, "_rep1.RData")))
male_topic <- model_output[[1]]

para_female <- topics2weights(female_data, ds_list, degree_freedom = df_P , topics = female_topic)
para_male <- topics2weights(male_data, ds_list, degree_freedom = df_P , topics = male_topic)

loadings_female_perds <- sapply(1:para_female$D, function(j)
  colMeans(para_female$unlist_zn[para_female$ds_list[[j]]$id,,drop = F]) ) %>% t
loadings_male_perds <-  sapply(1:para_male$D, function(j)
  colMeans(para_male$unlist_zn[para_male$ds_list[[j]]$id,,drop = F])) %>% t

# order the topics; adjust based on the values
chosen_female <- c()
corr_female <- c()
for(topic_id in order_ds){
  topic_ids <- sapply(1:K, function(j) cor(loadings_per_ds[, topic_id], loadings_female_perds[,j]))
  print(topic_ids)
  chosen_female <- c(chosen_female, which.max(topic_ids))
  corr_female <- c(corr_female, max(topic_ids))
}
chosen_female[4] <- 3
mean(corr_female) # 0.788

chosen_male <- c()
corr_male <- c()
for(topic_id in order_ds){
  topic_ids <- sapply(1:K, function(j) cor(loadings_per_ds[, topic_id], loadings_male_perds[,j]))
  print(topic_ids)
  chosen_male <- c(chosen_male, which.max(topic_ids))
  corr_male <- c(corr_male, max(topic_ids))
}
chosen_male[9] <- 2
mean(corr_male) # 0.773

loadings_sex_diff <- matrix(NA, nrow = para$D, ncol = 2*para$K)
for(i in 1:length(order_ds)){
  loadings_sex_diff[,2*i -1] <- loadings_female_perds[,chosen_female[i]]
  loadings_sex_diff[,2*i] <- loadings_male_perds[,chosen_male[i]]
}

longData<-melt(loadings_sex_diff) %>%
  mutate(Var1 =ds.system$diag_icd10[Var1])

ds.groups <- ds.system %>%
  select(diag_icd10, exclude_name) %>%
  mutate(Var1 = diag_icd10, group_range = (as.integer(exclude_name) %% 2)) %>%
  mutate(group_range = ifelse(is.na(group_range), 1, group_range)) %>%
  mutate(group_range = as.factor(group_range))
longData <- longData %>%
  left_join(ds.groups, by = c("Var1")) %>%
  mutate(Var1 = factor(Var1, levels = rev(ds.system$diag_icd10)))
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
df_segments <- data.frame(x = 0:10 * 2 + 0.5, xend =  0:10 * 2+ 0.5, y = rep(0,11), yend = rep(para$D,11))
plt <- plt + geom_segment(data=df_segments, aes(x,y,xend=xend, yend=yend), size=.5, alpha = 0.3, inherit.aes=F)
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/topic_sex_difference.png",plt, width = 6, height = 10)

##################################
# save age topics for each gender
K <- 10
df_P <- 5
DIR <- "~/Desktop/comorbidity/Results/"
load(paste0(DIR,paste0("female_best_output_AgeLDA_RunNumber5K", K,"_P",df_P, "_rep1.RData")))

para$D <- model_output[[3]]
pal_age <- colorRampPalette(c(blue, red, orange, purple, green))
pal_age_vector <- pal_age(para$D)
thre_pick <- 10/para$D
phe_phecode <- read.csv("info_phenotype_phecode.csv")
ds_list <- para$list_above500occu %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") )
icd10 <- ds_list$phenotype
# icd10 <- rep("", length(ds_list$phenotype))
# using the larger results file: ordering are changed in "model_output"
betas <- model_output[[1]]
for(topic_id in 1:K){
  trajs <- betas[30:79,,topic_id] # trajectories
  # plot_title <- paste0("Topic: ", topic_id)
  plot_title <- paste0("")
  plt <- plot_age_topics(icd10, trajs, pal_age_vector, plot_title, start_age = 30)
  ggsave(paste0("~/Desktop/comorbidity/figures/topics/female_K",K,"P",df_P,"age_topics",topic_id,".png"),
         plt, width = 8, height = 8)
}


##################################
# topic sparsity histogram
##################################
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
loadings_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[para$ds_list[[j]]$id,])) %>% t
top_1_value <- sapply(1:dim(loadings_per_ds)[1], function(x) sort(loadings_per_ds[x,], decreasing = T)) %>% t
df_box <- melt(top_1_value) %>%
  rename(topic_order = Var2, weights = value, record_number = Var1) %>%
  mutate(topic_order = factor(topic_order, levels = 1:10))

a <- which(para$list_above500occu$occ <= 2000)
df_box$record_number[which(df_box$record_number %in% a)] <- "less_2000"
a <- which(para$list_above500occu$occ > 2000 & para$list_above500occu$occ  <= 5000)
df_box$record_number[which(df_box$record_number %in% a)] <- "2000_to_5000"
a <- which(para$list_above500occu$occ > 5000)
df_box$record_number[which(df_box$record_number %in% a)] <- "more_5000"
df_box <- df_box %>%
  mutate(record_number = factor(record_number, levels = c("less_2000", "2000_to_5000", "more_5000")))

plt <- ggplot(data=df_box) +
  geom_boxplot(mapping=aes(x=topic_order, y=weights, fill = record_number),  width=0.8, alpha = 0.6, outlier.shape = NA) +
  # geom_jitter(mapping=aes(x=topic_order, y=weights),width=0.2, size = 0.1, alpha=0.4) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=c("less_2000" = blue, "2000_to_5000" = red, "more_5000" = green)) +
  labs(x = "Topic value rank", y = "Disease topic weight") +
  theme_bw(base_size = 20) +
  theme(legend.position = "None",panel.background=element_blank()) +
  geom_hline(aes(yintercept = 1), color = "red", linetype = "dashed")
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/disease_sparsity_incidence_number.png", plt, width = 5, height = 5)
# df_hist <- data.frame(top_topic = top_1_value[,1],  second_topic = top_1_value[,2]) %>%
#   pivot_longer(cols = everything(), names_to = "type", values_to = "weights")
# plt <- ggplot(filter(df_hist, type == "top_topic"), aes(x=weights, fill=type)) +
#   geom_histogram( aes(y = (..count..)/sum(..count..)), color=red, alpha=0.6, center = 0.025, binwidth = 0.05) +
#   scale_y_continuous(labels = scales::percent) +
#   scale_x_continuous(limits = c(0,1),breaks = c(0, 0.2,0.4,0.6,0.8, 1.0)) +
#   theme_bw(base_size = 20) +
#   labs(x = "Largest topic value", y = "Density") +
#   theme(legend.position = "none")
plt_all <- ggplot(data=df_box) +
  geom_boxplot(mapping=aes(x=topic_order, y=weights, fill = red), color = red,   width=0.5, alpha = 0.6, outlier.shape = NA) +
  # geom_jitter(mapping=aes(x=topic_order, y=weights),width=0.2, size = 0.1, alpha=0.4) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Topic value rank", y = "Disease topic weight") +
  theme_bw(base_size = 20) +
  theme(legend.position = "None",panel.background=element_blank()) +
  geom_hline(aes(yintercept = 1), color = "red", linetype = "dashed")
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/fig3_disease_sparsity.png", plt_all, width = 5, height = 5)
print(paste0("proportion of max topic value > 0.95: ", mean(top_1_value[,1] > 0.95) ) )


patient_loadings <- sweep((para$alpha_z-1), 1, rowSums(para$alpha_z - 1), FUN="/")
top_1_value <- sapply(1:dim(patient_loadings)[1], function(x) sort(patient_loadings[x,], decreasing = T)) %>% t
df_box <- melt(top_1_value) %>%
  rename(topic_order = Var2, weights = value, record_number = Var1) %>%
  mutate(topic_order = factor(topic_order, levels = 1:10))
a <- which(para$Ns <= 5)
df_box$record_number[which(df_box$record_number %in% a)] <- "less_5"
a <- which(para$Ns > 5 & para$Ns <= 10)
df_box$record_number[which(df_box$record_number %in% a)] <- "5_to_10"
a <- which(para$Ns >  10)
df_box$record_number[which(df_box$record_number %in% a)] <- "more_10"
df_box <- df_box %>%
  mutate(record_number = factor(record_number, levels = c("less_5", "5_to_10", "more_10")))

plt <- ggplot(data=df_box) +
  geom_boxplot(mapping=aes(x=topic_order, y=weights, fill = red), color = red,  width=0.5, alpha = 0.6, outlier.shape = NA) +
  # geom_jitter(mapping=aes(x=topic_order, y=weights),width=0.2, size = 0.1, alpha=0.4) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Topic value rank", y = "Patient topic weight") +
  theme_bw(base_size = 20) +
  theme(legend.position = "None",panel.background=element_blank()) +
  geom_hline(aes(yintercept = 1), color = "red", linetype = "dashed")
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/fig3_patient_sparsity.png", plt, width = 5, height = 5)
# plt <- ggplot(filter(df_hist, type == "top_2_value"), aes(x=weights, fill=type)) +
#   geom_histogram( aes(y = (..count..)/sum(..count..)), color=red, alpha=0.6, center = 0.025, binwidth = 0.05) +
#   scale_y_continuous(labels = scales::percent) +
#   scale_x_continuous(limits = c(0,1),breaks = c(0, 0.2,0.4,0.6,0.8, 1.0)) +
#   theme_bw(base_size = 20) +
#   labs(x = "Largest topic value", y = "Density") +
#   theme(legend.position = "none")

plt <- ggplot(data=df_box) +
  geom_boxplot(mapping=aes(x=topic_order, y=weights, fill = record_number), width=0.8, alpha = 0.5, outlier.shape = NA) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=c("less_5" = blue, "5_to_10" = red, "more_10" = green)) +
  labs(x = "Topic value rank", y = "Patient topic weight") +
  theme_bw(base_size = 20) +
  theme(legend.position = "None",panel.background=element_blank()) +
  geom_hline(aes(yintercept = 1), color = "red", linetype = "dashed")
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/Sfig_patient_sparsity_record_number.png", plt, width = 5, height = 5)
plt <- ggplot(data=df_box) +
  geom_boxplot(mapping=aes(x=topic_order, y=weights, fill = record_number), width=0.5, alpha = 0.6, outlier.shape = NA) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=c("less_5" = blue, "5_to_10" = red, "more_10" = green)) +
  labs(x = "Topic value rank", y = "Patient topic weight") +
  theme_bw(base_size = 20) +
  theme(legend.position = "right",panel.background=element_blank()) +
  geom_hline(aes(yintercept = 1), color = "red", linetype = "dashed")
legend_plt <- cowplot::get_legend(plt)
grid.newpage()
legend_plt <- grid.draw(legend_plt)
print(paste0("proportion of max topic value > 0.95: ", mean(top_1_value[,1] > 0.95) ) )
##################################
# topic coverage of diseases correlation matrix
##################################
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
# we only need para$ds_list for these cooccurrence analysis
ds_eid_list <- para$unlist_Ds_id %>%
  group_by(Ds_id) %>%
  group_split()

cmb <- expand.grid(i=1:para$D, j=1:para$D)
cos_overlap <- function(idx, X){
  A = X[[idx[1]]]$eid
  B = X[[idx[2]]]$eid
  return( length(intersect(A,B))/(sqrt(length(A))*sqrt(length(B))) )
}
matrix_coocurre <- matrix(apply(cmb,1,function(x) cos_overlap(x, ds_eid_list)),para$D,para$D)

longData<-melt(matrix_coocurre) %>%
  filter(Var1 != Var2)
longData<-longData[longData$value >= 0.05 ,]
plt <- ggplot(longData, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low="grey90", high="red") +
  scale_x_continuous(breaks =cumsum(c(10,33, 22, 17, 23,43, 24, 55, 29, 28, 13, 50, 2)), labels = NULL) +
  scale_y_continuous(breaks =cumsum(c(10,33, 22, 17, 23,43, 24, 55, 29, 28, 13, 50, 2)), labels = NULL) +
  labs(x=NULL, y=NULL, title="Disease correlation matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
ggsave("../figures/disease_correlation.png", plt, width = 11, height = 10)

ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
phe_phecode <- read.csv("info_phenotype_phecode.csv")
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") )


# divide phecodes into subgroups

# types: A&B (1-10) infections; C&D (11-43) Neoplasm; E (44-65) metabolism; F&G (66 - 82) neurological;
# H (83-105) sensory; I (106-148) cardiovascular; J (149 - 172) respiratory; K (173 - 227) digestive;
# N (228 - 256) urinary&malegenetal; N (257 - 284) femalegenetal; L (285-297) skin; M (298-347) bone; 348-349 death&sepsis
# labels = c("infections", "neoplasm", "metabolism","neurological",
#            "sensory", "cardiovascular", "respiratory", "digestive",
#            "urinary&malegenetal", "femalegenetal", "skin", "bone", "death&sepsis")
# diseases_type <- c(rep("infections", 10),rep("neoplasm",33), rep("metabolism", 22),  rep("neurological", 17),
#                    rep("sensory", 23), rep("cardiovascular" ,43), rep("respiratory",24 ), rep("digestive", 55),
#                    rep("urinary&malegenetal", 29),rep("femalegenetal", 28),rep("skin", 13), rep("bone", 50), rep("death&sepsis", 2))

# plotting topic coverage Matrix

assignment_per_ds <- as.factor(sapply(1:para$D, function(j) which.max(colMeans(para$unlist_zn[para$ds_list[[j]]$id,])) ))

topic_data <- melt(matrix_coocurre) %>%
  filter(assignment_per_ds[Var2] == assignment_per_ds[Var1]) %>%
  mutate(topic = assignment_per_ds[Var2])
library(colorBlindness)
plt <- ggplot(topic_data, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=topic)) +
  scale_fill_manual(values= PairedColor12Steps[c(1:2,5:12)],na.value = "white") +
  scale_x_continuous(breaks =cumsum(c(10,33, 22, 17, 23,43, 24, 55, 29, 28, 13, 50, 2)), labels = NULL) +
  scale_y_continuous(breaks =cumsum(c(10,33, 22, 17, 23,43, 24, 55, 29, 28, 13, 50, 2)), labels = NULL) +
  labs(x=NULL, y=NULL, title="Topic coverage matrix") +
  theme_bw() + theme(axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
ggsave("../figures/topic_disease_distribution.png", plt, width = 11, height = 10)


###############################################
# save normalised topics by total incidence rate
# using heatmap
################################################
source("plotting_functions.R")
# first find the best rep
DIR <- "~/Desktop/comorbidity/Results/"
pt <- paste0("^rec2CVB0_model_output_PheCode_age_dependent_K", K_chosen,"_P",P_chosen, "_rep*")
temp <- list.files(paste(DIR, sep=""), pattern=pt)
lb_rep <- data_frame(df_P = as.integer(), lower_bound = as.numeric())
for(rep_id in 1:length(temp)){
  load(paste0(DIR,temp[rep_id]))
  cvrg_lb <-  model_output[[2]] %>%
    filter(!is.na(Lower_bound)) %>%
    slice_tail %>%
    pull(2)
  lb_rep <- lb_rep %>%
    add_row(df_P = P_chosen, lower_bound = cvrg_lb)
}
rep_id <- order(lb_rep$lower_bound,decreasing = T)[1]
################## note: the order of temp is not rep_id from 1-10 (1, 10, 2 ..)
load(paste0(DIR,temp[rep_id]))

para$D <- model_output[[3]]
para$list_above500occu <- read.csv("listAbove1000include_deaths_PheCode.csv")
phe_phecode <- read.csv("info_phenotype_phecode.csv")
phecodeList <- para$list_above500occu %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode"))

# normalise model output
# compute the baseline risk across whole population; will be used for normalization purpose
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")

pop_rate <- rec_data %>%
  group_by(diag_icd10) %>%
  summarise(counts = n()) %>%
  ungroup() %>%
  mutate(incidence_rate = counts/sum(counts))
pop_rate <- para$list_above500occu %>%
  left_join(pop_rate, by = "diag_icd10")

pal_age <- colorRampPalette(c(blue, red, orange, purple, green))
pal_age_vector <- pal_age(para$D)
thre_pick <- 5
phecodeList$phenotype <- as.character(phecodeList$phenotype)
phecodeList$phenotype[349]<- "death"
for(topic_id in 1:K){
  trajs <- sapply(1:para$D, function(j) model_output[[1]][30:80,j,topic_id]/pop_rate$incidence_rate[j]) # trajectories
  select_ds <- which(colMeans(trajs) > thre_pick)
  trajs <- trajs[,select_ds]
  phenotypes_selected <- phecodeList$phenotype[select_ds]
  longData<-melt(trajs)
  longData<-longData[longData$value!=0,, drop = F] %>%
    mutate(age = Var1 + 30, disease = phenotypes_selected[Var2])

  plt <- ggplot(longData, aes(x = age, y = disease)) +
    geom_raster(aes(fill=value)) +
    scale_fill_gradient(low="grey90", high="red") +
    labs(x="Age (years)", y="Disease", title= paste("topic: ",topic_id)) +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))
  ggsave(paste0("~/Desktop/comorbidity/figures/topics/rep", rep_id, "K",K,"P",df_P,"normalised_age_topics",topic_id,".png"),
         plt, width = 10, height = 10)
}

##################################
# Step 2+: visualise all individual topic distribution using UMAP!
##################################

load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
phe_phecode <- read.csv("info_phenotype_phecode.csv")
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") )

###########################################################
# need to add small amount of noise to avoid numeric error
# umap_patients.R
###########################################################
library(umap)
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
set.seed(19940110)
patient_loadings <- sweep((para$alpha_z - 1), 1, rowSums(para$alpha_z -1), FUN="/") # using mode as the individual loading prediction
onehot_anchors <- matrix(0, nrow = 10, ncol = 10)
onehot_anchors[cbind(1:10, 1:10)] <- 1
pertubation_rate <- 0.005
perturbed_data <- rbind(patient_loadings + matrix( pertubation_rate * rnorm(para$K * para$M), nrow = para$M) , onehot_anchors)
umap.fit <- umap(perturbed_data)
save(umap.fit, file = paste0("umap.fit", pertubation_rate, ".RData") )

load(paste0("umap.fit", pertubation_rate, ".RData"))
# plot some distributions
all_data <- read.csv(file="../UKBB_interim_files/ukb4775.csv", header=TRUE, sep=",")
ethnic_subgroup <- all_data %>%
  select(eid, X21000.0.0) %>%
  filter(!is.na(X21000.0.0)) %>%
  mutate(fid = eid) %>%
  mutate(ethnicity = if_else(X21000.0.0== 1001 | X21000.0.0== 1002 | X21000.0.0== 1003, "white", NULL)) %>%
  mutate(ethnicity = if_else(X21000.0.0== 3001 | X21000.0.0== 3002 | X21000.0.0== 3003, "SouthAsian", ethnicity)) %>%
  mutate(ethnicity = if_else(X21000.0.0== 4001 | X21000.0.0== 4002 | X21000.0.0== 4003, "black", ethnicity)) %>%
  mutate(ethnicity = if_else(X21000.0.0== 5, "Chinese", ethnicity)) %>%
  select(eid, ethnicity)

birthyear <- read.csv("~/Desktop/genetics_longitudinal_data/longitudinal_data/Year_of_birth.csv") %>%
  rename(sex = X31.0.0, BMI = X23104.0.0) %>%
  select(eid, sex, BMI)

df_umap_patient <- data.frame(x = umap.fit$layout[1:para$M,1], y = umap.fit$layout[1:para$M,2],
                      eid = para$eid) %>%
  left_join(ethnic_subgroup, by = "eid") %>%
  left_join(birthyear,  by = "eid")
df_umap_anchor <- data.frame(x = umap.fit$layout[(para$M+1):(para$M+10),1], y = umap.fit$layout[(para$M+1):(para$M+10),2],
                              label = sapply(1:10, function(x) paste0("topic: ",x)))

# plot ethnicity distribution
df_minority <- df_umap_patient %>%
  filter(ethnicity != "white")
plt <- ggplot(df_umap_patient) +
  geom_point(aes(x = x, y =y), color = grey, alpha = 0.1, size = 0.1) +
  geom_point(data = df_minority, aes(x = x, y =y, color = ethnicity), alpha = 0.5, size = 1) +
  scale_colour_manual(name="Self-reported ethnicity",values=c("SouthAsian" = red, "black" = grey, "Chinese" = grey)) +
  geom_label(data= df_umap_anchor, aes(x = x, y =y, label = label), size = 5)
ggsave("~/Desktop/comorbidity/figures/ethnicity_distribution.png", width = 10, height = 10,plt)

# plot sex distribution
df_umap_patient <- df_umap_patient %>%
  mutate(sex = if_else(sex == 1, blue, red))
plt <- ggplot(df_umap_patient) +
  geom_point(aes(x = x, y =y), color = grey, alpha = 0.1, size = 0.1) +
  geom_point(aes(x = x, y =y, color = sex), alpha = 0.1, size = 0.1) +
  geom_label(data= df_umap_anchor, aes(x = x, y =y, label = label), size = 5)
ggsave("~/Desktop/comorbidity/figures/sex_distribution.png", width = 10, height = 10,plt)

# plot BMI distribution
df_BMI <- df_umap_patient %>%
  filter(BMI > 45)
plt <- ggplot(df_umap_patient) +
  geom_point(aes(x = x, y =y), color = grey, alpha = 0.5, size = 0.1) +
  geom_point(data = df_BMI, aes(x = x, y =y, color = red), alpha = 1, size = 0.3) +
  geom_label(data= df_umap_anchor, aes(x = x, y =y, label = label), size = 5)

ggsave("~/Desktop/comorbidity/figures/BMI_distribution.png", width = 10, height = 10,plt)

################################
# all disease distribution
################################
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
label_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
df_umap_anchor <- data.frame(x = umap.fit$layout[(para$M+1):(para$M+10),1], y = umap.fit$layout[(para$M+1):(para$M+10),2],
                             label = sapply(1:10, function(x) paste0("topic: ", label_name[x])))

phe_phecode <- read.csv("info_phenotype_phecode.csv")
all.ds.system <- rec_data %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") )
all.ds.system <- all.ds.system %>%
  group_by(eid, exclude_name) %>%
  summarise(fill_value = n())
all.ds.system <- all.ds.system %>%
  filter(exclude_name %in% c("digestive", "circulatory system", "genitourinary", "musculoskeletal", "neoplasms",  "endocrine/metabolic", "respiratory",  "sense organs") ) %>%
  select(eid, exclude_name, fill_value) %>%
  pivot_wider(names_from = exclude_name, values_from = fill_value, values_fill = list(fill_value = 0))

df_umap_ds.system <- data.frame(x = umap.fit$layout[1:para$M,1], y = umap.fit$layout[1:para$M,2],
                              eid = para$eid) %>%
  left_join(all.ds.system, by = "eid" ) %>%
  replace_na(list(digestive=0, respiratory=0, genitourinary=0, neoplasms=0, `circulatory system`=0, `endocrine/metabolic`=0, `sense organs` =0,musculoskeletal = 0))

df_umap_ds.system <- df_umap_ds.system %>%
  mutate(digestive = pmin(digestive, 10))
plt <- ggplot(df_umap_ds.system, aes(x = x, y =y)) +
  geom_point(aes(color = digestive), alpha = 0.5, size = 0.1) +
  scale_color_gradient2(low= grey,  high=red) +
  geom_label(data= df_umap_anchor, aes(x = x, y =y, label = label), size = 5)
ggsave("~/Desktop/comorbidity/figures/digestive_distribution.png", width = 10, height = 10,plt)

df_umap_ds.system <- df_umap_ds.system %>%
  mutate(`circulatory system` = pmin(`circulatory system`, 10))
plt <- ggplot(df_umap_ds.system, aes(x = x, y =y)) +
  geom_point(aes(color = `circulatory system`), alpha = 0.5, size = 0.1) +
  scale_color_gradient2(low= grey,  high=red) +
  geom_label(data= df_umap_anchor, aes(x = x, y =y, label = label), size = 5)
ggsave("~/Desktop/comorbidity/figures/circulatory_system_distribution.png", width = 10, height = 10,plt)

df_umap_ds.system <- df_umap_ds.system %>%
  mutate(genitourinary = pmin(genitourinary, 10))
plt <- ggplot(df_umap_ds.system, aes(x = x, y =y)) +
  geom_point(aes(color = genitourinary), alpha = 0.5, size = 0.1) +
  scale_color_gradient2(low= grey,  high=red) +
  geom_label(data= df_umap_anchor, aes(x = x, y =y, label = label), size = 5)
ggsave("~/Desktop/comorbidity/figures/genitourinary_distribution.png", width = 10, height = 10,plt)

df_umap_ds.system <- df_umap_ds.system %>%
  mutate(respiratory = pmin(respiratory, 10))
plt <- ggplot(df_umap_ds.system, aes(x = x, y =y)) +
  geom_point(aes(color = respiratory), alpha = 0.5, size = 0.1) +
  scale_color_gradient2(low= grey,  high=red) +
  geom_label(data= df_umap_anchor, aes(x = x, y =y, label = label), size = 5)
ggsave("~/Desktop/comorbidity/figures/respiratory_distribution.png", width = 10, height = 10,plt)

df_umap_ds.system <- df_umap_ds.system %>%
  mutate(neoplasms = pmin(neoplasms, 10))
plt <- ggplot(df_umap_ds.system, aes(x = x, y =y)) +
  geom_point(aes(color = neoplasms), alpha = 0.5, size = 0.1) +
  scale_color_gradient2(low= grey,  high=red) +
  geom_label(data= df_umap_anchor, aes(x = x, y =y, label = label), size = 5)
ggsave("~/Desktop/comorbidity/figures/neoplasms_distribution.png", width = 10, height = 10,plt)

df_umap_ds.system <- df_umap_ds.system %>%
  mutate(musculoskeletal = pmin(musculoskeletal, 10))
plt <- ggplot(df_umap_ds.system, aes(x = x, y =y)) +
  geom_point(aes(color = musculoskeletal), alpha = 0.5, size = 0.1) +
  scale_color_gradient2(low= grey,  high=red) +
  geom_label(data= df_umap_anchor, aes(x = x, y =y, label = label), size = 5)
ggsave("~/Desktop/comorbidity/figures/musculoskeletal_distribution.png", width = 10, height = 10,plt)

interest_sysmtes <- c("digestive","respiratory","genitourinary","neoplasms",
  "circulatory system", "endocrine/metabolic", "musculoskeletal", "sense organs")
df.all.system <- df_umap_ds.system %>%
  rowwise() %>%
  mutate(system = which.max(c(digestive,respiratory,genitourinary,neoplasms,
                            `circulatory system`, `endocrine/metabolic`, musculoskeletal, `sense organs`))) %>%
  mutate(system = ifelse(sum(c(digestive,respiratory,genitourinary,neoplasms,
                                      `circulatory system`, `endocrine/metabolic`, musculoskeletal, `sense organs`)) == 0,
                          NA, interest_sysmtes[system]) ) %>%
  mutate(system = factor(system))

library(colorBlindness)
displayAvailablePalette(color="white")

plt <- ggplot(df.all.system) +
  geom_point(aes(x = x, y =y, color = system), alpha = 0.3, size = 0.01, shape = 20) +
  scale_color_manual(values= PairedColor12Steps[c(1:2, 6, 8:12)],na.value = "white") +
  geom_label(data= df_umap_anchor, aes(x = x, y =y, label = label), size = 5) +
  guides(colour = guide_legend(override.aes = list(size = 10, alpha = 1, shape = 20)))
ggsave("~/Desktop/comorbidity/figures/all_system_distribution.png", width = 10, height = 10,plt)
#######################################
# step 3: using topics to nominate disease pairs
#######################################

rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")

pca_rslt <- list()
lda_rslt <- list()
for(reps in 1:10){
  load(paste0("../Results/training_best_output_AgeLDA_RunNumber5K10_P5_rep", reps, ".RData"))
  testing_eid <- model_output[[7]]
  all_eid <- rec_data %>%
    group_by(eid) %>%
    summarise()
  training_eid <- all_eid %>%
    anti_join(testing_eid, by = "eid")
  training_data <- training_eid %>%
    left_join(rec_data, by = "eid")
  testing_data <- testing_eid %>%
    left_join(rec_data, by = "eid")

  # recover the disease assignments from the topics
  # re-estimate all the loadings from the topic
  training_para <- topics2weights(training_data, ds_list, 5, model_output[[1]])

  # perform PCA on the training set
  ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
  code2id <- function(x){
    return( match(x, ds_list$diag_icd10))
  }

  # here I am rounding the disease time to year for computation efficiency
  Ds_matrix <- training_data %>%
    arrange(eid)  %>%
    mutate(Ds_id = code2id(diag_icd10), ds_state = 1) %>%
    select(eid, Ds_id,ds_state) %>%
    pivot_wider(names_from = Ds_id, values_from = ds_state, values_fill = list(ds_state = 0))
  scaled_ds_matrix <- Ds_matrix %>%
    select(-eid) %>%
    scale()
  # reorder the columns so it matches other diseases
  idx_oder <- scaled_ds_matrix %>% colnames() %>% as.numeric %>% order()
  scaled_ds_matrix <- scaled_ds_matrix[,idx_oder]
  PCA_training <- scaled_ds_matrix %>%
    prcomp()
  data_pca <- Ds_matrix %>%
    select(eid) %>%
    cbind(PCA_training[["x"]])
  # using absolute value of PCs
  PCA_assigment <- sapply(1:para$D, function(j) which.max( abs(PCA_training$rotation[j,1:10]) ) )

  #######################################################################################
  # comorbidity proposal without subtypes -- one disease is only assigned to one subtype
  #######################################################################################
  # using the larger results file: ordering are changed in "model_output"
  betas <- training_para$pi_beta_basis
  topic_assign <- c()
  disease_idx <- c()
  mean_age <- c()
  num_cases <- c()
  for(j in 1:training_para$D){
    mean_age <- c(mean_age, mean(training_para$ds_list[[j]]$age_diag))
    topic_assign <- c(topic_assign, which.max( apply(training_para$unlist_zn[training_para$ds_list[[j]]$id,], 2, mean) ))
    disease_idx <- c(disease_idx, training_para$list_above500occu$diag_icd10[j])
    num_cases <- c(num_cases, dim(training_para$unlist_zn[training_para$ds_list[[j]]$id,])[1] )
  }
  common_disease_within_topics <- data.frame(disease = disease_idx, topic = topic_assign,
                                             mean_age = mean_age, case_number = num_cases, PCA_assigment = PCA_assigment)


  # propose top 100 diseases
  proposed_set <- common_disease_within_topics %>%
    arrange(desc(case_number)) %>%
    slice(1:100)

  # initialise the testing data parameter
  phe_phecode <- read.csv("info_phenotype_phecode.csv")
  ds_list <- training_para$list_above500occu %>%
    left_join(phe_phecode, by = c("diag_icd10" = "phecode") )
  testing_data <- testing_eid %>%
    left_join(rec_data, by = "eid")
  para_testing <- topic_init_age(testing_data, ds_list, 10, 5)

  comorbidity_testing_set <- list()
  pca_testing_set <- list()
  for(comorbidity_number in 2:5){
    print(comorbidity_number)
    for(pca_flag in 0:1){
      # for each comorbidity_number M , focusing on set of individuals who have at least M diseases
      high_record_indiv <- which(sapply(1:para_testing$M, function(s) dim(para_testing$w[[s]])[1] >= comorbidity_number))
      num_indiv <- length(high_record_indiv)
      num_each_disease <- bind_rows(para_testing$w[high_record_indiv])
      rate_each_disease <- sapply(1:para_testing$D, function(j) sum(num_each_disease$Ds_id == j)/num_indiv)
      OR_comorbidity <- list()
      for(topic_id in 1:para_testing$K){
        if(pca_flag){
          select_ds <- proposed_set %>%
            filter(PCA_assigment == topic_id)
        }else{
          select_ds <- proposed_set %>%
            filter(topic == topic_id)
        }
        # sometimes the topic does not have enough disease for proposal
        if(dim(select_ds)[1] >= comorbidity_number){
          # get all combinations
          disease_set <- combn(select_ds$disease, comorbidity_number)
          # using the index to facilitate computation
          disease_idx_set <- matrix(match(disease_set, para_testing$list_above500occu$diag_icd10), nrow = dim(disease_set)[1], ncol = dim(disease_set)[2])
          cases_comorbidity <- sapply(1:dim(disease_set)[2],
                                      function(idx) Reduce(intersect, lapply(disease_idx_set[,idx],
                                                                             function(x) para_testing$unlist_Ds_id[para_testing$ds_list[[x]]$id, ]$eid)) %>% length)
          odds_ratio <- (cases_comorbidity/num_indiv) / sapply(1:dim(disease_set)[2],
                                                               function(idx) (prod(rate_each_disease[disease_idx_set[,idx]])))
          comb_OR <- data.frame(cases_comorbidity = cases_comorbidity, odds_ratio = odds_ratio, diseases = t(apply(disease_set, 2, sort)))
          comb_OR$topic <-  topic_id
          comb_OR$comorbidity_number <- comorbidity_number
          OR_comorbidity[[topic_id]] <- comb_OR
        }

      }
      OR_comorbidity <- bind_rows(OR_comorbidity)
      if(pca_flag){
        pca_testing_set[[comorbidity_number]] <- OR_comorbidity
      }else{
        comorbidity_testing_set[[comorbidity_number]] <- OR_comorbidity
      }
    }
  }
  pca_rslt[[reps]] <- pca_testing_set
  lda_rslt[[reps]] <- comorbidity_testing_set
}
save(pca_rslt, file = paste0("../Results/","PCA_comorbidities.RData"))
save(lda_rslt, file = paste0("../Results/","ageLDA_comorbidities.RData"))

################################
# plot the odds ratio of nominated comorbidity
################################
load(paste0("../Results/","PCA_comorbidities.RData"))

pca_testing_set <- pca_testing_set %>%
  bind_rows() %>%
  mutate(method_topic = "PCA")
df_boxplot <- comorbidity_testing_set %>%
  bind_rows() %>%
  mutate(method_topic = "ATM") %>%
  bind_rows(pca_testing_set) %>%
  mutate(topic = factor(topic), comorbidity_number = factor(comorbidity_number), method_topic = factor(method_topic))

plt <- ggplot(data=df_boxplot,aes(x=comorbidity_number, y=log(odds_ratio), fill = method_topic)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width=.8,position=position_dodge(width=0.85)) +
  geom_point(position=position_jitterdodge(jitter.width = .02, dodge.width=0.85), size = 0.3, alpha=0.6) +
  # geom_jitter(width=0.1, size = 0.2, alpha=0.4) +
  labs(x = "Number of diseases in the comorbidity", y = "Prediction log Odds Ratio") +
  scale_fill_manual(name = "Models",values=cbPalette[2:5]) +
  theme(panel.background=element_blank())
ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/OR_ATM_pca.png"), plt, width = 10, height = 4)

################################
# save interesting comorbidities
################################
OR_comorbidity <- comorbidity_testing_set[[3]]
phe_phecode <- read.csv("info_phenotype_phecode.csv")
phecodeList <- OR_comorbidity %>%
  left_join(phe_phecode, by = c("diseases.1" = "phecode")) %>%
  left_join(phe_phecode, by = c("diseases.2" = "phecode")) %>%
  left_join(phe_phecode, by = c("diseases.3" = "phecode")) %>%
  filter(exclude_name.x != exclude_name.y, exclude_name.x != exclude_name, exclude_name != exclude_name.y ) %>%
  rename(phenotype.1 = phenotype.x, phenotype.2 = phenotype.y, phenotype.3 = phenotype) %>%
  select(cases_comorbidity, odds_ratio, diseases.1, diseases.2, diseases.3, topic, phenotype.1, phenotype.2, phenotype.3) %>%
  arrange(desc(odds_ratio))
phecodeList %>%
  write.csv("../paper_writing/Interesting_3_disease_comorbidity.csv", row.names = F)

#############################################
# predicted ordered comorbidity
#############################################

######################################
# compare disease topic assignments of LDA v.s. ageLDA
######################################
# basic LDA
load("../Results/Run_BaselinLDA_Phecode_K10_rep6.RData")
loadings_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[para$ds_list[[j]]$id,])) %>% t
hist(apply(loadings_per_ds, 1, max), breaks = 50)
longData<-melt(loadings_per_ds)
ggplot(longData, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Topic", y="Diseases", title="ageLDA") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
heterogeneity_nonage_test <- matrix(NA, nrow = para$D, ncol = 4)
loadings_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[para$ds_list[[j]]$id,])) %>% t

for(k_num in 2:5){
  Null_index <- which(sapply(1:para$D, function(j) sum(tail(sort(loadings_per_ds[j,]), (k_num-1) ))) > 0.9)
  for(j in 1:para$D){
    print(paste0("disease id: ", j))
    perm_sz <- 1000 # how many permutation samples to collect
    stats_target <- matrix(0, nrow = 2, ncol = 1)
    permutate_stats <- matrix(0, nrow = 2, ncol = perm_sz)

    # the ageLDA results
    loadings <- para$unlist_zn[para$ds_list[[j]]$id,]
    # for testing under null one topic:
    target_sz <- dim(loadings)[1]
    try({
      stats_target[1] <- rs_kmean(loadings, k_num) # get the target test statistic
      stats_target[2] <- rs_kmean(loadings, k_num - 1)})
    for(sp in 1:perm_sz){
      loadingsNull <- para$unlist_zn[para$ds_list[[sample(Null_index,1)]]$id,]
      tot_sz <- dim(loadingsNull)[1]
      try({ # if no possibility to cluster (all loadings equal to each other), just assume tot.betweenss to be 0
      permutate_stats[1, sp] <- rs_kmean(loadingsNull[sample(1:tot_sz, target_sz, replace = T), ], k_num)
      permutate_stats[2, sp] <- rs_kmean(loadingsNull[sample(1:tot_sz, target_sz, replace = T), ], k_num - 1)
      })
    }
    heterogeneity_nonage_test[j,k_num - 1] <- (1 + sum((stats_target[1] - stats_target[2]) <
                                                      (permutate_stats[1, ] - permutate_stats[2, ]) ) )/perm_sz
  }

}

basicLDA_ds_list <- ds_list
basicLDA_ds_list$p_2_subtype <- -log10(heterogeneity_nonage_test[,1])
basicLDA_ds_list$p_3_subtype <- -log10(heterogeneity_nonage_test[,2])
basicLDA_ds_list$p_4_subtype <- -log10(heterogeneity_nonage_test[,3])
basicLDA_ds_list$p_5_subtype <- -log10(heterogeneity_nonage_test[,4])
write.csv(basicLDA_ds_list, "test_basicLDA_heterogeneity_pvalues.csv")
basicLDA_ds_list <- basicLDA_ds_list %>%
  arrange(p_subtype)
basicLDA_ds_list$p_sim = sort(- log( runif(para$D)) )
ggplot(data = basicLDA_ds_list) +
  geom_point(aes(x = p_sim, y = p_subtype), color = grey, size = 2) +
  geom_abline(slope = 1) +
  geom_label_repel(aes(x = p_sim, y = p_subtype, label=ifelse(p_subtype > 4,as.character(phenotype),'')), max.overlaps = 50)
plt <- ggplot(data = basicLDA_ds_list) +
  geom_point(aes(x = p_sim, y = p_subtype), color = grey, size = 2) +
  geom_abline(slope = 1) +
  geom_label_repel(aes(x = p_sim, y = p_subtype, label=ifelse(p_subtype > 3,as.character(phenotype),'')),max.overlaps = 50)
ggsave("../figures/subtype_LDA.png", plt, width = 11, height = 10)

# compare power
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
phe_phecode <- read.csv("info_phenotype_phecode.csv")
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode"))
ds_list$p_subtype <- -log(heterogeneity_age_test[,1])
ds_list$p_basicLDA <- -log(heterogeneity_nonage_test[,1])
ggplot(data = ds_list) +
  geom_point(aes(x = p_basicLDA, y = p_subtype), color = grey, size = 2) +
  geom_abline(slope = 1)
# testing three cases: (1) the p-value calibration when applied to a single diseases
# (2) combine two random disease together and test the power
# (3) varying the number of each subgroup and test the power


# first results: Umap could plot all of the incidences
subtype_disease <- read.table("top50_hererogeneous_disease.txt")
for(id in 1:length(subtype_disease$V1)){
  ds_id <- as.character(subtype_disease$V1[id])
  j <- match(ds_id, para$list_above500occu$diag_icd10)
  loadings <- para$unlist_zn[para$ds_list[[j]]$id,]
  kmfit <- kmeans(loadings, 2)
  new_umap <- umap(rbind(assignment_per_ds, kmfit$centers))
  newdf_umap <- data.frame(x = new_umap$layout[,1], y = new_umap$layout[,2],
                        disease_type = as.factor(c(diseases_type, diseases_type[j],diseases_type[j])),
                        phecode = c(ds_list$diag_icd10, ds_id,  ds_id), phenotype = c(ds_list$phenotype, "subtp1","subtp2"),
                        subtypes = c(rep(grey, 349), blue, red))
  # plot histogram for each group
  plt1 <- ggplot(newdf_umap) +
    geom_point(aes(x = x, y =y, color = disease_type), size = 3) +
    labs(title=paste0("Disease: ", phecodeList$phenotype[j]))
  subtypes <-  c(rep(grey, 349), blue, red)
  plt2 <- ggplot(newdf_umap) +
    geom_point(aes(x = x, y =y ), color = subtypes, size = 3) +
    labs(title=paste0("Disease: ", phecodeList$phenotype[j]))

  ggsave(paste0("../figures/","umap_grouping",ds_id,".png"), plt1, width = 10, height = 10)
  ggsave(paste0("../figures/","umap_subtype",ds_id,".png"), plt2, width = 10, height = 10)

}


