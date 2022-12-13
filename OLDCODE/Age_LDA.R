library(topicmodels)
library(dplyr)
library(ggplot2)
require(survival)
library(stringr)
library(tidyverse)
library(gtools)
library(maxLik)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# set the commonly used color
blue <- cbPalette[6]
red <- cbPalette[7]
green <- cbPalette[4]
orange <- cbPalette[2]
grey <- cbPalette[1]
yellow <- cbPalette[5]
purple <- cbPalette[8]
skyblue <- cbPalette[3]

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
  # we will lose a few patient as some of them don't have birthyear info
  merge(birthyear, by="eid") %>% 
  mutate(birth_year=as.Date(paste(X34.0.0,X52.0.0,"01", sep="-"))) %>% 
  mutate(age_diag = difftime(epistart, birth_year, units = "days")/365) %>%
  filter(!is.na(age_diag))

# find all of the ICD-10 that is above 500 occurence (> 0.1%)
ds_occ_thre <- 500 # 500
list_above500occu <- new_data %>%
  select(eid, diag_icd10, age_diag) %>%
  filter(str_detect(diag_icd10, "^[A-N]")) %>%
  group_by(eid, diag_icd10) %>% 
  slice(1) %>%
  group_by(diag_icd10) %>%
  summarise(occ = n()) %>% 
  filter(occ >= ds_occ_thre)

# only keep the first incidence of disease 
first_incidence_age <- new_data %>%
  select(eid, diag_icd10, age_diag) %>%
  filter(diag_icd10 %in% list_above500occu$diag_icd10) %>%
  group_by(eid, diag_icd10) %>%
  slice(1) %>% # just pick the first record which since the record are in order; we use the the code below in case some records are not in order.
  # filter(n() == 1 | age_diag == min(age_diag) ) %>% %>%
  ungroup() # 1879,452 rows

ptm <- proc.time()
first_incidence_age <- new_data %>%
  select(eid, diag_icd10, age_diag) %>%
  filter(diag_icd10 %in% list_above500occu$diag_icd10) %>%
  group_by(eid, diag_icd10) %>% 
  # slice(1) %>% # just pick the first record and usually it is ranked with age
  filter(n() == 1 | age_diag == min(age_diag) ) %>% # this row is highly optimized, a lot faster the slice_min ### don't change
  slice(1) %>% # avoid the ties in min
  ungroup() 
print(proc.time() - ptm) # 166.086s to run


# create another disease record where each individual has to have at least two records
individual_two_records <- first_incidence_age %>%
  group_by(eid) %>% 
  filter(n() > 1)
##################################
# save the data for future use
##################################
# A2N is the range of codes considered here
write.csv(first_incidence_age, paste0("DiseaseAbove",ds_occ_thre,"occur_ICDA2N.csv"), row.names = F)

write.csv(list_above500occu, paste0("listAbove",ds_occ_thre,"_ICDA2N.csv"), row.names = F)

# save the data with death event
death_age <- read.csv("death_age.csv") %>%  # include death as one of the event
  mutate(eid = as.character(eid))

first_incidence_age %>%
  mutate(age_diag = as.double(age_diag)) %>%
  bind_rows(death_age) %>%
  arrange(eid) %>% # important step!
  write.csv(paste0("DiseaseAbove",ds_occ_thre,"occur_include_death_ICDA2N.csv"), row.names = F)

list_above500occu %>% 
  add_row(diag_icd10 = "Death", occ = dim(death_age)[1]) %>%
  write.csv(paste0("listAbove",ds_occ_thre,"include_deaths_ICDA2N.csv"), row.names = F)

# also save the data for individuals with at least one records 
write.csv(individual_two_records, paste0("rec2subjectAbove",ds_occ_thre,"occur_ICDA2N.csv"), row.names = F)

# only include death events to those who have two disease records
eid_list <- individual_two_records %>% 
  group_by(eid) %>%
  summarise() 

death_age <- read.csv("death_age.csv") %>%  # include death as one of the event
  mutate(eid = as.character(eid)) %>%
  semi_join(eid_list, by = "eid")

individual_two_records %>%
  mutate(age_diag = as.double(age_diag)) %>%
  bind_rows(death_age) %>%
  arrange(eid) %>% # important step!
  write.csv(paste0("rec2subjectAbove",ds_occ_thre,"occur_include_death_ICDA2N.csv"), row.names = F)

####################################
# data exploration 
####################################
first_incidence_age <- read.csv(paste0("DiseaseAbove500occur.csv"))
list_above500occu <- read.csv(paste0("listAbove500.csv"))

# plot the number distribution of indiviudal diseases
df_number_records <- first_incidence_age %>%
  group_by(eid) %>%
  summarise(n())

df_simu_pois <- data.frame(num_records = floor(1+rexp(370545, 1/(-0.5+5.87))))
ggplot(df_number_records) + 
  geom_histogram(aes(x = `n()`, fill = "true"), alpha = 1, binwidth = 1) + 
  geom_histogram(data = df_simu_pois, aes(x = num_records, fill = "simulated"), alpha = .5,, binwidth = 1) + 
  lims(x = c(0,40)) +
  scale_fill_manual(values = c("true" = grey, "simulated" = red)) + 
  theme(legend.position = c(.8,.8),panel.background=element_blank()) 

# print the disease distribution over age 
onset_by_year <- first_incidence_age %>% 
  mutate(age_diag = floor(age_diag)) %>%
  group_by(age_diag) %>%
  summarise(record_per_age = n())

ggplot(data = onset_by_year) + 
  geom_line(aes(x = age_diag, y = record_per_age))

# plot the record diagnosis age range (last diag_age - first diag_age) distribution within each indiviual
age_range <- first_incidence_age %>% 
  group_by(eid) %>%
  summarise(age_range = max(age_diag) - min(age_diag))

plt <- ggplot(age_range) + 
  geom_histogram(aes(x = age_range, fill = "record span"), alpha = 1, binwidth = 1) +
  scale_fill_manual(values = c("record span" = red)) + 
  theme(legend.position = c(.8,.8),panel.background=element_blank())+ 
  labs(x="Records span for each individual (years)", y="Number", title="Record span for each individual")
ggsave("~/Desktop/comorbidity/figures/rec_age_span.png",plt,width = 4, height =4)

##########################
# load the packages & data
##########################
source("topic_functions.R")
rec_data <- read.csv("DiseaseAbove500occur_include_death_ICDA2N.csv")
ds_list <- read.csv("listAbove500include_deaths_ICDA2N.csv")
topic_num <- 10
degree_free_num <- 4
para <- topic_init_age(rec_data, ds_list, topic_num, degree_free_num)

############################
# start optimization
############################
# set the number of update 
para$max_itr <- 500
para$alpha <- rep(1.1, para$K)
para$lb <- data.frame("Iteration" = 0,"Lower_bound" = comp_lda_lb(para))
para$itr_beta <- 1
para$itr_check_lb <- 1
para$itr_save <- 100
para$tol <- 10^(-7)
for(itr in 1:para$max_itr){
  print(paste0("Interation: ",itr))
  for(itr_inside in 1:para$itr_beta){ # in practice we perform quick steps a few times before move on.
    para <- comp_E_zn(para)
    para <- comp_E_lntheta(para)
  }
  para <- fast_update_age_depend_lda(para) 
  # para <- update_alpha(para) # we use non-informative alpha
  para$lb[nrow(para$lb) + 1,] <- c(itr, comp_lda_lb(para))
  if(itr %% para$itr_check_lb  ==0){
    curr_lb <- pull(filter(para$lb, Iteration == itr), Lower_bound) 
    prev_lb <- pull(filter(para$lb, Iteration == (itr -para$itr_check_lb )), Lower_bound) 
    print(paste0("Current Lower bound ", curr_lb, " at iteration: ",itr))
    try({
      if(is.finite((curr_lb - prev_lb)) & (curr_lb - prev_lb)/abs(prev_lb) < para$tol ){
        print(paste0("Optimization converged at step ", itr))
        break
      }
    })
  }
  if(itr %% para$itr_save ==0){
    save(para, file = paste0("Results/","Run_age_dependent_K",para$K,"_P",para$P,"_rep",para$rep_ID,"_intermediate",itr, ".RData"))
  }
}

# save the parameter
save(para, file = paste0("~/Desktop/comorbidity/Results/","Run_age_dependent_",Sys.Date(),".RData"))

load(file = paste0("~/Desktop/comorbidity/Results/Run_age_dependent_2021-02-03.RData"))
topics <- data.frame(ds_id = para$list_above500occu$diag_icd10, topics = para$beta)
df_topics <- list()
for(i in 1:para$K){
df_topics[[i]]  <- topics %>% filter(get(paste0("topics.",i)) > 10/para$D) %>%
                      pull(ds_id)
}

topics <- data.frame(ds_id = para$list_above500occu$diag_icd10, topics = para$beta)

longData<-melt(para$beta)
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

# lb plotting for different models
constant_effect_at_50 <- c(-10742047, -10734306, -10721206,-10719147, -10721083, -10772395)
age_effect_at_30 <- c(-10548999, -10497078, -10564816, -10579397, -10568686,-10659908)
df_log_lb <- data_frame(number_topics = 5:10, constant_effect_at_50, age_effect_at_30)
plt <- ggplot(df_log_lb) + 
  geom_line(aes(x = number_topics, y = constant_effect_at_50, color = "lowerbound without age effect"), linetype = "dashed") + 
  geom_line(aes(x = number_topics, y = age_effect_at_30, color = "lower bound fitted using age effect"), linetype = "dashed") + 
  geom_point(aes(x = number_topics, y = constant_effect_at_50, color = "lowerbound without age effect"), size = 1) + 
  geom_point(aes(x = number_topics, y = age_effect_at_30, color = "lower bound fitted using age effect"), size = 1) +
  scale_color_manual(name="Model Type", values=c("lowerbound without age effect" = red, "lower bound fitted using age effect"  = green))
ggsave("~/Desktop/comorbidity/figures/lower_bound_comparison.png",plt,width = 6, height =4)  

###################################################
# simulate a topic data set to test the inference
###################################################
source("topic_functions.R")
sample_sz <- 20000 # 20000
topic_number <- 5 # 5
disease_number <- 50
degree_freedom <- 3
para <- simulate_age_topic_data(sample_sz, topic_number, disease_number, degree_freedom)
############################
# run on the simulated data
############################
# set the number of update 
para$max_itr <- 500
para$alpha <- rep(1, para$K)
para$lb <- data.frame("Iteration" = 0,"Lower_bound" = comp_lda_lb(para))
para$itr_beta <- 1
para$itr_check_lb <- 1
para$tol <- 10^(-6)
cvb_tag <- T
cvb0_tag <- T
for(itr in 1:para$max_itr){
  print(paste0("Interation: ",itr))
  for(itr_inside in 1:para$itr_beta){ # in practice we perform quick steps a few times before move on.
    if(cvb_tag){
      if(cvb0_tag){
        para <- CVB0_E_zn(para)
      }else{
        para <- CVB_E_zn(para)
      }
    }else{
      para <- comp_E_zn(para)
      para <- comp_E_lntheta(para)
    }
  }
  para <- fast_update_age_depend_lda(para)
  # para <- update_alpha(para) # we use non-informative alpha
  if(cvb_tag){
    para$lb[nrow(para$lb) + 1,] <- c(itr, CVB_lb(para))
  }else{
    para$lb[nrow(para$lb) + 1,] <- c(itr, comp_lda_lb(para))
  }
  if(itr %% para$itr_check_lb ==0){
    curr_lb <- pull(filter(para$lb, Iteration == itr), Lower_bound) 
    prev_lb <- pull(filter(para$lb, Iteration == (itr -para$itr_check_lb )), Lower_bound) 
    print(paste0("Current Lower bound ", curr_lb, " at iteration: ",itr))
    try({
      if(is.finite((curr_lb - prev_lb)) & (curr_lb - prev_lb)/abs(prev_lb) < para$tol ){
        print(paste0("Optimization converged at step ", itr))
        break
      }
    })
  }
}
# compare the convergence pattern of CVB0 and CVB; studies has suggested CVB0 might converge faster
if(cvb0_tag){
  lb_cvb0 <- para$lb
}else{
  lb_cvb2 <- para$lb
}
iter_keep <- (min(dim(lb_cvb0)[1], dim(lb_cvb2)[1]))
df_cvb_lb_compare <- data.frame(lb_cvb0[1:iter_keep,], CVB2_lB = lb_cvb2$Lower_bound[1:iter_keep]) %>%
  filter(!is.infinite(Lower_bound))
ggplot(data = df_cvb_lb_compare, aes(x = Iteration)) +
  geom_line(aes(y = Lower_bound, color = "CVB0")) + 
  geom_line(aes(y = CVB2_lB, color = "CVB")) + 
  scale_color_manual(values=c("CVB" = red, "CVB0"  = blue))

if(cvb_tag){
  para_cvb <- para
}else{
  para_mean_field <- para
}

# get quantile of data
trim_90 <- para_cvb$beta[which(para_cvb$beta < quantile(para_cvb$beta, .95) & 
                      para_cvb$beta > quantile(para_cvb$beta, .05))]
var(trim_90)

# compare lower bound
iter_keep <- (min(dim(para_mean_field$lb)[1], dim(para_cvb$lb)[1]))
df_lb_compare <- data.frame(para_mean_field$lb[1:iter_keep,], CVB_lB = para_cvb$lb$Lower_bound[1:iter_keep]) %>%
  filter(!is.infinite(Lower_bound))
ggplot(data = df_lb_compare, aes(x = Iteration)) +
  geom_line(aes(y = Lower_bound), color = blue) + 
  geom_line(aes(y = CVB_lB), color = red)
  
# save the estimated profiles
library(reshape2)
library(ggplot2)
plt <- list()
for(i in 1:para$K){
  if(cvb_tag){
    longData<-melt(para_cvb$pi_beta_basis[30:81,,i])
  }else{
    longData<-melt(para_mean_field$pi_beta_basis[30:81,,i])
  }
  plt[[i]] <- ggplot(longData, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill=value)) + 
    scale_fill_gradient(low="grey90", high="red") +
    labs(x="Disease", y="Age", title=paste0("Topic ",i)) +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))
  # if(cvb_tag){
  #   ggsave(paste0("~/Desktop/comorbidity/figures/Topics/cvb_infer_simulate_age_K",para$K,"_topics",i,".png"),
  #          plt[[i]], width = 4, height = 6)
  # }else{
  #   ggsave(paste0("~/Desktop/comorbidity/figures/Topics/mean_field_infer_simulate_age_K",para$K,"_topics",i,".png"),
  #          plt[[i]], width = 4, height = 6)
  # }
  
}
# save the true profiles for comparison
for(i in 1:para$K){
  # plot the true topics for comparison
  longData<-melt(para_cvb$true_beta[30:81,,i])
  
  plt <- ggplot(longData, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill=value)) + 
    scale_fill_gradient(low="grey90", high="red") +
    labs(x="Disease", y="Age", title=paste0("Topic ",i)) +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))
  ggsave(paste0("~/Desktop/comorbidity/figures/Topics/true_simulate_age_K",para$K,"_topics",i,".png"), 
         plt, width = 4, height = 6)
}

# ordering
rank_beta <- sapply(1:para$K, function(i) 
  which.min(sapply(1:para$K, function(i_inf) 
    sum( (para$pi_beta_basis[30:81,,i_inf] - true_beta[30:81,,i] )^2 ) ) )) 
ordered_beta <- para$pi_beta_basis[30:81, ,rank_beta]

# compare individual topic composition
topic_id <- 2
df_topics <- data.frame(age=30:81, true_beta =  true_beta[30:81,,topic_id], inferred_topics = ordered_beta[,,topic_id])


plt <- ggplot(data = df_topics, aes(x = age)) + 
  geom_line(aes(y = inferred_topics.1), color = blue, linetype = "solid", alpha = .5, size = 1) + 
  geom_line(aes(y = true_beta.1), color = blue, linetype = "dashed",size = 1) + 
  geom_line(aes(y = inferred_topics.2), color = red, linetype = "solid", alpha = .5, size = 1) + 
  geom_line(aes(y = true_beta.2), color = red, linetype = "dashed",size = 1) + 
  geom_line(aes(y = inferred_topics.3), color = green, linetype = "solid", alpha = .5, size = 1) + 
  geom_line(aes(y = true_beta.3), color = green, linetype = "dashed",size = 1) + 
  geom_line(aes(y = inferred_topics.4), color = yellow, linetype = "solid", alpha = .5, size = 1) + 
  geom_line(aes(y = true_beta.4), color = yellow, linetype = "dashed",size = 1) + 
  geom_line(aes(y = inferred_topics.5), color = purple, linetype = "solid", alpha = .5, size = 1) + 
  geom_line(aes(y = true_beta.5), color = purple, linetype = "dashed",size = 1) + 
  geom_line(aes(y = inferred_topics.6), color = skyblue, linetype = "solid", alpha = .5, size = 1) + 
  geom_line(aes(y = true_beta.6), color = skyblue, linetype = "dashed",size = 1) 
  

# compare the topic composition 
para$fit_theta <- t(apply(para$alpha_z, 1, function(x) x/sum(x)))[,rank_beta]

df_theta_compare <- data.frame(true_theta =para$theta, fit_theta = para$fit_theta)

ggplot(sample_n(df_theta_compare,size = 1000)) + 
  geom_point(aes(x = true_theta.1, y =fit_theta.1), color = red, alpha = .5) +
  geom_point(aes(x = true_theta.3, y =fit_theta.3), color = blue, alpha = .5) + 
  geom_abline(intercept = 0, slope = 1) 


#####################################
# Analyse data that have converged
#####################################
load("~/Desktop/comorbidity/Results/Run_age_dependent_K5_P3rep1.RData")
topic_id <- 5
# focus on 30-81 years old
dominant_ds_id <- sapply(1:para$D, function(j) max(para$pi_beta_basis[30:81,j,topic_id]) > 20/para$D)
para$list_above500occu$diag_icd10[dominant_ds_id]
filtered_topics <- para$pi_beta_basis[30:81,dominant_ds_id,topic_id]

longData<-melt(filtered_topics)
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

for(topic_id in 1:para$K){
  dominant_ds_id <- sort(sapply(1:para$D, function(j) mean(para$pi_beta_basis[30:81,j,topic_id]) ), index.return = T, decreasing = T)$ix[1:30]
  selected_ds_id <- sort(sapply(dominant_ds_id, function(j) range(para$pi_beta_basis[40:70,j,topic_id])[2] - range(para$pi_beta_basis[30:81,j,topic_id])[1]), 
                         index.return = T, decreasing = T)$ix[1:6]
  
  df_topics <- data.frame(age=30:81, inferred_topics =  para$pi_beta_basis[30:81,selected_ds_id,topic_id])
  
  legend <- data.frame( ds_label = as.character(para$list_above500occu$diag_icd10[selected_ds_id]),
                        x_pos = rep(30, length(selected_ds_id)) ,
                        y_pos = unlist(df_topics[1,2:(1+length(selected_ds_id))]) )
  plt <- ggplot(data = df_topics, aes(x = age)) + 
    geom_line(aes(y = inferred_topics.1), color = blue, linetype = "solid", alpha = .5, size = 1) + 
    geom_line(aes(y = inferred_topics.2), color = red, linetype = "solid", alpha = .5, size = 1) + 
    geom_line(aes(y = inferred_topics.3), color = green, linetype = "solid", alpha = .5, size = 1) + 
    geom_line(aes(y = inferred_topics.4), color = yellow, linetype = "solid", alpha = .5, size = 1) + 
    geom_line(aes(y = inferred_topics.5), color = purple, linetype = "solid", alpha = .5, size = 1) + 
    geom_line(aes(y = inferred_topics.6), color = skyblue, linetype = "solid", alpha = .5, size = 1) + 
    geom_text(data = legend, aes(x = x_pos, y = y_pos, label = ds_label, hjust = -1, vjust = 1))
  ggsave(paste0("~/Desktop/comorbidity/figures/age_topics",topic_id,".png"), 
         plt, width = 6, height = 4)
}


##########################
# speed comparison dplyr
##########################
ptm <- proc.time()
first_incidence_age <- new_data %>%
  select(eid, diag_icd10, age_diag) %>%
  filter(diag_icd10 %in% list_above500occu$diag_icd10) %>%
  slice(1:100000) %>%
  group_by(eid, diag_icd10) %>% 
  # slice(1) %>% # just pick the first record and usually it is ranked with age
  slice_min(age_diag, n = 1) %>%
  ungroup() 
print(proc.time() - ptm) 

ptm <- proc.time()
first_incidence_age <- new_data %>%
  select(eid, diag_icd10, age_diag) %>%
  filter(diag_icd10 %in% list_above500occu$diag_icd10) %>%
  slice(1:1000000) %>%
  group_by(eid, diag_icd10) %>% 
  # slice(1) %>% # just pick the first record and usually it is ranked with age
  filter(n() == 1 | age_diag == min(age_diag) ) %>%
  ungroup() 
print(proc.time() - ptm) # this is perhaps doable 1million ~ 27 sec  

ptm <- proc.time()
first_incidence_age <- new_data %>%
  select(eid, diag_icd10, age_diag) %>%
  filter(diag_icd10 %in% list_above500occu$diag_icd10) %>%
  slice(1:100000) %>%
  group_by(eid, diag_icd10) %>% 
  # slice(1) %>% # just pick the first record and usually it is ranked with age
  filter(row_number() == 1) %>%
  ungroup() 
print(proc.time() - ptm) 

ptm <- proc.time()
first_incidence_age <- new_data %>%
  select(eid, diag_icd10, age_diag) %>%
  filter(diag_icd10 %in% list_above500occu$diag_icd10) %>%
  slice(1:100000) %>%
  group_by(eid, diag_icd10) %>% 
  # slice(1) %>% # just pick the first record and usually it is ranked with age
  slice(1) %>%
  ungroup() 
print(proc.time() - ptm) 




# testing 
para$exp_age_beta_basis <- array( sapply(1:para$K, 
                                         function(i) exp(para$age_basis %*% para$beta[,,i])),
                                  dim=c(dim(para$age_basis)[1], para$D, para$K))
identical(sapply(1:para$K, function(i) exp(para$age_basis %*% para$beta[,,i]), simplify = F)[[1]], 
          para$exp_age_beta_basis[,,1])

ptm <- proc.time()
for(i in 1:1){
  sapply(1:para$M, function(s) data.matrix(para$w[[s]]) ) 
}
print(proc.time() - ptm)

# comparison for lapply of para$w with first compute big matrix then split matrix
ptm <- proc.time()
for(i in 1:1){
  para$zeta <- lapply(para$w,
                      function(w) para$sum_exp_age_beta_basis[w$age_diag,,drop=F]) 
}
print(proc.time() - ptm)

ptm <- proc.time()
for(i in 1:1){
  para$zeta_full <- para$sum_exp_age_beta_basis[para$unlist_Ds_id$age_diag,,drop=F]
  para$zeta <- lapply(para$patient_lst, function(x) para$zeta_full[x,,drop=F]) 
}
print(proc.time() - ptm)
##########################################################################################
# optimizing para$beta is the most time consuming step. Each step will be evaluated here 
##########################################################################################
# test the time cost of compute two terms separately

ptm <- proc.time()
for(tt in 1:6800){
  term1 <- para$unlist_zn[para$ds_list[[j]]$id, i] %*%  (para$basis_phi[para$ds_list[[j]]$id, ] %*% beta_ij )
  term2 <- para$z_zeta[para$ds_list[[j]]$id, i] %*%  exp(para$basis_phi[para$ds_list[[j]]$id, ] %*% beta_ij )
}
print(proc.time() - ptm)

ptm <- proc.time()
for(tt in 1:6800){
  term1 <- para$unlist_zn[para$ds_list[[j]]$id, i] %*%  (para$basis_phi[para$ds_list[[j]]$id, ] %*% beta_ij ) -
    para$z_zeta[para$ds_list[[j]]$id, i] %*%  exp(para$basis_phi[para$ds_list[[j]]$id, ] %*% beta_ij )
}
print(proc.time() - ptm)

# testing to find the constraint of optim
para$fun_beta <- function(beta_ij,i,j){
  term <- para$unlist_zn[para$ds_list[[j]]$id, i] %*%  (para$basis_phi[para$ds_list[[j]]$id, ] %*% beta_ij ) - 
    para$z_zeta[para$ds_list[[j]]$id, i] %*%  exp(para$basis_phi[para$ds_list[[j]]$id, ] %*% beta_ij )
  return(term)
}
para$beta_grad <- function(beta_ij,i,j){
  term <-   t(para$basis_phi[para$ds_list[[j]]$id, ]) %*%
    ( para$unlist_zn[para$ds_list[[j]]$id, i] -
        (para$z_zeta[para$ds_list[[j]]$id, i] *exp(para$basis_phi[para$ds_list[[j]]$id, ] %*% beta_ij) ) ) 
  return(term)
}

i <- 3
j <- 3
itr_max <- 1000
sz_stp <- 10^(-5)
ptm <- proc.time()
for(tt in 1:1){
  a <- rnorm(para$P)
  lb_lst <- rep(NA, itr_max)
  gd_lst <- rep(NA, itr_max)
  for(itr in 1:itr_max){
    # print(paste0("iteration at: ", itr, "lower bound at: ", para$fun_beta(para$beta[[i]][,j],i,j)) )
    prev <- a
    a <- a + sz_stp * para$beta_grad(a,i,j)
    if(sum(abs(para$beta_grad(a,i,j))) < 10^(-3)){
      print("converged")
      est_gd <- a
      break
    }
    lb_lst[itr] <- para$fun_beta(a, i,j)
    gd_lst[itr] <- mean(para$beta_grad(a, i,j))
  }
  est_gd <- a
  # print(a)
}
print(proc.time() - ptm)
plot(lb_lst)

ptm <- proc.time()
for(tt in 1:100){
  a <- rnorm(para$P)
  opt <- optim(par = a,
             fn = function(x) fun_age_beta(x,i,j, para), 
             gr = function(x) grad_age_beta(x,i,j, para), method ="BFGS",
             control = list(fnscale = -1) )$par
  # print(opt$par)
  # print(para$fun_beta(opt$par,i,j))
  #print(para$fun_beta(opt$par,i,j))
}
print(proc.time() - ptm) # BFGS is much faster

plot(para$age_basis %*% opt$par)
plot(para$age_basis %*% est_gd)


para$fun_beta <- function(beta_ij,i_fn,j_fn){
  para$unlist_zn[para$ds_list[[j_fn]]$id, i_fn] %*%  (para$basis_phi[para$ds_list[[j_fn]]$id, ] %*% beta_ij ) - 
    para$z_zeta[para$ds_list[[j_fn]]$id, i_fn] %*%  exp(para$basis_phi[para$ds_list[[j_fn]]$id, ] %*% beta_ij )
}
para$beta_grad <- function(beta_ij,i_gd,j_gd){
   t(para$basis_phi[para$ds_list[[j_gd]]$id, ]) %*%
    ( para$unlist_zn[para$ds_list[[j_gd]]$id, i_gd] -
        (para$z_zeta[para$ds_list[[j_gd]]$id, i_gd] *exp(para$basis_phi[para$ds_list[[j_gd]]$id, ] %*% beta_ij) ) ) 
}
ptm <- proc.time()
for(tt in 1:1){
  para$beta <- array(rnorm(para$P * para$D,sd =0.1), dim = c(para$P, para$D, para$K))
  para$beta <- sapply(1:para$K, function(i) sapply(1:para$D, 
                                                   function(j)
                                                     optim(par = para$beta[,j,i,drop=F],
                                                           fn = function(x) fun_age_beta(x,i,j, para), 
                                                           gr = function(x) grad_age_beta(x,i,j, para), method ="BFGS",
                                                           control = list(fnscale = -1) )$par), simplify = F)

}
print(proc.time() - ptm) 
para$beta <- lapply(1:para$K, function(i) sapply(1:para$D, 
                                                 function(j)
                                                   optim(par = para$beta[,j,i],
                                                         fn = function(x) fun_age_beta(x,i,j, para), 
                                                         gr = function(x) grad_age_beta(x,i,j, para), method ="BFGS",
                                                         control = list(fnscale = -1) )$par ))


ptm <- proc.time()
for(tt in 1:1){
  para$beta <- array(rnorm(para$P * para$D,sd =0.1), dim = c(para$P, para$D, para$K))
  para$beta <- lapply(1:para$K, function(i) sapply(1:para$D, 
                                                   function(j)
                                                     optim(par = para$beta[[i]][,j],
                                                           fn = para$fun_beta, i_fn = i, j_fn = j,
                                                           gr = para$beta_grad, method ="BFGS",i_gd = i, j_gd = j,
                                                           control = list(fnscale = -1))$par
  ) )
}
print(proc.time() - ptm) 

# test if the lowerbound is actually increased by zeta
for(tt in 1:10){
  # first compute a full matrix then split: this is complicated but is much more efficient than directly lapply of para$w
  
  # compute M-step for beta: basic case, direct maximize the upper bound
  para$unlist_zn <- do.call(rbind, para$E_zn)
  
  para$beta <- array(rnorm(para$P * para$D,sd =0.1), dim = c(para$P, para$D, para$K))
  # I need to update zeta between each j; which will make this really slow....
  for(j in 1:para$D){
    # first compute the whole age_beta_basis, for each topic it is T-by-D
    para$exp_age_beta_basis <- array( sapply(1:para$K, 
                                             function(i) exp(para$age_basis %*% para$beta[,,i])),
                                      dim=c(dim(para$age_basis)[1], para$D, para$K))
    para$sum_exp_age_beta_basis <- apply(para$exp_age_beta_basis, c(1,3), sum)
    # this is the basis of softmax function
    para$pi_beta_basis <- apply(para$exp_age_beta_basis, c(1,3), function(x) x/sum(x)) %>% 
      aperm(perm = c(2,1,3))
    para$beta_w_full <- apply(para$pi_beta_basis, 3, function(x) 
      x[as.matrix(select(para$unlist_Ds_id, age_diag, Ds_id))]) 
    para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])
    
    # term3 <- sum(sapply(1:para$M, function(s) sum(para$E_zn[[s]] * log(para$beta_w[[s]]) ) ) )
    # 
    # term_mimic_3 <- sapply(1:para$K, function(i)
    #   sapply(1:para$D, function(j) t(para$basis_phi[para$ds_list[[j]]$id, ] %*% para$beta[, j, i] -
    #                                    (rowSums(exp(para$basis_phi[para$ds_list[[j]]$id, ]) %*% para$beta[, , i]))/
    #                                    para$zeta_full[para$ds_list[[j]]$id, i] -
    #                                    log(para$zeta_full[para$ds_list[[j]]$id, i]) + 1) %*% 
    #            para$unlist_zn[para$ds_list[[j]]$id, i] ) ) %>% sum
    
    # compute the variational parameter zeta nsi
    para$zeta_full <- para$sum_exp_age_beta_basis[para$unlist_Ds_id$age_diag,,drop=F]
    # para$zeta <- lapply(para$patient_lst, function(x) para$zeta_full[x,,drop=F]) 
    para$z_zeta <- para$unlist_zn/para$zeta_full  # this term helps speed up
    
    # debug
    ################################
    # need to implement an lower bound with zeta: maybe in the function file
    print(j)
    print(lb_beta_zeta(para))
    
    para$beta[,j,] <- sapply(1:para$K, function(i) optim(par = para$beta[,j,i],
                                                    fn = function(x) fun_age_beta(x,i,j, para), 
                                                    gr = function(x) grad_age_beta(x,i,j, para), method ="BFGS",
                                                    control = list(fnscale = -1) )$par )
    # debug
    print(lb_beta_zeta(para))
  }
  
  # compute the beta_w: exponential divided by sum
  # first compute the whole age_beta_basis, for each topic it is T-by-D
  para$exp_age_beta_basis <- array( sapply(1:para$K, 
                                           function(i) exp(para$age_basis %*% para$beta[,,i])),
                                    dim=c(dim(para$age_basis)[1], para$D, para$K))
  para$sum_exp_age_beta_basis <- apply(para$exp_age_beta_basis, c(1,3), sum)
  # this is the basis of softmax function
  para$pi_beta_basis <- apply(para$exp_age_beta_basis, c(1,3), function(x) x/sum(x)) %>% 
    aperm(perm = c(2,1,3))
  
  # update beta_w: list of Ns-by-K 
  para$beta_w_full <- apply(para$pi_beta_basis, 3, function(x) 
    x[as.matrix(select(para$unlist_Ds_id, age_diag, Ds_id))]) 
  para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])
  
}

# comparison of parallel over j and joint sapply 
ptm <- proc.time()
for(tt in 1:1){
  para$beta <- array(rnorm(para$P * para$D,sd =0.1), dim = c(para$P, para$D, para$K))
  para$beta <- array(sapply(1:para$K, function(i) sapply(1:para$D, function(j) 
    optim(par = para$beta[,j,i],
          fn = function(x) fun_age_beta(x,i,j, para), 
          gr = function(x) grad_age_beta(x,i,j, para), method ="BFGS",
          control = list(fnscale = -1) )$par )), 
    dim = c(para$P, para$D, para$K))
}
print(proc.time() - ptm) 

ptm <- proc.time()
for(tt in 1:1){
  para$beta <- array(rnorm(para$P * para$D,sd =0.1), dim = c(para$P, para$D, para$K))
  for(j in 1:para$D){
    para$beta[,j,] <- sapply(1:para$K, function(i)
      optim(par = para$beta[,j,i],
            fn = function(x) fun_age_beta(x,i,j, para),
            gr = function(x) grad_age_beta(x,i,j, para), method ="BFGS",
            control = list(fnscale = -1) )$par )
  }
  
}
print(proc.time() - ptm) # not much difference; might use for loop since it could save memory


source("topic_functions.R")
# check zeta lower bound with respect to each step
for(itr in 1:para$max_itr){
  print(paste0("Interation: ",itr))
  para <- comp_E_zn(para)
  para <- comp_E_lntheta(para)
  print(paste0(itr," before_beta ", lb_beta_zeta(para)))
  para <- update_age_depend_lda(para)
  print(paste0(itr," after_beta ", lb_beta_zeta(para)))
}

# maximize the expensive step after a few alternation of the cheap steps
source("topic_functions.R")
para$max_itr <- 500
para$lb <- data.frame("Iteration" = 0,"Lower_bound" = comp_lda_lb(para))
para$alpha <- rep(1, para$K)
para$itr_beta <- 5 
for(itr in 1:para$max_itr){
  print(paste0("Interation: ",itr))
  for(inner_itr in 1:3){
    para <- comp_E_zn(para)
    para <- comp_E_lntheta(para)
    print(paste0("inner iteration:", inner_itr," before_beta ", comp_lda_lb(para)))
  }
  para <- update_age_depend_lda(para)
  print(paste0(itr," after_beta ", comp_lda_lb(para)))
  # para <- update_alpha(para) # we use non-informative alpha
  if(itr %% 5 ==0){
    para$lb[nrow(para$lb) + 1,] <- c(itr, comp_lda_lb(para))
    print(paste0("Current Lower bound ", pull(filter(para$lb, Iteration == itr), Lower_bound), " at iteration: ",itr))
    if(abs( pull(filter(para$lb, Iteration == itr), Lower_bound) - 
            pull(filter(para$lb, Iteration == (itr -5)), Lower_bound) ) < .1 ){
      print(paste0("Optimization converged at step ", itr))
      break
    }
  }
}

#################################################
# optimize the gradient and likelihood function
# 2021-02-01
#################################################
j_gd <- j_fn <- which.max(sapply(1:para$D, function(j) length(para$ds_list[[j]]$id)))
i_gd <- i_fn <- 1
beta_ij <- para$beta[,j_gd,i_gd]
####################################################################################
# test the gradient starting point: original gradient code 
ptm <- proc.time()
for(tt in 1:100){
  t(para$basis_phi[para$ds_list[[j_gd]]$id, ]) %*% para$unlist_zn[para$ds_list[[j_gd]]$id, i_gd] -
    t(para$basis_phi) %*%   (para$z_zeta[, i_gd] *exp(para$basis_phi %*% beta_ij) )
}
print(proc.time() - ptm) # test function 

ptm <- proc.time()
for(tt in 1:100){
  crossprod(para$basis_phi[para$ds_list[[j_gd]]$id, ], para$unlist_zn[para$ds_list[[j_gd]]$id, i_gd]) -
    crossprod(para$basis_phi,  (para$z_zeta[, i_gd] *exp(para$basis_phi %*% beta_ij) ) ) 
}
print(proc.time() - ptm) # test function use crossprod # don't use t() for big matrix


ptm <- proc.time()
para$phi_z_zeta <- lapply(1:para$K, function(i) para$basis_phi*para$z_zeta[,i])
for(tt in 1:100){
  crossprod(para$basis_phi[para$ds_list[[j_gd]]$id, ], para$unlist_zn[para$ds_list[[j_gd]]$id, i_gd]) -
    crossprod(para$phi_z_zeta[[i_gd]], exp(para$basis_phi %*% beta_ij) ) 
}
print(proc.time() - ptm) # save para$phi_z_zeta provide another speedup
identical(crossprod(para$basis_phi,  (para$z_zeta[, i_gd] *exp(para$basis_phi %*% beta_ij) ) ),
          crossprod(para$phi_z_zeta[[i_gd]], exp(para$basis_phi %*% beta_ij) ) )


identical(as.vector(crossprod(para$basis_phi,  (para$z_zeta[, i_gd] *exp(para$basis_phi %*% beta_ij) ) )),
          phi_z_zeta_exp_phi_beta(x = para$basis_phi,beta = beta_ij, phi_z_zeta = para$phi_z_zeta[[i_gd]]) )

# test Rcpp on gradient function
para$phi_z_zeta <- lapply(1:para$K, function(i) para$basis_phi*para$z_zeta[,i])
ptm <- proc.time()
for(tt in 1:100){
  crossprod(para$basis_phi[para$ds_list[[j_gd]]$id, ], para$unlist_zn[para$ds_list[[j_gd]]$id, i_gd]) -
    phi_z_zeta_exp_phi_beta(x = para$basis_phi,beta = beta_ij, phi_z_zeta = para$phi_z_zeta[[i_gd]]) 
}
print(proc.time() - ptm) # this is about 1/2 of the original code


##########################################
# optimize likelihood function
# starting point: original gradient code 
j_gd <- j_fn <- which.max(sapply(1:para$D, function(j) length(para$ds_list[[j]]$id)))
i_gd <- i_fn <- 1
beta_ij <- para$beta[,j_gd,i_gd]

ptm <- proc.time()
for(tt in 1:100){
  para$unlist_zn[para$ds_list[[j_fn]]$id, i_fn] %*%  (para$basis_phi[para$ds_list[[j_fn]]$id, ] %*% beta_ij ) -
    para$z_zeta[, i_fn] %*%  exp(para$basis_phi %*% beta_ij)
}
print(proc.time() - ptm) 

identical(para$z_zeta[, i_fn] %*%  exp(para$basis_phi %*% beta_ij), 
          z_zeta_exp_phi_beta(para$basis_phi, beta_ij, para$z_zeta[, i_fn]))
# Rcpp of likelihood
ptm <- proc.time()
for(tt in 1:100){
  para$unlist_zn[para$ds_list[[j_fn]]$id, i_fn] %*%  (para$basis_phi[para$ds_list[[j_fn]]$id, ] %*% beta_ij ) - 
    z_zeta_exp_phi_beta(para$basis_phi, beta_ij, para$z_zeta[, i_fn])
}
print(proc.time() - ptm) # about 1/3 improvement
####################################################################
# save the terms phi z term for each disease: should be a column vector
####################################################################
para$phi_z <- lapply(1:para$D, function(j) 
  crossprod(para$basis_phi[para$ds_list[[j]]$id, ], para$unlist_zn[para$ds_list[[j]]$id, ]) ) 

# likelihood: starting point
ptm <- proc.time()
for(tt in 1:100){
  para$unlist_zn[para$ds_list[[j_fn]]$id, i_fn] %*%  (para$basis_phi[para$ds_list[[j_fn]]$id, ] %*% beta_ij ) - 
    z_zeta_exp_phi_beta(para$basis_phi, beta_ij, para$z_zeta[, i_gd])
}
print(proc.time() - ptm) # rcpp about 1/3 improvement

# Rcpp of likelihood
ptm <- proc.time()
for(tt in 1:100){
  para$phi_z[[j_fn]][,i_fn] %*% beta_ij - 
    z_zeta_exp_phi_beta(para$basis_phi, beta_ij, para$z_zeta[, i_fn])
}
print(proc.time() - ptm) # another 1/3 improvement

# gradient: starting point
ptm <- proc.time()
for(tt in 1:100){
  crossprod(para$basis_phi[para$ds_list[[j_gd]]$id, ], para$unlist_zn[para$ds_list[[j_gd]]$id, i_gd]) -
    phi_z_zeta_exp_phi_beta(x = para$basis_phi,beta = beta_ij, phi_z_zeta = para$phi_z_zeta[[i_gd]]) 
}
print(proc.time() - ptm) # this is about 1/2 of the original code

ptm <- proc.time()
for(tt in 1:100){
  para$phi_z[[j_fn]][,i_fn] -
    phi_z_zeta_exp_phi_beta(x = para$basis_phi,beta = beta_ij, phi_z_zeta = para$phi_z_zeta[[i_gd]]) 
}
print(proc.time() - ptm) # save another one seconds

##################################
# playing with RcppArmadillo
sourceCpp("armadillo_topic_functions.cpp")
ptm <- proc.time()
for(tt in 1:100){
  para$phi_z[[j_fn]][,i_fn] %*% beta_ij - 
    arma_z_zeta_exp_phi_beta(para$basis_phi, beta_ij, para$z_zeta[, i_gd])
}
print(proc.time() - ptm) # slower... 

##################################
# write the whole gradient and likelihood in Rcpp
##################################
j_gd <- j_fn <- j <- which.max(sapply(1:para$D, function(j) length(para$ds_list[[j]]$id)))
i_gd <- i_fn <- i <- 1
beta_ij <- para$beta[,j,i]

# likelihood
ptm <- proc.time()
for(tt in 1:1000){
  fun_age_beta(para$beta[,j,i],i,j, para)
}
print(proc.time() - ptm)

ptm <- proc.time()
for(tt in 1:1000){
  fun_age_betaC(para$basis_phi, para$beta[,j,i], para$z_zeta[, i], para$phi_z[[j]][,i])
}
print(proc.time() - ptm) # no real difference




################################
# test on the whole optim
################################
source("topic_functions.R")

j <- which.max(sapply(1:para$D, function(j) length(para$ds_list[[j]]$id)))
i <- 1

ptm <- proc.time()
para$beta <- array(rnorm(para$P * para$D,sd =0.1), dim = c(para$P, para$D, para$K))
for(tt in 1:1){
  for(j in 1:para$D){
    print(j)
    para$beta[,j,i] <- optim(par = para$beta[,j,i],
                             fn = function(x) fun_age_beta(x,i,j, para),
                             gr = function(x) grad_age_beta(x,i,j, para), method ="BFGS",
                             control = list(fnscale = -1) )$par
  }
}
print(proc.time() - ptm) 

# any parallel for using sapply?
ptm <- proc.time()
para$beta <- array(rnorm(para$P * para$D,sd =0.1), dim = c(para$P, para$D, para$K))
for(tt in 1:1){
  sapply(1:para$D, function(j) 
    optim(par = para$beta[,j,i],
          fn = function(x) fun_age_beta(x,i,j, para), 
          gr = function(x) grad_age_beta(x,i,j, para), method ="BFGS",
          control = list(fnscale = -1) )$par )
}
print(proc.time() - ptm) # about 10-20s faster: now it takes about 15min for each iteration of beta 

# any chance for earlier steps to converge quicker? reltol = 10^(-4) # doesn't seem to make d difference
ptm <- proc.time()
para$reltol <- 10^(-8)
para$beta <- array(rnorm(para$P * para$D * para$K,sd =0.1), dim = c(para$P, para$D, para$K))
for(tt in 1:1){
    for(j in 1:para$D){
      print(j)
      para$beta[,j,i] <- optim(par = para$beta[,j,i],
                               fn = function(x) fun_age_beta(x,i,j, para),
                               gr = function(x) grad_age_beta(x,i,j, para), method ="BFGS",
                               control = list(fnscale = -1, reltol = para$reltol) )$par
      print(para$beta[,j,i])
    }
}
print(proc.time() - ptm) 

# test the for loop on the simlated dataset 
ptm <- proc.time()
para$beta <- array(rnorm(para$P * para$D * para$K,sd =0.1), dim = c(para$P, para$D, para$K))
for(tt in 1:10){
  para$beta <- array(rnorm(para$P * para$D * para$K,sd =0.1), dim = c(para$P, para$D, para$K))
  for(j in 1:para$D){
    for(i in 1:para$K){
      para$beta[,j,i] <- optim(par = para$beta[,j,i],
                               fn = function(x) fun_age_beta(x,i,j, para),
                               gr = function(x) grad_age_beta(x,i,j, para), method ="BFGS",
                               control = list(fnscale = -1) )$par
    }
  }
}
print(proc.time() - ptm) 

ptm <- proc.time()
for(tt in 1:10){
  para$beta <- array(rnorm(para$P * para$D * para$K,sd =0.1), dim = c(para$P, para$D, para$K))
  para$beta <- array(sapply(1:para$K, function(i) sapply(1:para$D, function(j) 
      optim(par = para$beta[,j,i],
            fn = function(x) fun_age_beta(x,i,j, para),
            gr = function(x) grad_age_beta(x,i,j, para), method ="BFGS",
            control = list(fnscale = -1) )$par )),
      dim = c(para$P, para$D, para$K))
}
print(proc.time() - ptm) 

###########################################################
# one more way to speed up the algorithms significantly
###########################################################
# making use of the discrete age within the gradient exponential!
j_gd <- j_fn <- which.max(sapply(1:para$D, function(j) length(para$ds_list[[j]]$id)))
i_gd <- i_fn <- 1
beta_ij <- para$beta[,j_gd,i_gd]
ptm <- proc.time()
for(tt in 1:100){
  para$phi_z[[j_fn]][,i_fn] -
    phi_z_zeta_exp_phi_beta(x = para$basis_phi,beta = beta_ij, phi_z_zeta = para$phi_z_zeta[[i_gd]]) 
}
print(proc.time() - ptm) # save another one seconds

para$zeta_full <- para$sum_exp_age_beta_basis[para$unlist_Ds_id$age_diag,,drop=F]
# para$zeta <- lapply(para$patient_lst, function(x) para$zeta_full[x,,drop=F]) # we actually only need zeta as a full list
para$z_zeta <- para$unlist_zn/para$zeta_full # this term helps speed up
######## add following line to functions ##############
para$z_zeta_sum_by_age <- sapply(min(para$unlist_Ds_id$age_diag):max(para$unlist_Ds_id$age_diag), 
                                 function(y) colSums(para$z_zeta[which(para$unlist_Ds_id$age_diag == y),,drop=F]) ) %>% t

para$phi_z_zeta <- lapply(1:para$K, function(i) para$basis_phi*para$z_zeta[,i]) # help to speed up
para$phi_z <- lapply(1:para$D, function(j) 
  crossprod(para$basis_phi[para$ds_list[[j]]$id, ,drop = F], para$unlist_zn[para$ds_list[[j]]$id, ]) ) 

a1 <- z_zeta_exp_phi_beta(para$basis_phi, beta_ij, para$z_zeta[, i_fn])
b1 <- para$z_zeta_sum_by_age[,i_fn] %*% exp(para$age_basis[min(para$unlist_Ds_id$age_diag):max(para$unlist_Ds_id$age_diag), ] %*% beta_ij)
identical(a1,b1[1,1])

ptm <- proc.time()
para$basis_phi <- para$age_basis[para$unlist_Ds_id$age_diag,,drop = F]
for(tt in 1:1000){
    phi_z_zeta_exp_phi_beta(x = para$basis_phi,beta = beta_ij, phi_z_zeta = para$phi_z_zeta[[i_gd]]) 
}
print(proc.time() - ptm) # save another one seconds
# improve the speed by a lot!!!!
ptm <- proc.time()
para$age_basis_discrte <- para$age_basis[min(para$unlist_Ds_id$age_diag):max(para$unlist_Ds_id$age_diag), ]
para$basis_phi <- para$age_basis[para$unlist_Ds_id$age_diag,,drop = F]
for(tt in 1:1000){
    para$z_zeta_sum_by_age[,i_fn] %*% exp(para$age_basis_discrte %*% beta_ij)
}
print(proc.time() - ptm) # minimal time spent!


a <- phi_z_zeta_exp_phi_beta(x = para$basis_phi,beta = beta_ij, phi_z_zeta = para$phi_z_zeta[[i_gd]])
para$age_basis_discrte <- para$age_basis[min(para$unlist_Ds_id$age_diag):max(para$unlist_Ds_id$age_diag), ]
b <- crossprod(para$age_basis_discrte, para$z_zeta_sum_by_age[,i_fn] * exp(para$age_basis_discrte %*% beta_ij) )  
identical(a,b)

ptm <- proc.time()
para$age_basis_discrte <- para$age_basis[min(para$unlist_Ds_id$age_diag):max(para$unlist_Ds_id$age_diag), ]
para$basis_phi <- para$age_basis[para$unlist_Ds_id$age_diag,,drop = F]
for(tt in 1:1000){
  crossprod(para$age_basis_discrte, para$z_zeta_sum_by_age[,i_fn] * exp(para$age_basis_discrte %*% beta_ij) ) 
}
print(proc.time() - ptm) # minimal time spent!

# compare the full update beta
ptm <- proc.time()
for(tt in 1:10){
  fast_update_age_depend_lda(para)
}
print(proc.time() - ptm) # 67s spend! doable (before was 10h+)

# testing if accessing the parameters took a lot of time
fun_age_beta_fast <- function(beta_ij,i_fn,j_fn, para){
  para$phi_z[[j_fn]][,i_fn] %*% beta_ij - 
    para$z_zeta_sum_by_age[,i_fn] %*% exp(para$age_basis_discrte %*% beta_ij)
}

ptm <- proc.time()
for(tt in 1:10000){
  fun_age_beta_fast(beta_ij,i_fn,j_fn, para)
}
print(proc.time() - ptm) 
grad_age_beta_fast <- function(beta_ij,i_gd,j_gd, para){
  para$phi_z[[j_gd]][,i_gd] -
    crossprod(para$age_basis_discrte, para$z_zeta_sum_by_age[,i_fn] * exp(para$age_basis_discrte %*% beta_ij) )  
}
ptm <- proc.time()
for(tt in 1:10000){
  grad_age_beta_fast(beta_ij,i_gd,j_gd, para)
}
print(proc.time() - ptm) 

# by not access the para, we save about half of the time
phi_z <- para$phi_z[[j_fn]][,i_fn]
z_zeta <- para$z_zeta_sum_by_age[,i_fn]
age_bases <- para$age_basis_discrte
fun_age_beta_fast <- function(beta_ij,phi_z,z_zeta,age_bases){
  phi_z %*% beta_ij - 
    z_zeta %*% exp(age_bases %*% beta_ij)
}
ptm <- proc.time()
for(tt in 1:10000){
  fun_age_beta_fast(beta_ij,phi_z,z_zeta,age_bases)
}
print(proc.time() - ptm)

grad_age_beta_fast <- function(beta_ij,phi_z,z_zeta,age_bases){
  phi_z -
    crossprod(age_bases, z_zeta * exp(age_bases %*% beta_ij) )  
}
ptm <- proc.time()
for(tt in 1:10000){
  grad_age_beta_fast(beta_ij,phi_z,z_zeta,age_bases)
}
print(proc.time() - ptm)

# compare the full update beta
ptm <- proc.time()
for(tt in 1:10){
  fast_update_age_depend_lda(para)
}
print(proc.time() - ptm) # 67s! much faster than before (before was 10h+)


###########################################
# debugging the algorithms: why NAs occurs?
###########################################
example_file <- "../Results/Run_2rec_A2N_age_dependent_K11_P5_rep2.RData"
load(example_file) # it is caused by inf/inf in para$exp_age_beta_basis 

example_training_file <- "../Results/training_2rec_age_K5_P3_rep3.RData"
load(example_training_file)
