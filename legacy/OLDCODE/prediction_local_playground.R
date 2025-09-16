# tree LDA
library(topicmodels)
library(dplyr)
library(ggplot2)
require(survival)
library(stringr)
library(tidyverse)
library(gtools)
library(maxLik)
library(reshape2)

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

############################################
# plot accuracy over topics for half_half case
############################################
source("topic_functions.R")
rep_number <- 10
P <- 3
df_predict_lik_P_K <- data_frame(percentile = as.character(),df_K = as.integer(), specificity = as.numeric())
for(K in 5:20){
  lb_lst <- list()
    for(rep_id in 1:rep_number){
      try({load(paste0("../Results/prediction_age_K",K,"_P",P,"_rep",rep_id, ".RData"))
        df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
          add_row(percentile = "1_percentile", df_K = K, specificity = predictions[[1]]$mean_rank_top1)
        df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
          add_row(percentile = "2_percentile", df_K = K, specificity = predictions[[1]]$mean_rank_top2)
        df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
          add_row(percentile = "5_percentile", df_K = K, specificity = predictions[[1]]$mean_rank_top5)
        df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
          add_row(percentile = "10_percentile", df_K = K, specificity = predictions[[1]]$mean_rank_top10)
      })
    }
}
plt <- list()
# plot a box plot for selecting best topic number and degree of freedom
for(types in 1:4){
  percentile_range <- c("1_percentile", "2_percentile", "5_percentile", "10_percentile")[types]
  df_boxplot <- df_predict_lik_P_K %>%
    filter(df_K  < 21, df_K > 4) %>% # help visualise one or two Ks
    mutate(df_K = as.factor(df_K)) %>%
    filter(percentile == percentile_range)
  intcpt <- c(0.1395,0.2030,0.3558,0.4932)[types]
  plt[[types]] <- ggplot(data=df_boxplot,aes(x=df_K, y=specificity)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, width=.8,position=position_dodge(width=0.85)) +
    # geom_point(position=position_jitterdodge(jitter.width = .02, dodge.width=0.85), size = 0.3, alpha=0.6) +
    # geom_jitter(width=0.1, size = 0.2, alpha=0.4) + 
    labs(x = "Number of topics", y = "sensitivity") + 
    scale_fill_manual(values=cbPalette[2:7]) + 
    theme(panel.background=element_blank()) +
    geom_hline(yintercept = intcpt, linetype = "dashed", color = red)# add three dashed lines
  
}

# plot different number of records
source("topic_functions.R")
rep_number <- 10
P <- 3
df_predict_lik_P_K <- data_frame(percentile = as.character(),df_K = as.integer(), specificity = as.numeric(), record_number = as.integer())
for(K in c(5,10,15,20)){
  lb_lst <- list()
  for(rep_id in 1:rep_number){
    try({load(paste0("../Results/prediction_age_K",K,"_P",P,"_rep",rep_id, ".RData"))
      for(rec_num in 2:10){
        df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
          add_row(percentile = "1_percentile", df_K = K, specificity = predictions[[rec_num]]$mean_rank_top1,record_number = rec_num)
        df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
          add_row(percentile = "2_percentile", df_K = K, specificity = predictions[[rec_num]]$mean_rank_top2,record_number = rec_num)
        df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
          add_row(percentile = "5_percentile", df_K = K, specificity = predictions[[rec_num]]$mean_rank_top5,record_number = rec_num)
        df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
          add_row(percentile = "10_percentile", df_K = K, specificity = predictions[[rec_num]]$mean_rank_top10,record_number = rec_num)
      }
    })
  }
}
plt <- list()
# plot a box plot for selecting best topic number and degree of freedom
for(types in 1:4){
  percentile_range <- c("1_percentile", "2_percentile", "5_percentile", "10_percentile")[types]
  df_boxplot <- df_predict_lik_P_K %>%
    filter(df_K  < 21, df_K > 4) %>% # help visualise one or two Ks
    mutate(df_K = as.factor(df_K), record_number  = as.factor(record_number)) %>%
    filter(percentile == percentile_range)
  intcpt <- c(0.1395,0.2030,0.3558,0.4932)[types]
  plt[[types]] <- ggplot(data=df_boxplot,aes(x=record_number, y=specificity, fill=df_K)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, width=.8,position=position_dodge(width=0.85)) +
    geom_point(position=position_jitterdodge(jitter.width = .02, dodge.width=0.85), size = 0.3, alpha=0.6) +
    # geom_jitter(width=0.1, size = 0.2, alpha=0.4) + 
    labs(x = "Number of records", y = "sensitivity") + 
    scale_fill_manual(values=cbPalette[1:4]) + 
    theme(panel.background=element_blank()) +
    geom_hline(yintercept = intcpt, linetype = "dashed", color = red)# add three dashed lines
  
}


source("topic_functions.R")
rep_number <- 10
P <- 3
df_predict_lik_P_K <- data_frame(percentile = as.character(),df_K = as.integer(), specificity = as.numeric())
for(K in 5:20){
  lb_lst <- list()
  for(rep_id in 1:rep_number){
    try({load(paste0("../Results/prediction_age_K",K,"_P",P,"_rep",rep_id, ".RData"))
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "1_percentile", df_K = K, specificity = predictions[[1]]$mean_rank_top1)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "2_percentile", df_K = K, specificity = predictions[[1]]$mean_rank_top2)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "5_percentile", df_K = K, specificity = predictions[[1]]$mean_rank_top5)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "10_percentile", df_K = K, specificity = predictions[[1]]$mean_rank_top10)
    })
  }
}
plt <- list()
# plot a box plot for selecting best topic number and degree of freedom
for(types in 1:4){
  percentile_range <- c("1_percentile", "2_percentile", "5_percentile", "10_percentile")[types]
  df_boxplot <- df_predict_lik_P_K %>%
    filter(df_K  < 21, df_K > 4) %>% # help visualise one or two Ks
    mutate(df_K = as.factor(df_K)) %>%
    filter(percentile == percentile_range)
  intcpt <- c(0.1395,0.2030,0.3558,0.4932)[types]
  plt[[types]] <- ggplot(data=df_boxplot,aes(x=df_K, y=specificity)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, width=.8,position=position_dodge(width=0.85)) +
    # geom_point(position=position_jitterdodge(jitter.width = .02, dodge.width=0.85), size = 0.3, alpha=0.6) +
    # geom_jitter(width=0.1, size = 0.2, alpha=0.4) + 
    labs(x = "Number of topics", y = "sensitivity") + 
    scale_fill_manual(values=cbPalette[2:7]) + 
    theme(panel.background=element_blank()) +
    geom_hline(yintercept = intcpt, linetype = "dashed", color = red)# add three dashed lines
  
}
#####################################
# plot accuracy for each nth records! (one by one strategy)
#####################################
source("topic_functions.R")
rep_number <- 10
P <- 3
df_predict_lik_P_K <- data_frame(percentile = as.character(),df_K = as.integer(), specificity = as.numeric(), record_number = as.integer())
for(K in c(5,10,15,20)){
  lb_lst <- list()
  for(rep_id in 1:rep_number){
    try({
      load(paste0("../Results/prediction_onebyone_age_K",K,"_P",P,"_rep",rep_id, ".RData"))
      for(rec_num in 2:10){
        df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
          add_row(percentile = "1_percentile", df_K = K, specificity = prediction_onebyone[[rec_num]]$mean_rank_top1,record_number = rec_num)
        df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
          add_row(percentile = "2_percentile", df_K = K, specificity = prediction_onebyone[[rec_num]]$mean_rank_top2,record_number = rec_num)
        df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
          add_row(percentile = "5_percentile", df_K = K, specificity = prediction_onebyone[[rec_num]]$mean_rank_top5,record_number = rec_num)
        df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
          add_row(percentile = "10_percentile", df_K = K, specificity = prediction_onebyone[[rec_num]]$mean_rank_top10,record_number = rec_num)
      }
    })
  }
}
plt <- list()
# plot a box plot for selecting best topic number and degree of freedom
for(types in 1:4){
  percentile_range <- c("1_percentile", "2_percentile", "5_percentile", "10_percentile")[types]
  df_boxplot <- df_predict_lik_P_K %>%
    filter(df_K  < 21, df_K > 4) %>% # help visualise one or two Ks
    mutate(df_K = as.factor(df_K), record_number  = as.factor(record_number)) %>%
    filter(percentile == percentile_range)
  intcpt <- c(0.1395,0.2030,0.3558,0.4932)[types]
  plt[[types]] <- ggplot(data=df_boxplot,aes(x=record_number, y=specificity, fill=df_K)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, width=.8,position=position_dodge(width=0.85)) +
    geom_point(position=position_jitterdodge(jitter.width = .02, dodge.width=0.85), size = 0.3, alpha=0.6) +
    # geom_jitter(width=0.1, size = 0.2, alpha=0.4) + 
    labs(x = "nth records", y = "sensitivity") + 
    scale_fill_manual(values=cbPalette[1:4]) + 
    theme(panel.background=element_blank()) +
    geom_hline(yintercept = intcpt, linetype = "dashed", color = red)# add three dashed lines
  
}
########################################################
# plot the prediction accuracy over parameter setting
########################################################
source("topic_functions.R")
rep_number <- 10
maxP <- 7
df_predict_lik_P_K <- data_frame(df_P = as.integer(),df_K = as.integer(), likelihood = as.numeric())
for(K in 5:20){
  lb_lst <- list()
  for(df_P in 2:maxP){
    for(rep_id in 1:rep_number){
      try({load(paste0("../Results/prediction_onebyone_age_K",K,"_P",df_P,"_rep",rep_id, ".RData"))
        df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
          add_row(df_P = df_P, df_K = K, likelihood = prediction_onebyone_rslt[[1]][[1]])
      })
    }
  }
}

# plot a box plot for selecting best topic number and degree of freedom
df_boxplot <- df_predict_lik_P_K %>%
  filter(df_K  < 21, df_K > 4) %>% # help visualise one or two Ks
  mutate(df_P = as.character(df_P), df_K = as.factor(df_K)) 
plt <- ggplot(data=df_boxplot,aes(x=df_K, y=likelihood, fill = df_P)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width=.8,position=position_dodge(width=0.85)) +
  geom_point(position=position_jitterdodge(jitter.width = .02, dodge.width=0.85), size = 0.3, alpha=0.6) +
  # geom_jitter(width=0.1, size = 0.2, alpha=0.4) + 
  labs(x = "Number of topics", y = "sensitivity") + 
  scale_fill_manual(values=cbPalette[2:7]) + 
  theme(panel.background=element_blank()) 
ggsave(paste0("~/Desktop/comorbidity/figures/prediction_onebyone_topic_number_degree_freedom.png"), plt, width = 10, height = 4)


#######################################
# compare with and without age models
#######################################
source("topic_functions.R")
rep_number <- 10
P <- 5
df_predict_lik_P_K <- data_frame(percentile = as.character(),type = as.character(),df_K = as.integer(), specificity = as.numeric(), record_number = as.integer())
K <- 10
for(rep_id in 1:rep_number){
  try({load(paste0("../Results/prediction_onebyone_age_K",K,"_P",P,"_rep",rep_id, ".RData"))
    prediction_onebyone <- prediction_onebyone_rslt
    load(paste0("../Results/prediction_onebyone_baseline_K",K,"_rep",rep_id, ".RData"))
    for(rec_num in 2:10){
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "1_percentile", type = "age", df_K = K, specificity = prediction_onebyone[[rec_num]]$mean_rank_top1,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "2_percentile", type = "age",df_K = K, specificity = prediction_onebyone[[rec_num]]$mean_rank_top2,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "5_percentile", type = "age",df_K = K, specificity = prediction_onebyone[[rec_num]]$mean_rank_top5,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "10_percentile", type = "age",df_K = K, specificity = prediction_onebyone[[rec_num]]$mean_rank_top10,record_number = rec_num)
      # baseline LDA
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "1_percentile", type = "baseline", df_K = K, specificity = prediction_onebyone_rslt[[rec_num]]$mean_rank_top1,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "2_percentile", type = "baseline", df_K = K, specificity = prediction_onebyone_rslt[[rec_num]]$mean_rank_top2,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "5_percentile", type = "baseline", df_K = K, specificity = prediction_onebyone_rslt[[rec_num]]$mean_rank_top5,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "10_percentile", type = "baseline", df_K = K, specificity = prediction_onebyone_rslt[[rec_num]]$mean_rank_top10,record_number = rec_num)
    }
  })
}

# sumarise the mean and se
df_predict_se <- df_predict_lik_P_K %>%
  group_by(percentile, type, record_number) %>%
  summarise( se_specificity = sd(specificity)/sqrt(n())) %>%
  ungroup()

df_predict_mean <- df_predict_lik_P_K %>%
  group_by(percentile, type, record_number) %>%
  summarise( mean_specificity = mean(specificity)) %>%
  ungroup()

df_plot <- df_predict_se %>%
  left_join(df_predict_mean, by = c("percentile", "type", "record_number")) 

plt <- list()
# plot a box plot for selecting best topic number and degree of freedom
for(types in 1:4){
  percentile_range <- c("1_percentile", "2_percentile", "5_percentile", "10_percentile")[types]
  intcpt <- c(0.1395,0.2030,0.3558,0.4932)[types]
  df_predict_age <- df_plot %>% filter(percentile == percentile_range,type == "age")
  df_predict_baseline <- df_plot %>% filter(percentile == percentile_range,type == "baseline")
  plt[[types]] <- ggplot(df_predict_age, aes(x = record_number)) +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.position = "none",panel.background=element_blank(),plot.title = element_text(size = 10, face = "bold"))  + 
    labs(title = paste0("Percentile:", percentile_range), x = "nth records", y = "sensitivity") + 
    geom_line(aes(y = mean_specificity), linetype = "dashed", color = red) + 
    geom_pointrange(aes(y = mean_specificity, ymin = mean_specificity - 1.96*se_specificity, ymax = mean_specificity + 1.96*se_specificity), size = .5, color = red) +
    geom_line(data = df_predict_baseline ,aes( x = record_number, y = mean_specificity), linetype = "dashed", color = blue) + 
    geom_pointrange(data = df_predict_baseline, aes(x = record_number, y = mean_specificity, ymin = mean_specificity - 1.96*se_specificity, ymax = mean_specificity + 1.96*se_specificity), size = .5, color = blue) + 
    theme(panel.background=element_blank()) +
    geom_hline(yintercept = intcpt, linetype = "dashed", color = red)# add three dashed lines
  ggsave(paste0("~/Desktop/comorbidity/figures/agelda_lda_comparison_",percentile_range,"_K",K,"P",P,".png"), 
         plt[[types]], width = 4, height = 4)
}

################################################
# plotting the Risk ratio (the risk for those who did get the specified diseases v.s. those who get other disease) risk as manhatten plot
################################################
source("topic_functions.R")
source("plotting_functions.R")
ds_list <- read.csv("listAbove500include_deaths_ICDA2N.csv")
K <- 10
P <- 3
rep_number <- 10
DIR <- "~/Desktop/comorbidity/Results/"
# pt <- paste0("^rec2CVB0_model_output_A2N_age_dependent_K", K,"_P", P, "_rep*")
pt <- paste0("^OR_each_disease_K", K,"_P", P, "_rep*")
temp <- list.files(paste(DIR, sep=""), pattern=pt)
ORs <- list()
for(rep_id in 1:length(temp)){
  load(paste0(DIR,temp[[rep_id]]))
  ORs[[rep_id]]  <- prediction_OR_rslt[[1]][[1]]
}
ORs <- do.call(cbind,ORs)
logORs <- log(ORs)
mean_ORs <- rowMeans(logORs)
se_ORs <- apply(logORs, 1, function(x) sd(x)/sqrt(dim(logORs)[2]))
df_mht_plt <- data.frame(disease = ds_list$diag_icd10, occ = ds_list$occ, Mean = mean_ORs, SE = se_ORs)
ggplot(df_mht_plt, aes(x = 1:length(disease))) + 
  geom_pointrange(aes(y = Mean, ymin = Mean + 1.96 * SE, ymax = Mean - 1.96 * SE), size = .2) +
  labs(x = "disease", y= "log risk ratio")

df_filter_disease <- df_mht_plt %>%
  filter(Mean > .5, occ > 5000)

ggplot(df_filter_disease, aes(x = disease)) + 
  geom_pointrange(aes(y = Mean, ymin = Mean + 1.96 * SE, ymax = Mean - 1.96 * SE), size = .2) 

# take out the distributions
rep_id <- 1
load(paste0(DIR,temp[[rep_id]]))
loadings_each_disease <- prediction_OR_rslt[[2]][[1]]
loadings_each_disease[which(ds_list$diag_icd10 %in% df_filter_disease$disease),]
topics <- prediction_OR_rslt[[3]]

# plot the average risk profile for individuals with this code
ds_num <- dim(ds_list)[1]
pal_age <- colorRampPalette(c(blue, red, orange, purple, green))
pal_age_vector <- pal_age(ds_num)
thre_pick <- 10/ds_num
icd10 <- ds_list$diag_icd10
for(ds_code in df_filter_disease$disease){
  loadings_ds <- loadings_each_disease[which(ds_list$diag_icd10 == ds_code),]
  age_profile_ds <- sapply(1:dim(topics)[1], function(t) topics[t,,] %*% loadings_ds) %>% t
  
  trajs <- age_profile_ds[30:80,]  # trajectories
  plt <- plot_age_topics(icd10, trajs, thre_pick, pal_age_vector, start_age = 30)
  ggsave(paste0("~/Desktop/comorbidity/figures/topics/",ds_code, "risk_profile.png"), 
         plt, width = 4, height = 4)
}

###################################################################################
# compute OR per s.d. of log Risk ratio, using logistic regression for each disease
# also get a p-value for each disease
###################################################################################
# following previous discussion over odds ratio, we formulate a different strategy here
# load example data for coding 
rec_data <- read.csv("rec2subjectAbove500occur_include_death_ICDA2N.csv")
ds_list <- read.csv("listAbove500include_deaths_ICDA2N.csv")
load(file = paste0("../Results/training_2rec_age_K5_P3_rep3.RData"))
all_eid <- rec_data %>%
  group_by(eid) %>%
  summarise() 
training_eid <- data.frame(eid = para$eid)
testing_eid <- all_eid %>% 
  anti_join(training_eid, by = "eid")

testing_data <- testing_eid %>%  
  left_join(rec_data, by = "eid")

max_predict <- 11 # predict the records until which step?
para_training <- para

prediction_OR_rslt <- prediction_OR_onebyone(testing_data, ds_list, para_training, max_predict) 

# plot the per-sd logOR for diseases
source("topic_functions.R")
source("plotting_functions.R")
ds_list <- read.csv("listAbove500include_deaths_ICDA2N.csv")
K <- 10
P <- 3
rep_number <- 10
DIR <- "~/Desktop/comorbidity/Results/"
# pt <- paste0("^rec2CVB0_model_output_A2N_age_dependent_K", K,"_P", P, "_rep*")
pt <- paste0("^OR_each_disease_K", K,"_P", P, "_rep*")
temp <- list.files(paste(DIR, sep=""), pattern=pt)
logORs <- list()
p_values <- list()
for(rep_id in 1:length(temp)){
  load(paste0(DIR,temp[[rep_id]]))
  logORs[[rep_id]]  <- prediction_OR_rslt[[4]][,1]
  p_values[[rep_id]] <- prediction_OR_rslt[[4]][,4]
}
# per-sd logRR
logORs <- do.call(cbind,logORs)
mean_ORs <- rowMeans(logORs)
se_ORs <- apply(logORs, 1, function(x) sd(x)/sqrt(dim(logORs)[2]))

# p-values
p_values <- do.call(cbind,p_values)
mean_p <- rowMeans(p_values)
se_p <- apply(p_values, 1, function(x) sd(x)/sqrt(dim(p_values)[2]))
df_mht_plt <- data.frame(disease = ds_list$diag_icd10, occ = ds_list$occ, mean_ORs, se_ORs, mean_p, se_p)
ggplot(df_mht_plt, aes(x = 1:length(disease))) + 
  geom_pointrange(aes(y = mean_ORs, ymin = mean_ORs + 1.96 * se_ORs, ymax = mean_ORs - 1.96 * se_ORs), size = .2) +
  labs(x = "disease", y= "log OR per-sd log risk ratio")

ggplot(df_mht_plt, aes(x = 1:length(disease))) + 
  geom_point(aes(y = -log(mean_p)), size = .2) +
  labs(x = "disease", y= "p-value")

##########################################
# plot disease by disease cooccurrence rate
# here we are interested in the pattern where disease tend to cooccure
##########################################
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
topic_num <- 5
degree_free_num <- 2
para <- topic_init_age(rec_data, ds_list, topic_num, degree_free_num)
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

longData<-melt(matrix_coocurre)
longData<-longData[longData$value > 0.1 ,]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))


# this function compute the average difference for two diseases
age_diff <- function(idx, X){
  A <- X[[idx[1]]]$eid
  B <- X[[idx[2]]]$eid
  ids <- intersect(A,B)
  if(length(ids) > 50){ # only care about disease with more than 50 couccurence 
    A_location <- sapply(ids, function(x) which(A == x))
    B_location <- sapply(ids, function(x) which(B == x))
    return( mean(X[[idx[1]]]$age_diag[A_location] - X[[idx[2]]]$age_diag[B_location]))
  }else{
    return(0)
  }
}
matrix_year_diff <- matrix(apply(cmb,1,function(x) age_diff(x, ds_eid_list)),para$D,para$D)

# plot the age difference for those with high coocurrence rate
idx_long <- melt(matrix_coocurre)
longData<- melt(matrix_year_diff)
longData<-longData[idx_long$value > 0.1 ,]

# what this plot did was 
plt <- ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradient2(low="blue", mid = "grey", high="red") +
  labs(x="Diseases", y="Diseases", title="Average age gap (for disease pairs that have high comorbidities (cosine similarity > 0.1))") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
ggsave(paste0("~/Desktop/comorbidity/figures/comorbidity_age_ordering.png"), plt, width = 10, height = 10)

# save all these pair of diseases
disease_pair <- melt(matrix_year_diff) %>%
  rename(age_gap = value) %>%
  left_join(melt(matrix_coocurre), by = c("Var1", "Var2")) %>%
  rename(comorbidity_rate = value) %>%
  filter(age_gap > 2, comorbidity_rate > 0.1) %>%
  mutate(Disease1 = para$list_above500occu$diag_icd10[Var1], Disease2 = para$list_above500occu$diag_icd10[Var2]) %>%
  select(-Var1,-Var2)
write.csv(disease_pair, file = "~/Desktop/comorbidity/figures/comorbidity_pair.csv", row.names = F, col.names = T)

# create a figure for connection
source("plotting_functions.R")
age_gap_min <- 3
comorbidity_rate_thre <- .1
comorbidity_graph <- create_net(matrix_year_diff, matrix_coocurre, age_gap_min,comorbidity_rate_thre)

# create a location map by the age of the diseases:
# step 1: choose the disease with the most number of childrens
start <- edges %>% 
  group_by(source) %>%
  summarise(number_children = n()) %>%
  filter(number_children == max(number_children)) %>%
  pull(1)
  
# step 2: recursively find all childrens and chidren of children, position on x determined by age
l_x <- rep(NA, vcount(comorbidity_graph))
l_y <- rep(NA, vcount(comorbidity_graph))
l_x[which(filter_nodes$diag_icd10 == start)] <- 0
l_y[which(filter_nodes$diag_icd10 == start)] <- sum(!is.na(l_x))

# children will recursively search all tree, both parent and children
rslt <- recurs_children(node = start, l_x, l_y, edges,filter_nodes)
l_x <- rslt[[1]]
l_y <- rslt[[2]]
# step 3: recursively find all parents of a node
while(anyNA(l_x)){
  empty <- filter_nodes$diag_icd10[which(is.na(l_x))]
  nextnd <- edges %>% 
    group_by(source) %>%
    summarise(number_children = n()) %>%
    filter(source %in% empty) %>%
    filter(number_children == max(number_children)) %>%
    slice(1) %>%
    pull(1)
  l_x[which(filter_nodes$diag_icd10 == nextnd)] <- 0
  l_y[which(filter_nodes$diag_icd10 == nextnd)] <- sum(!is.na(l_x))
  
  rslt <- recurs_children(node = nextnd, l_x, l_y, edges,filter_nodes)
  l_x <- rslt[[1]]
  l_y <- rslt[[2]]
}

l <- cbind(l_x, l_y)
pdf(file = paste0("~/Desktop/comorbidity/figures/comorbidities_full_network.pdf"), width = 8, height = 8)
plot(comorbidity_graph, layout = l, edge.arrow.size = .2, vertex.label.color="black",vertex.label.cex=.5, vertex.label.dist=0.2, vertex.label.font = 2, edge.curved = 0)
dev.off()
#######################################################
# additional analysis 1: get the subgraph of one disease
# also perform analysis over only one node, start (ignore the others)
start <- "K227"

# step 2: recursively find all childrens and chidren of children, position on x determined by age
l_x <- rep(NA, vcount(comorbidity_graph))
l_y <- rep(NA, vcount(comorbidity_graph))
l_x[which(filter_nodes$diag_icd10 == start)] <- 0
l_y[which(filter_nodes$diag_icd10 == start)] <- sum(!is.na(l_x))

# children will recursively search all tree, both parent and children
rslt <- recurs_children(node = start, l_x, l_y, edges,filter_nodes)
l_x <- rslt[[1]]
l_y <- rslt[[2]]
# if started with a endnode, search backwards
rslt <- recurs_parent(node = start, l_x, l_y, edges,filter_nodes)
l_x <- rslt[[1]]
l_y <- rslt[[2]]

keep_node <- which(!is.na(l_y))
sub_g <- induced_subgraph(comorbidity_graph,keep_node)
l <- cbind(l_x[keep_node], l_y[keep_node])
dev.off()
pdf(file = paste0("~/Desktop/comorbidity/figures/comorbidities_network",start,".pdf"), width = 8, height = 8)
plot(sub_g, layout = l,  edge.arrow.size = .2, vertex.label.color="black",vertex.label.cex=.5, vertex.label.dist=.2, vertex.label.font = 2, edge.curved = 0)
dev.off()

#######################################################
# additional analysis 2: separate edges connected to one node
# only keep edges of one node 
disease_list <- melt(matrix_year_diff) %>%
  rename(age_gap = value) %>%
  left_join(melt(matrix_coocurre), by = c("Var1", "Var2")) %>%
  rename(comorbidity_rate = value) %>%
  filter(age_gap > age_gap_min | age_gap < -age_gap_min, comorbidity_rate > comorbidity_rate_thre) %>%  # to include both order
  mutate(Disease = para$list_above500occu$diag_icd10[Var1]) %>%
  group_by(Disease) %>%
  summarise(number_connection = n()) %>%
  arrange(-number_connection)

for(dsid in 1:10){
  disease <- disease_list$Disease[dsid]
  
  age_gap_min <- 2
  comorbidity_rate_thre <- .1
  disease_pair <- melt(matrix_year_diff) %>%
    rename(age_gap = value) %>%
    left_join(melt(matrix_coocurre), by = c("Var1", "Var2")) %>%
    rename(comorbidity_rate = value) %>%
    filter(age_gap > age_gap_min, comorbidity_rate > comorbidity_rate_thre) %>%
    mutate(Disease1 = para$list_above500occu$diag_icd10[Var1], Disease2 = para$list_above500occu$diag_icd10[Var2]) %>%
    select(-Var1,-Var2) %>%
    filter(Disease1 == disease | Disease2 == disease)
  
  pal_tree <- colorRampPalette(c(red, grey))
  pal_tree_vector <- rev(pal_tree(100))
  
  edges <- disease_pair  %>% # remove edges that have same source and target
    rename(source = Disease2, target = Disease1) %>%
    select(source, target, age_gap, comorbidity_rate)
  
  filter_nodes <- data.frame(diag_icd10 = unique(c(as.character(edges$source), as.character(edges$target)) ) )  %>%
    left_join(para$list_above500occu, by = "diag_icd10")
  
  comorbidity_graph <- graph_from_data_frame(d=edges, vertices=filter_nodes, directed=T) 
  E(comorbidity_graph)$arrow.mode <- 2
  E(comorbidity_graph)$width <- edges$comorbidity_rate*10
  E(comorbidity_graph)$color <- pal_tree_vector[floor(edges$comorbidity_rate/max(edges$comorbidity_rate) * 100)]
  V(comorbidity_graph)$size <- sqrt(V(comorbidity_graph)$occ)/20
  V(comorbidity_graph)$color <- pal_tree_vector[pmin(floor(V(comorbidity_graph)$occ/max(filter_nodes$occ) * 100) + 1, 100)]
  V(comorbidity_graph)$frame.color <- NA
  
  l_x <- rep(NA, vcount(comorbidity_graph))
  l_y <- rep(NA, vcount(comorbidity_graph))
  l_x[which(filter_nodes$diag_icd10 == disease)] <- 0
  l_y[which(filter_nodes$diag_icd10 == disease)] <- sum(!is.na(l_x))
  
  # children will recursively search all tree, both parent and children
  rslt <- recurs_children(node = disease, l_x, l_y, edges,filter_nodes)
  l_x <- rslt[[1]]
  l_y <- rslt[[2]]
  # if started with a endnode, search backwards
  rslt <- recurs_parent(node = disease, l_x, l_y, edges,filter_nodes)
  l_x <- rslt[[1]]
  l_y <- rslt[[2]]
  l <- cbind(l_x, l_y)
  pdf(file = paste0("~/Desktop/comorbidity/figures/Topics/pair_comorbidities",disease,".pdf"), width = 8, height = 8)
  plot(comorbidity_graph, layout = l,  edge.arrow.size = .2, vertex.label.color="black",vertex.label.cex=.5, vertex.label.dist=.2, vertex.label.font = 2, edge.curved = 0)
  dev.off()
}

################################################
# do the age sex baseline model 
################################################
source("topic_functions.R")
require("nnet")
require("mlogit")
require("mnlogit")
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")

ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")

all_eid <- rec_data %>%
  group_by(eid) %>%
  summarise() 

birthyear <- read.csv(file = "~/Desktop/genetics_longitudinal_data/longitudinal_data/Year_of_birth.csv") %>%
  select( eid, X31.0.0) %>%
  rename(sex = X31.0.0)

code2id <- function(x){
  return( match(x, ds_list$diag_icd10))
}
rec_data <- rec_data %>% 
  left_join(birthyear, by = "eid") %>%
  mutate(diag_icd10 = code2id(diag_icd10))
  
write.csv(rec_data,"~/Desktop/comorbidity/Multi-morbidity_biobank/age_sex_disease.csv", row.names = F)
# create 10 testing/training division
for(div in 1:10){
  training_eid <- all_eid %>% 
    sample_frac(size  = .8)
  
  testing_eid <- all_eid %>% 
    anti_join(training_eid, by = "eid")
  
  training_data <- training_eid %>%  
    left_join(rec_data, by = "eid") %>%
    mutate(age_diag = age_diag/max(age_diag))
  
  write.csv(training_data,paste0("~/Desktop/comorbidity/Multi-morbidity_biobank/age_sex_data/training",div,"_age_sex_Phecode.csv"), row.names = F)
  
  testing_data <- testing_eid %>%  
    left_join(rec_data, by = "eid") %>%
    mutate(age_diag = age_diag/max(age_diag))
  
  write.csv(testing_data,paste0("~/Desktop/comorbidity/Multi-morbidity_biobank/age_sex_data/testing",div,"_age_sex_Phecode.csv"), row.names = F)
}



###########################################################################
# using pytorch for this task! The data is too huge for batch training!!
###########################################################################
# here we load the data from pytorch for the results here
for(div in 1:10){
  print(paste0("division id: ", div))
  max_predict <- 11
  testing_data <- read.csv(paste0("~/Desktop/comorbidity/Multi-morbidity_biobank/age_sex_data/testing",div,"_age_sex_Phecode.csv"))
  # load the weights learned from python and perform the analysis
  # first two columns are the weights and the last column is the bias
  weights <- read.csv(paste0("~/Desktop/comorbidity/Multi-morbidity_biobank/age_sex_data/weights",div,"_age_sex_Phecode.csv"), header = F)
  testing_data <- testing_data %>%
    group_by(eid) %>%
    arrange(age_diag, .by_group = TRUE)
  
  prediction_onebyone_agesex <- list()
  
  for(predict_num in 1:max_predict){
    testing_eid_included <- testing_data %>% # set incidence number threshold for inclusion in testing set
      group_by(eid) %>%
      summarise(n()) %>%
      filter(`n()` >= predict_num) %>%
      select(eid)
    
    testing_data_included <- testing_eid_included %>%
      left_join(testing_data, by = "eid")
    
    # for the last iteration, just predict all left
    if(predict_num == max_predict){
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        slice(predict_num:n()) %>%
        ungroup
    }else{
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        slice(predict_num) %>%
        ungroup
    }
    # exclude those disease that have happened
    estimating_test_set <- testing_data_included %>%
      group_by(eid) %>%
      slice(1:(predict_num-1)) 
    estimate_eid <- estimating_test_set %>%
      group_by(eid) %>%
      summarise(n())
    predicting_test_set <- predicting_test_set %>%
      mutate(eid = match(eid, estimate_eid$eid)) 
    
    estimating_test_set <- estimating_test_set %>%
      group_split()
    
    # compute the average rank: for both the prediction and those used in estimation
    beta_rank_all <- sapply(1:dim(predicting_test_set)[1], function(x) rank( - (weights$V3 + weights$V1 * predicting_test_set$age_diag[x] + 
                                                                                 weights$V2 * predicting_test_set$sex[x]))[c(predicting_test_set$diag_icd10[x],
                                                                                                                      estimating_test_set[[predicting_test_set$eid[x]]]$diag_icd10) + 1] )  # add one is to fix the index difference between python and R
  
    
    # this is a complicated step: we want to remove the ICD-10 codes that have been used in estimation 
    beta_rank <- beta_rank_all[1, ] - sapply(1:dim(beta_rank_all)[2], function(n) 
      sum(beta_rank_all[1,n] > beta_rank_all[2:dim(beta_rank_all)[1],n]))
    
    
    # compute the specificity: first step remove the number of diseases used in estimation
    total_number <- dim(weights)[1]
    predictions <- list()
    predictions$mean_rank <- mean(beta_rank)/total_number
    predictions$mean_rank_top1 <- mean( beta_rank <= total_number/100) 
    predictions$mean_rank_top2 <- mean( beta_rank <= total_number/50) 
    predictions$mean_rank_top5 <- mean( beta_rank <= total_number/20) 
    predictions$mean_rank_top10 <- mean( beta_rank <= total_number/10) 
    
    prediction_onebyone_agesex[[predict_num]] <- predictions
  }
  save(prediction_onebyone_agesex, file = paste0("../Results/prediction_onebyone_age_sex_Phecode",div, ".RData"))
}

###############################################
# making prediction purely based on frequency
# plot prediction accurracy with all different profiles
###############################################
# order the diseases by the frequency
rec_data <- read.csv("DiseaseAbove1000occur_include_death_PheCode.csv")
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
ds_list$rank <- rank(-ds_list$occ)
for(div in 1:10){
  print(paste0("division id: ", div))
  max_predict <- 11
  testing_data <- read.csv(paste0("~/Desktop/comorbidity/Multi-morbidity_biobank/age_sex_data/testing",div,"_age_sex_Phecode.csv"))
  # load the weights learned from python and perform the analysis
  # first two columns are the weights and the last column is the bias
  testing_data <- testing_data %>%
    group_by(eid) %>%
    arrange(age_diag, .by_group = TRUE)
  
  prediction_onebyone_freq <- list()
  
  for(predict_num in 1:max_predict){
    testing_eid_included <- testing_data %>% # set incidence number threshold for inclusion in testing set
      group_by(eid) %>%
      summarise(n()) %>%
      filter(`n()` >= predict_num) %>%
      select(eid)
    
    testing_data_included <- testing_eid_included %>%
      left_join(testing_data, by = "eid")
    
    # for the last iteration, just predict all left
    if(predict_num == max_predict){
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        slice(predict_num:n()) %>%
        ungroup
    }else{
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        slice(predict_num) %>%
        ungroup
    }
    # exclude those disease that have happened
    estimating_test_set <- testing_data_included %>%
      group_by(eid) %>%
      slice(1:(predict_num-1)) 
    estimate_eid <- estimating_test_set %>%
      group_by(eid) %>%
      summarise(n())
    predicting_test_set <- predicting_test_set %>%
      mutate(eid = match(eid, estimate_eid$eid)) 
    
    estimating_test_set <- estimating_test_set %>%
      group_split()
    
    # compute the average rank: for both the prediction and those used in estimation
    beta_rank_all <- sapply(1:dim(predicting_test_set)[1], function(x) ds_list$rank[c(predicting_test_set$diag_icd10[x],
                                                                          estimating_test_set[[predicting_test_set$eid[x]]]$diag_icd10)] )  # add one is to fix the index difference between python and R
    
    
    # this is a complicated step: we want to remove the ICD-10 codes that have been used in estimation 
    beta_rank <- beta_rank_all[1, ] - sapply(1:dim(beta_rank_all)[2], function(n) 
      sum(beta_rank_all[1,n] > beta_rank_all[2:dim(beta_rank_all)[1],n]))
    
    # compute the specificity: first step remove the number of diseases used in estimation
    total_number <- dim(ds_list)[1]
    predictions <- list()
    predictions$mean_rank <- mean(beta_rank)/total_number
    predictions$mean_rank_top1 <- mean( beta_rank <= total_number/100) 
    predictions$mean_rank_top2 <- mean( beta_rank <= total_number/50) 
    predictions$mean_rank_top5 <- mean( beta_rank <= total_number/20) 
    predictions$mean_rank_top10 <- mean( beta_rank <= total_number/10) 
    
    prediction_onebyone_freq[[predict_num]] <- predictions
  }
  save(prediction_onebyone_freq, file = paste0("../Results/prediction_onebyone_freq_div",div, ".RData"))
}


# plot the prediction power with two baseline: one with age and sex; one with simple LDA
source("topic_functions.R")
rep_number <- 10
P <- 5
df_predict_lik_P_K <- data_frame(percentile = as.character(),type = as.character(),df_K = as.integer(), specificity = as.numeric(), record_number = as.integer())
K <- 10
for(rep_id in 1:rep_number){
  try({load(paste0("../Results/prediction_onebyone_age_K",K,"_P",P,"_rep",rep_id, ".RData"))
    prediction_onebyone <- prediction_onebyone_rslt
    load(paste0("../Results/prediction_onebyone_baseline_K",K,"_rep",rep_id, ".RData"))
    load(paste0("../Results/prediction_onebyone_age_sex_Phecode",rep_id, ".RData"))
    load(paste0("../Results/prediction_onebyone_freq_div",rep_id, ".RData"))
    for(rec_num in 2:10){
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "1_percentile", type = "age", df_K = K, specificity = prediction_onebyone[[rec_num]]$mean_rank_top1,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "2_percentile", type = "age",df_K = K, specificity = prediction_onebyone[[rec_num]]$mean_rank_top2,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "5_percentile", type = "age",df_K = K, specificity = prediction_onebyone[[rec_num]]$mean_rank_top5,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "10_percentile", type = "age",df_K = K, specificity = prediction_onebyone[[rec_num]]$mean_rank_top10,record_number = rec_num)
      # baseline LDA
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "1_percentile", type = "LDAbaseline", df_K = K, specificity = prediction_onebyone_rslt[[rec_num]]$mean_rank_top1,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "2_percentile", type = "LDAbaseline", df_K = K, specificity = prediction_onebyone_rslt[[rec_num]]$mean_rank_top2,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "5_percentile", type = "LDAbaseline", df_K = K, specificity = prediction_onebyone_rslt[[rec_num]]$mean_rank_top5,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "10_percentile", type = "LDAbaseline", df_K = K, specificity = prediction_onebyone_rslt[[rec_num]]$mean_rank_top10,record_number = rec_num)
      # baseline softmax regression
      df_predict_lik_P_K  <- df_predict_lik_P_K %>%
        add_row(percentile = "1_percentile", type = "agesex_baseline", df_K = K, specificity = prediction_onebyone_agesex[[rec_num]]$mean_rank_top1,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>%
        add_row(percentile = "2_percentile", type = "agesex_baseline", df_K = K, specificity = prediction_onebyone_agesex[[rec_num]]$mean_rank_top2,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>%
        add_row(percentile = "5_percentile", type = "agesex_baseline", df_K = K, specificity = prediction_onebyone_agesex[[rec_num]]$mean_rank_top5,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>%
        add_row(percentile = "10_percentile", type = "agesex_baseline", df_K = K, specificity = prediction_onebyone_agesex[[rec_num]]$mean_rank_top10,record_number = rec_num)
      # baseline freqency
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "1_percentile", type = "freq_baseline", df_K = K, specificity = prediction_onebyone_freq[[rec_num]]$mean_rank_top1,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "2_percentile", type = "freq_baseline", df_K = K, specificity = prediction_onebyone_freq[[rec_num]]$mean_rank_top2,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "5_percentile", type = "freq_baseline", df_K = K, specificity = prediction_onebyone_freq[[rec_num]]$mean_rank_top5,record_number = rec_num)
      df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
        add_row(percentile = "10_percentile", type = "freq_baseline", df_K = K, specificity = prediction_onebyone_freq[[rec_num]]$mean_rank_top10,record_number = rec_num)
    }
  })
}

# sumarise the mean and se
df_predict_se <- df_predict_lik_P_K %>%
  group_by(percentile, type, record_number) %>%
  summarise( se_specificity = sd(specificity)/sqrt(n())) %>%
  ungroup()

df_predict_mean <- df_predict_lik_P_K %>%
  group_by(percentile, type, record_number) %>%
  summarise( mean_specificity = mean(specificity)) %>%
  ungroup()

df_plot <- df_predict_se %>%
  left_join(df_predict_mean, by = c("percentile", "type", "record_number")) 

plt <- list()
# plot a box plot for selecting best topic number and degree of freedom
for(types in 1:4){
  percentile_range <- c("1_percentile", "2_percentile", "5_percentile", "10_percentile")[types]
  # intcpt <- c(0.1395,0.2030,0.3558,0.4932)[types]
  df_predict_age <- df_plot %>% filter(percentile == percentile_range,type == "age")
  df_predict_baseline <- df_plot %>% filter(percentile == percentile_range,type == "LDAbaseline")
  df_agesex_baseline <- df_plot %>% filter(percentile == percentile_range,type == "agesex_baseline")
  df_freq_baseline <- df_plot %>% filter(percentile == percentile_range,type == "freq_baseline")
  plt[[types]] <- ggplot(df_predict_age, aes(x = record_number)) +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.position = "none",panel.background=element_blank(),plot.title = element_text(size = 10, face = "bold"))  + 
    labs(title = paste0("Percentile:", percentile_range), x = "nth records", y = "sensitivity") + 
    geom_line(aes(y = mean_specificity), linetype = "dashed", color = red) + 
    geom_pointrange(aes(y = mean_specificity, ymin = mean_specificity - 1.96*se_specificity, ymax = mean_specificity + 1.96*se_specificity), size = .5, color = red) +
    geom_line(data = df_predict_baseline ,aes( x = record_number, y = mean_specificity), linetype = "dashed", color = blue) + 
    geom_pointrange(data = df_predict_baseline, aes(x = record_number, y = mean_specificity, ymin = mean_specificity - 1.96*se_specificity, ymax = mean_specificity + 1.96*se_specificity), size = .5, color = blue) + 
    geom_line(data = df_agesex_baseline ,aes( x = record_number, y = mean_specificity), linetype = "dashed", color = orange) + 
    geom_pointrange(data = df_agesex_baseline, aes(x = record_number, y = mean_specificity, ymin = mean_specificity - 1.96*se_specificity, ymax = mean_specificity + 1.96*se_specificity), size = .5, color = orange) + 
    theme(panel.background=element_blank()) +
    geom_line(data = df_freq_baseline ,aes( x = record_number, y = mean_specificity), linetype = "dashed", color = grey) + 
    geom_pointrange(data = df_freq_baseline, aes(x = record_number, y = mean_specificity, ymin = mean_specificity - 1.96*se_specificity, ymax = mean_specificity + 1.96*se_specificity), size = .5, color = grey) 
    # geom_hline(yintercept = intcpt, linetype = "dashed", color = red)# add three dashed lines
  ggsave(paste0("~/Desktop/comorbidity/figures/agelda_baseline_comparison_",percentile_range,"_K",K,"P",P,".png"), 
         plt[[types]], width = 4, height = 4)
}

####################################################################################
# compute recall rate as a function of percentile of diseases selected
####################################################################################
source("topic_functions.R")
source("plotting_functions.R")
ds_list <- read.csv("listAbove500include_deaths_ICDA2N.csv")
K <- 12
P <- 7
rep_number <- 10
DIR <- "~/Desktop/comorbidity/Results/"
# use ageLDA prediction rank data
pt <- paste0("^rank_of_target_disease", K,"_P", P, "_rep*")
temp <- list.files(paste(DIR, sep=""), pattern=pt)
rank_prediction <- list()
for(rep_id in 1:length(temp)){
  load(paste0(DIR,temp[[rep_id]]))
  rank_prediction[[rep_id]] <- sapply(1:dim(ds_list)[1], function(x) mean(collect_prediction_ranks <= x))
}
recall_ageLDA <- Reduce("+", rank_prediction)/length(rank_prediction)
# use LDA prediction rank data
pt <- paste0("^baselineLDA_rank_of_target_disease", K,"_P_rep*")
temp <- list.files(paste(DIR, sep=""), pattern=pt)
rank_prediction <- list()
for(rep_id in 1:length(temp)){
  load(paste0(DIR,temp[[rep_id]]))
  rank_prediction[[rep_id]] <- sapply(1:dim(ds_list)[1], function(x) mean(collect_prediction_ranks <= x))
}
recall_baselineLDA <- Reduce("+", rank_prediction)/length(rank_prediction)

# get the recall for the other two methods
# here we load the data from pytorch for the results here
ds_list <- read.csv("listAbove500include_deaths_ICDA2N.csv")
ds_list$rank <- rank(-ds_list$occ)
recall_freq <- list()
recall_agesex <- list()
for(div in 1:10){
  print(paste0("division id: ", div))
  max_predict <- 11
  testing_data <- read.csv(paste0("~/Desktop/comorbidity/Multi-morbidity_biobank/age_sex_data/testing",div,"_age_sex_disease.csv"))
  # load the weights learned from python and perform the analysis
  # first two columns are the weights and the last column is the bias
  weights <- read.csv(paste0("~/Desktop/comorbidity/Multi-morbidity_biobank/age_sex_data/weights",div,"_age_sex_predict.csv"), header = F)
  testing_data <- testing_data %>%
    group_by(eid) %>%
    arrange(age_diag, .by_group = TRUE)
  
  collect_beta_rank_agesex <- c()
  collect_beta_rank_freq <- c()
  
  for(predict_num in 1:max_predict){
    testing_eid_included <- testing_data %>% # set incidence number threshold for inclusion in testing set
      group_by(eid) %>%
      summarise(n()) %>%
      filter(`n()` >= predict_num) %>%
      select(eid)
    
    testing_data_included <- testing_eid_included %>%
      left_join(testing_data, by = "eid")
    
    # for the last iteration, just predict all left
    if(predict_num == max_predict){
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        slice(predict_num:n()) %>%
        ungroup
    }else{
      predicting_test_set <- testing_data_included %>%
        group_by(eid) %>%
        slice(predict_num) %>%
        ungroup
    }
    # exclude those disease that have happened
    estimating_test_set <- testing_data_included %>%
      group_by(eid) %>%
      slice(1:(predict_num-1)) 
    estimate_eid <- estimating_test_set %>%
      group_by(eid) %>%
      summarise(n())
    predicting_test_set <- predicting_test_set %>%
      mutate(eid = match(eid, estimate_eid$eid)) 
    
    estimating_test_set <- estimating_test_set %>%
      group_split()
    
    # compute the average rank: for both the prediction and those used in estimation
    beta_rank_all <- sapply(1:dim(predicting_test_set)[1], function(x) rank( - (weights$V3 + weights$V1 * predicting_test_set$age_diag[x] + 
                                                                                  weights$V2 * predicting_test_set$sex[x]))[c(predicting_test_set$diag_icd10[x],
                                                                                                                             estimating_test_set[[predicting_test_set$eid[x]]]$diag_icd10) + 1] )  # add one is to fix the index difference between python and R
    
    # this is a complicated step: we want to remove the ICD-10 codes that have been used in estimation 
    beta_rank <- beta_rank_all[1, ] - sapply(1:dim(beta_rank_all)[2], function(n) 
      sum(beta_rank_all[1,n] > beta_rank_all[2:dim(beta_rank_all)[1],n]))
    
    collect_beta_rank_agesex <- c(collect_beta_rank_agesex, beta_rank) 
    # frequency data
    beta_rank_all <- sapply(1:dim(predicting_test_set)[1], function(x) ds_list$rank[c(predicting_test_set$diag_icd10[x],
                                                                                      estimating_test_set[[predicting_test_set$eid[x]]]$diag_icd10)] )  # add one is to fix the index difference between python and R
    
    
    # this is a complicated step: we want to remove the ICD-10 codes that have been used in estimation 
    beta_rank <- beta_rank_all[1, ] - sapply(1:dim(beta_rank_all)[2], function(n) 
      sum(beta_rank_all[1,n] > beta_rank_all[2:dim(beta_rank_all)[1],n]))
    collect_beta_rank_freq <- c(collect_beta_rank_freq, beta_rank) 
  }
  recall_agesex[[div]] <- sapply(1:dim(ds_list)[1], function(x) mean(collect_beta_rank_agesex <= x))
  recall_freq[[div]] <- sapply(1:dim(ds_list)[1], function(x) mean(collect_beta_rank_freq <= x))
}
recall_agesex <- Reduce("+", recall_agesex)/length(recall_agesex)
recall_freq <- Reduce("+", recall_freq)/length(recall_freq)

risk_percentile <- 1:dim(ds_list)[1]/dim(ds_list)[1]
df_recall <- data.frame(recall_ageLDA, recall_baselineLDA, risk_percentile, recall_agesex, recall_freq)

plt <- ggplot(df_recall) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  theme(legend.position = "none",plot.title = element_text(size = 10, face = "bold"))  + 
  labs(title = "Recall v.s. Percentile of predicted risk selected", x = "Risk percentile", y = "Recall") + 
  geom_line(aes(x= risk_percentile, y= recall_baselineLDA), color = blue) + 
  geom_line(aes(x= risk_percentile, y= recall_ageLDA), color = red) + 
  geom_line(aes(x= risk_percentile, y= recall_agesex), color = orange) + 
  geom_line(aes(x= risk_percentile, y= recall_freq), color = grey) + 
  geom_abline(intercept = 0, slope = 1, color = grey)  

ggsave(paste0("~/Desktop/comorbidity/figures/recall_curve_comparison_baselines.png"), plt, width = 6, height = 4)




################################################################
# prediction fore each disease using PheRS
################################################################
source("topic_functions.R")
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode"))

survive_age <- read.csv("survive_age.csv")
load(file = paste0("../Results/training_2rec_age_K10_P5_rep2.RData"))
all_eid <- rec_data %>%
  group_by(eid) %>%
  summarise() 
training_eid <- data.frame(eid = para$eid)
testing_eid <- all_eid %>% 
  anti_join(training_eid, by = "eid")

testing_data <- testing_eid %>%  
  left_join(rec_data, by = "eid")


PheRS_rslt <- prediction_PheRS_by_disease(testing_data, ds_list, para) 

LASSO_rslt <- LASSO_predict(rec_data, para)

rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode"))

survive_age <- read.csv("survive_age.csv")

load(file = paste0("../Results/training_2rec_age_K10_P5_rep2.RData"))

# here, use the loadings of learned topics to train another predictor. 
# first order the incidences by age
all_eid <- rec_data %>%
  group_by(eid) %>%
  summarise() 
####################################################
# to do: for better performance: each loading (variable) should be normalised.
####################################################
topic_loadings <- sweep(para$alpha_z, 1, rowSums(para$alpha_z), FUN="/")
# topic_loadings <- para$alpha_z
training_eid <- data.frame(eid = para$eid, topic_loadings)
training_data <- training_eid %>%  
  left_join(rec_data, by = "eid")

testing_eid <- all_eid %>% 
  anti_join(training_eid, by = "eid")

testing_data <- testing_eid %>%  
  left_join(rec_data, by = "eid") %>%
  group_by(eid) %>%
  arrange(age_diag, .by_group = TRUE)
para_training <- para

#####################################
######################################
#######################################
# we will be comparing three methods at this point:
# methods 1: logistic regression with individual topic loadings + age
# methods 2: logistic regression with loadings of prior diseases + age
# methods 3: logistic regression with topic loadings + age interaction
# methods 3: logistic regression with loadings of prior diseases + age interaction
AUC_methods <- matrix(nrow = para_training$D, ncol = 6)
# for(j in 10:para$D){
for(j in 1:para$D){
  # j <- 18
  ds_id <- para$list_above500occu$diag_icd10[j]
  #####################
  # create training data
  #####################
  cases_training_eid <- training_data %>% # select all the cases
    filter(diag_icd10 == ds_id) %>%
    rename(target_age = age_diag) %>%
    select(eid, target_age) 
  training_cases_data <- cases_training_eid %>% # select all disease that has a diagnosis that are earlier than target!
    left_join(training_data, by = "eid") %>% 
    filter(age_diag < target_age) %>% 
    group_by(eid) %>%
    arrange(desc(age_diag), .by_group = T) %>%
    group_by(eid) %>%
    slice(1)
  
  training_control_eid <- training_data %>% 
    anti_join(cases_training_eid, by = "eid") %>%
    group_by(eid) %>%
    slice(1) %>%
    select(eid)
  
  training_control_eid$target_age <- sample(training_cases_data$target_age, dim(training_control_eid)[1], replace = T)
  training_control_data <- training_control_eid %>% # select all disease that has a diagnosis that are earlier than target!
    left_join(training_data, by = "eid") %>% 
    filter(age_diag < target_age) %>%
    group_by(eid) %>%
    arrange(desc(age_diag), .by_group = T) %>%
    group_by(eid) %>%
    slice(1)
  
  training_set <- bind_rows(training_cases_data, training_control_data) %>%
    select(-target_age, -diag_icd10)
  
  training_cases_eid <- training_cases_data %>%
    group_by(eid) %>%
    slice_tail(n = 1) %>%
    mutate(outcome = 1) %>%
    select(eid, outcome, age_diag)
  training_control_eid <- training_control_data %>%
    group_by(eid) %>%
    slice_tail(n = 1) %>%
    mutate(outcome = 0) %>%
    select(eid, outcome, age_diag)
  training_eid <- bind_rows(training_cases_eid, training_control_eid)
  
  # create the non_normalised z_n data set
  target_age <- bind_rows(training_cases_data, training_control_data) %>%
    select(eid, target_age) %>%
    group_by(eid) %>%
    slice_tail(n = 1) 
  age_filter <- data.frame(eid = para$eid) %>% # create an age filter that are of the order in para$eid
    left_join(target_age, by = "eid")
  # use a helper function to compute the sum of loadings
  helper_age_filtering <- function(s, para, age_filter){
    ids <- which(para$w[[s]]$age_diag <  age_filter$target_age[s])
    if(!length(ids)){
      return(rep(0, para$K))
    }else{
      # return(colMeans(para$E_zn[[s]][ids, , drop =F]))
      return(colSums(para$E_zn[[s]][ids, , drop =F]))
    }
  }
  loading_sums <- sapply(1:para$M, function(s) helper_age_filtering(s, para, age_filter)) %>% t
  loading_sums <- data_frame(eid = para$eid, loadings = loading_sums)
  # loading_sums <- data_frame(eid = para$eid, loadings = para$alpha_z - para$alpha)
  loading_sums <- training_eid %>%
    left_join(loading_sums, by = "eid")
  
  # method 1 normal regression outperformed lasso? 
  x_train <- as.matrix(training_set[,2:12])
  y_train <- training_eid$outcome
  glmfitting <- glm(y_train ~ x_train, family = "binomial")
  auc(y_train, predict(glmfitting))
  
  cvfit_testing <- cv.glmnet(x_train, y_train, alpha = 1, family = "binomial")
  auc(y_train, predict(cvfit_testing, newx = x_train, s = "lambda.min"))
  
  # Method 2: now perform the logistic regression over the loadingsum
  loadingsums_train <- as.matrix(loading_sums[,3:4])
  loadingsums_y_train <- loading_sums$outcome
  loadingsumfitting <- glm(loadingsums_y_train ~ loadingsums_train, family = "binomial")
  auc(loadingsums_y_train, predict(loadingsumfitting))
  
  # method 3: incorporate age interaction
  age_threshold <- 1*(training_set$age_diag < median(training_set$age_diag))
  age_itr <- sapply(2:11, function(x) as.matrix(training_set[,x]) * age_threshold)
  x_train <- as.matrix(cbind(training_set[,2:12], age_itr))
  y_train <- training_eid$outcome
  age_itr_fitting <- glm(y_train ~ x_train, family = "binomial")
  auc(y_train, predict(age_itr_fitting))
  
  normalised_age <- (training_set$age_diag - mean(training_set$age_diag))/sd(training_set$age_diag)
  age_itr <- sapply(2:11, function(x) as.matrix(training_set[,x]) * normalised_age)
  x_train <- as.matrix(cbind(training_set[,2:11], normalised_age, age_itr))
  age_itr_ridge <- cv.glmnet(x_train, y_train, alpha = 1, family = "binomial")
  auc(y_train, predict(age_itr_ridge, newx = x_train, s = "lambda.min"))
  
  # methods 4: incoporate age interaction with loadingsums
  age_itr <- sapply(1:10, function(x) as.matrix(loading_sums[,4])[,x] * training_set$age_diag)
  x_train <- cbind(as.matrix(loading_sums[,3:4]), age_itr)
  y_train <- loading_sums$outcome
  age_itr2_fitting <- glm(y_train ~ x_train, family = "binomial")
  auc(loadingsums_y_train, predict(age_itr2_fitting))
  
  #####################
  # create testing data
  #####################
  cases_testing_eid <- testing_data %>% # select all the cases
    filter(diag_icd10 == ds_id) %>%
    rename(target_age = age_diag) %>%
    select(-diag_icd10)
  
  cases_testing_data <- cases_testing_eid %>% # select all disease that has a diagnosis that are earlier than target!
    left_join(testing_data, by = "eid") %>% 
    filter(age_diag < target_age)
  
  # exclude cases that have no prior diagnosis 
  # (i.e. target disease happen to be the first one; no way to predict!)
  targe_age_distribution <- cases_testing_data %>%
    group_by(eid) %>%
    slice(1) 
  control_testing_eid <- testing_data %>% 
    anti_join(cases_testing_eid, by = "eid") %>%
    group_by(eid) %>%
    slice(1) %>%
    select(eid)
  control_testing_eid$target_age <- sample(targe_age_distribution$target_age, dim(control_testing_eid)[1], replace = T)
  
  control_testing_data <- control_testing_eid %>% # select all disease that has a diagnosis that are earlier than target!
    left_join(testing_data, by = "eid") %>% 
    filter(age_diag < target_age)
  
  # save the case and control outcomes; save the age of the last event before target for computing area
  case_testing_eid <- cases_testing_data %>%
    group_by(eid) %>%
    arrange(desc(age_diag), .by_group = T) %>%
    slice(1) %>%
    mutate(outcome = 1) %>%
    select(eid, outcome, age_diag)
  control_testing_eid <- control_testing_data %>%
    group_by(eid) %>%
    arrange(desc(age_diag), .by_group = T) %>%
    slice(1) %>%
    mutate(outcome = 0) %>%
    select(eid, outcome, age_diag)
  outcomes <- bind_rows(case_testing_eid, control_testing_eid) 
  # below is how we estimate the area under curve as well as the 
  estimating_test_set <- bind_rows(cases_testing_data, control_testing_data) %>%
    select(- target_age) %>%
    ungroup()
  
  if(is.null(para_training$P)){ # using P to determine if age is included
    para_training$pi_beta_basis <- array( rep(para_training$beta, each = 81),dim =  c(81, para_training$D, para_training$K) )
    # for both treeLDA and baseline lda, we only need to initialise baseline case here
    para_testing <- topic_init_baseline(estimating_test_set, ds_list, para_training$K)
  }else{
    para_testing <- topic_init_age(estimating_test_set, ds_list, para_training$K, para_training$P)
  }
  # assigning beta to the testing set
  para_testing$beta_w_full <- apply(para_training$pi_beta_basis, 3, function(x) 
    x[as.matrix(select(para_testing$unlist_Ds_id, age_diag, Ds_id))]) 
  para_testing$beta_w <- lapply(para_testing$patient_lst, function(x) para_testing$beta_w_full[x,,drop=F])
  
  # updating Z_n
  para_testing$max_itr <- 10
  para_testing$alpha <- para_training$alpha 
  para_testing$lb <- data.frame("Iteration" = 0,"Lower_bound" = CVB_lb(para_testing))
  para_testing$tol <- 10^(-6)
  for(itr in 1:para_testing$max_itr){
    print(paste0("Interation: ",itr))
    para_testing <- CVB0_E_zn(para_testing) # we choose CVB0 as papers shown it could converge quicker
    para_testing$lb[nrow(para_testing$lb) + 1,] <- c(itr, CVB_lb(para_testing))
    curr_lb <- pull(filter(para_testing$lb, Iteration == itr), Lower_bound) 
    prev_lb <- pull(filter(para_testing$lb, Iteration == (itr - 1 )), Lower_bound) 
    print(paste0("Current Lower bound ", curr_lb, " at iteration: ",itr))
    try({
      if(is.finite((curr_lb - prev_lb)) & abs(curr_lb - prev_lb)/abs(prev_lb) < para_testing$tol ){
        print(paste0("Optimization converged at step ", itr))
        break
      }
    })
  }
  alpha_z <- para_testing$alpha_z
  alpha_z <- sweep(alpha_z, 1, rowSums(alpha_z), FUN="/")
  # normalise alpha_z by row to get estimation of theta
  loadings <- sweep(para_testing$alpha_z, 1, rowSums(para_testing$alpha_z), FUN="/")
  data_testing <- data.frame(eid = para_testing$eid, loadings) %>%
    left_join(outcomes, by = "eid")
  
  # method 1
  predict_x <- cbind(rep(1, dim(data_testing)[1]),as.matrix(data_testing[,c(2:10, 13)]))
  scores <- predict_x %*% glmfitting$coefficients[c(1:10,12)]
  AUC_methods[j,1] <- auc(data_testing$outcome, scores)
  
  x_predict <- as.matrix(data_testing[,c(2:11, 13)])
  auc(data_testing$outcome, predict(cvfit_testing, newx = x_predict, s = "lambda.min"))

  
  # Method 2: also compute it without normalisation
  data_alpha_z <- data.frame(eid = para_testing$eid, alpha_z ) %>%
    left_join(outcomes, by = "eid")
  predict_x <- cbind(rep(1, dim(data_alpha_z)[1]),as.matrix(data_alpha_z[,c(13, 2:11)]))
  scores <- predict_x %*% loadingsumfitting$coefficients
  AUC_methods[j,2] <- auc(data_alpha_z$outcome, scores)
  
  # Method 3: compute the with interaction
  age_threshold <- 1*(data_testing$age_diag < median(training_set$age_diag))
  age_itr <- sapply(2:11, function(x) as.matrix(data_testing[,x]) * age_threshold)
  predict_x <- cbind(rep(1, dim(data_testing)[1]),as.matrix(data_testing[,c(2:10, 13)]), as.matrix(age_itr)) 
  scores <- predict_x %*% age_itr_fitting$coefficients[c(1:10,12:22)] # jump through one of them
  AUC_methods[j,3] <- auc(data_testing$outcome, scores) 
  
  normalised_age <- (data_testing$age_diag - mean(training_set$age_diag))/sd(training_set$age_diag)
  age_itr <- sapply(2:11, function(x) as.matrix(data_testing[,x]) * normalised_age)
  predict_x <- as.matrix(cbind(data_testing[,2:11], normalised_age, age_itr))
  auc(data_testing$outcome, predict(age_itr_ridge, newx = predict_x, s = "lambda.min"))
  
  # Method 4: age interaction with loading sum
  age_itr <- sapply(2:11, function(x) as.matrix(data_alpha_z[,x]) * data_alpha_z[,13])
  predict_x <- cbind(rep(1, dim(data_alpha_z)[1]),as.matrix(data_alpha_z[,c(13, 2:11)]),age_itr)
  scores <- predict_x %*% age_itr2_fitting$coefficients
  AUC_methods[j,4] <- auc(data_alpha_z$outcome, scores)
  
  # age_itr <- sapply(2:11, function(x) as.matrix(data_testing[,x]) * age_threshold  ) # data_testing$age_diag)
  # predict_x <- cbind(as.matrix(data_testing[,c(2:11, 13)]), as.matrix(age_itr))
  # scores <- predict(age_itr_ridge, newx = predict_x, s = "lambda.min")
  # auc(data_testing$outcome, scores)
  
  print(paste0(j, " ", ds_id, " AUC loadings ", AUC_methods[j,1], "; Non_normalised ", AUC_methods[j,2], "; age interaction ", AUC_methods[j,3],
               "; age interaction 2", AUC_methods[j,4]))
}

plot(x = df_auc$prediction_cumulative_loadings, y = df_auc$preiction_loadings)
abline(coef = c(0,1))
df_auc <- data.frame(dieases = ds_list$diag_icd10, preiction_loadings = AUC_methods[,1],
                     prediction_cumulative_loadings = AUC_methods[,2], phenotype = ds_list$phenotype)
write.csv(x = df_auc, file = "AUC_loading_prediction.csv", row.names = F)
ggplot(data = df_auc) + 
  geom_point(aes(x = prediction_cumulative_loadings, y = preiction_loadings), color = grey) + 
  geom_abline(slope = 1)
AUC_loadings <- read.csv("AUC_loading_prediction.csv")
a <- load("../Results/LASSO_prediction10_P3_rep9.RData")
df_auc$LASSO <- AUC_per_ds
ggplot(data = df_auc) + 
  geom_point(aes(x = LASSO, y = preiction_loadings), color = grey) + 
  geom_abline(slope = 1)
###############################################
# plot single disease prediction accuracy
###############################################
source("topic_functions.R")
source("plotting_functions.R")
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode"))
K <- 10
P <- 5
DIR <- "~/Desktop/comorbidity/Results/"
pt <- paste0("^single_diseae_prediction", K,"_P", P, "_rep*")
temp <- list.files(paste(DIR, sep=""), pattern=pt)
AUCs <- list()
p_values <- list()
logORs <- list()
for(rep_id in 1:length(temp)){
  load(paste0(DIR,temp[[rep_id]]))
  AUCs[[rep_id]]  <- PheRS_rslt[[2]][,1]
  # flip the wrong ones
  AUCs[[rep_id]][which(PheRS_rslt[[1]][,1] < 0 )] <- 1 - AUCs[[rep_id]][which(PheRS_rslt[[1]][,1] < 0 )]
  p_values[[rep_id]] <- PheRS_rslt[[1]][,4]
  logORs[[rep_id]] <- PheRS_rslt[[1]][,1]
}
# per-sd logRR
AUCs <- do.call(cbind,AUCs)
mean_AUCs <- rowMeans(AUCs)
se_AUCs <- apply(AUCs, 1, function(x) sd(x)/sqrt(dim(AUCs)[2]))

# p-values
p_values <- do.call(cbind,p_values)
mean_p <- rowMeans(p_values)
se_p <- apply(p_values, 1, function(x) sd(x)/sqrt(dim(p_values)[2]))

# per-sd logRR
logORs <- do.call(cbind,logORs)
mean_ORs <- rowMeans(logORs)
se_ORs <- apply(logORs, 1, function(x) sd(x)/sqrt(dim(logORs)[2]))

df_mht_plt <- data.frame(disease = ds_list$diag_icd10, occ = ds_list$occ, mean_AUCs, se_AUCs, mean_p, se_p, mean_ORs, se_ORs)
ggplot(df_mht_plt, aes(x = 1:length(disease))) + 
  geom_pointrange(aes(y = mean_AUCs, ymin = mean_AUCs + 1.96 * se_AUCs, ymax = mean_AUCs - 1.96 * se_AUCs), size = .2) +
  labs(x = "disease", y= "AUC")

ggplot(df_mht_plt) + 
  geom_point(aes(x = mean_ORs,y = -log(mean_p)), size = .5) +
  labs(x = "effect size", y= "- log p-value") +
  lims(x = c(-1,1))

ds_list$phenotype[which(mean_AUCs >0.7)]



