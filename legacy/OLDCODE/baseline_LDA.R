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
  mutate(birth_year=as.Date(paste(X34.0.0,"01-01", sep="-"))) %>% 
  mutate(age_diag = difftime(epistart, birth_year, units = "days")/365) %>%
  filter(!is.na(age_diag))

ds_occ_thre <- 100 # > 0.5%
# find all of the ICD-10 that is above 500 occurence (> 0.1%)
list_above500occu <- new_data %>%
  select(eid, diag_icd10, age_diag) %>%
  filter(str_detect(diag_icd10, "^[A-Q]")) %>%
  group_by(eid, diag_icd10) %>% 
  slice(1) %>%
  group_by(diag_icd10) %>%
  summarise(occ = n()) %>% 
  filter(occ >= ds_occ_thre)

# only keep the first incidence of disease 
first_icidence_age <- new_data %>%
  select(eid, diag_icd10, age_diag) %>%
  filter(diag_icd10 %in% list_above500occu$diag_icd10) %>%
  group_by(eid, diag_icd10) %>% 
  # slice(1) %>% # just pick the first record and usually it is ranked with age
  filter(n() == 1 | age_diag == min(age_diag) ) %>% # this row is highly optimized, a lot faster the slice_min ### don't change
  slice(1) %>% # avoid the ties in min
  ungroup() 

matrix_data <- first_icidence_age %>%
  select(eid, diag_icd10) %>%
  mutate(eid = as.factor(eid))

##################################
# save the data for future use
##################################
write.csv(first_icidence_age, paste0("DiseaseAbove",ds_occ_thre,"occur.csv"), row.names = F)

write.csv(list_above500occu, paste0("listAbove",ds_occ_thre,".csv"), row.names = F)

####################################
# prepare data for LDA
####################################
first_icidence_age <- read.csv(paste0("DiseaseAbove500occur.csv"))
list_above500occu <- read.csv(paste0("listAbove500.csv"))

# plot the number distribution of indiviudal diseases
df_number_records <- first_icidence_age %>%
  group_by(eid) %>%
  summarise(n())

############################################
# plotting record number per individual
############################################
df_simu_pois <- data.frame(num_records = floor(1+rexp(370545, 1/(-0.5+5.87))))
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

##########################
# load the packages
##########################
source("topic_functions.R")
##############################

##########################
# process the true data
##########################
rec_data <- read.csv("rec2subjectAbove500occur_include_death_ICDA2N.csv")
ds_list <- read.csv("listAbove500include_deaths_ICDA2N.csv")
topic_num <- 10
para <- topic_init_baseline(rec_data, ds_list, topic_num)

############################
# start optimization: mean-field 
############################
# set the number of update 
para$max_itr <- 500
para$lb <- data.frame("Iteration" = 0,"Lower_bound" = comp_lda_lb(para))
para$alpha <- rep(1, para$K)
para$tol <- 10^(-7) # tolerance of lower bound step
para$itr_check_lb <- 1
for(itr in 1:para$max_itr){
  print(paste0("Interation: ",itr))
  para <- comp_E_zn(para)
  para <- comp_E_lntheta(para)
  para <- update_beta_basic_lda(para)
  # print(paste0("Current Lower bound at beta ", comp_lda_lb(para)))
  # para <- update_alpha(para) # we do not need to update alpha
  # print(paste0("Current Lower bound at alpha ", comp_lda_lb(para)))
  para$lb[nrow(para$lb) + 1,] <- c(itr, comp_lda_lb(para))
  if(itr %% para$itr_check_lb  ==0){
    print(paste0("Current Lower bound ", pull(filter(para$lb, Iteration == itr), Lower_bound), " at iteration: ",itr))
    curr_lb <- pull(filter(para$lb, Iteration == itr), Lower_bound) 
    prev_lb <- pull(filter(para$lb, Iteration == (itr -para$itr_check_lb )), Lower_bound) 
    if(is.finite((curr_lb - prev_lb)) & (curr_lb - prev_lb)/abs(prev_lb) < para$tol ){
      print(paste0("Optimization converged at step ", itr))
      break
    }
  }
}


# save the parameter
save(para, file = paste0("~/Desktop/comorbidity/Results/","Run_",Sys.Date(),".RData"))

load(file = paste0("~/Desktop/comorbidity/Results/Run_2021-02-10.RData"))

# plot topics
topics <- data.frame(ds_id = para$list_above500occu$diag_icd10, topics = para$beta)
df_topics <- list()
for(i in 1:para$K){
df_topics[[i]]  <- topics %>% filter(get(paste0("topics.",i)) > 5/para$D) %>%
                      pull(ds_id)
}

topics <- data.frame(ds_id = para$list_above500occu$diag_icd10, topics = para$beta)

library(reshape2)
library(ggplot2)
longData<-melt(para$beta)
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

# visualize the topic proportion: only select those with more than 20 records
idx_high <- which(rowSums(para$alpha_z - 1) > 20)
posti_theta <- para$alpha_z/rowSums(para$alpha_z)

longData<-melt(posti_theta[idx_high[1:40],])
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

# compute the correlation of topic proportion with basic information in birthyear file
birthyear <- read.csv("~/Desktop/genetics_longitudinal_data/longitudinal_data/Year_of_birth.csv") %>%
  select(-X) %>%
  rename(year_birth = X34.0.0, sex = X31.0.0, BMI = X23104.0.0)
df_topic_proportion <- data.frame(eid = para$eid, theta = posti_theta) %>% 
  left_join(birthyear, by="eid")

cor_year_birth <- rcorr(posti_theta, df_topic_proportion$year_birth)
cor_year_birth <- cor_year_birth$r[11,1:10]

cor_sex <- rcorr(posti_theta, df_topic_proportion$sex)
cor_sex <- cor_sex$r[11,1:10]

cor_BMI <- rcorr(posti_theta, df_topic_proportion$BMI)
cor_BMI <- cor_BMI$r[11,1:10]

cor_df <- rbind(cor_year_birth, cor_sex, cor_BMI)^2
colnames(cor_df) <- 1:10
longData<-melt(cor_df)
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))


###################################################
# simulate a topic data set to test the inference
###################################################
para <- list()
# set the data size
para$M <- 20000 # tarting with 10000 subjects
para$Ns <- floor(rexp(para$M, rate = 1/(-0.5+5.87)) + 1)
para$K <- 10 # start with 3 component
para$D <- 50 # at the moment, only simulate 6 diseases

# initiate the variables
para$alpha <- rgamma(para$K, shape = 50, rate = 50)
print(paste0("simulated alpha: ", para$alpha))
# each column is a topic; D*K matrix
para$beta <- matrix(c(0.4,0.4, 0,0,0,0,0,0,0,0,
                      0,0.4,0.4,0,0,0,0,0,0,0,
                      0,0,0.4,0.4,0,0,0,0,0,0) + 0.2/para$D, nrow = para$D, ncol = para$K)

para$theta <- rdirichlet(para$M, para$alpha)
para$zn <- list()
para$w <-list()
for(s in 1:para$M){
  para$zn[[s]] <- sample( 1:para$K, para$Ns[s], replace=TRUE, prob=para$theta[s,])
  para$w[[s]] <- data.frame(Ds_id = sapply(1:para$Ns[s], function(n) sample( 1:para$D, size= 1, replace=TRUE, 
                                                                             prob=para$beta[,para$zn[[s]][n]]))) %>%
    mutate(eid = s)
}
# build dataframes for computation
para$unlist_Ds_id <- bind_rows(para$w)
para$ds_list <- para$unlist_Ds_id %>%
  mutate(id = row_number()) %>%
  group_by(Ds_id) %>%
  group_split()

para$patient_lst <- para$unlist_Ds_id %>%
  mutate(id = row_number()) %>%
  select(eid, id) %>%
  group_by(eid) %>%
  group_split(keep = F) %>%
  lapply(pull)

# this beta_w parameter save the beta for each word: it is a list of M elements and each contain a K*Ns matrix
# initiate beta
para$eta <- rgamma(para$D,shape = 100, rate = 100)
# each column is a topic; D*K matrix
para$beta <- t(rdirichlet(para$K, para$eta))
para$beta_w <- list()
para$beta_w <- lapply(para$w, function(w) para$beta[w$Ds_id,,drop=FALSE] )
para$beta_w_full <- do.call(rbind, para$beta_w)
# Matrix of M*K
para$E_lntheta <- t(sapply(rgamma(para$M,shape = 100, rate = 100), function(x) x*(digamma(para$alpha) - digamma(sum(para$alpha))) ))

# update E_zn: list of M; each element is matrix of Ns*K
para$E_zn <- list()
para <- comp_E_zn(para)

############################
# run on the simulatd data
############################
# set the number of update 
para$max_itr <- 500
para$lb <- data.frame("Iteration" = 0,"Lower_bound" = comp_lda_lb(para))
para$alpha <- rep(1.01, para$K)
para$tol <- 10^(-6) # tolerance of lower bound step
para$itr_check_lb <- 1
for(itr in 1:para$max_itr){
  print(paste0("Interation: ",itr))
  para <- comp_E_zn(para)
  para <- comp_E_lntheta(para)
  para <- update_beta_basic_lda(para)
  # print(paste0("Current Lower bound at beta ", comp_lda_lb(para)))
  # para <- update_alpha(para) # we do not need to update alpha
  # print(paste0("Current Lower bound at alpha ", comp_lda_lb(para)))
  para$lb[nrow(para$lb) + 1,] <- c(itr, comp_lda_lb(para))
  if(itr %% para$itr_check_lb  ==0){
    print(paste0("Current Lower bound ", pull(filter(para$lb, Iteration == itr), Lower_bound), " at iteration: ",itr))
    if( ( pull(filter(para$lb, Iteration == itr), Lower_bound) - 
            pull(filter(para$lb, Iteration == (itr -para$itr_check_lb )), Lower_bound) )/
        abs(pull(filter(para$lb, Iteration == (itr -para$itr_check_lb )), Lower_bound)) < para$tol ){
      print(paste0("Optimization converged at step ", itr))
      break
    }
  }
}

library(reshape2)
library(ggplot2)

longData<-melt(para$beta)
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

longData<-melt(matrix(c(0.4,0.4, 0,0,0,0,0,0,0,0,
                        0,0,0.4,0.4,0,0,0,0,0,0,
                        0,0,0,0,0.4,0.4,0,0,0,0) + 0.2/para$D, nrow = para$D, ncol = para$K))
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

# visualize the theta and compare with the fitted theta
true_beta <- matrix(c(0.4,0.4, 0,0,0,0,0,0,0,0,
         0,0,0.4,0.4,0,0,0,0,0,0,
         0,0,0,0,0.4,0.4,0,0,0,0) + 0.2/para$D, nrow = para$D, ncol = para$K)

# ordering
rank_beta <- sapply(1:para$K, function(i) which.min(colSums( (para$beta - true_beta[,i] )^2 ) ) ) 
ordered_beta <- para$beta[, rank_beta]

# compare individual topic composition
para$fit_theta <- t(apply(para$alpha_z, 1, function(x) x/sum(x)))[,rank_beta]
 
df_theta_compare <- data.frame(true_theta =para$theta, fit_theta = para$fit_theta)
ggplot(df_theta_compare) + 
  geom_point(aes(x = true_theta.1, y =fit_theta.1), color = red, alpha = .5) +
  # geom_point(aes(x = true_theta.2, y =fit_theta.2), color = blue, alpha = .5) +
  geom_abline(intercept = 0, slope = 1) 

ggplot(sample_n(df_theta_compare,size = 10000)) + 
    geom_point(aes(x = true_theta.1, y =fit_theta.1), color = red, alpha = .5) +
    geom_point(aes(x = true_theta.3, y =fit_theta.3), color = blue, alpha = .5) +
    geom_abline(intercept = 0, slope = 1) 

  

#####################################
# compare with packages 
#####################################
lda_x <- sapply(1:para$M, function(s) 1*(1:para$D %in% para$w[[s]]$Ds_id)) %>% t
k <- 10
control_LDA_VEM <-
  list(estimate.alpha = F, alpha = 1, estimate.beta = TRUE,
       verbose = 0, prefix = tempfile(), save = 0, keep = 0,
       seed = as.integer(Sys.time()), nstart = 1, best = TRUE,
       var = list(iter.max = 500, tol = 10^-6),
       em = list(iter.max = 1000, tol = 10^-4),
       initialize = "random")
model_lda <- LDA(lda_x, k, method = "VEM", control = control_LDA_VEM)
lda_Result <- posterior(model_lda)
attributes(lda_Result)

library(reshape2)
library(ggplot2)

longData<-melt(t(lda_Result$terms))
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
# when estimating the alpha: package doesn't perform well. 
topics <- data.frame(ds_id = para$list_above500occu$diag_icd10, topics = t(lda_Result$terms))
df_topics <- list()
for(i in 1:para$K){
  df_topics[[i]]  <- topics %>% filter(get(paste0("topics.",i)) > 5/para$D) %>%
    pull(ds_id)
}
#################################################################
# testing for consistency over different hyper-parameter settings
#################################################################
# simulate dataset 
source("topic_functions.R")
para <- simulate_basic_topic_data(sample_sz = 20000, topic_number=8, disease_number=50, overlap = 2)
# plot true beta
longData<-melt(para$save_simulating_beta)
longData<-longData[longData$value!=0,]
ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

lda_x <- sapply(1:para$M, function(s) 1*(1:para$D %in% para$w[[s]]$Ds_id)) %>% t
k <- 8
topic_prior <- 1
control_LDA_VEM <-
  list(estimate.alpha = F, alpha = topic_prior, estimate.beta = TRUE,
       verbose = 0, prefix = tempfile(), save = 0, keep = 0,
       seed = as.integer(Sys.time()), nstart = 1, best = TRUE,
       var = list(iter.max = 500, tol = 10^-6),
       em = list(iter.max = 1000, tol = 10^-4),
       initialize = "random")
model_VEM <- LDA(lda_x, k, method = "VEM", control = control_LDA_VEM)
# plot true topics
# plot an example topics
lda_Result <- posterior(model_VEM)
attributes(lda_Result)
longData<-melt(t(lda_Result$terms))
longData<-longData[longData$value!=0,]
ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

# plot the true topics for comparison

#####################
# plot the lower-bound with respect to K number
#####################

##############
# try Gibbs: it simply works much better!
##############
control_LDA_Gibbs <-
  list( alpha = 1, estimate.beta = TRUE,
       verbose = 0, prefix = tempfile(), save = 0, keep = 0,
       seed = as.integer(Sys.time()), nstart = 1, best = TRUE,
       delta = 1,
       initialize = "random")

model_Gibbs <- LDA(lda_x, k, method = "Gibbs", control = control_LDA_Gibbs)
lda_Result <- posterior(model_Gibbs)
attributes(lda_Result)
longData<-melt(t(lda_Result$terms))
longData<-longData[longData$value!=0,]
ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

# first step: try put a prior for the VB


#############################
# working with the collapse VB
#############################
source("topic_functions.R")
para <- simulate_basic_topic_data(sample_sz = 20000, topic_number=8, disease_number=50, overlap = 2)
# plot true beta
longData<-melt(para$save_simulating_beta)
longData<-longData[longData$value!=0,]
ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

para$max_itr <- 500
para$lb <- data.frame("Iteration" = 0,"Lower_bound" = comp_lda_lb(para))
para$alpha <- rep(1, para$K)
para$tol <- 10^(-6) # tolerance of lower bound step
para$itr_check_lb <- 1
cvb_tag <- T
for(itr in 1:para$max_itr){
  print(paste0("Interation: ",itr))
  if(cvb_tag){
    para <- CVB0_E_zn(para)
  }else{
    para <- comp_E_zn(para)
    para <- comp_E_lntheta(para)
  }
  para <- update_beta_basic_lda(para)
  if(cvb_tag){
    para$lb[nrow(para$lb) + 1,] <- c(itr, CVB_lb(para))
  }else{
    para$lb[nrow(para$lb) + 1,] <- c(itr, comp_lda_lb(para))
  }
  if(itr %% para$itr_check_lb  ==0){
    print(paste0("Current Lower bound ", pull(filter(para$lb, Iteration == itr), Lower_bound), " at iteration: ",itr))
    if( ( abs(pull(filter(para$lb, Iteration == itr), Lower_bound) -
          pull(filter(para$lb, Iteration == (itr -para$itr_check_lb )), Lower_bound) ))/
        abs(pull(filter(para$lb, Iteration == (itr -para$itr_check_lb )), Lower_bound)) < para$tol ){
      print(paste0("Optimization converged at step ", itr))
      break
    }
  }
}
if(cvb_tag){
  para_cvb <- para
}else{
  para_mean_field <- para
}
longData<-melt(para$beta)
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") + 
  labs(x="letters", y="LETTERS", title="Matrix") + 
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))



###############################
# also compare with softImpute
###############################

ordered_lst <- list()
for(rep in 1:num_rep){
  model_lda <- model_lst[[rep]]
  lda_Result <- posterior(model_lda)
  ordered_lst[[rep]] <- lda_Result$terms
}
beta_for_clustering <- do.call(rbind, ordered_lst)
n <- nrow(beta_for_clustering) 
cmb <- expand.grid(i=1:n, j=1:n) 
cos.sim <- function(idx, X){
  A = X[idx[1],]
  B = X[idx[2],]
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}
C <- matrix(apply(cmb,1,function(x) cos.sim(x, beta_for_clustering)),n,n)
library(reshape2)
library(ggplot2)

longData<-melt(C)
longData<-longData[longData$value!=0,]
plt <- ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="topics", y="topics", title="cosine similarity") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))


# find the best way to approximate lnGamma
x0 <- c(seq(0.1,2,0.01))
lgma <- lgamma(x0)
approximate_lgama <- (x0-.5)*log(x0) - x0 + .5*log(2*pi) + 1/(12*x0) 
df_lgma_approx <- data_frame(x0, lgma, approximate_lgama)
ggplot(data = df_lgma_approx) +
  geom_line(aes(x=x0, y=lgma)) + 
  geom_line(aes(x=x0, y=approximate_lgama), linetype = "dashed")
####################################
# code speed timing for optimization
####################################
# time comsuming step 1: update E_zn
ptm <- proc.time()
for(i in 1:10){
  para$E_zn <-sapply(1:para$M, 
                     function(s)
                       (para$beta_w[[s]] %*% diag(exp(para$E_lntheta[s,]))  )/
                       (para$beta_w[[s]]) %*% t(exp(para$E_lntheta[rep(s,para$K),])),
                     simplify = FALSE)
}
print(proc.time() - ptm) # 16s

ptm <- proc.time()
for(i in 1:10){
  para$E_zn <- sapply(1:para$M, 
                      function(s) sapply(1:dim(para$w[[s]])[1] ,
                                         function(n) (exp(para$E_lntheta[s,]) * para$beta_w[[s]][n,])/
                                           sum(exp(para$E_lntheta[s,]) * para$beta_w[[s]][n,])), 
                      simplify = FALSE)
}
print(proc.time() - ptm) # 50s

# time comsuming step 2: update beta_w
ptm <- proc.time()
for(i in 1:10){
  para$beta_w <- lapply(para$w, function(w) as.matrix(para$beta[w$Ds_id,]) )
}
print(proc.time() - ptm) # as.matrix is pretty slow

ptm <- proc.time()
for(i in 1:10){
  para$beta_w <- lapply(para$w, function(w) para$beta[w$Ds_id,,drop=FALSE] )
}
print(proc.time() - ptm) # 43s 

ptm <- proc.time()
for(i in 1:10){
  para$beta_w <- sapply(1:para$M, function(s) para$beta[para$w[[s]]$Ds_id,,drop=FALSE] )
}
print(proc.time() - ptm) # 45s similar time

# compute speed between mapply and sapply
ptm <- proc.time()
for(i in 1:10){
  term3 <- Reduce('+', mapply(function(x,y) sum(x*y), para$E_zn, para$beta_w)) 
}
print(proc.time() - ptm) # 14s

ptm <- proc.time()
for(i in 1:10){
  term3 <- sum(sapply(1:para$M, function(s) sum(para$E_zn[[s]] * para$beta_w[[s]]) ) )
}
print(proc.time() - ptm) # 9s use sapply to avoid two big matrix here

# update E ln(theta) 
ptm <- proc.time()
for(i in 1:5){
  para$E_lntheta <- sapply(1:para$M,
                           function(s) digamma(para$alpha + colSums(para$E_zn[[s]])) - 
                             digamma(sum(para$alpha + colSums(para$E_zn[[s]]) ))) %>% t 
}
print(proc.time() - ptm) # 60s

ptm <- proc.time()
for(i in 1:5){
  para$alpha_z <- sapply(1:para$M, function(s) para$alpha + colSums(para$E_zn[[s]])) %>% t
  para$E_lntheta <- sapply(1:para$M,
                           function(s) digamma(para$alpha_z[s,]) - 
                             digamma(sum(para$alpha_z[s,]))) %>% t 
}
print(proc.time() - ptm) # 38s

# update beta
ptm <- proc.time()
for(i in 1:5){
  para$beta <- sapply(1:para$D, function(j) colSums(para$unlist_zn[para$unlist_Ds_id$Ds_id == j,]) ) %>% t 
}
print(proc.time() - ptm)  # 20s

ptm <- proc.time()
for(i in 1:5){
  beta_list <- para$unlist_Ds_id %>%
    mutate(id = row_number()) %>%
    group_by(Ds_id) %>%
    group_split()
  para$beta <- sapply(1:para$D, function(j) colSums(para$unlist_zn[beta_list[[j]]$id,]) ) %>% t 
}
print(proc.time() - ptm) # 2s: use dplyr whener possible!

ptm <- proc.time()
for(i in 1:5){
  term2 <- mclapply(1:para$M,
                           function(s) sum( para$E_zn[[s]] %*% para$E_lntheta[s,] ), mc.cores = 14) 
}
print(proc.time() - ptm) #  12s mclapply is not much quicker

ptm <- proc.time()
for(i in 1:5){
  term2 <- sapply(1:para$M,
                  function(s) sum( para$E_zn[[s]] %*% para$E_lntheta[s,] )) 
}
print(proc.time() - ptm) # 7.5s 

# testing other methods to avoid numeric NaN
ptm <- proc.time()
for(i in 1:10){
  term3 <- sum(para$unlist_zn * log(para$beta_w_full) ) 
}
print(proc.time() - ptm) # 8.9s

ptm <- proc.time()
for(i in 1:10){
  term3 <- sum( log(para$beta_w_full^para$unlist_zn) ) 
}
print(proc.time() - ptm) #  



