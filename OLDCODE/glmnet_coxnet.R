library(glmnet)
library(dplyr)
library(ggplot2)
require(survival)
library(stringr)
library(tidyverse)

survive_age <- read.csv("~/Desktop/genetics_longitudinal_data/longitudinal_data/survive_age.csv")
keep <- read.table(file="~/Desktop/genetics_longitudinal_data/longitudinal_data/keep.txt", header=FALSE, sep=" ") 
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


# find all of the ICD-10 that is above 50 occurence (> 0.1%)
list_above500occu <- new_data %>%
  select(eid, diag_icd10, age_diag) %>%
  filter(str_detect(diag_icd10, "^[A-Q]")) %>%
  group_by(eid, diag_icd10) %>% 
  slice(1) %>%
  group_by(diag_icd10) %>%
  summarise(occ = n()) %>% 
  filter(occ >= 500)

# only keep the first incidence of disease 
first_icidence_age <- new_data %>%
  select(eid, diag_icd10, age_diag) %>%
  filter(diag_icd10 %in% list_above50occu$diag_icd10) %>%
  group_by(eid, diag_icd10) %>% 
  slice(1) %>% # just pick the first record and usually it is ranked with age
  # summarise(age_diag = min(age_diag)) %>%
  ungroup() 

##################################
# save the data for future use
##################################
write.csv(first_icidence_age, "DiseaseAbove500occur.csv", row.names = F)
first_icidence_age <- read.csv(paste0("DiseaseAbove500occur.csv"))


##################################
# create the surv data 
##################################
coxnet_data <- spread(first_icidence_age, key = diag_icd10, value = age_diag, fill = 0) %>%
  semi_join(keep, by=c("eid" = "V1"))
  

icd_10 <- "I251"

predictor <- coxnet_data %>% 
  select( -UQ(icd_10), -eid) %>%
  as.matrix()
predictor[which(predictor > 0)] <- 1

survive_data <-  coxnet_data %>% 
  left_join(survive_age) %>%
  mutate(time = ifelse(I251 > 0, I251, survive_year)) %>%
  mutate(status = ifelse(I251 > 0, 1, 0)) 

fit <- glmnet(predictor, Surv(survive_data$time,survive_data$status), family = "cox", maxit = 1000)
plot(fit)

Coefficients <- coef(fit, s = 0.05)
Active.Index <- which(Coefficients != 0) 
Active.Coefficients <- Coefficients[Active.Index]



