library(dplyr)
library(softImpute)
library(ggplot2)
# import the UK_bb_dataset
covariates <- read.csv(file="~/Desktop/genetics_longitudinal_data/longitudinal_data/ukb4775.csv", header=TRUE, sep=",")
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

# apply the soft impute
matrix_data <- new_data %>%
  select(eid, diag_icd10) %>%
  group_by(eid, diag_icd10) %>% 
  slice(1) %>%
  ungroup() %>%
  mutate(eid = as.factor(eid))
# still keep the matrix data for future reference
numeric_data <- matrix_data %>%
  mutate(eid = as.numeric(eid), diag_icd10 = as.numeric(diag_icd10))
idx_patient <- numeric_data[[1]]
idx_disease <- numeric_data[[2]]
# also make a lookup table 
disease_label_map <- matrix_data %>%
  select(diag_icd10) %>%
  mutate( diag_icd10_num = as.numeric(diag_icd10)) %>%
  distinct()
# making an incomplete-class 
sparse_matrtix_patient_disease <- Incomplete(idx_patient, idx_disease, rep(1, length(idx_disease)))

# first get the uppper limit of constraint
l0 <- lambda0(sparse_matrtix_patient_disease)
fit1 <- softImpute(sparse_matrtix_patient_disease, rank.max = 10, lambda = 40, maxit = 200)

# recover the loadings and disease id
id_disease <- list()
loading_disease <- list()
plt <- lst()
for(colidx in 1:(length(fit1$d)-1)){
  id_disease[[colidx]] <- which(abs(fit1$v[,colidx]) > 0.1)
  loading_disease[[colidx]] <- fit1$v[which(abs(fit1$v[,colidx]) > 0.1),colidx]
}
for(colidx in 1:(length(fit1$d)-1)){
  df_plt <- data_frame(ds = id_disease[[colidx]], loading = loading_disease[[colidx]]) %>%
    left_join(disease_label_map, c("ds" = "diag_icd10_num")) %>%
    mutate(loading_dim = colidx)
  plt[[colidx]] <- ggplot(df_plt, aes(loading_dim, diag_icd10,fill = loading)) +
    geom_tile()
  ggsave(paste0("softImpute_component", colidx, ".png"),plt[[colidx]], width = 2, height = 5)
}

# run a lasso path 
fit_list <- list()
for(lbd in seq(from = 100, to = 20, by = -10)){
  print(paste0("lambda value: ",lbd))
  fit_list[[lbd]] <- softImpute(sparse_matrtix_patient_disease, rank.max = 16, lambda = lbd, maxit = 200)
}
save(fit_list, file = "fitting_softImpute_path.RData")

# plot the path
df_plt <- setNames(data.frame(matrix(0,ncol = 16, nrow = 10)), 
                   sapply(1:16, function(x) paste0("component", x)))
df_plt$lambda <- c(l0, seq(from = 100, to = 20, by = -10))

for(k in 1:16){
  df_plt[[k]] <- 0
}
for(i in 1:9){
  lbd <- seq(from = 100, to = 20, by = -10)[i]
  d <- fit_list[[lbd]]$d
  df_plt[(i+1), 1:length(d)] <- d
}
ggplot(df_plt) + 
  geom_line(aes(x = lambda, y = component1)) +
  geom_line(aes(x = lambda, y = component2))+
  geom_line(aes(x = lambda, y = component3))+
  geom_line(aes(x = lambda, y = component4))+
  geom_line(aes(x = lambda, y = component5))+
  geom_line(aes(x = lambda, y = component6)) +
  geom_line(aes(x = lambda, y = component7)) 
  

path_plt <- ggplot(df_plt)
for(k in 2:16){
  path_plt <- path_plt + geom_line(aes(x = lambda, y = paste0("component", k)))
}
path_plt <- path_plt + scale_y_continuous(trans = 'log2')

# analysis on females only 
keep_woman <- read.table(file="keep_women.txt", header=FALSE, sep=" ") %>%
  mutate(eid = as.character(V1))
woman_matrix_data <- matrix_data %>% 
  semi_join(keep_woman, by = "eid" )

woman_numeric_data <- woman_matrix_data %>%
  mutate(eid = as.numeric(eid), diag_icd10 = as.numeric(diag_icd10))
idx_woman <- woman_numeric_data[[1]]
idx_woman_disease <- woman_numeric_data[[2]]

woman_sparse_matrtix_patient_disease <- Incomplete(idx_woman, idx_woman_disease, rep(1, length(idx_woman_disease)))

# do the analysis on females
woman_l0 <- lambda0(woman_sparse_matrtix_patient_disease)
fit1_female <- softImpute(woman_sparse_matrtix_patient_disease, rank.max = 10, lambda = 20, maxit = 200)

id_disease <- list()
loading_disease <- list()
plt <- lst()
for(colidx in 1:(length(fit1_female$d)-1)){
  id_disease[[colidx]] <- which(abs(fit1_female$v[,colidx]) > 0.1)
  loading_disease[[colidx]] <- fit1_female$v[which(abs(fit1_female$v[,colidx]) > 0.1),colidx]
}
for(colidx in 1:(length(fit1_female$d)-1)){
  df_plt <- data_frame(ds = id_disease[[colidx]], loading = loading_disease[[colidx]]) %>%
    left_join(disease_label_map, c("ds" = "diag_icd10_num")) %>%
    mutate(loading_dim = colidx)
  plt[[colidx]] <- ggplot(df_plt, aes(loading_dim, diag_icd10,fill = loading)) +
    geom_tile()
  ggsave(paste0("female_softImpute_component", colidx, ".png"),plt[[colidx]], width = 2, height = 5)
}

#############################
# plot the age distribution 
# for difference disease
#############################


