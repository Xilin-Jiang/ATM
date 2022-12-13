library(Hmisc)
#####################################
# summary of individual information
#####################################
a <- read.csv(file="~/Desktop/comorbidity/UKBB_interim_files/ukb4775.csv", header=TRUE, sep=",")
# coding information available at file:///Users/xilin/Desktop/genetics_longitudinal_data/longitudinal_data/ukb4775.html 

# delete withdraw group
with_draw <- read.csv(file="~/Desktop/genetics_longitudinal_data/longitudinal_data/with_draw_2020Feb11.csv", 
                      header=FALSE) 
keep <- read.table(file="~/Desktop/genetics_longitudinal_data/longitudinal_data/keep.txt", 
                   header=FALSE, sep=" ") %>% 
  anti_join(with_draw, by="V1")


########################################
# summary of age, sex, BMI and death
########################################
a <- read.csv(file="ukb4775.csv", header=TRUE, sep=",")
birthyear <- a %>% 
  select(eid, X34.0.0, X31.0.0,X52.0.0, X23104.0.0, X40000.0.0) 
write.csv(birthyear, file = "Year_of_birth.csv")

# create a survive time 
# note here the birth year contain X52.0.0
survive_age <- birthyear %>% # 2018-02-14 is the end of the censoring
  mutate(eid = as.character(eid), X40000.0.0 = as.character(X40000.0.0)) %>%
  semi_join(keep["V1"], by = c("eid" = "V1")) %>%
  mutate(birth_year=as.Date(paste(X34.0.0, X52.0.0, "01", sep="-"))) %>%
  mutate(survive_date=as.Date(replace(X40000.0.0, X40000.0.0 == "", "2018-02-14"))) %>%
  mutate(survive_year = difftime(survive_date ,birth_year,units = "days")/365) %>%
  select(eid, survive_year) %>%
  filter(!is.na(survive_year))

write.csv(survive_age, paste("survive_age.csv", sep = ""), row.names = FALSE)

# compute death as a disease
death_age <- birthyear %>% # 2018-02-14 is the end of the censoring
  mutate(eid = as.character(eid), X40000.0.0 = as.character(X40000.0.0)) %>%
  filter(X40000.0.0 != "") %>%
  mutate(birth_year=as.Date(paste(X34.0.0, X52.0.0, "01", sep="-"))) %>%
  mutate(death_date=as.Date(X40000.0.0)) %>%
  mutate(age_diag = difftime(death_date ,birth_year,units = "days")/365) %>%
  mutate(diag_icd10 = "Death") %>%
  select(eid,diag_icd10, age_diag) %>%
  filter(!is.na(age_diag))
write.csv(death_age, paste("death_age.csv", sep = ""), row.names = FALSE)
########################################################
# HES data: combining primary and secondary
########################################################

# primary diagnosis: the reasons that the patient came to the hospital, contains information of the visit and record id is unique
hes <- read.table(file = '~/Desktop/genetics_longitudinal_data/longitudinal_data/hesin_130418.tsv',quote = '', sep = '\t', header = TRUE, fill=TRUE)

# secondary diagnosis, use the "record_id" column to link back to the first file for the visit information
hes_diag <- read.table(file = '~/Desktop/genetics_longitudinal_data/longitudinal_data/hesin_diag10_130418.tsv',quote = '', sep = '\t', header = TRUE, fill=TRUE) %>%
  mutate(record_id = as.factor(record_id))


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

########################################################
# Genetic principle component
########################################################
# get the principle genetics components
fam <- read.table(file="/Users/xilin/PycharmProjects/longitudinal_EHR/ukb1062_cal_chr1_v2_s488366.fam", header=FALSE, sep=" ")
principle_com <- read.table(file="~/Desktop/genetics_longitudinal_data/longitudinal_data/ukb_sqc_v2.txt", header=TRUE, sep=" ") 
data <- cbind(fam, principle_com)  
data %>%  
  select(V1, V2, V5,  contains("PC",ignore.case=FALSE)) %>%
  write.table("covariates_for_asso_over_time.txt", sep="\t", col.names = FALSE, row.names = FALSE)

########################################################
# Genetic data for 3000 LD-independent variants
########################################################

# extraction of useful snp and encoding alleles as 0/1/2 is performed on the cluster to adhere the data privacy requirements.
# /apps/well/plink/1.90b2n/plink --bed /well/mcvean/xilin/UK_biobank/all_chromosome.bed --bim 
# /well/mcvean/xilin/UK_biobank/all_chromosome.bim --fam /well/mcvean/xilin/UK_biobank/all_chromosome.fam  
# --extract ${1}/${1}_snp_list.txt -keep /users/mcvean/xilin/xilin/UK_biobank/keep.txt --recode A  --out ${1}/${1}

# LD filtering by hand
library(Hmisc)
library(reshape2)
df <- load("~/Desktop/comorbidity/UKBB_interim_files/snp_data_291019.rdata")

pt <- "*.raw"
temp <- list.files(paste("Survival/", sep=""), pattern=pt)

ds_list <- lapply(temp, function(x) strsplit(x, "[.]")[[1]][1])
for(id in 1:length(ds_list)){
  ds_id <- ds_list[[id]]
  print(ds_id)
  genetics <- read.table(paste0("Survival/",ds_id,".raw"), header=T, check.names=FALSE) %>%
    select(-IID, -PAT, -MAT, -SEX, -PHENOTYPE)
  res <- rcorr(as.matrix(select(genetics, - FID)), type = "pearson")
  longData <- melt(res[[1]])
  longData<-longData[longData$value!=0,]
  longData<-longData[longData$value!=1,]
  cor_data <- longData %>% filter(abs(value) > 0.2)
  snp_lst <- names(genetics)[2:length(genetics)]
  rm_snp <- c()
  for(snp in snp_lst){
    cor <- cor_data %>%
      filter(Var1 == snp| Var2 == snp)
    if(dim(cor)[1] > 0){
      rm_snp <- c(rm_snp,snp)
      cor_data <- cor_data %>%
        filter(!Var1 == snp & !Var2 == snp)
    }
  }
  rm_snp <- sapply(rm_snp, function(x) strsplit(x, "_")[[1]][1])
  pos <- which(t$coding == ds_id)
  t[pos,3] <- t[pos,3] - length(rm_snp)
  tables[[pos]] <- tables[[pos]] %>%
    filter(! SNP %in% rm_snp)
}
save(t,tables,file = "SNP_data_LD_rm_20200521.rdata")
load("SNP_data_LD_rm_20200521.rdata")

