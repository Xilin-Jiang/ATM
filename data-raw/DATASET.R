## code to prepare `UKB_HES_10topics` dataset goes here
# the data is the one that has highest ELBO among random initializations
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
usethis::use_data(UKB_HES_10topics, overwrite = TRUE)


# code to prepare `HES_age_example` data set
rec_data <- read.csv("~/Desktop/comorbidity/Multi-morbidity_biobank/rec2subjectAbove1000occur_include_death_PheCode.csv")
sample_eid <- rec_data %>%
  group_by(eid) %>%
  summarise(per_patient_diag = n()) %>%
  filter(per_patient_diag > 10)
HES_age_example <- rec_data %>%
  filter(eid %in% sample_eid$eid) %>%
  group_by(eid) %>%
  sample_frac(0.5) %>%
  ungroup()
usethis::use_data(HES_age_example, overwrite = TRUE)
