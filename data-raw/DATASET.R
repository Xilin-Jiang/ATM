## code to import some common functions from dplyr and stats
usethis::use_pipe(export = TRUE)
usethis::use_import_from("dplyr", c("anti_join", "arrange", "group_by",
                                    "bind_rows", "filter", "left_join",
                                    "mutate", "n", "pull", "rename",
                                    "select", "slice", "summarise", "ungroup",
                                    "data_frame"))
usethis::use_import_from("grDevices", "colorRampPalette")
usethis::use_import_from("stats", c("binomial", "coef", "glm", "lm", "optim",
           "p.adjust", "pnorm", "quantile", "rexp", "rgamma", "rnorm",
           "runif", "sd"))
usethis::use_import_from("utils", "object.size")
usethis::use_import_from("ggplot2", c("aes", "theme_bw", "labs", "geom_line",
                                      "aes_string", "ggplot"))


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
load(paste0(DIR, "Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData"))
# UKB_HES_10topics <- model_output[[1]]
UKB_HES_10topics <- para$pi_beta_basis
usethis::use_data(UKB_HES_10topics, overwrite = TRUE)

# code to prepare `UKB_349_disease`
UKB_349_disease <- read.csv("~/Desktop/comorbidity/Multi-morbidity_biobank/listAbove1000include_deaths_PheCode.csv")
usethis::use_data(UKB_349_disease, overwrite = TRUE)

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

# code to create `phecode_icd10cm`, `phecode_icd10`, and `disease_info_phecode_icd10`; mapping between ICD2Phecode
# mapping all ICD-10 to phecodes; need to first deal with the non-mapping (one icd10 to different phecode)
DIR <- "~/Desktop/comorbidity/Multi-morbidity_biobank/"
non_one2one_map <- read.csv(paste0(DIR, "Phecode_map_v1_2_icd10cm_beta.csv")) %>%
  group_by(phecode) %>%
  summarise(occ = n())

phecode_icd10cm <- read.csv(paste0(DIR, "Phecode_map_v1_2_icd10cm_beta.csv")) %>%
  left_join(non_one2one_map, by = "phecode") %>%
  group_by(icd10cm) %>%
  arrange( desc(occ), .by_group = T ) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(ICD10 = sub("[.]","",icd10cm)) %>% # remove the dots
  select(ICD10, phecode,exclude_range, exclude_name)
usethis::use_data(phecode_icd10cm, overwrite = TRUE)

# only use the first 4 letters of the ICD10
short_icd10cm <- phecode_icd10cm %>%
  mutate(ICD10 = substring(ICD10, 1,4)) %>%
  left_join(non_one2one_map, by = "phecode") %>%
  group_by(ICD10) %>%
  arrange( desc(occ), .by_group = T, ) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  rename(parent_phecode = phecode)
usethis::use_data(short_icd10cm, overwrite = TRUE)

phecode_icd10 <- read.csv(paste0(DIR,"phecode_icd10.csv")) %>%
  left_join(non_one2one_map, by = c("PheCode" = "phecode")) %>%
  group_by(ICD10) %>%
  arrange( desc(occ), .by_group = T, ) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(ICD10 = sub("[.]","",ICD10)) %>% # remove the dots
  select(ICD10, PheCode, Excl..Phecodes, Excl..Phenotypes)
usethis::use_data(phecode_icd10, overwrite = TRUE)

# save a phecode classifications for future plotting
disease_info_phecode_icd10 <- read.csv(paste0(DIR, "Phecode_map_v1_2_icd10cm_beta.csv")) %>%
  group_by(phecode) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(ICD10 = sub("[.]","",icd10cm)) %>% # remove the dots
  select(ICD10, phecode, phecode_str, exclude_range, exclude_name) %>%
  rename(phenotype = phecode_str)
usethis::use_data(disease_info_phecode_icd10, overwrite = TRUE)

# save an icd10 rec_data as an example
# code to prepare `HES_icd10_example` data set
rec_data <- read.csv("~/Desktop/comorbidity/Multi-morbidity_biobank/rec2subjectAbove500occur_ICDA2N.csv")
sample_eid <- rec_data %>%
  group_by(eid) %>%
  summarise(per_patient_diag = n()) %>%
  filter(per_patient_diag > 10)
HES_icd10_example <- rec_data %>%
  filter(eid %in% sample_eid$eid) %>%
  group_by(eid) %>%
  sample_frac(0.5) %>%
  ungroup()
usethis::use_data(HES_icd10_example, overwrite = TRUE)


