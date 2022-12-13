###############################
# analysis of subtypes: for figure 4 of the paper
###############################
source("topic_functions.R")
source("plotting_functions.R")

##################################
# supplementary figures and tables about subtype property
# h2g v.s. subtype size, age
# subtype distribution (size)
##################################
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
loadings_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[para$ds_list[[j]]$id,])) %>% t
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds.system <- para$list_above500occu %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") ) 
# combine the last four systems to other:congenital anomalies, symptoms, injuries & poisonings, <NA>  
ds.system <- ds.system %>%
  mutate(exclude_name = if_else(exclude_name %in% c("congenital anomalies", "symptoms", "injuries & poisonings","", NA), "others", exclude_name))
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

subtypes_list <- read.table( "subtype_disease_topic_list.txt", sep="\t", header = F)
subtype_diseases_list <- subtypes_list %>%
  group_by(V1) %>%
  slice(1)
all_ldsc_seg <- read.csv("h2g_imputed.csv", header =  T)
library(colorBlindness)
displayAvailablePalette(color="white")
subtype_info <- list()
for(idx in 1:dim(subtype_diseases_list)[1]){
  ds_id <- subtype_diseases_list$V1[idx]
  j <- match(as.numeric(ds_id), para$list_above500occu$diag_icd10)
  cases_eid <- list()
  for(topic_id in 1:para$K){
    topic_specific_id <- which(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max) == topic_id)
    topic_specific_set <- para$unlist_Ds_id[para$ds_list[[j]]$id,][topic_specific_id,]
    
    cases_eid[[topic_id]] <- topic_specific_set %>%
      mutate(topic = topic_id) %>%
      mutate(Ds_id = para$list_above500occu[Ds_id, 1])
  }
  cases_eid <- bind_rows(cases_eid)
  
  subtype_info[[idx]] <- cases_eid %>% 
    group_by(topic) %>%
    summarise(age = mean(age_diag), subtype_count = n()) %>%
    mutate(subtype_percentage = subtype_count/sum(subtype_count), disease = ds_id)

}
subtype_info <- bind_rows(subtype_info)
subtype_count <- subtype_info %>% 
  mutate(topic = factor(topic, levels = order_ds)) %>%
  left_join(select(ds.system, diag_icd10, phenotype), by=c("disease" = "diag_icd10")) %>%
  arrange(topic) %>%
  pivot_wider(id_cols = c("disease", "phenotype"), names_from = "topic", values_from = "subtype_count") %>%
  replace(is.na(.), 0) %>%
  arrange(disease)
topic_ordered_namte <- c("NRI", "CER","SRD", "CVD", "UGI", "LGI", "FGND",  "MGND", "MDS", "ARP")
names(subtype_count) <- c("Phecode", "Description",topic_ordered_namte)
subtype_count %>% 
  write.csv("~/Desktop/comorbidity/paper_writing/supplementary_files/Subtype_counts.csv", row.names = F)

# also same a table of average age for each subtype. 
subtype_age <- subtype_info %>% 
  mutate(topic = factor(topic, levels = order_ds), age = round(age, 2)) %>%
  left_join(select(ds.system, diag_icd10, phenotype), by=c("disease" = "diag_icd10")) %>%
  arrange(topic) %>%
  pivot_wider(id_cols = c("disease", "phenotype"), names_from = "topic", values_from = "age") %>%
  arrange(disease)
topic_ordered_namte <- c("NRI", "CER","SRD", "CVD", "UGI", "LGI", "FGND",  "MGND", "MDS", "ARP")
names(subtype_age) <- c("Phecode", "Description",topic_ordered_namte)
subtype_age %>% 
  write.csv("~/Desktop/comorbidity/paper_writing/supplementary_files/Subtype_age.csv", row.names = F)

# save a distribution of topoic loading figure
subtype_pct <- subtype_info %>% 
  mutate(topic = factor(topic, levels = order_ds)) %>%
  left_join(select(ds.system, diag_icd10, phenotype), by=c("disease" = "diag_icd10")) %>%
  arrange(topic) %>%
  pivot_wider(id_cols = c("disease", "phenotype"), names_from = "topic", values_from = "subtype_percentage") %>%
  replace(is.na(.), 0) %>%
  select(-disease, - phenotype) %>%
  as.matrix()

top_value <- sapply(1:dim(subtype_pct)[1], function(x) sort(subtype_pct[x,], decreasing = T)) %>% t
colnames(top_value) <- 1:10
df_box <- melt(top_value) %>%
  rename(topic_order = Var2, weights = value, record_number = Var1) %>%
  mutate(topic_order = factor(topic_order, levels = 1:10))
plt <- ggplot(data=df_box) +
  geom_boxplot(mapping=aes(x=topic_order, y=weights), fill = red, width=0.8, alpha = 0.5, outlier.shape = NA) +
  scale_y_continuous(labels = scales::percent) + 
  labs(x = "Topic value rank", y = "Subtype proportion") + 
  theme_bw(base_size = 20) +
  theme(legend.position = "None",panel.background=element_blank()) + 
  geom_hline(aes(yintercept = 1), color = "red", linetype = "dashed")
ggsave("~/Desktop/comorbidity/paper_writing/supplementary_files/Sfig_subtype_topic_weight_dist.png", plt, width = 5, height = 5)

# plot h2g (subtype v.s. all case) with respect to sample size and age-at-onset 
sub_h2g <- all_ldsc_seg %>%
  filter(topic != "all")
all_h2g <- all_ldsc_seg %>%
  filter(topic == "all") %>%
  rename(all_age = age, all_h2g = h2g, all_seh2g=seh2g, all_N = N) %>%
  select( - topic, - z_score)
plot_h2g <- sub_h2g %>%
  left_join(all_h2g, by = "disease") %>%
  mutate(age_diff = age - all_age, h2g_diff = h2g - all_h2g, se_diff = sqrt(seh2g^2 + all_seh2g^2))

plt <- ggplot(plot_h2g) +
  geom_pointrange(aes(x = age_diff, y = h2g_diff, ymin =h2g_diff -1.96*se_diff, ymax = h2g_diff + 1.96*se_diff), color = red) + 
  labs(x = "Age difference (subtype vs. all cases)", y = "Heritability difference") + 
  theme_bw(base_size = 20) +
  theme(legend.position = "None",panel.background=element_blank()) 
ggsave("~/Desktop/comorbidity/paper_writing/supplementary_files/h2g_vs_age.png", plt, width = 6, height = 5)
plt <- ggplot(plot_h2g) +
  geom_pointrange(aes(x = N, y = h2g, ymin =h2g -1.96*seh2g, ymax = h2g + 1.96*seh2g), color = red) + 
  labs(x = "Sample size", y = "Heritability") + 
  theme_bw(base_size = 20) +
  theme(legend.position = "None",panel.background=element_blank()) 
ggsave("~/Desktop/comorbidity/paper_writing/supplementary_files/h2g_vs_sampe_size.png", plt, width = 6, height = 5)

###################################################### 
# step 1: filter the diseases / subtypes for each topics
###################################################### 
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
thre_pick <- 500
thre_confident <- .5
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
    topic_specific_id <- which(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max) == topic_id & apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, max) > thre_confident)
    if(length(topic_specific_id) > thre_pick){
      number_case <- c(number_case, length(topic_specific_id))
      mean_age <- c(mean_age, mean(para$ds_list[[j]][topic_specific_id, ]$age_diag))
      topic_assign <- c(topic_assign, topic_id)
      disease_idx <- c(disease_idx, para$list_above500occu$diag_icd10[j])
    }
  }
}
common_disease_within_topics <- data.frame(disease = disease_idx, topic = topic_assign, mean_age = mean_age)
subtypes_list <- common_disease_within_topics %>%
  group_by(disease) %>%
  summarise(topic_number = n()) %>%
  filter(topic_number > 1)
subtypes_data <- common_disease_within_topics %>% 
  filter(disease %in% subtypes_list$disease) %>%
  left_join(select(ds_list, diag_icd10, phenotype), by = c("disease" = "diag_icd10"))
  
write.table(subtypes_data, "subtype_disease_topic_list.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)
# write.table(subtypes_data$disease, "subtype_disease_list.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)


###########################################
# Figure 4: subtypes heatmap of the larger figure
###########################################

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

# load subtype list
subtypes_list <- read.table( "subtype_disease_topic_list.txt", sep="\t", header = F)

df_subplot <- longData %>%
  filter( diag_icd10 %in% subtypes_list$V1) 
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
df_segments <- data.frame(x = 0:10 * 2 + 0.5, xend =  0:10 * 2+ 0.5, y = rep(0.5,11), yend = rep(length(unique(df_subplot$Var1)) + 0.5,11)) 
plt <- plt + geom_segment(data=df_segments, aes(x,y,xend=xend, yend=yend), size=.5, alpha = 0.3, inherit.aes=F)
ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/subtype_overview.png"), plt, width = 5, height = 10)

categories <- ds.system$exclude_name
categories[346] <- categories[347] # fixing one error in Phecode coding
longData <- melt(as.matrix(categories)) %>%
  mutate(Var1 = factor(ds.system$diag_icd10[Var1], levels = rev(ds.system$diag_icd10))) %>%
  filter(Var1 %in% as.character(subtypes_list$V1) )
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
ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/pallete_subtypes.png") ,plt_pallete, width = 1, height = 10)
subtypes_list <- subtypes_list %>%
  mutate(V1 = as.character(V1)) %>%
  group_by(V1) %>%
  slice(1)
longData %>% left_join(subtypes_list, by = c("Var1" = "V1") ) %>%
  select(V4) %>%
  write.table(paste0("~/Desktop/comorbidity/Multi-morbidity_biobank/diseases_subtypes.txt"), col.names = FALSE, row.names = FALSE, quote = F)

longData %>% left_join(subtypes_list, by = c("Var1" = "V1") ) %>%
  select(Var1) %>%
  write.table(paste0("~/Desktop/comorbidity/Multi-morbidity_biobank/top_hererogeneous_disease.txt"), col.names = FALSE, row.names = FALSE, quote = F)


#######################################################
# plotting disease incidence topic assignment over age
#######################################################
subtypes_list <- read.table( "subtype_disease_topic_list.txt", sep="\t", header = F)
subtype_diseases_list <- subtypes_list %>%
  group_by(V1) %>%
  slice(1)

library(colorBlindness)
displayAvailablePalette(color="white")

for(idx in 1:dim(subtype_diseases_list)[1]){
  ds_id <- subtype_diseases_list$V1[idx]
  j <- match(as.numeric(ds_id), para$list_above500occu$diag_icd10)
  cases_eid <- list()
  for(topic_id in 1:para$K){
    topic_specific_id <- which(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max) == topic_id)
    topic_specific_set <- para$unlist_Ds_id[para$ds_list[[j]]$id,][topic_specific_id,]
    
    cases_eid[[topic_id]] <- topic_specific_set %>%
      mutate(topic = topic_id) %>%
      mutate(Ds_id = para$list_above500occu[Ds_id, 1])
  }
  cases_eid <- bind_rows(cases_eid)
  
  colors_topic <- c(cbPalette, paletteMartin[c(8,11)])
  names(colors_topic) <- as.character(1:10)
  # plot histogram for each group
  plt <- ggplot(cases_eid,aes(x = age_diag, fill = factor(topic, levels = as.character(order_ds)), 
                              colour = factor(topic, levels = as.character(order_ds)))) +
    geom_histogram(alpha = 0.5, position = "stack", binwidth = 1) + 
    scale_fill_manual(values = colors_topic) + 
    scale_color_manual(values =colors_topic) +
    theme_bw(base_size = 15) + 
    theme(legend.position = "None") + 
    labs(x="Age (years)", y="incidence count", title=paste0("Disease: ", subtype_diseases_list$V4[idx])) 
  ggsave(paste0("../figures/","subtypes_age_distribution_",ds_id,".png"), plt, width = 6, height = 4)
}

# plot the color legend
color_data <-  matrix(rev(order_ds), ncol = 1) %>%
  melt %>%
  mutate(value = as.character(value))
plt_palette <- ggplot() + 
  geom_tile(data = color_data, aes(x = Var2, y = Var1, fill=value), alpha = 0.5) + 
  scale_fill_manual(values = colors_topic) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/subtypes_topic_color.png"), plt_palette, width = 1, height = 8)


###################################################################
# comorbidity composition with subtypes (i.e. one Phecode could be considered as multiple diseases)
# --- this will be very useful when computing the odds ratio with ordering
###################################################################

####################################################################
# compute PRS score with subtype categories
####################################################################
DIR <- "~/Desktop/comorbidity/association_results/"
ds_list <- read.table("top_hererogeneous_disease.txt")
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
# playing with the code
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("V1" = "phecode"))

colors_topic <- c(cbPalette, paletteMartin[c(8,11)])
names(colors_topic) <- as.character(1:10)

regression_md <- list()
p_values <- matrix(nrow = length(ds_list$V1), ncol = 1)
data_percentile <- list()
for(idx in 1:length(ds_list$V1)){
  ds_id <- as.character(ds_list$V1[idx])
  # check the PRS data are indeed correct
  ds_eid <- rec_data %>% 
    filter(diag_icd10  == ds_id) %>%
    select(eid)
  PRS_profile <- read.table(paste0(DIR, ds_id, ".profile"), header =T)
  
  # get the disease of different groups
  j <- match(as.numeric(ds_id), para$list_above500occu$diag_icd10)
  cases_eid <- list()
  for(topic_id in 1:para$K){
    topic_specific_id <- which(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max) == topic_id)
    topic_specific_set <- para$unlist_Ds_id[para$ds_list[[j]]$id,][topic_specific_id,]
    if(dim(topic_specific_set)[1] > 500 & dim(topic_specific_set)[1]/dim(ds_eid)[1] > 0.05){
      cases_eid[[topic_id]] <- topic_specific_set %>%
        mutate(topic = topic_id) %>%
        mutate(Ds_id = para$list_above500occu[Ds_id, 1])
    }
  }
  cases_eid <- bind_rows(cases_eid)
  topic_list <- unique(cases_eid$topic)
  
  # plot percentile PRS with incidence rate
  percentile <- PRS_profile %>% mutate(percentiles = ntile(SCORE,100))
  prevalence <- rep(NA, 100)
  prevalence_df <- list()
  for(k in topic_list){
    eid_group <- cases_eid %>%
      filter(topic == k)
    prevalence_subgroup <- rep(NA, 100)
    for(pct in 1:100){
      pct_eid <- percentile %>%
        filter(percentiles == pct) %>%
        pull(1)
      # compute log(Relative risk)
      prevalence_subgroup[pct] <- 100 * sum(pct_eid %in% eid_group$eid )/length(eid_group$eid)
    }
    prevalence_df[[k]] <- data.frame(percentile = 1:100/100, topic = k, prevalence_subgroup = prevalence_subgroup)
  }

  prevalence_df <- bind_rows(prevalence_df)
  plt <- ggplot(prevalence_df) +
    geom_point(aes(x = percentile, y = prevalence_subgroup, color = as.character(topic))) +
    scale_color_manual(values =colors_topic) + 
    scale_x_continuous(labels = scales::percent) + 
    theme_bw(base_size = 15) + 
    theme(legend.position = "None") + 
    labs(x="PRS percentile", y="Relative risk", title=paste0("Disease: ", ds_list$phenotype[idx])) 
  ggsave(paste0("../figures/","subtypes_PRS_",ds_id,".png"), plt, width = 6, height = 4)
  
  # regression analysis
  regression_df <- cases_eid %>% 
    left_join(select(PRS_profile, FID, SCORE), by = c("eid" = "FID")) %>%
    mutate(topic = factor(topic))

  try({
    regression_md[[idx]] <- lm(SCORE ~ topic, data = regression_df)
    p_values[idx] <- anova(regression_md[[idx]])$`Pr(>F)`[1]
    print(paste(ds_list$phenotype[idx], " has a p-value ", p_values[idx], " and ", length(topic_list), " clusters"))
  })
}
ds_list$p_true <- -log10(p_values)
ds_list <- ds_list %>%
  arrange(p_true)
ds_list$p_sim = sort(- log( runif(length(p_values))) )
ggplot(data = ds_list) + 
  geom_point(aes(x = p_sim, y = p_true), color = grey, size = 2) + 
  geom_abline(slope = 1) + 
  geom_label_repel(aes(x = p_sim, y = p_true, label=ifelse(p_true> 3,as.character(phenotype),''))) 

ds_list$phenotype[which(ds_list$p_values > 2)] 

#########################################
# example PRS percentile 
#########################################
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
DIR <- "~/Desktop/comorbidity/association_results/"
ds_list <- read.table("top_hererogeneous_disease.txt")
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
# playing with the code
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("V1" = "phecode"))

colors_topic <- c(cbPalette, paletteMartin[c(8,11)])
names(colors_topic) <- as.character(1:10)

for(ds_id in c("401.1", "272.11", "250.2")){
  idx <- which(ds_list$V1 == ds_id)
  
  if(ds_id == "272.11"){
    topic_list <- c(5,7)
  }else{
    topic_list <- c(5,6,7)
  }
  
  # check the PRS data are indeed correct
  ds_eid <- rec_data %>% 
    filter(diag_icd10  == ds_id) %>%
    select(eid)
  PRS_profile <- read.table(paste0(DIR, ds_id, ".profile"), header =T)
  
  # get the disease of different groups
  j <- match(as.numeric(ds_id), para$list_above500occu$diag_icd10)
  cases_eid <- list()
  for(topic_id in 1:para$K){
    topic_specific_id <- which(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max) == topic_id)
    topic_specific_set <- para$unlist_Ds_id[para$ds_list[[j]]$id,][topic_specific_id,]
    if(dim(topic_specific_set)[1] > 500 & 
       dim(topic_specific_set)[1]/dim(ds_eid)[1] > 0.05 & 
       (topic_id %in% topic_list)){
      cases_eid[[topic_id]] <- topic_specific_set %>%
        mutate(topic = topic_id) %>%
        mutate(Ds_id = para$list_above500occu[Ds_id, 1])
    }
  }
  cases_eid <- bind_rows(cases_eid)
  topic_list <- unique(cases_eid$topic)
  
  ###################### test
  # rate <- c()
  # for(i in 1:100){
  #   fid <- percentile %>% 
  #     filter(percentiles == i) %>%
  #     pull(FID)
  #   rate[i] <- mean(fid %in% ds_eid$eid)
  # }
  ###################### test
  
  # plot percentile PRS with incidence rate
  tiles_num <- 50
  percentile <- PRS_profile %>% mutate(percentiles = ntile(SCORE,tiles_num))
  prevalence <- rep(NA, tiles_num)
  prevalence_df <- list()
  for(k in topic_list){
    eid_group <- cases_eid %>%
      filter(topic == k)
    prevalence_subgroup <- rep(NA, tiles_num)
    for(pct in 1:tiles_num){
      pct_eid <- percentile %>%
        filter(percentiles == pct) %>%
        pull(1)
      # compute log(Relative risk)
      prevalence_subgroup[pct] <- tiles_num * sum(pct_eid %in% eid_group$eid )/length(eid_group$eid)
    }
    prevalence_df[[k]] <- data.frame(percentile = 1:tiles_num/tiles_num, topic = k, prevalence_subgroup = prevalence_subgroup)
  }
  
  prevalence_df <- bind_rows(prevalence_df)
  plt <- ggplot(prevalence_df) +
    geom_point(aes(x = percentile, y = prevalence_subgroup, color = as.character(topic))) +
    scale_color_manual(values =colors_topic) + 
    scale_x_continuous(labels = scales::percent) + 
    theme_bw(base_size = 15) + 
    theme(legend.position = "None") + 
    labs(x="PRS percentile", y="Relative risk", title=paste0("Disease: ", ds_list$phenotype[idx])) 
  ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/","subtypesExample_PRS_",ds_id,".png"), plt, width = 6, height = 4)
  
}

######################################
# perform regression over the patient loadings
######################################
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
subtypes_list <- read.table( "subtype_disease_topic_list.txt", sep="\t", header = F)
subtype_diseases_list <- subtypes_list %>%
  group_by(V1) %>%
  slice(1)
# step 1: normalise both loadings and score -- otherwise the unit doesn't make sense

cc_dir <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/CC_gwas/"
patient_loadings <- sweep((para$alpha_z - 1), 1, rowSums(para$alpha_z -1), FUN="/") 
DIR <- "~/Desktop/comorbidity/association_results/"
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
ds_list <- read.table("top_hererogeneous_disease.txt")
# playing with the code
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("V1" = "phecode"))
coefficients <- matrix(nrow = length(ds_list$V1), ncol = 10)
p_values <- matrix(nrow = length(ds_list$V1), ncol = 10)
for(idx in 1:length(ds_list$V1)){
  ds_id <- as.character(ds_list$V1[idx])
  # check the PRS data are indeed correct
  ds_eid <- rec_data %>% 
    filter(diag_icd10  == ds_id) %>%
    select(eid)
  PRS_profile <- read.table(paste0(DIR, ds_id, ".profile"), header =T) %>%
    mutate(SCORE = SCORE/sd(SCORE)) %>%
    semi_join(ds_eid, c("FID" = "eid"))
  # get the disease of different groups
  j <- match(as.numeric(ds_id), para$list_above500occu$diag_icd10)
  subtype <- read.table(paste0(cc_dir, ds_id, "topic_id_subtp.txt"))
  topic_list <- subtype$V3
  
  # regression analysis
  try({
    for(topic_id in topic_list){
      df_umap_patient <- data.frame(eid = para$eid, loadings = patient_loadings[,topic_id])  %>%
        mutate(loadings = loadings/sd(loadings))
      regression_df <- PRS_profile %>% 
        left_join(df_umap_patient, by = c("FID" = "eid" )) 
      regression_md <- lm(SCORE ~ loadings, data = regression_df)
      coefficients[idx, topic_id] <- summary(regression_md)$coefficients[2,1]
      p_values[idx, topic_id] <- summary(regression_md)$coefficients[2,4]
    }
  })
  print(paste(ds_list$phenotype[idx], " done"))

}

pasted_rslt <- matrix(mapply(function(x,y) paste0(as.character(x), " (P = ", as.character(y), ")"), 
                             round(coefficients,digits = 4), round(p_values, digits = 4)), dim(ds_list)[1])
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
colnames(pasted_rslt) <- topic_name
df_subtype_regression <- data.frame(Phecode = ds_list$V1, Phenotype = ds_list$phenotype, pasted_rslt) 
write.csv(df_subtype_regression, paste0("~/Desktop/comorbidity/paper_writing/Production_figures/subtype_regression.csv" ))
subtype_regression <- list(ds_list, p_values, coefficients)
save(subtype_regression, file = paste0("~/Desktop/comorbidity/paper_writing/Production_figures/subtype_regression.Rdata" ))
ds_list <- subtype_regression[[1]]
p_values <- subtype_regression[[2]]
p_values <- bind_cols(ds_list$phenotype, data.frame(p_values))
colnames(p_values) <-  c("phenotype","MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
p_values <- p_values %>%
  pivot_longer(cols = !phenotype, names_to = "topic", values_to = "p_value") %>%
  filter(!is.na(p_value))
coefficients <- subtype_regression[[3]]
coefficients <- bind_cols(ds_list$phenotype, data.frame(coefficients))
colnames(coefficients) <- c("phenotype","MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
coefficients <- coefficients %>%
  pivot_longer(cols = !phenotype, names_to = "topic", values_to = "coefficients") %>%
  filter(!is.na(coefficients))
regression_PRS <- p_values %>%
  left_join(coefficients, by = c("phenotype", "topic")) 
regression_PRS$FDR <- p.adjust(regression_PRS$p_value)
regression_PRS %>% select(phenotype, topic, coefficients, p_value, FDR) %>%
  filter(FDR < 0.1) %>%
  write.csv("~/Desktop/comorbidity/paper_writing/Production_figures/PRS_regression.csv", row.names = F)

# Fig 4 sub figure
# first step: add one analysis of the ARP topic for hypertension
load("~/Desktop/comorbidity/paper_writing/Production_figures/subtype_regression.Rdata" )
p_values <- subtype_regression[[2]]
coefficients <- subtype_regression[[3]]
idx <- 18
topic_id <- 2
ds_id <- as.character(ds_list$V1[idx])
# check the PRS data are indeed correct
ds_eid <- rec_data %>% 
  filter(diag_icd10  == ds_id) %>%
  select(eid)
PRS_profile <- read.table(paste0(DIR, ds_id, ".profile"), header =T) %>%
  mutate(SCORE = SCORE/sd(SCORE)) %>%
  semi_join(ds_eid, c("FID" = "eid"))
df_umap_patient <- data.frame(eid = para$eid, loadings = patient_loadings[,topic_id])  %>%
  mutate(loadings = loadings/sd(loadings))
regression_df <- PRS_profile %>% 
  left_join(df_umap_patient, by = c("FID" = "eid" )) 
regression_md <- lm(SCORE ~ loadings, data = regression_df)
coefficients[idx, topic_id] <- summary(regression_md)$coefficients[2,1]
p_values[idx, topic_id] <- summary(regression_md)$coefficients[2,4]

idx_plot <- which(ds_list$V1 %in% c(272.11, 250.2, 401.1, 278.1))
p_plot <- p_values[idx_plot, ]
colnames(p_plot) <- topic_name
rownames(p_plot) <- c("Type 2 diabetes", "Hypercholesterolemia", "Obesity", "Essential hypertension")
p_plot <- melt(p_plot[,c(2,5,6,7)]) %>%
  rename(P = value)
coef_plot <- coefficients[idx_plot, ]
colnames(coef_plot) <- topic_name
rownames(coef_plot) <- c("Type 2 diabetes", "Hypercholesterolemia", "Obesity", "Essential hypertension")
coef_plot <- melt(coef_plot[,c(2,5,6,7)])%>%
  rename(coefficient = value)

df_plot <- coef_plot %>%
  left_join(p_plot, by = c("Var1", "Var2")) %>%
  mutate(P_levels = -log10(P) ) 
# use color that are different
paletteMartin[[4]]
Brown2Blue12Steps[[10]]
plt <- ggplot(df_plot) +
  geom_point(aes(x = Var2, y = Var1, fill = coefficient, size = P_levels), shape = 21) +
  scale_fill_gradient2(low = "#32E3FF", mid = "white", high = "#ff6db6") + 
  scale_size_area(name = expression(-log[10](P)), max_size = 20) + 
  labs(x = "", y = "") +
  theme_bw(base_size = 20)+
  theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
  )
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/Subtype_liability_interaction.png", plt, width = 6, height = 4.5)
plt <- ggplot(df_plot) +
  geom_point(aes(x = Var2, y = Var1, fill = coefficient, size = P_levels), shape = 21) +
  scale_fill_gradient2(low = "#32E3FF", mid = "white", high = "#ff6db6") + 
  scale_size_area(name = expression(-log[10](P)), max_size = 40) + 
  labs(x = "", y = "") +
  theme_bw(base_size = 20)
legend_plt <- cowplot::get_legend(plt)
grid.newpage()
legend_plt <- grid.draw(legend_plt)


###################
# perform Fst
###################
subtypes_list <- read.table( "subtype_disease_topic_list.txt", sep="\t", header = F)
subtype_diseases_list <- subtypes_list %>%
  group_by(V1) %>%
  slice(1)
subtype_fst <- read.csv("Association_analysis/subtypes_Fst.csv") %>%
  left_join(select(subtype_diseases_list,V1, V4 ) , by = c("disease" = "V1")) %>%
  rename(phenotype = V4)

# plot Fst P-value
subtype_fst$q_value <- p.adjust(subtype_fst$p_fst, method = "fdr")
Fst_p <- subtype_fst %>%
  filter(!is.na(subtype_fst$p_fst)) %>%
  mutate(p_fst = -log10(p_fst)) %>%
  arrange(p_fst) 
Fst_p$sim_p <- sort(-log10(runif(dim(Fst_p)[1])))
plt <- ggplot(Fst_p, aes(x= sim_p, y = p_fst, label=phenotype)) +
  geom_point(aes( color = Fst), size = 3) + 
  scale_color_gradient(low = "grey80",  high=red) + 
  geom_abline(slope = 1) + 
  geom_label_repel(label=ifelse(Fst_p$Fst > 5 * 10^(-5),as.character(Fst_p$phenotype),''), max.overlaps = 15, 
                   nudge_y = ifelse(Fst_p$p_fst > 2.5, -.5, 0),
                   nudge_x = ifelse(Fst_p$p_fst > 2.5, -.2, 0), size = 5) + # do a bit of nudging 
  labs(x="Expected -log10P", y="Observed -log10P") + 
  theme_bw(base_size = 15) +
  theme(legend.position = "None")
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/Fst_subtypes.png", plt, width = 8, height = 8)
plt <- ggplot(Fst_p, aes(x= sim_p, y = p_fst, label=phenotype)) +
  geom_point(aes( color = Fst), size = 3) + 
  scale_color_gradient(low = "grey80",  high=red) + 
  geom_abline(slope = 1) + 
  geom_label_repel(label=ifelse(Fst_p$p_fst > 2,as.character(Fst_p$phenotype),''), max.overlaps = 20) +
  labs(x="Expected -log10P", y="Observed -log10P") + 
  theme_bw(base_size = 15) 
legend_plt <- cowplot::get_legend(plt)
grid.newpage()
legend_plt <- grid.draw(legend_plt)
print(paste0("proportion of max topic value > 0.95: ", mean(top_1_value[,1] > 0.95) ) )
sum(p.adjust(Fst_p$p_fst, method = "fdr") < 0.1)

##########################################
# Fst with topic matched controls
##########################################
subtypes_list <- read.table( "subtype_disease_topic_list.txt", sep="\t", header = F)
subtype_diseases_list <- subtypes_list %>%
  group_by(V1) %>%
  slice(1)
subtype_fst <- read.csv("Association_analysis/subtypes_Fst_matched_topic.csv") %>%
  left_join(select(subtype_diseases_list,V1, V4 ) , by = c("disease" = "V1")) %>%
  rename(phenotype = V4)

# plot Fst P-value
subtype_fst$q_value <- p.adjust(subtype_fst$p_fst, method = "fdr")
Fst_p <- subtype_fst %>%
  filter(!is.na(subtype_fst$p_fst)) %>%
  mutate(p_fst = -log10(p_fst)) %>%
  arrange(p_fst) 
Fst_p$sim_p <- sort(-log10(runif(dim(Fst_p)[1])))
plt <- ggplot(Fst_p, aes(x= sim_p, y = p_fst, label=phenotype)) +
  geom_point(aes( color = Fst), size = 3) + 
  scale_color_gradient(low = "grey80",  high=red) + 
  geom_abline(slope = 1) + 
  geom_label_repel(label=ifelse(Fst_p$Fst > 5 * 10^(-5),as.character(Fst_p$phenotype),''), max.overlaps = 15, 
                   nudge_y = ifelse(Fst_p$p_fst > 2.5, -.5, 0),
                   nudge_x = ifelse(Fst_p$p_fst > 2.5, -.2, 0), size = 5) + # do a bit of nudging 
  labs(x="Expected -log10P", y="Observed -log10P") + 
  theme_bw(base_size = 15) +
  theme(legend.position = "None")
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/Fst_topic_matched_subtypes.png", plt, width = 8, height = 8)
plt <- ggplot(Fst_p, aes(x= sim_p, y = p_fst, label=phenotype)) +
  geom_point(aes( color = Fst), size = 3) + 
  scale_color_gradient(low = "grey80",  high=red) + 
  geom_abline(slope = 1) + 
  geom_label_repel(label=ifelse(Fst_p$p_fst > 2,as.character(Fst_p$phenotype),''), max.overlaps = 20) +
  labs(x="Expected -log10P", y="Observed -log10P") + 
  theme_bw(base_size = 15) 
legend_plt <- cowplot::get_legend(plt)
grid.newpage()
legend_plt <- grid.draw(legend_plt)
print(paste0("proportion of max topic value > 0.95: ", mean(top_1_value[,1] > 0.95) ) )
sum(Fst_p$q_value < 0.1)

subtype_fst %>%
  filter(!is.na(subtype_fst$p_fst)) %>%
  rename(P = p_fst, Phecode = disease, Disease = phenotype) %>%
  write.csv("~/Desktop/comorbidity/paper_writing/supplementary_files/fst.csv", row.names = F)

# count number of subtypes involved
significant_ds <- Fst_p %>% filter(q_value < 0.1) %>% pull(disease)  
not_included_ds <- c(716.90, 401.1, 272.11)
subtypes_list %>%
  filter(!(V1 %in% not_included_ds)) %>%
  dim()
subtypes_list %>%
  filter((V1 %in% significant_ds)) %>%
  dim()
ds_id <- 250.2
permutation_Fst <- read.csv(paste("Association_analysis/",ds_id,"_control_Fst.csv",sep=""))
targetFst <- Fst_p$Fst[which(Fst_p$disease == ds_id)]
Fstmean <- mean(targetFst - permutation_Fst$weighted_Fst)
Fstse <- sd(targetFst - permutation_Fst$weighted_Fst)/sqrt(length(permutation_Fst$weighted_Fst))
###########################################
# CC-gwas
###########################################
library(qqman)
cc_dir <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/CC_gwas/"
ds_list <- read.table("top_hererogeneous_disease.txt")
# playing with the code
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("V1" = "phecode"))

colors_topic <- c(cbPalette, paletteMartin[c(8,11)])
names(colors_topic) <- as.character(1:10)
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
# maybe save hits using the data frame
hits <- matrix(nrow = length(ds_list$V1), ncol = 1)
for(idx in 1:length(ds_list$V1)){
  ds_id <- ds_list$V1[idx]
  subtype <- read.table(paste0(cc_dir, ds_id, "topic_id_subtp.txt"))
  ds_name <- ds_list$phenotype[idx]
  if(dim(subtype)[1] > 1){
    try({
      for(ccgwas in 1:(dim(subtype)[1])){
        assoc.linear <- read.table(paste0(cc_dir, ds_id, "linear_ccgwas_", ccgwas, ".assoc.linear"), header = T)
        plt <- plt_manhattan(assoc.linear, paste0(ds_name, ": topic ", topic_name[subtype$V3[ccgwas]]))
        ggsave(paste0("../figures/GxTopic/","ccgwas_",ds_id,"_topic", topic_name[subtype$V3[ccgwas]],".png"), plt, width = 10, height = 4)
      }
    })
  }
}

# plot specific disease for main figure
ds_id <- "272.11"
ccgwas <- 3
idx <- 11
subtype <- read.table(paste0(cc_dir, ds_id, "topic_id_subtp.txt"))
ds_name <- ds_list$phenotype[idx]
assoc.linear <- read.table(paste0(cc_dir, ds_id, "linear_ccgwas_", ccgwas, ".assoc.linear"), header = T)

don <- assoc.linear %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(assoc.linear, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)
axisdf <- don %>% 
  group_by(CHR) %>% 
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
plt <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
  scale_color_manual(values = rep(c(blue, grey), 22 )) +
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0.1, 0) ) +     # remove space between plot area and x axis
  geom_hline(aes(yintercept = 7.3), color = "red", linetype = "dashed") +
  labs(x = "Chromosome number", y = expression(-log[10](P)), title = paste0(ds_name, ": ", topic_name[subtype$V3[ccgwas]], " topic liability")) +
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


ggsave(paste0("../paper_writing/Production_figures/","ccgwas_",ds_id,"_topic", topic_name[subtype$V3[ccgwas]],".png"), plt, width = 10, height = 4)

# qq plot 
df_qq.logistic <- assoc.linear %>%
  mutate(logP = -log10(P)) %>%
  arrange(logP)
df_qq.logistic$sim_logP <- sort( - log10(0.1 * runif(dim(df_qq.logistic)[1])))
ggplot(df_qq.logistic) + 
  geom_point(aes(x = sim_logP, y = logP)) + 
  geom_abline(slope = 1) 
  
#############################################
# effect size comparison: CaseCase v.s. Case-control
#############################################
cc_dir <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/CC_gwas/"
# step 1: extract all the SNPs
ds_id <- "272.11"
assoc_CaseControl <- read.table(paste0(cc_dir, ds_id, ".assoc.logistic"), header = T)
clump_CaseControl <- read.table(paste0(cc_dir, ds_id, "_CaseControl.clumped"), header = T) %>%
  select(SNP, CHR, BP) 


subtype <- read.table(paste0(cc_dir, ds_id, "topic_id_subtp.txt"))
CaseCase_id <- 1
assoc_CaseCase <- read.table(paste0(cc_dir, ds_id, "no_P_threshold_linear_ccgwas_", CaseCase_id, ".assoc.linear"), header = T) 

joint_set <- clump_CaseControl %>%
  group_by(CHR) %>%
  arrange(BP, .by_group = TRUE) %>%
  left_join(select(assoc_CaseCase, SNP, A1, BETA), by = "SNP") %>%
  left_join(select(assoc_CaseControl, SNP, A1, OR), by = c("SNP", "A1")) %>%
  mutate(logOR = log(OR))
plt <- ggplot(joint_set) +
  geom_point(aes(x = BETA, y = logOR), alpha = 0.5) +
  scale_y_continuous(limits = c(-0.4,0.4)) + 
  scale_x_continuous(limits = c(-0.02,0.02)) +  
  labs(x = "Effect size on topic liability", y = paste0("Effect size on disease " , expression(log(OR)))) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
ggsave(paste0("../paper_writing/Production_figures/","Effect_comparison",ds_id,"_topic", topic_name[subtype$V3[CaseCase_id]],".png"), plt, width = 4, height = 4)


# comparing effect with two topics
cc_dir <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/CC_gwas/"
# step 1: extract all the SNPs
ds_id <- "250.2" # 272.11
clump_CaseControl <- read.table(paste0(cc_dir, ds_id, "_CaseControl.clumped"), header = T) %>%
  select(SNP, CHR, BP) 

subtype <- read.table(paste0(cc_dir, ds_id, "topic_id_subtp.txt"))
CaseCase_id <- 1
assoc_CaseCase <- read.table(paste0(cc_dir, ds_id, "no_P_threshold_linear_ccgwas_", CaseCase_id, ".assoc.linear"), header = T) 
CaseCase_id2 <- 3
assoc_CaseCase2 <- read.table(paste0(cc_dir, ds_id, "no_P_threshold_linear_ccgwas_", CaseCase_id2, ".assoc.linear"), header = T) 
joint_set <- clump_CaseControl %>%
  group_by(CHR) %>%
  arrange(BP, .by_group = TRUE) %>%
  left_join(select(assoc_CaseCase, SNP, A1, BETA, P), by = "SNP") %>%
  left_join(select(assoc_CaseCase2, SNP, A1, BETA, P), by = c("SNP", "A1")) # %>%
  # filter(P.x < 0.05 | P.y < 0.05)

assoc_CaseControl <- read.table(paste0(cc_dir, ds_id, ".assoc.logistic"), header = T)
background_SNP <- assoc_CaseControl %>%
  filter(P > 0.01) %>%
  sample_n(10000) %>%
  select(SNP, CHR, BP) %>%
  left_join(select(assoc_CaseCase, SNP, A1, BETA, P), by = "SNP") %>%
  left_join(select(assoc_CaseCase2, SNP, A1, BETA, P), by = c("SNP", "A1")) # %>%
  # filter(P.x < 0.05 | P.y < 0.05)

plt <- ggplot(background_SNP) +
  geom_point(aes(x = BETA.x, y = BETA.y),  color = grey, alpha = 0.1) +
  scale_y_continuous(limits = c(-0.02,0.02)) + 
  scale_x_continuous(limits = c(-0.02,0.02)) +  
  labs(x = paste0("Effect size on topic ", topic_name[subtype$V3[CaseCase_id]], " liability"), 
       y = paste0("Effect size on topic ", topic_name[subtype$V3[CaseCase_id2]], " liability")) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) + 
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed") +
  geom_vline(aes(xintercept = 0), color = "red", linetype = "dashed") + 
  geom_point(data = joint_set, aes(x = BETA.x, y = BETA.y),color = red, alpha = 0.8) 

ggsave(paste0("../paper_writing/Production_figures/","Effect_comparison",ds_id,"_topic", topic_name[subtype$V3[CaseCase_id]],"_",topic_name[subtype$V3[CaseCase_id2]], ".png"), plt, width = 4, height = 4)

  


# clump_CaseCase <- read.table(paste0(cc_dir, ds_id, "_CaseCase_Topic", CaseCase_id, ".clumped"), header = T) %>%
#   select(SNP, CHR, BP)
# joint_set <- bind_rows(clump_CaseControl, clump_CaseCase) %>%
#   group_by(CHR) %>%
#   arrange(BP, .by_group = TRUE) %>%
#   left_join(select(assoc_CaseCase, SNP, BETA), by = "SNP") %>%
#   left_join(select(assoc_CaseControl, SNP, OR), by = "SNP") %>%
#   mutate(logOR = log(OR))

#############################################
# GxTopic interaction
#############################################
# compute qq plot for all interaction analysis
cc_dir <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/CC_gwas/"
ds_list <- read.table("top_hererogeneous_disease.txt")
# playing with the code
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("V1" = "phecode"))
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
hits <- read.table(paste0(cc_dir, "272.11", ".assoc.logistic"), header = T)
for(idx in 1:length(ds_list$V1)){
  ds_id <- ds_list$V1[idx]
  subtype <- read.table(paste0(cc_dir, ds_id, "topic_id_subtp.txt"))
  ds_name <- ds_list$phenotype[idx]
  if(dim(subtype)[1] > 1){
    try({
      for(ccgwas in 1:(dim(subtype)[1])){
        gxtopic <- read.table(paste0(cc_dir, ds_id, "GxTopic_", ccgwas, ".qassoc.gxe"), header = T) %>%
          mutate(logP = -log10(P_GXE)) %>%
          filter(!is.na(logP)) %>%
          left_join(select(hits, SNP, BP), by = "SNP" )
        df_qq.logistic <- gxtopic %>%
          arrange(logP)
        df_qq.logistic$sim_logP <- sort( - log10( runif(dim(df_qq.logistic)[1])))
        plt <- ggplot(df_qq.logistic) + 
          geom_point(aes(x = sim_logP, y = logP), size = 0.5, alpha = 0.5) + 
          geom_abline(slope = 1, linetype = "dashed", color = red) +
          theme(legend.position = "none",panel.background=element_blank()) + 
          xlab(expression("Observed" * -log[10](P)) ) + ylab(expression("Expected" * -log[10](P))) +
          ggtitle(paste0(ds_name, " x topic: ", topic_name[subtype$V3[ccgwas]]))
        
        ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/", ds_id, "_topic", 
                      topic_name[subtype$V3[ccgwas]],"qq_GxTopic.png"), width = 4, height = 4,plt)
      }
      })
    }
}

# GxTopic for BMI
subtype <- read.table(paste0(cc_dir,"278.1topic_id_subtp.txt"))
for(ccgwas in 1:3){
  gxtopic <- read.table(paste0(cc_dir, "BMI_allcov_GxTopic_", ccgwas, ".assoc.linear"), header = T) %>%
    filter(TEST == "ADDxCOV1") %>%
    mutate(logP = -log10(P)) %>%
    filter(!is.na(logP)) %>%
    left_join(select(hits, SNP, BP), by = "SNP" )
  df_qq.logistic <- gxtopic %>%
    arrange(logP)
  df_qq.logistic$sim_logP <- sort( - log10( runif(dim(df_qq.logistic)[1])))
  plt <- ggplot(df_qq.logistic) + 
    geom_point(aes(x = sim_logP, y = logP), size = 0.5, alpha = 0.5) + 
    geom_abline(slope = 1, linetype = "dashed", color = red) +
    theme(legend.position = "none",panel.background=element_blank()) + 
    xlab(expression("Observed" * -log[10](P)) ) + ylab(expression("Expected" * -log[10](P))) +
    ggtitle(paste0("BMI x topic: ", topic_name[subtype$V3[ccgwas]]))
  
  ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/BMI_topic", 
                topic_name[subtype$V3[ccgwas]],"qq_GxTopic.png"), width = 4, height = 4,plt)
}
# also save the effect on BMI controlled the topic loading
subtype <- read.table(paste0(cc_dir,"278.1topic_id_subtp.txt"))
for(ccgwas in 1:3){
  gxtopic <- read.table(paste0(cc_dir, "BMI_allcov_GxTopic_", ccgwas, ".assoc.linear"), header = T) %>%
    filter(TEST == "ADD") %>%
    mutate(logP = -log10(P)) %>%
    filter(!is.na(logP)) %>%
    left_join(select(hits, SNP, BP), by = "SNP" )
  df_qq.logistic <- gxtopic %>%
    arrange(logP)
  df_qq.logistic$sim_logP <- sort( - log10( runif(dim(df_qq.logistic)[1])))
  plt <- ggplot(df_qq.logistic) + 
    geom_point(aes(x = sim_logP, y = logP), size = 0.5, alpha = 0.5) + 
    geom_abline(slope = 1, linetype = "dashed", color = red) +
    theme(legend.position = "none",panel.background=element_blank()) + 
    xlab(expression("Observed" * -log[10](P)) ) + ylab(expression("Expected" * -log[10](P))) +
    ggtitle(paste0("BMI controled by topic: ", topic_name[subtype$V3[ccgwas]]))
  
  ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/BMI_topic", 
                topic_name[subtype$V3[ccgwas]],"qq_G.png"), width = 4, height = 4,plt)
}


###############################################
# plot GxTopic specific example: hypercholesterolemia
###############################################
cc_dir <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/CC_gwas/"
# baseline
ds_id <- "272.11" # 272.11
hits <- read.table(paste0(cc_dir, ds_id, ".assoc.logistic"), header = T)
qqplot(-log10(runif(dim(hits)[1])) , -log10(hits$P))
abline(0,1)

# model: Disease = G + Topic + Topic*G (not very netural)
ds_id <- "272.11" # 272.11
CaseCase_id <- 1
gxtopic <- read.table(paste0(cc_dir, ds_id, "GxTopic_", CaseCase_id, ".assoc.logistic"), header = T) %>%
  filter(TEST == "USER_1DF", STAT != -9)
gxtopic %>% 
  filter(P < 10^(-5))
qqplot(-log10(runif(dim(gxtopic)[1])) , -log10(gxtopic$P))
abline(0,1)

# model: Topic = G + Disease + Disease*G (stratefied by disaease)
ds_id <- "272.11" # 272.11
CaseCase_id <- 3
gxtopic <- read.table(paste0(cc_dir, ds_id, "GxTopic_", CaseCase_id, ".qassoc.gxe"), header = T) %>%
  mutate(logP = -log10(P_GXE)) %>%
  filter(!is.na(logP)) %>%
  left_join(select(hits, SNP, BP), by = "SNP" )
df_qq.logistic <- gxtopic %>%
  arrange(logP)
df_qq.logistic$sim_logP <- sort( - log10( runif(dim(df_qq.logistic)[1])))
plt <- ggplot(df_qq.logistic) + 
  geom_point(aes(x = sim_logP, y = logP), size = 0.5, alpha = 0.5) + 
  geom_abline(slope = 1, linetype = "dashed", color = red) +
  theme(legend.position = "none",panel.background=element_blank()) + 
  xlab(expression("Observed" * -log[10](P)) ) + ylab(expression("Expected" * -log[10](P)))

ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/", ds_id, "_topic", topic_name[subtype$V3[CaseCase_id]],"qq_GxTopic.png"), width = 4, height = 4,plt)

# plot a manhattan plot for all the hits
ds_id <- "272.11"
ccgwas <- 3
idx <- 11
subtype <- read.table(paste0(cc_dir, ds_id, "topic_id_subtp.txt"))
ds_name <- "Hypercholesterolemia"

don <- gxtopic %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gxtopic, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)
axisdf <- don %>% 
  group_by(CHR) %>% 
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
plt <- ggplot(don, aes(x=BPcum, y=-log10(P_GXE))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
  scale_color_manual(values = rep(c(blue, grey), 22 )) +
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0.1, 0) ) +     # remove space between plot area and x axis
  geom_hline(aes(yintercept = 7.3), color = "red", linetype = "dashed") +
  labs(x = "Chromosome number", y = expression(-log[10](P)), title = paste0(ds_name, " x ", topic_name[subtype$V3[ccgwas]], " topic liability")) +
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


ggsave(paste0("../paper_writing/Production_figures/","GxTopic_",ds_id,"_topic", topic_name[subtype$V3[ccgwas]],".png"), plt, width = 10, height = 4)

hits <- read.table(paste0(cc_dir, ds_id, ".assoc.logistic"), header = T)
don <- hits %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(hits, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)
axisdf <- don %>% 
  group_by(CHR) %>% 
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
plt <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
  scale_color_manual(values = rep(c(blue, grey), 22 )) +
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0.1, 0) ) +     # remove space between plot area and x axis
  geom_hline(aes(yintercept = 7.3), color = "red", linetype = "dashed") +
  labs(x = "Chromosome number", y = expression(-log[10](P)), title = paste0(ds_name,  " GWAS")) +
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


ggsave(paste0("../paper_writing/Production_figures/","GWAS_",ds_id,".png"), plt, width = 10, height = 4)


# small illustration figure
N <- 500
topic_causal_state <- (runif(N) > 0.5) * 1
topic_lia <- (1 - topic_causal_state) + rnorm(N) * 0.2/(topic_causal_state + 1)
topic_eff <- 1 + 0.2*topic_causal_state + rnorm(N) * 0.2 + 0.1 * topic_lia + 0 * topic_causal_state
# plot the genetic effect v.s. topic liability
df_liability <- data_frame(snp_label = as.character(topic_causal_state), liability = topic_lia, effect_size = topic_eff)
plt <- ggplot(df_liability, aes(x = liability, y = effect_size)) + 
  geom_point(aes(color = snp_label), size = 2) +
  theme(legend.position = "none",panel.background=element_blank(), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  xlab("") + ylab("") + 
  scale_colour_manual(name="Genotype",values=c("0" = grey, "1" = red)) 
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/effect_topic.png", width = 4, height = 4,plt)


# quality check: qq plot of GxTopic for all SNPs in CVD topic 
c_dir <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/CC_gwas/"
subtype <- read.table(paste0(cc_dir, ds_id, "topic_id_subtp.txt"))
# baseline
ds_id <- "272.11" # 272.11
hits <- read.table(paste0(cc_dir, ds_id, ".assoc.logistic"), header = T)
# exclude the hits that are nominally significant for hypercholesterolemia
exclude_hits <- hits %>% 
  filter(P < 0.05)
topic_id <- 7
topic_gwas <- read.table(paste0("../association_results/topic", topic_id, "_loading_K10_rep10.assoc.linear"), header = T) %>%
  anti_join(exclude_hits, by = "SNP") %>%
  filter(P < 0.001)

ds_id <- "272.11" # 272.11
CaseCase_id <- 3
gxtopic <- read.table(paste0(cc_dir, ds_id, "GxTopic_", CaseCase_id, ".qassoc.gxe"), header = T) %>%
  semi_join(topic_gwas, by = "SNP") %>%
  mutate(logP = -log10(P_GXE)) %>%
  filter(!is.na(logP)) 
df_qq.logistic <- gxtopic %>%
  arrange(logP)
df_qq.logistic$sim_logP <- sort( - log10( runif(dim(df_qq.logistic)[1])))
plt <- ggplot(df_qq.logistic) + 
  geom_point(aes(x = sim_logP, y = logP), size = 0.5, alpha = 0.5) + 
  geom_abline(slope = 1, linetype = "dashed", color = red) +
  theme(legend.position = "none",panel.background=element_blank()) + 
  xlab(expression("Observed" * -log[10](P)) ) + ylab(expression("Expected" * -log[10](P)))

ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/topic_only", topic_name[subtype$V3[CaseCase_id]],"qq_GxTopic.png"), width = 4, height = 4,plt)




############################################
# compare PCA, LDA, ageLDA
############################################
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
code2id <- function(x){
  return( match(x, ds_list$diag_icd10))
}

# here I am rounding the disease time to year for computation efficiency
Ds_matrix <- rec_data %>%
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
PCA_ukbb <- scaled_ds_matrix %>%
  prcomp()
data_pca <- Ds_matrix %>%
  select(eid) %>%
  cbind(PCA_ukbb[["x"]])
PCA_results <- list(data_pca, PCA_ukbb)
save(PCA_results, file = paste0("../Results/","UKBB_phecodesPCA.RData"))


variance_pc <- PCA_ukbb$sdev^2/sum(PCA_ukbb$sdev^2)
plot(variance_pc[1:20])
longData<-melt(PCA_ukbb$rotation[,1:10])
longData<-longData[longData$value!=0,]
plt <- ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient2(low=blue, mid = "white", high=red) +
  labs(x="PC", y="Diseases", title=paste0("Top PCs")) +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
ggsave(paste0("../figures/","UKBB_PCA_top10.png"), plt, width = 6, height = 10)

# LDA
# first find the best rep
DIR <- "~/Desktop/comorbidity/Results/"
K_chosen <- 10
pt <- paste0("^BaselinLDA_model_output_PheCode_K", K_chosen,"_P_rep*")
temp <- list.files(paste(DIR, sep=""), pattern=pt)
lb_rep <- data_frame(lower_bound = as.numeric())
for(rep_id in 1:length(temp)){
  load(paste0(DIR,paste0("BaselinLDA_model_output_PheCode_K", K_chosen,"_P_rep",rep_id, ".Rdata")))
  cvrg_lb <-  model_output[[2]] %>% 
    filter(!is.na(Lower_bound)) %>%
    slice_tail %>% 
    pull(2)
  lb_rep <- lb_rep %>% 
    add_row(lower_bound = cvrg_lb)
}
rep_id <- order(lb_rep$lower_bound,decreasing = T)[1]
print(temp[rep_id])
load(paste0(DIR,"Run_BaselinLDA_Phecode_K10_rep6.Rdata"))
longData<-melt(para$beta)
longData<-longData[longData$value!=0,]
plt <- ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient2(low=blue, mid = "white", high=red) +
  labs(x="topics", y="Diseases", title=paste0("LDA")) +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
ggsave(paste0("../figures/","UKBB_LDA.png"), plt, width = 6, height = 10)

# basicLDA coocure
assignment_per_ds <- as.factor(sapply(1:para$D, function(j) which.max(colMeans(para$unlist_zn[para$ds_list[[j]]$id,])) ))
matrix_coocurre <- matrix(NA,para$D,para$D)
topic_data <- melt(matrix_coocurre) %>%
  filter(assignment_per_ds[Var2] == assignment_per_ds[Var1]) %>%
  mutate(topic = assignment_per_ds[Var2])
library(colorBlindness)
plt <- ggplot(topic_data, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=topic)) + 
  scale_fill_manual(values= PairedColor12Steps[c(1:2,5:12)],na.value = "white") +
  labs(x="Diseases", y="Diseases", title="basicLDA") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11)) 
# scale_x_discrete(breaks=c("5","20","55", "75", "90", "125", "160", "200", "240", "270", "290", "320","348"),labels=unique(diseases_type)) 
ggsave("../figures/BasicLDAtopic_disease_distribution.png", plt, width = 11, height = 10)

########################################
# extract s.d. change of PRS in subtypes
########################################
ds_id <- 495 # 250.2
j <- match(as.numeric(ds_id), para$list_above500occu$diag_icd10)
ds_eid <- rec_data %>% 
  filter(diag_icd10  == ds_id) %>%
  select(eid)
PRS_profile <- read.table(paste0(DIR, ds_id, ".profile"), header =T) %>%
  mutate(SCORE = SCORE/sd(SCORE)) %>%
  semi_join(ds_eid, c("FID" = "eid"))
topic_id <- 10 # 7
topic_specific_id <- which(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max) == topic_id)
topic_specific_set <- para$unlist_Ds_id[para$ds_list[[j]]$id,][topic_specific_id,]

cases_eid <- topic_specific_set %>%
  mutate(topic = 1) %>%
  mutate(Ds_id = para$list_above500occu[Ds_id, 1]) %>%
  left_join(PRS_profile, by = c("eid" = "FID"))

# get the set that are not the target topic
topic_specific_id <- which(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max) != topic_id)
topic_specific_set <- para$unlist_Ds_id[para$ds_list[[j]]$id,][topic_specific_id,]
controls_eid <- topic_specific_set %>%
  mutate(topic = 0) %>%
  mutate(Ds_id = para$list_above500occu[Ds_id, 1]) %>%
  left_join(PRS_profile, by = c("eid" = "FID"))

data_all <- bind_rows(cases_eid, controls_eid)
md_rg <- glm(data = data_all, formula = topic ~ SCORE, family = binomial)
summary(md_rg)
mean(cases_eid$SCORE, na.rm = T) -  mean(controls_eid$SCORE, na.rm = T)

##########################################
# disease heterogeneity analysis
##########################################

#######################################
# k-means permutation test for diseases heterogeneity
######################################
# plot the distribution of largest loading size
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
loadings_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[para$ds_list[[j]]$id,])) %>% t
hist(apply(loadings_per_ds, 1, max), breaks = 50)
which(apply(loadings_per_ds, 1, max) > 0.99)


# ageLDA
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
loadings_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[para$ds_list[[j]]$id,])) %>% t
order_ds <- order(apply(loadings_per_ds, 1, which.max))
longData<-melt(loadings_per_ds[order_ds, ]) 
ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Topic", y="Diseases", title="ageLDA") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
samples <- para$unlist_zn[sample(1:dim(para$unlist_zn)[1],size = 10000), ]
order_incidences <- order(apply(samples, 1, which.max))
longData<-melt(samples[order_incidences, ])
ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Topic", y="Diseases", title="ageLDA") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
# plot examples T2D & T1D
j <- 48 # T2D: 48; obesity: 61; hypertension:110 
order_incidences <- order(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max))
longData<-melt(para$unlist_zn[para$ds_list[[j]]$id,][order_incidences, ]) 
ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Topic", y="Diseases", title="ageLDA") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

heterogeneity_age_test <- matrix(NA, nrow = para$D, ncol = 1)
loadings_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[para$ds_list[[j]]$id,])) %>% t

for(k_num in 2){ # the initial test is only for more than 1 subtypes
  Null_index <- which(sapply(1:para$D, function(j) sum(tail(sort(loadings_per_ds[j,]), (k_num-1) ))) > 0.95)
  for(j in 1:para$D){
    print(paste0("disease id: ", j))
    perm_sz <- 1000 # how many permutation samples to collect
    stats_target <- matrix(0, nrow = 2, ncol = 1) 
    permutate_stats <- matrix(0, nrow = 2, ncol = perm_sz) 
    
    # the ageLDA results
    loadings <- para$unlist_zn[para$ds_list[[j]]$id,]
    # for testing under null one topic:
    target_sz <- dim(loadings)[1]
    stats_target[1] <- rs_kmean(loadings, k_num) # get the target test statistic
    stats_target[2] <- rs_kmean(loadings, k_num - 1)
    for(sp in 1:perm_sz){
      loadingsNull <- para$unlist_zn[para$ds_list[[sample(Null_index,1)]]$id,]
      tot_sz <- dim(loadingsNull)[1]
      permutate_stats[1, sp] <- rs_kmean(loadingsNull[sample(1:tot_sz, target_sz, replace = T), ], k_num)
      permutate_stats[2, sp] <- rs_kmean(loadingsNull[sample(1:tot_sz, target_sz, replace = T), ], k_num - 1)
    }
    heterogeneity_age_test[j,k_num - 1] <- (1 + sum((stats_target[1] - stats_target[2]) < 
                                                      (permutate_stats[1, ] - permutate_stats[2, ]) ) )/perm_sz
  }
  
}

ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode"))
ds_list$p_2_subtype <- -log10(heterogeneity_age_test[,1])
# ds_list$p_3_subtype <- -log10(heterogeneity_age_test[,2])
# ds_list$p_4_subtype <- -log10(heterogeneity_age_test[,3])
# ds_list$p_5_subtype <- -log10(heterogeneity_age_test[,4])
write.csv(ds_list, "test_disease_heterogeneity_pvalues.csv")

ds_list <- ds_list %>%
  arrange(p_2_subtype)
ds_list$p_sim = sort(- log10( runif(para$D)) )
plt <- ggplot(data = ds_list) + 
  geom_point(aes(x = p_sim, y = p_2_subtype), color = grey, size = 2) + 
  geom_abline(slope = 1) + 
  geom_label_repel(aes(x = p_sim, y = p_2_subtype, label=ifelse(p_2_subtype > 2,as.character(phenotype),'')), max.overlaps = 50) +
  labs(x="Expected -Log10(P)", y="Actual -Log10(P)", title="ageLDA disease subtypes")
ggsave("../figures/subtype_ageLDA.png", plt, width = 11, height = 10)

# using the FDR to progressively find the best number of K for each disease
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
loadings_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[para$ds_list[[j]]$id,])) %>% t
ds_list <- read.csv("test_disease_heterogeneity_pvalues.csv")

Null_set <- list()
subtypes_ds <- list()
Null_set[[1]] <- which(sapply(1:para$D, function(j) sum(tail(sort(loadings_per_ds[j,]), 1 ))) > 0.95)
subtypes_ds[[1]] <- ds_list
subtypes_ds[[1]] <- subtypes_ds[[1]] %>% 
  filter(! X %in% Null_set[[1]])
subtypes_ds[[1]]$Q_value <- p.adjust(10^(- subtypes_ds[[1]]$p_2_subtype), method = "fdr")
# using the matrix below to save p value for testing for K > 2
pvalues_progressive_K <- matrix(NA, nrow = para$D, ncol = 4)
for(k_num in 3:4){
  # we will only look at the set that has Q < 0.1
  significant_set_from_previous_K <- subtypes_ds[[k_num - 2]] %>%
    filter(Q_value < 0.2) 
  Null_index <- significant_set_from_previous_K$X[sapply(significant_set_from_previous_K$X, function(j) sum(tail(sort(loadings_per_ds[j,]), (k_num-1) ))) > (1 - 1/(10 * 2^(k_num - 1))) ]
  testing_set <- significant_set_from_previous_K$X[! significant_set_from_previous_K$X %in% Null_index]
  # testing_set <- significant_set_from_previous_K$X
  for(j in testing_set){
    print(paste0("disease id: ", j))
    perm_sz <- 1000 # how many permutation samples to collect
    stats_target <- matrix(0, nrow = 2, ncol = 1) 
    permutate_stats <- matrix(0, nrow = 2, ncol = perm_sz) 
    
    # the ageLDA results
    loadings <- para$unlist_zn[para$ds_list[[j]]$id,]
    # for testing under null one topic:
    target_sz <- dim(loadings)[1]
    stats_target[1] <- rs_kmean(loadings, k_num) # get the target test statistic
    stats_target[2] <- rs_kmean(loadings, k_num - 1)
    for(sp in 1:perm_sz){
      loadingsNull <- para$unlist_zn[para$ds_list[[sample(Null_index,1)]]$id,]
      tot_sz <- dim(loadingsNull)[1]
      permutate_stats[1, sp] <- rs_kmean(loadingsNull[sample(1:tot_sz, target_sz, replace = T), ], k_num)
      permutate_stats[2, sp] <- rs_kmean(loadingsNull[sample(1:tot_sz, target_sz, replace = T), ], k_num - 1)
    }
    pvalues_progressive_K[j,k_num - 1] <- (1 + sum((stats_target[1] - stats_target[2]) < 
                                                     (permutate_stats[1, ] - permutate_stats[2, ]) ) )/perm_sz
  }
  Null_set[[k_num - 1]] <- Null_index
  subtypes_ds[[k_num - 1]] <- significant_set_from_previous_K %>%
    filter(X %in% testing_set)
  p_values <- pvalues_progressive_K[which(!is.na(pvalues_progressive_K[,k_num - 1])),k_num - 1] 
  subtypes_ds[[k_num - 1]]$pvalues <- p_values
  subtypes_ds[[k_num - 1]]$Q_value <- p.adjust(p_values, method = "fdr")
}

type_3 <- subtypes_ds[[2]] %>% 
  filter(Q_value < 0.1) %>%
  pull(X)
type_4 <- subtypes_ds[[3]] %>% 
  filter(Q_value < 0.1) %>%
  pull(X)
subtype_diseases_list <- subtypes_ds[[1]] %>% 
  mutate(cluster_number = if_else(Q_value < 0.01, 2, NULL)) %>%
  mutate(cluster_number = if_else(X %in% type_3, 3, cluster_number)) %>% 
  mutate(cluster_number = if_else(X %in% type_4, 4, cluster_number)) %>%
  filter(!is.na(cluster_number)) %>%
  select(X, diag_icd10, occ, ICD10, phenotype, Q_value, p_2_subtype,cluster_number) %>%
  left_join(select(subtypes_ds[[2]], X, pvalues, Q_value), by = c("X") ) %>%
  rename(disease_idx = X, Q_value.1 = Q_value.x, P_value.1 = p_2_subtype, Q_value.2 = Q_value.y, P_value.2 = pvalues)

# save which are the diseases that we are interested in performing heritability analysis
write.csv(subtype_diseases_list, file = "./subtype_disesea_list.csv", row.names = F)
write.table(subtype_diseases_list$diag_icd10, "top_hererogeneous_disease.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)

