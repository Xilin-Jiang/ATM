source("simulation_functions.R")
############################
# figure 2: disease comorbidity heterogeneity identification
############################
# first plot a small schematic
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
ds_id <- 174.11
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
# cases_eid <- cases_eid%>%
#   mutate(age_diag = age_diag + rnorm(dim(cases_eid)[1], sd = 1)) %>%
#   mutate(age_diag = age_diag + rnorm(dim(cases_eid)[1], sd = 1)) %>%
#   mutate(age_diag = age_diag + rnorm(dim(cases_eid)[1], sd = 1))

colors_topic <- c(cbPalette, paletteMartin[c(8,11)])
colors_topic[3] <- "#D55E00"
colors_topic[4] <- "#0072B2"
names(colors_topic) <- as.character(1:10)
cases_eid_little_diff <- cases_eid %>%
  mutate(age_diag = ifelse(topic == 3, age_diag - 3, age_diag + 3))
# plot histogram for each group
plt <- ggplot(cases_eid_little_diff,aes(x = age_diag, fill = factor(topic, levels = as.character(order_ds)), 
                                   colour = factor(topic, levels = as.character(order_ds)))) +
  stat_density(alpha = 0.5, position = "identity", adjust = 2) + # adjust refers to the bandwidth of the kernel
  scale_fill_manual(values = colors_topic) + 
  scale_color_manual(values =colors_topic) +
  theme_bw(base_size = 20) + 
  theme(legend.position = "None") + 
  labs(x="Age", y="Incidence count", title="") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank())
ggsave(paste0("../paper_writing/Production_figures/","simulation_age_distribution_5ydiff.png"), plt, width = 7, height = 4)
# moderate diff
# plot histogram for each group 
plt <- ggplot(cases_eid,aes(x = age_diag, fill = factor(topic, levels = as.character(order_ds)), 
                                        colour = factor(topic, levels = as.character(order_ds)))) +
  geom_density(alpha = 0.5, position = "identity", adjust = 2) + 
  scale_fill_manual(values = colors_topic) + 
  scale_color_manual(values =colors_topic) +
  theme_bw(base_size = 20) + 
  theme(legend.position = "None") + 
  labs(x="Age", y="Incidence count", title="") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank())
ggsave(paste0("../paper_writing/Production_figures/","simulation_age_distribution_10ydiff.png"), plt, width = 7, height = 4)
# with larger difference
cases_eid_larger_diff <- cases_eid %>% 
  mutate(age_diag = ifelse(topic == 3, age_diag + 3, age_diag - 3))
plt <- ggplot(cases_eid_larger_diff,aes(x = age_diag, fill = factor(topic, levels = as.character(order_ds)), 
                            colour = factor(topic, levels = as.character(order_ds)))) +
  geom_density(alpha = 0.5, position = "identity", adjust = 2) + 
  scale_fill_manual(values = colors_topic) + 
  scale_color_manual(values =colors_topic) +
  theme_bw(base_size = 20) + 
  theme(legend.position = "None") + 
  labs(x="Age", y="Incidence count", title="") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank())
ggsave(paste0("../paper_writing/Production_figures/","simulation_age_distribution_20ydiff.png"), plt, width = 7, height = 4)

# another example for supplementary simulation figure
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
ds_id <- 591
j <- match(as.numeric(ds_id), para$list_above500occu$diag_icd10)
cases_eid <- list()
for(topic_id in 1:para$K){
  topic_specific_id <- which(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max) == topic_id)
  topic_specific_set <- para$unlist_Ds_id[para$ds_list[[j]]$id,][topic_specific_id,]
  
  cases_eid[[topic_id]] <- topic_specific_set %>%
    mutate(topic = topic_id) %>%
    mutate(Ds_id = para$list_above500occu[Ds_id, 1])
}
cases_eid <- bind_rows(cases_eid) %>%
  filter(topic %in% c(6,4))

colors_topic <- c(cbPalette, paletteMartin[c(8,11)])
colors_topic[6] <- "#D55E00"
colors_topic[4] <- "#0072B2"
names(colors_topic) <- as.character(1:10)
# plot histogram for each group
plt <- ggplot(cases_eid,aes(x = age_diag, fill = factor(topic, levels = as.character(order_ds)), 
                            colour = factor(topic, levels = as.character(order_ds)))) +
  geom_histogram(alpha = 0.5, position = "stack", binwidth = 1) + 
  scale_fill_manual(values = colors_topic) + 
  scale_color_manual(values =colors_topic) +
  theme_bw(base_size = 20) + 
  theme(legend.position = "None") + 
  labs(x="Age (years)", y="incidence count", title="") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank())
ggsave(paste0("../paper_writing/Production_figures/","simulation_age_distribution_",ds_id,".png"), plt, width = 7, height = 4)

#########################################################
# 2022-07-28 plot schematic using the simulated data
#########################################################
para_orginal_sim <- simulate_age_topic_data(10000, 2, 20, 3)

# need to make sure each diseaes only occur for once
rec_data <- para_orginal_sim$unlist_Ds_id %>% 
  rename(diag_icd10 = Ds_id) %>%
  select(eid, diag_icd10, age_diag) %>%
  group_by(eid, diag_icd10) %>% 
  arrange(age_diag, .by_group = T) %>% # keep only the first records of repeated diagnosis
  slice(1) %>%
  ungroup() %>%
  mutate(rowid = 1:n())

std_age_ds <-rec_data %>%
  group_by(diag_icd10) %>%
  summarise(std_age_ds = sd(age_diag)) 

# compute the sd of age-at-onset for the population

# fig 1: 20 year difference
plot_sche <- rec_data %>%
  filter(diag_icd10 == 1 | diag_icd10 == 11)
plt <- ggplot(plot_sche, aes(x = age_diag, fill = factor(diag_icd10, levels = c(1,11)), colour = factor(diag_icd10, levels = c(1,11)))) +
  geom_histogram(alpha = 0.5, position = "stack", binwidth = 1) + 
  scale_fill_manual(values = c(blue, red)) + 
  scale_color_manual(values = c(blue, red)) +
  theme_bw(base_size = 20) + 
  theme(legend.position = "None") + 
  labs(x="Age (years)", y="incidence count", title="") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank())
ggsave(paste0("../paper_writing/Production_figures/simulation_age_20y_diff.png"), plt, width = 7, height = 4)
# fig 2: 10 year difference
plot_sche <- rec_data %>%
  filter(diag_icd10 == 1 | diag_icd10 == 11) %>%
  mutate(age_diag = ifelse(diag_icd10 == 11, age_diag - 10, age_diag))
plt <- ggplot(plot_sche, aes(x = age_diag, fill = factor(diag_icd10, levels = c(1,11)), colour = factor(diag_icd10, levels = c(1,11)))) +
  geom_histogram(alpha = 0.5, position = "stack", binwidth = 1) + 
  scale_fill_manual(values = c(blue, red)) + 
  scale_color_manual(values = c(blue, red)) +
  theme_bw(base_size = 20) + 
  theme(legend.position = "None") + 
  labs(x="Age (years)", y="incidence count", title="") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank())
ggsave(paste0("../paper_writing/Production_figures/simulation_age_10y_diff.png"), plt, width = 7, height = 4)
# fig : 0 year difference
plot_sche <- rec_data %>%
  filter(diag_icd10 == 1 | diag_icd10 == 10)
plt <- ggplot(plot_sche, aes(x = age_diag, fill = factor(diag_icd10, levels = c(1,10)), colour = factor(diag_icd10, levels = c(1,11)))) +
  geom_histogram(alpha = 0.5, position = "stack", binwidth = 1) + 
  scale_fill_manual(values = c(blue, red)) + 
  scale_color_manual(values = c(blue, red)) +
  theme_bw(base_size = 20) + 
  theme(legend.position = "None") + 
  labs(x="Age (years)", y="incidence count", title="") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank())
ggsave(paste0("../paper_writing/Production_figures/simulation_age_0y_diff.png"), plt, width = 7, height = 4)

#######################################################
# subtype case 1: plot the subtype  recall, precision, accuracy, and AUC
# sample size ratio
#######################################################
recall_age <- NULL
precision_age <- NULL
auc_age <- NULL
# save the results where the late onset group is larger
recall_age_rev <- NULL
precision_age_rev <- NULL
auc_age_rev <- NULL

recall_lda <- NULL
precision_lda <- NULL
auc_lda <- NULL
# save the results where the late onset group is larger
recall_lda_rev <- NULL
precision_lda_rev <- NULL
auc_lda_rev <- NULL

pt <- paste0("^subtype_tests_simulation_type1_rep*")
temp <- list.files(paste("simulation_results", sep=""), pattern=pt)
for(repid in 1:length(temp)){
  load(file = paste0("simulation_results/", temp[repid]))
  
  recall_age <- rbind(recall_age, subtype_tests[[2]][2,])
  precision_age <- rbind(precision_age, subtype_tests[[3]][2,])
  auc_age <- rbind(auc_age, subtype_tests[[4]][2,])
  recall_age_rev <- rbind(recall_age_rev, subtype_tests[[2]][3,])
  precision_age_rev <- rbind(precision_age_rev, subtype_tests[[3]][3,])
  auc_age_rev <- rbind(auc_age_rev, subtype_tests[[4]][3,])
  
  recall_lda <- rbind(recall_lda, subtype_tests[[6]][2,])
  precision_lda <- rbind(precision_lda, subtype_tests[[7]][2,])
  auc_lda <- rbind(auc_lda, subtype_tests[[8]][2,])
  recall_lda_rev <- rbind(recall_lda_rev, subtype_tests[[6]][3,])
  precision_lda_rev <- rbind(precision_lda_rev, subtype_tests[[7]][3,])
  auc_lda_rev <- rbind(auc_lda_rev, subtype_tests[[8]][3,])
}
# recall_age <- cbind(recall_age, recall_age_rev[, ncol(recall_age_rev):1])
# precision_age <- cbind(precision_age, precision_age_rev[, ncol(precision_age_rev):1])
# auc_age <- cbind(auc_age, auc_age_rev[, ncol(auc_age_rev):1])
# 
# recall_lda <- cbind(recall_lda, recall_lda_rev[, ncol(recall_lda_rev):1])
# precision_lda <- cbind(precision_lda, precision_lda_rev[, ncol(precision_lda_rev):1])
# auc_lda <- cbind(auc_lda, auc_lda_rev[, ncol(auc_lda_rev):1])

plot_df <- data.frame(mix_ratio = c( 0.1 * (0:9)/( 0.1 * (0:9) + 1) ), # , 1/( 0.1 * (9:0) + 1) ),
                      mean_recall_age = colMeans(recall_age, na.rm = T),
                      mean_precision_age = colMeans(precision_age, na.rm = T),
                      mean_auc_age =colMeans(auc_age, na.rm = T),
                      se_recall_age = apply(recall_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_precision_age = apply(precision_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_auc_age = apply(auc_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      mean_recall_lda = colMeans(recall_lda, na.rm = T),
                      mean_precision_lda = colMeans(precision_lda, na.rm = T),
                      mean_auc_lda = colMeans(auc_lda, na.rm = T),
                      se_recall_lda = apply(recall_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_precision_lda = apply(precision_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_auc_lda = apply(auc_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))) )
df_20y_diff <- plot_df
plt_recall <- ggplot(plot_df,aes(x = mix_ratio)) +
  geom_line(aes(y = mean_recall_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_recall_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_recall_age, ymin=mean_recall_age-1.96*se_recall_age, ymax = mean_recall_age + 1.96*se_recall_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_recall_lda, ymin=mean_recall_lda-1.96*se_recall_lda, ymax = mean_recall_lda + 1.96*se_recall_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Subtype sample size ratio", y = "Recall:  identify the minor subtype") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,0.5)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","recall_simulation.png"), plt_recall, width = 4, height = 4)

plt_precision <- ggplot(plot_df,aes(x = mix_ratio)) +
  geom_line(aes(y = mean_precision_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_precision_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_precision_age, ymin=mean_precision_age-1.96*se_precision_age, ymax = mean_precision_age + 1.96*se_precision_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_precision_lda, ymin=mean_precision_lda-1.96*se_precision_lda, ymax = mean_precision_lda + 1.96*se_precision_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Subtype sample size ratio", y = "Precision: identify the minor subtype") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,0.5)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","precision_simulation.png"), plt_precision, width = 4, height = 4)

plt_auc <- ggplot(plot_df,aes(x = mix_ratio)) +
  geom_line(aes(y = mean_auc_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_auc_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_auc_age, ymin=mean_auc_age-1.96*se_auc_age, ymax = mean_auc_age + 1.96*se_auc_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_auc_lda, ymin=mean_auc_lda-1.96*se_auc_lda, ymax = mean_auc_lda + 1.96*se_auc_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Subtype sample size ratio", y = "AUPRC") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,0.5)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","auc_simulation.png"), plt_auc, width = 4, height = 4)

# plot AUPRC for classifying the second disease
plot_second_ds_20y <- data.frame(mix_ratio = c( 0.1 * (0:9)/( 0.1 * (0:9) + 1) ), # , 1/( 0.1 * (9:0) + 1) ),
                      mean_auc_age =colMeans(auc_age_rev, na.rm = T),
                      se_auc_age = apply(auc_age_rev, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      mean_auc_lda = colMeans(auc_lda, na.rm = T),
                      se_auc_lda = apply(auc_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))) ) %>%
  filter(mix_ratio != 0)
plt_auc <- ggplot(plot_second_ds_20y,aes(x = mix_ratio)) +
  geom_line(aes(y = mean_auc_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_auc_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_auc_age, ymin=mean_auc_age-1.96*se_auc_age, ymax = mean_auc_age + 1.96*se_auc_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_auc_lda, ymin=mean_auc_lda-1.96*se_auc_lda, ymax = mean_auc_lda + 1.96*se_auc_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Subtype sample size ratio", y = "AUPRC") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,0.5)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","auc_20y_2ndDs.png"), plt_auc, width = 4, height = 4)


# second panel figures for precision and recall, when disease have the same age distribution
recall_age <- NULL
precision_age <- NULL
auc_age <- NULL

recall_lda <- NULL
precision_lda <- NULL
auc_lda <- NULL
pt <- paste0("^subtype_tests_simulation_type1_rep*")
temp <- list.files(paste("simulation_results", sep=""), pattern=pt)
for(repid in 1:length(temp)){
  load(file = paste0("simulation_results/", temp[repid]))
  
  recall_age <- rbind(recall_age, subtype_tests[[2]][1,])
  precision_age <- rbind(precision_age, subtype_tests[[3]][1,])
  auc_age <- rbind(auc_age, subtype_tests[[4]][1,])
  
  recall_lda <- rbind(recall_lda, subtype_tests[[6]][1,])
  precision_lda <- rbind(precision_lda, subtype_tests[[7]][1,])
  auc_lda <- rbind(auc_lda, subtype_tests[[8]][1,])
}
plot_df <- data.frame(mix_ratio = 0.1 * (0:9)/( 0.1 * (0:9) + 1),
                      mean_recall_age = colMeans(recall_age, na.rm = T),
                      mean_precision_age = colMeans(precision_age, na.rm = T),
                      mean_auc_age =colMeans(auc_age, na.rm = T),
                      se_recall_age = apply(recall_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_precision_age = apply(precision_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_auc_age = apply(auc_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      mean_recall_lda = colMeans(recall_lda, na.rm = T),
                      mean_precision_lda = colMeans(precision_lda, na.rm = T),
                      mean_auc_lda = colMeans(auc_lda, na.rm = T),
                      se_recall_lda = apply(recall_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_precision_lda = apply(precision_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_auc_lda = apply(auc_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))) )

plt_recall <- ggplot(plot_df,aes(x = mix_ratio)) +
  geom_line(aes(y = mean_recall_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_recall_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_recall_age, ymin=mean_recall_age-1.96*se_recall_age, ymax = mean_recall_age + 1.96*se_recall_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_recall_lda, ymin=mean_recall_lda-1.96*se_recall_lda, ymax = mean_recall_lda + 1.96*se_recall_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Subtype sample size ratio", y = "Recall:  identify the minor subtype") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,0.5)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","recall_same_age_dist.png"), plt_recall, width = 4, height = 4)

plt_precision <- ggplot(plot_df,aes(x = mix_ratio)) +
  geom_line(aes(y = mean_precision_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_precision_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_precision_age, ymin=mean_precision_age-1.96*se_precision_age, ymax = mean_precision_age + 1.96*se_precision_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_precision_lda, ymin=mean_precision_lda-1.96*se_precision_lda, ymax = mean_precision_lda + 1.96*se_precision_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Subtype sample size ratio", y = "Precision: identify the minor subtype") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,0.5)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","precision_same_age_dist.png"), plt_precision, width = 4, height = 4)

plt_auc <- ggplot(plot_df,aes(x = mix_ratio)) +
  geom_line(aes(y = mean_auc_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_auc_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_auc_age, ymin=mean_auc_age-1.96*se_auc_age, ymax = mean_auc_age + 1.96*se_auc_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_auc_lda, ymin=mean_auc_lda-1.96*se_auc_lda, ymax = mean_auc_lda + 1.96*se_auc_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Subtype sample size ratio", y = "AUPRC") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,0.5)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","auc_same_age_dist.png"), plt_auc, width = 4, height = 4)


# middle panel where there are moderate age difference
recall_age <- NULL
precision_age <- NULL
auc_age <- NULL
# save the results where the late onset group is larger
recall_age_rev <- NULL
precision_age_rev <- NULL
auc_age_rev <- NULL

recall_lda <- NULL
precision_lda <- NULL
auc_lda <- NULL
# save the results where the late onset group is larger
recall_lda_rev <- NULL
precision_lda_rev <- NULL
auc_lda_rev <- NULL

pt <- paste0("^subtype_tests_simulation_type5_rep*")
temp <- list.files(paste("simulation_results", sep=""), pattern=pt)
for(repid in 1:length(temp)){
  load(file = paste0("simulation_results/", temp[repid]))
  
  recall_age <- rbind(recall_age, subtype_tests[[2]][2,])
  precision_age <- rbind(precision_age, subtype_tests[[3]][2,])
  auc_age <- rbind(auc_age, subtype_tests[[4]][2,])
  recall_age_rev <- rbind(recall_age_rev, subtype_tests[[2]][3,])
  precision_age_rev <- rbind(precision_age_rev, subtype_tests[[3]][3,])
  auc_age_rev <- rbind(auc_age_rev, subtype_tests[[4]][3,])
  
  recall_lda <- rbind(recall_lda, subtype_tests[[6]][2,])
  precision_lda <- rbind(precision_lda, subtype_tests[[7]][2,])
  auc_lda <- rbind(auc_lda, subtype_tests[[8]][2,])
  recall_lda_rev <- rbind(recall_lda_rev, subtype_tests[[6]][3,])
  precision_lda_rev <- rbind(precision_lda_rev, subtype_tests[[7]][3,])
  auc_lda_rev <- rbind(auc_lda_rev, subtype_tests[[8]][3,])
}
# recall_age <- cbind(recall_age, recall_age_rev[, ncol(recall_age_rev):1])
# precision_age <- cbind(precision_age, precision_age_rev[, ncol(precision_age_rev):1])
# auc_age <- cbind(auc_age, auc_age_rev[, ncol(auc_age_rev):1])
# 
# recall_lda <- cbind(recall_lda, recall_lda_rev[, ncol(recall_lda_rev):1])
# precision_lda <- cbind(precision_lda, precision_lda_rev[, ncol(precision_lda_rev):1])
# auc_lda <- cbind(auc_lda, auc_lda_rev[, ncol(auc_lda_rev):1])

plot_df <- data.frame(mix_ratio = c( 0.1 * (0:9)/( 0.1 * (0:9) + 1)), # 1/( 0.1 * (9:0) + 1) ),
                      mean_recall_age = colMeans(recall_age, na.rm = T),
                      mean_precision_age = colMeans(precision_age, na.rm = T),
                      mean_auc_age =colMeans(auc_age, na.rm = T),
                      se_recall_age = apply(recall_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_precision_age = apply(precision_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_auc_age = apply(auc_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      mean_recall_lda = colMeans(recall_lda, na.rm = T),
                      mean_precision_lda = colMeans(precision_lda, na.rm = T),
                      mean_auc_lda = colMeans(auc_lda, na.rm = T),
                      se_recall_lda = apply(recall_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_precision_lda = apply(precision_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_auc_lda = apply(auc_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))) )
df_10y_diff <- plot_df
plt_recall <- ggplot(plot_df,aes(x = mix_ratio)) +
  geom_line(aes(y = mean_recall_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_recall_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_recall_age, ymin=mean_recall_age-1.96*se_recall_age, ymax = mean_recall_age + 1.96*se_recall_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_recall_lda, ymin=mean_recall_lda-1.96*se_recall_lda, ymax = mean_recall_lda + 1.96*se_recall_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Subtype sample size ratio", y = "Recall:  identify the minor subtype") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,0.5)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","recall_10y_simulation.png"), plt_recall, width = 4, height = 4)

plt_precision <- ggplot(plot_df,aes(x = mix_ratio)) +
  geom_line(aes(y = mean_precision_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_precision_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_precision_age, ymin=mean_precision_age-1.96*se_precision_age, ymax = mean_precision_age + 1.96*se_precision_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_precision_lda, ymin=mean_precision_lda-1.96*se_precision_lda, ymax = mean_precision_lda + 1.96*se_precision_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Subtype sample size ratio", y = "Precision: identify the minor subtype") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,0.5)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","precision_10y_simulation.png"), plt_precision, width = 4, height = 4)

plt_auc <- ggplot(plot_df,aes(x = mix_ratio)) +
  geom_line(aes(y = mean_auc_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_auc_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_auc_age, ymin=mean_auc_age-1.96*se_auc_age, ymax = mean_auc_age + 1.96*se_auc_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_auc_lda, ymin=mean_auc_lda-1.96*se_auc_lda, ymax = mean_auc_lda + 1.96*se_auc_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Subtype sample size ratio", y = "AUPRC") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,0.5)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","auc_10y_simulation.png"), plt_auc, width = 4, height = 4)

# second disease
plot_second_ds_10y <- data.frame(mix_ratio = c( 0.1 * (0:9)/( 0.1 * (0:9) + 1) ), # , 1/( 0.1 * (9:0) + 1) ),
                                 mean_auc_age =colMeans(auc_age_rev, na.rm = T),
                                 se_auc_age = apply(auc_age_rev, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                                 mean_auc_lda = colMeans(auc_lda, na.rm = T),
                                 se_auc_lda = apply(auc_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))) ) %>%
  filter(mix_ratio != 0)
plt_auc <- ggplot(plot_second_ds_10y,aes(x = mix_ratio)) +
  geom_line(aes(y = mean_auc_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_auc_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_auc_age, ymin=mean_auc_age-1.96*se_auc_age, ymax = mean_auc_age + 1.96*se_auc_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_auc_lda, ymin=mean_auc_lda-1.96*se_auc_lda, ymax = mean_auc_lda + 1.96*se_auc_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Subtype sample size ratio", y = "AUPRC") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,0.5)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","auc_10y_2ndDs.png"), plt_auc, width = 4, height = 4)

# lower panel where there are small age difference
recall_age <- NULL
precision_age <- NULL
auc_age <- NULL
# save the results where the late onset group is larger
recall_age_rev <- NULL
precision_age_rev <- NULL
auc_age_rev <- NULL

recall_lda <- NULL
precision_lda <- NULL
auc_lda <- NULL
# save the results where the late onset group is larger
recall_lda_rev <- NULL
precision_lda_rev <- NULL
auc_lda_rev <- NULL

pt <- paste0("^subtype_tests_simulation_type6_rep*")
temp <- list.files(paste("simulation_results", sep=""), pattern=pt)
for(repid in 1:length(temp)){
  load(file = paste0("simulation_results/", temp[repid]))
  
  recall_age <- rbind(recall_age, subtype_tests[[2]][2,])
  precision_age <- rbind(precision_age, subtype_tests[[3]][2,])
  auc_age <- rbind(auc_age, subtype_tests[[4]][2,])
  recall_age_rev <- rbind(recall_age_rev, subtype_tests[[2]][3,])
  precision_age_rev <- rbind(precision_age_rev, subtype_tests[[3]][3,])
  auc_age_rev <- rbind(auc_age_rev, subtype_tests[[4]][3,])
  
  recall_lda <- rbind(recall_lda, subtype_tests[[6]][2,])
  precision_lda <- rbind(precision_lda, subtype_tests[[7]][2,])
  auc_lda <- rbind(auc_lda, subtype_tests[[8]][2,])
  recall_lda_rev <- rbind(recall_lda_rev, subtype_tests[[6]][3,])
  precision_lda_rev <- rbind(precision_lda_rev, subtype_tests[[7]][3,])
  auc_lda_rev <- rbind(auc_lda_rev, subtype_tests[[8]][3,])
}
# recall_age <- cbind(recall_age, recall_age_rev[, ncol(recall_age_rev):1])
# precision_age <- cbind(precision_age, precision_age_rev[, ncol(precision_age_rev):1])
# auc_age <- cbind(auc_age, auc_age_rev[, ncol(auc_age_rev):1])
# 
# recall_lda <- cbind(recall_lda, recall_lda_rev[, ncol(recall_lda_rev):1])
# precision_lda <- cbind(precision_lda, precision_lda_rev[, ncol(precision_lda_rev):1])
# auc_lda <- cbind(auc_lda, auc_lda_rev[, ncol(auc_lda_rev):1])

plot_df <- data.frame(mix_ratio = c( 0.1 * (0:9)/( 0.1 * (0:9) + 1)), # 1/( 0.1 * (9:0) + 1) ),
                      mean_recall_age = colMeans(recall_age, na.rm = T),
                      mean_precision_age = colMeans(precision_age, na.rm = T),
                      mean_auc_age =colMeans(auc_age, na.rm = T),
                      se_recall_age = apply(recall_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_precision_age = apply(precision_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_auc_age = apply(auc_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      mean_recall_lda = colMeans(recall_lda, na.rm = T),
                      mean_precision_lda = colMeans(precision_lda, na.rm = T),
                      mean_auc_lda = colMeans(auc_lda, na.rm = T),
                      se_recall_lda = apply(recall_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_precision_lda = apply(precision_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_auc_lda = apply(auc_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))) )
df_5y_diff <- plot_df
plt_recall <- ggplot(plot_df,aes(x = mix_ratio)) +
  geom_line(aes(y = mean_recall_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_recall_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_recall_age, ymin=mean_recall_age-1.96*se_recall_age, ymax = mean_recall_age + 1.96*se_recall_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_recall_lda, ymin=mean_recall_lda-1.96*se_recall_lda, ymax = mean_recall_lda + 1.96*se_recall_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Subtype sample size ratio", y = "Recall:  identify the minor subtype") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,0.5)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","recall_5y_simulation.png"), plt_recall, width = 4, height = 4)

plt_precision <- ggplot(plot_df,aes(x = mix_ratio)) +
  geom_line(aes(y = mean_precision_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_precision_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_precision_age, ymin=mean_precision_age-1.96*se_precision_age, ymax = mean_precision_age + 1.96*se_precision_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_precision_lda, ymin=mean_precision_lda-1.96*se_precision_lda, ymax = mean_precision_lda + 1.96*se_precision_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Subtype sample size ratio", y = "Precision: identify the minor subtype") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,0.5)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","precision_5y_simulation.png"), plt_precision, width = 4, height = 4)

plt_auc <- ggplot(plot_df,aes(x = mix_ratio)) +
  geom_line(aes(y = mean_auc_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_auc_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_auc_age, ymin=mean_auc_age-1.96*se_auc_age, ymax = mean_auc_age + 1.96*se_auc_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_auc_lda, ymin=mean_auc_lda-1.96*se_auc_lda, ymax = mean_auc_lda + 1.96*se_auc_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Subtype sample size ratio", y = "AUPRC") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,0.5)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","auc_5y_simulation.png"), plt_auc, width = 4, height = 4)
# second disease
plot_second_ds_5y <- data.frame(mix_ratio = c( 0.1 * (0:9)/( 0.1 * (0:9) + 1) ), # , 1/( 0.1 * (9:0) + 1) ),
                                 mean_auc_age =colMeans(auc_age_rev, na.rm = T),
                                 se_auc_age = apply(auc_age_rev, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                                 mean_auc_lda = colMeans(auc_lda, na.rm = T),
                                 se_auc_lda = apply(auc_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))) ) %>%
  filter(mix_ratio != 0)
plt_auc <- ggplot(plot_second_ds_5y,aes(x = mix_ratio)) +
  geom_line(aes(y = mean_auc_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_auc_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_auc_age, ymin=mean_auc_age-1.96*se_auc_age, ymax = mean_auc_age + 1.96*se_auc_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_auc_lda, ymin=mean_auc_lda-1.96*se_auc_lda, ymax = mean_auc_lda + 1.96*se_auc_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Subtype sample size ratio", y = "AUPRC") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,0.5)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","auc_5y_2ndDs.png"), plt_auc, width = 4, height = 4)

######################################################
# save Supplementary Table for main figure 2
######################################################
df_5y_diff <- df_5y_diff %>% 
  mutate(diff5 = paste0(scales::percent(mean_auc_age, accuracy = 0.01), " (SE = ", scales::percent(se_auc_age, accuracy = 0.01), ")")) %>%
  mutate(lda_result = paste0(scales::percent(mean_auc_lda, accuracy = 0.01), " (SE = ", scales::percent(se_auc_lda, accuracy = 0.01), ")")) %>%
  select(mix_ratio, diff5, lda_result) 

df_10y_diff <- df_10y_diff %>% 
  mutate(diff10= paste0(scales::percent(mean_auc_age, accuracy = 0.01), " (SE = ", scales::percent(se_auc_age, accuracy = 0.01), ")")) %>%
  select(mix_ratio, diff10) 

df_20y_diff <- df_20y_diff %>% 
  mutate(diff20= paste0(scales::percent(mean_auc_age, accuracy = 0.01), " (SE = ", scales::percent(se_auc_age, accuracy = 0.01), ")")) %>%
  select(mix_ratio, diff20) 

df_save <- df_5y_diff %>%
  left_join(df_10y_diff, by = "mix_ratio") %>%
  left_join(df_20y_diff, by = "mix_ratio") %>%
  filter(mix_ratio != 0) %>%
  mutate(mix_ratio = scales::percent(mix_ratio, accuracy = 0.01)) %>%
  select(mix_ratio, diff20, diff10, diff5, lda_result) 

names(df_save) <- c("Sample size ratio", "AUPRC: 20 year difference (ATM)", 
                    "AUPRC: 10 year difference (ATM)", "AUPRC: 5 year difference (ATM)",
                    "AUPRC: 20/10/5 year difference (LDA)")
write.csv(df_save, "~/Desktop/comorbidity/paper_writing/supplementary_files/Figure2DataTable.csv", row.names = F)

#######################################################
# subtype case 2: plot the subtype  AUPRC
# population size
#######################################################
recall_age <- NULL
precision_age <- NULL
auc_age <- NULL

recall_lda <- NULL
precision_lda <- NULL
auc_lda <- NULL
pt <- paste0("^subtype_tests_simulation_type2_rep*")
temp <- list.files(paste("simulation_results", sep=""), pattern=pt)
for(repid in 1:length(temp)){
  load(file = paste0("simulation_results/", temp[repid]))
  
  recall_age <- rbind(recall_age, subtype_tests[[2]][2,])
  precision_age <- rbind(precision_age, subtype_tests[[3]][2,])
  auc_age <- rbind(auc_age, subtype_tests[[4]][2,])
  
  recall_lda <- rbind(recall_lda, subtype_tests[[6]][2,])
  precision_lda <- rbind(precision_lda, subtype_tests[[7]][2,])
  auc_lda <- rbind(auc_lda, subtype_tests[[8]][2,])
}

plot_df <- data.frame(pop_sz = c(200, 500, 1000, 2000, 5000, 10000),
                      mean_recall_age = colMeans(recall_age, na.rm = T),
                      mean_precision_age = colMeans(precision_age, na.rm = T),
                      mean_auc_age =colMeans(auc_age, na.rm = T),
                      se_recall_age = apply(recall_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_precision_age = apply(precision_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_auc_age = apply(auc_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      mean_recall_lda = colMeans(recall_lda, na.rm = T),
                      mean_precision_lda = colMeans(precision_lda, na.rm = T),
                      mean_auc_lda = colMeans(auc_lda, na.rm = T),
                      se_recall_lda = apply(recall_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_precision_lda = apply(precision_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_auc_lda = apply(auc_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))) )

plt_auc <- ggplot(plot_df,aes(x = pop_sz)) +
  geom_line(aes(y = mean_auc_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_auc_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_auc_age, ymin=mean_auc_age-1.96*se_auc_age, ymax = mean_auc_age + 1.96*se_auc_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_auc_lda, ymin=mean_auc_lda-1.96*se_auc_lda, ymax = mean_auc_lda + 1.96*se_auc_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Population size", y = "AUPRC") + 
  # scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","auc_population_sz.png"), plt_auc, width = 4, height = 4)

#######################################################
# subtype case 3: plot the subtype  AUPRC
# number of distinct diseases
#######################################################
recall_age <- NULL
precision_age <- NULL
auc_age <- NULL

recall_lda <- NULL
precision_lda <- NULL
auc_lda <- NULL
pt <- paste0("^subtype_tests_simulation_type3_rep*")
temp <- list.files(paste("simulation_results", sep=""), pattern=pt)
for(repid in 1:length(temp)){
  load(file = paste0("simulation_results/", temp[repid]))
  
  recall_age <- rbind(recall_age, subtype_tests[[2]][2,])
  precision_age <- rbind(precision_age, subtype_tests[[3]][2,])
  auc_age <- rbind(auc_age, subtype_tests[[4]][2,])
  
  recall_lda <- rbind(recall_lda, subtype_tests[[6]][2,])
  precision_lda <- rbind(precision_lda, subtype_tests[[7]][2,])
  auc_lda <- rbind(auc_lda, subtype_tests[[8]][2,])
}

plot_df <- data.frame(disease_number = 10 * (2:6),
                      mean_recall_age = colMeans(recall_age, na.rm = T),
                      mean_precision_age = colMeans(precision_age, na.rm = T),
                      mean_auc_age =colMeans(auc_age, na.rm = T),
                      se_recall_age = apply(recall_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_precision_age = apply(precision_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_auc_age = apply(auc_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      mean_recall_lda = colMeans(recall_lda, na.rm = T),
                      mean_precision_lda = colMeans(precision_lda, na.rm = T),
                      mean_auc_lda = colMeans(auc_lda, na.rm = T),
                      se_recall_lda = apply(recall_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_precision_lda = apply(precision_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_auc_lda = apply(auc_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))) )

plt_auc <- ggplot(plot_df,aes(x = disease_number)) +
  geom_line(aes(y = mean_auc_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_auc_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_auc_age, ymin=mean_auc_age-1.96*se_auc_age, ymax = mean_auc_age + 1.96*se_auc_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_auc_lda, ymin=mean_auc_lda-1.96*se_auc_lda, ymax = mean_auc_lda + 1.96*se_auc_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Number of diseases", y = "AUPRC") + 
  # scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","auc_disease_number.png"), plt_auc, width = 4, height = 4)

#######################################################
# subtype case 4: plot the subtype  AUPRC
# number of incidences for each individual
#######################################################
recall_age <- NULL
precision_age <- NULL
auc_age <- NULL

recall_lda <- NULL
precision_lda <- NULL
auc_lda <- NULL
pt <- paste0("^subtype_tests_simulation_type4_rep*")
temp <- list.files(paste("simulation_results", sep=""), pattern=pt)
for(repid in 1:length(temp)){
  load(file = paste0("simulation_results/", temp[repid]))
  
  recall_age <- rbind(recall_age, subtype_tests[[2]][2,])
  precision_age <- rbind(precision_age, subtype_tests[[3]][2,])
  auc_age <- rbind(auc_age, subtype_tests[[4]][2,])
  
  recall_lda <- rbind(recall_lda, subtype_tests[[6]][2,])
  precision_lda <- rbind(precision_lda, subtype_tests[[7]][2,])
  auc_lda <- rbind(auc_lda, subtype_tests[[8]][2,])
}

plot_df <- data.frame(disease_per_indiv = c(2,4,6,8,10),
                      mean_recall_age = colMeans(recall_age, na.rm = T),
                      mean_precision_age = colMeans(precision_age, na.rm = T),
                      mean_auc_age =colMeans(auc_age, na.rm = T),
                      se_recall_age = apply(recall_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_precision_age = apply(precision_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_auc_age = apply(auc_age, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      mean_recall_lda = colMeans(recall_lda, na.rm = T),
                      mean_precision_lda = colMeans(precision_lda, na.rm = T),
                      mean_auc_lda = colMeans(auc_lda, na.rm = T),
                      se_recall_lda = apply(recall_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_precision_lda = apply(precision_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))),
                      se_auc_lda = apply(auc_lda, 2, function(x) sd(x, na.rm = T)/sqrt(length(x))) )

plt_auc <- ggplot(plot_df,aes(x = disease_per_indiv)) +
  geom_line(aes(y = mean_auc_age), color = red, linetype = "dashed") + 
  geom_line(aes(y =  mean_auc_lda), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = mean_auc_age, ymin=mean_auc_age-1.96*se_auc_age, ymax = mean_auc_age + 1.96*se_auc_age, color = "Age dependent modelling"), size = .3) + 
  geom_pointrange(aes(y = mean_auc_lda, ymin=mean_auc_lda-1.96*se_auc_lda, ymax = mean_auc_lda + 1.96*se_auc_lda, color = "Latent Dirichlet Allocation"),size = .3) + 
  scale_color_manual(values = c("Age dependent modelling" = red, "Latent Dirichlet Allocation" = green)) + 
  labs(x ="Average number of disease per individual", y = "AUPRC") + 
  # scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","auc_ds_per_idv.png"), plt_auc, width = 4, height = 4)



#####################################################
# topic loading accuracy 
#####################################################
v_id <- 1
rep_number <- 20
df_inf_accuracy <- data.frame(sample_sz = as.numeric(), weight_cor = as.numeric(), loading_sim = as.numeric())
for(s_id in 1:5){
  print(s_id)
  for(rep_id in 1:rep_number){
    load(paste0("~/Desktop/comorbidity/Multi-morbidity_biobank/simulation_results/nongen_simulation_v",v_id,"s",s_id,"rep",rep_id, ".RData"))
    
    truth <- nongenetic_tests[[1]]
    ageLDA <- nongenetic_tests[[2]]
    basicLDA <- nongenetic_tests[[3]]
    pcaResults<- nongenetic_tests[[4]]
    para_age <- nongenetic_tests[[5]]
    
    para_inf_age <- nongenetic_tests[[5]]
    para_sim <- nongenetic_tests[[7]]
    para_inf_lda <- nongenetic_tests[[8]]
    
    true_loading <- para_sim$true_beta[30:81,,]
    true_weights <- para_sim$theta
    
    infer_loading <- para_inf_age$pi_beta_basis[30:81,,]
    infer_weights <- sweep((para_inf_age$alpha_z ), 1, rowSums(para_inf_age$alpha_z), FUN="/") 
    
    # order the topics based on the correlation matching
    correlation <- cor(true_weights, infer_weights)
    order_topic <- c()
    for(i in 1:dim(true_weights)[2]){
      if(length(order_topic) > 0){
        # it is a bit complicated but I need to exclude previous chosen topic only matching the remaining
        relative_pos <- which.max(cor(true_weights[,i], infer_weights[, -order_topic]))
        order_topic <- c(order_topic, c(1:dim(true_weights)[2])[-order_topic][relative_pos])
      }else{
        order_topic <- c(order_topic, which.max(cor(true_weights[,i], infer_weights)))
      }
      
    }
    if(length(unique(order_topic)) == para_sim$K){
      infer_loading <- infer_loading[,,order_topic]
      infer_weights <- infer_weights[,order_topic]
      
      weight_pearson_r <- diag(cor(true_weights, infer_weights))
      
      cosine.sim <- function(A,B) sum(A*B)/sqrt(sum(A^2)*sum(B^2))
      loading_cos <- sapply(1:length(order_topic), function(x) cosine.sim(true_loading[,,x], infer_loading[,,x]) ) 
      
      df_inf_accuracy <- df_inf_accuracy %>%
        add_row(sample_sz = c(500, 2000, 5000,  10000, 20000)[s_id], weight_cor = mean(weight_pearson_r), loading_sim = mean(loading_cos))
      
    }else{
      print("poor fitting")
    }
  }
}

df_plot <- df_inf_accuracy %>%
  group_by(sample_sz) %>%
  summarise(mean_weight_cor = mean(weight_cor), se_weight_cor = sd(weight_cor)/sqrt(n()), 
            mean_loading_sim = mean(loading_sim), se_loading_sim = sd(loading_sim)/sqrt(n()))
plt_topic_weight_accuracy <- ggplot(df_plot,aes(x = sample_sz)) +
  geom_line(aes(y = mean_weight_cor), color = red, linetype = "dashed") + 
  geom_pointrange(aes(y = mean_weight_cor, ymin = mean_weight_cor - 1.96 * se_weight_cor, ymax = mean_weight_cor + 1.96 * se_weight_cor), color = red, size = 1) + 
  labs(x ="Simulated sample size", y = "Pearson correlation",) + 
  scale_y_continuous(labels = scales::percent, limits = c(0.4,0.8)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","topic_weight_sample_sz.png"), plt_topic_weight_accuracy, width = 4, height = 4)

plt_loading_accuracy <- ggplot(df_plot,aes(x = sample_sz)) +
  geom_line(aes(y = mean_loading_sim), color = red, linetype = "dashed") + 
  geom_pointrange(aes(y = mean_loading_sim, ymin = mean_loading_sim - 1.96 * se_loading_sim, ymax = mean_loading_sim + 1.96 * se_loading_sim), color = red, size = 1) + 
  labs(x ="Simulated sample size", y = "Cosine similarity",) + 
  scale_y_continuous(labels = scales::percent, limits = c(0.6,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","topic_loading_sample_sz.png"), plt_loading_accuracy, width = 4, height = 4)

# varying number of disease per individual 
v_id <- 3
rep_number <- 20
df_inf_accuracy <- data.frame(rec_per_indiv = as.numeric(), weight_cor = as.numeric(), loading_sim = as.numeric())
for(s_id in 1:5){
  print(s_id)
  for(rep_id in 1:rep_number){
    load(paste0("~/Desktop/comorbidity/Multi-morbidity_biobank/simulation_results/nongen_simulation_v",v_id,"s",s_id,"rep",rep_id, ".RData"))
    
    truth <- nongenetic_tests[[1]]
    ageLDA <- nongenetic_tests[[2]]
    basicLDA <- nongenetic_tests[[3]]
    pcaResults<- nongenetic_tests[[4]]
    para_age <- nongenetic_tests[[5]]
    
    para_inf_age <- nongenetic_tests[[5]]
    para_sim <- nongenetic_tests[[7]]
    para_inf_lda <- nongenetic_tests[[8]]
    
    true_loading <- para_sim$true_beta[30:81,,]
    true_weights <- para_sim$theta
    
    infer_loading <- para_inf_age$pi_beta_basis[30:81,,]
    infer_weights <- sweep((para_inf_age$alpha_z ), 1, rowSums(para_inf_age$alpha_z), FUN="/") 
    
    # order_topic <- apply(cor(true_weights, infer_weights), 1, which.max) # this one is matching with replacement
    # order the topics based on the correlation matching
    correlation <- cor(true_weights, infer_weights)
    order_topic <- c()
    for(i in 1:dim(true_weights)[2]){
      if(length(order_topic) > 0){
        # it is a bit complicated but I need to exclude previous chosen topic only matching the remaining
        relative_pos <- which.max(cor(true_weights[,i], infer_weights[, -order_topic]))
        order_topic <- c(order_topic, c(1:dim(true_weights)[2])[-order_topic][relative_pos])
      }else{
        order_topic <- c(order_topic, which.max(cor(true_weights[,i], infer_weights)))
        }
      
    }
    
    if(length(unique(order_topic)) == para_sim$K){
      infer_loading <- infer_loading[,,order_topic]
      infer_weights <- infer_weights[,order_topic]
      
      weight_pearson_r <- diag(cor(true_weights, infer_weights))
      
      cosine.sim <- function(A,B) sum(A*B)/sqrt(sum(A^2)*sum(B^2))
      loading_cos <- sapply(1:length(order_topic), function(x) cosine.sim(true_loading[,,x], infer_loading[,,x]) ) 
      
      df_inf_accuracy <- df_inf_accuracy %>%
        add_row(rec_per_indiv = c(2,4,6,8,10)[s_id], weight_cor = mean(weight_pearson_r), loading_sim = mean(loading_cos))
      
    }else{
      print("poor fitting")
    }
  }
}

df_plot <- df_inf_accuracy %>%
  group_by(rec_per_indiv) %>%
  summarise(mean_weight_cor = mean(weight_cor), se_weight_cor = sd(weight_cor)/sqrt(n()), 
            mean_loading_sim = mean(loading_sim), se_loading_sim = sd(loading_sim)/sqrt(n()))
plt_topic_weight_accuracy <- ggplot(df_plot,aes(x = rec_per_indiv)) +
  geom_line(aes(y = mean_weight_cor), color = red, linetype = "dashed") + 
  geom_pointrange(aes(y = mean_weight_cor, ymin = mean_weight_cor - 1.96 * se_weight_cor, ymax = mean_weight_cor + 1.96 * se_weight_cor), color = red, size = 1) + 
  labs(x ="Disease per individual", y = "Pearson correlation",) + 
  scale_y_continuous(labels = scales::percent, limits = c(0.4,0.8)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","topic_weight_rec_per_indiv.png"), plt_topic_weight_accuracy, width = 4, height = 4)

plt_loading_accuracy <- ggplot(df_plot,aes(x = rec_per_indiv)) +
  geom_line(aes(y = mean_loading_sim), color = red, linetype = "dashed") + 
  geom_pointrange(aes(y = mean_loading_sim, ymin = mean_loading_sim - 1.96 * se_loading_sim, ymax = mean_loading_sim + 1.96 * se_loading_sim), color = red, size = 1) + 
  labs(x ="Disease per individuale", y = "Cosine similarity",) + 
  scale_y_continuous(labels = scales::percent, limits = c(0.6,1)) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","topic_loading_rec_per_indiv.png"), plt_loading_accuracy, width = 4, height = 4)




#######################################################
# grouping accuracy
#######################################################

# grouping accuracy: similar between PCA, lda, ageLDA
v_id <- 2
s_id <- 2
load(paste0("~/Desktop/comorbidity/Multi-morbidity_biobank/simulation_results/nongen_simulation_v",v_id,"s",s_id,".RData"))

truth <- nongenetic_tests[[1]]
ageLDA <- nongenetic_tests[[2]]
basicLDA <- nongenetic_tests[[3]]
pcaResults<- nongenetic_tests[[4]]
para <- nongenetic_tests[[5]]
non_noise <- 1:sum(truth != 0)
grouping <- expand.grid(non_noise,non_noise)
grouping <- grouping %>%
  mutate(truth_same_group = (truth[Var1] == truth[Var2]),
         ageLDA =  (ageLDA[Var1] == ageLDA[Var2]),
         basicLDA =  (basicLDA[Var1] == basicLDA[Var2]),
         pcaResults =  (pcaResults[Var1] == pcaResults[Var2]))
mean(grouping$truth_same_group == grouping$ageLDA)
mean(grouping$truth_same_group == grouping$basicLDA)
mean(grouping$truth_same_group == grouping$pcaResults)

longData<-melt(para$true_beta[30:81,,1]) 
ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Diseases", y="Diseases", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

longData<-melt(para$true_beta[30:81,,2]) 
ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Diseases", y="Diseases", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

# plot prediction odds ratio for different topic number
v_id <- 2
rep_number <- 20
df_predict_OR <- data_frame(num_topic = s_id,OR_top5 = as.numeric(), ELBO = as.numeric(), model = as.character())
for(rep_id in 1:rep_number){
  for(s_id in 1:14){
  try({
    load(paste0("~/Desktop/comorbidity/Multi-morbidity_biobank/simulation_results/nongen_simulation_v",v_id,"s",s_id,"rep",rep_id, ".RData"))
    prediction_onebyone_rslt <- nongenetic_tests[[6]]
    para <- nongenetic_tests[[5]]
    df_predict_OR  <- df_predict_OR %>% 
      add_row(num_topic = (2:14)[s_id],
              OR_top5 = (prediction_onebyone_rslt[[1]][[4]]/(1 - prediction_onebyone_rslt[[1]][[4]]))/(0.05), 
              ELBO = para$lb[[2]][dim(para$lb)[1]], 
              model = "ATM"
      )
  })
  }
}
df_plot <- df_predict_OR %>%
  group_by(num_topic) %>%
  summarise(mean_predict_OR = mean(OR_top5), se_predict_OR = sd(OR_top5)/sqrt(n()), 
            mean_ELBO = mean(ELBO), se_ELBO = sd(ELBO)/sqrt(n()))
plt_topic_number <- ggplot(df_plot,aes(x = num_topic)) +
  geom_line(aes(y = mean_predict_OR), color = red, linetype = "dashed") + 
  geom_pointrange(aes(y = mean_predict_OR, ymin = mean_predict_OR - 1.96 * se_predict_OR, ymax = mean_predict_OR + 1.96 * se_predict_OR), color = red, size = 1) + 
  labs(x ="Number of topics used in model", y = "Prediction Odds Ratio",) + 
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","predictOR_simulation.png"), plt_topic_number, width = 4, height = 4)

# plot lower bound for differnt topic number
plt_topic_number <- ggplot(df_plot,aes(x = num_topic)) +
  geom_line(aes(y = mean_ELBO), color = red, linetype = "dashed") + 
  geom_pointrange(aes(y = mean_ELBO, ymin = mean_ELBO - 1.96 * se_ELBO, ymax = mean_ELBO + 1.96 * se_ELBO), color = red, size = 1) + 
  labs(x ="Number of topics", y = "Prediction Odds Ratio",) + 
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","topic_ELBO_simulation.png"), plt_topic_number, width = 4, height = 4)

longData<-melt(para$pi_beta_basis[31:80,,3])
longData<-longData[longData$value!=0,]
ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient2(low=blue, mid = "white", high=red) +
  labs(x="topics", y="Diseases", title=paste0("LDA")) +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

# grouping accuracy for sample size
v_id <- 1
rep_number <- 20
df_grouping_accuracy <- data_frame(sample_sz = as.numeric(), grouping_accuracy = as.numeric())
for(s_id in 1:4){
  for(rep_id in 1:rep_number){
    load(paste0("~/Desktop/comorbidity/Multi-morbidity_biobank/simulation_results/nongen_simulation_v",v_id,"s",s_id,"rep",rep_id, ".RData"))
    
    truth <- nongenetic_tests[[1]]
    ageLDA <- nongenetic_tests[[2]]
    basicLDA <- nongenetic_tests[[3]]
    pcaResults<- nongenetic_tests[[4]]
    para <- nongenetic_tests[[5]]
    non_noise <- 1:sum(truth != 0)
    grouping <- expand.grid(non_noise,non_noise)
    grouping <- grouping %>%
      mutate(truth_same_group = (truth[Var1] == truth[Var2]),
             ageLDA =  (ageLDA[Var1] == ageLDA[Var2]),
             basicLDA =  (basicLDA[Var1] == basicLDA[Var2]),
             pcaResults =  (pcaResults[Var1] == pcaResults[Var2]))
    df_grouping_accuracy  <- df_grouping_accuracy %>% 
      add_row(sample_sz = c(500, 2000, 5000,  10000)[s_id],
              grouping_accuracy = mean(grouping$truth_same_group == grouping$ageLDA))
  }
}
df_plot <- df_grouping_accuracy %>%
  group_by(sample_sz) %>%
  summarise(mean_grouping_accuracy = mean(grouping_accuracy), se_grouping_accuracy = sd(grouping_accuracy)/sqrt(n()))

plt_sample_sz <- ggplot(df_plot,aes(x = sample_sz)) +
  geom_line(aes(y = mean_grouping_accuracy), color = red, linetype = "dashed") + 
  geom_pointrange(aes(y = mean_grouping_accuracy, 
                      ymin = mean_grouping_accuracy - 1.96 * se_grouping_accuracy, 
                      ymax = mean_grouping_accuracy + 1.96 * se_grouping_accuracy), color = red, size = 1) + 
  scale_y_continuous(labels = scales::percent, limits = c(0.8,1)) +
  labs(x ="Simulated population size", y = "Grouping accuracy for disease pairs") + 
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","grouping_sim_pop_size.png"), plt_sample_sz, width = 4, height = 4)

# grouping accuracy for number of cases per-person 
v_id <- 3
rep_number <- 20
df_grouping_accuracy <- data_frame(incidence_per_indiv = as.numeric(), grouping_accuracy = as.numeric())
for(s_id in 1:4){
  for(rep_id in 1:rep_number){
    load(paste0("~/Desktop/comorbidity/Multi-morbidity_biobank/simulation_results/nongen_simulation_v",v_id,"s",s_id,"rep",rep_id, ".RData"))
    
    truth <- nongenetic_tests[[1]]
    ageLDA <- nongenetic_tests[[2]]
    basicLDA <- nongenetic_tests[[3]]
    pcaResults<- nongenetic_tests[[4]]
    para <- nongenetic_tests[[5]]
    non_noise <- 1:sum(truth != 0)
    grouping <- expand.grid(non_noise,non_noise)
    grouping <- grouping %>%
      mutate(truth_same_group = (truth[Var1] == truth[Var2]),
             ageLDA =  (ageLDA[Var1] == ageLDA[Var2]),
             basicLDA =  (basicLDA[Var1] == basicLDA[Var2]),
             pcaResults =  (pcaResults[Var1] == pcaResults[Var2]))
    df_grouping_accuracy  <- df_grouping_accuracy %>% 
      add_row(incidence_per_indiv = c(2,4,6,8)[s_id],
              grouping_accuracy = mean(grouping$truth_same_group == grouping$ageLDA))
  }
}

df_plot <- df_grouping_accuracy %>%
  group_by(incidence_per_indiv) %>%
  summarise(mean_grouping_accuracy = mean(grouping_accuracy), se_grouping_accuracy = sd(grouping_accuracy)/sqrt(n()))

plt_per_sample_rec <- ggplot(df_plot,aes(x = incidence_per_indiv)) +
  geom_line(aes(y = mean_grouping_accuracy), color = red, linetype = "dashed") + 
  geom_pointrange(aes(y = mean_grouping_accuracy, 
                      ymin = mean_grouping_accuracy - 1.96 * se_grouping_accuracy, 
                      ymax = mean_grouping_accuracy + 1.96 * se_grouping_accuracy), color = red, size = 1) + 
  scale_y_continuous(labels = scales::percent, limits = c(0.6,1)) +
  labs(x ="Simulated average incidence per individual", y = "Grouping accuracy for disease pairs") + 
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","grouping_sim_per_sample_rec.png"), plt_per_sample_rec, width = 4, height = 4)



######################################
# simulation with genetics 
######################################
# preparation: extracting per-SNP -> topic effect; per-SNP disease effect; Topic 
topic_gwas <- read.table(file = 'Topic_GWAS/topic1_loading_K10_rep10.assoc.linear',header = TRUE) 
topic_filter <- topic_gwas %>%
  filter(P < 5 * 10^(-8)) 
mean(topic_filter$BETA) 
topic_filter %>% arrange(desc(abs(BETA))) # use 0.01 as effect sizes
don <- topic_filter %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(topic_filter, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)
axisdf <- don %>% 
  group_by(CHR) %>% 
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
ggplot(don, aes(x=BPcum, y=-log10(P))) +
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
# for topic: use 20 causal SNP of 0.04 effect size ~ making the h2g of topic ~1-2%

# Disease SNP size
disease_gwas <- read.table(file = 'CC_gwas/272.11.assoc.logistic',header = TRUE) 
disease_filter <- disease_gwas %>%
  filter(P < 5 * 10^(-8)) %>%
  mutate(logOR = log(OR))
disease_filter %>% arrange(desc(logOR)) # use .15 as effect sizes
don <- disease_filter %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(disease_filter, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)
axisdf <- don %>% 
  group_by(CHR) %>% 
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
ggplot(don, aes(x=BPcum, y=-log10(P))) +
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
# for disease, use 20 SNP each of size 0.15 
md <- glm(disease_matrix[,(sim_para$disease_number + 1)] ~ genetics_population[,1:20], family = binomial)
null_md <- glm(disease_matrix[,(sim_para$disease_number + 1)] ~matrix(1, nrow = 10000, ncol = 20), family = binomial)
1-logLik(md)/logLik(null_md)
# have to use larger effect size as this effect is too small for h2g (maybe use 0.3)


# disease -> topic effect
rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
patient_loadings <- sweep((para$alpha_z - 1), 1, rowSums(para$alpha_z -1), FUN="/")
ds_data <- rec_data %>% 
  filter(diag_icd10 == 272.11)
regression_df <- data_frame(eid = para$eid, loadings = patient_loadings[,5]) %>%
  left_join(ds_data, by="eid") %>%
  mutate(ds = 1-is.na(diag_icd10))
regression_md <- lm(loadings ~ ds, data = regression_df)
summary(regression_md) # use 0.1 for the disease -> topic effect
regression_md <- glm(ds ~ loadings, family = binomial, data = regression_df)
summary(regression_md) # use 2 for the topic -> disease effect


# interaction effect size
disease_gwas <- read.table(file = 'CC_gwas/BMI_allcov_GxTopic_2.assoc.linear',header = TRUE) 
disease_filter <- disease_gwas %>%
  filter(P < 5 * 10^(-8), TEST == "ADDxCOV1") 
disease_filter %>% arrange(desc(BETA)) # use 1 as interaction sizes

##############################
# toy model with 1 SNP
#############################
N <- 10000
genetics <-  as.matrix(generate_genetics(0.3, N))
maptan <- function(x)  (atan(x)/(pi/2) + 1)/2
# 0. SNP -> T -> T
snp2T <- 2
t2d <- 2
theta <- rnorm(N) + genetics * snp2T 
disease <- rnorm(N) +  theta * t2d 
gxT <- lm(disease ~ theta * genetics)
summary(gxT)
gxD <- lm(theta ~ disease * genetics)
summary(gxD)


# 1. SNP->T->D
snp2T <- 0.2 # combined with norm correspond to ~0.05 regression effect
t2d <- 2
theta <- pnorm(rnorm(N) + genetics * snp2T - mean(genetics * snp2T)) 
# disease <- 1 * (runif(N) < 1/(1+exp(-(rnorm(N, mean = -1) +  theta * t2d ) ) ) )
# disease <- 1 * (runif(N) < maptan(rnorm(N, mean = -1) +  theta * t2d ) )
disease <- 1 * (runif(N) < pnorm(rnorm(N) +  theta * t2d ) )
p_gxT <- test_gxT(theta, disease, genetics)
p_gxD <- test_gxD(theta, disease, genetics)
gxT <- glm(disease ~ theta * genetics, family = binomial)
summary(gxT)
gxD <- lm(theta ~ disease * genetics)
summary(gxD)



# 2. SNP -> D -> T
snp2D <- 0.2
d2t <- 1
# disease <- 1 * (runif(N) < 1/(1+exp(-(rnorm(N) + genetics * snp2D))) )
# disease <- 1 * (runif(N) < maptan(rnorm(N) + genetics * snp2D)) 
# disease <- 1 * (runif(N) < pnorm(rnorm(N) + genetics * snp2D)) 
lia_ds <- rnorm(N) + genetics * snp2D
disease <- 1 * (quantile(lia_ds, 0.8) < lia_ds) 
theta<- pnorm(rnorm(N) + d2t * disease)
p_gxT <- test_gxT(theta, disease, genetics)
p_gxD <- test_gxD(theta, disease, genetics)
gxT <- glm(disease ~ theta * genetics, family = binomial)
summary(gxT)
gxD <- lm(theta ~ disease * genetics)
summary(gxD)



# 3. SNP * T -> D
t2d <- 2
theta <- runif(N) 
# disease <- 1 * (runif(N) < 1/(1+exp(-(rnorm(N, mean = -1) +  theta * t2d ) ) ) )
# disease <- 1 * (runif(N) < pnorm(rnorm(N) +  theta * genetics * t2d ) )
lia_ds <- rnorm(N) +  theta * genetics * t2d # + theta * t2d * 10
disease <- 1 * (quantile(lia_ds, 0.8) < lia_ds) 
p_gxT <- test_gxT(theta, disease, genetics)
p_gxD <- test_gxD(theta, disease, genetics)
gxT <- glm(disease ~ theta * genetics, family = binomial)
summary(gxT)
gxD <- lm(theta ~ disease * genetics)
summary(gxD)
###############################################
# GxTopic Power changes relative to SNP main effect 
###############################################
N <- 10000
t2d <- 2
snp2D <- 0.2
p_gxT <- matrix(nrow = 10, ncol = 100)
p_gxT_T2 <- matrix(nrow = 10, ncol = 100)
for(id in 1:10){
  scale <- ((1:10)/10)[id]
  for(rep in 1:100){
    theta <- runif(N) 
    genetics <-  as.matrix(generate_genetics(0.3, N))
    snp2D <- 0.2 * scale
    lia_ds <- rnorm(N) +  snp2D * t2d * theta * genetics + theta * t2d  + genetics * snp2D
    disease <- 1 * (quantile(lia_ds, 0.8) < lia_ds) 
    p_gxT[id, rep] <- test_gxT(theta, disease, genetics)
    p_gxT_T2[id, rep] <- test_gxT_T2(theta, disease, genetics)
  }
}
pw_thr <- 0.05
gxT_mean <- rowMeans(p_gxT < pw_thr,na.rm = T)
gxT_se <- sqrt(gxT_mean*(1-gxT_mean)/dim(p_gxT)[2])
# non-linear
gT2_mean <- rowMeans(p_gxT_T2 < pw_thr,na.rm = T)
gT2_se <- sqrt(gT2_mean*(1-gT2_mean)/dim(p_gxT_T2)[2])

df_1 <- data.frame(SNP_effect = t2d*snp2D*((1:10)/10),  gxT_mean= gxT_mean, gxT_se = gxT_se, gT2_mean=gT2_mean, gT2_se = gT2_se)
plt_gxT <- ggplot(df_1,aes(x = SNP_effect)) +
  geom_line(aes(y = gxT_mean), color = red, linetype = "dashed") + 
  geom_pointrange(aes(y = gxT_mean, ymin=gxT_mean-1.96*gxT_se, ymax = gxT_mean + 1.96*gxT_se, color = "G x Topic model"), size = .3) + 
  geom_line(aes(y = gT2_mean), color = blue, linetype = "dashed") + 
  geom_pointrange(aes(y = gT2_mean, ymin=gT2_mean-1.96*gT2_se, ymax = gT2_mean + 1.96*gT2_se, color = "Topic^2 + G x Topic model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "Topic^2 + G x Topic model" = blue)) + 
  labs(x ="Interaction effect", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","GxT_main_effect.png"), plt_gxT, width = 4, height = 4)

t2d <- 2
snp2D <- 0.2
p_gxT <- matrix(nrow = 10, ncol = 100)
p_gxT_T2 <- matrix(nrow = 10, ncol = 100)
for(id in 1:10){
  variant_number <- c(2 * (1:10))[id]
  for(rep in 1:(floor(100/variant_number) + 1)){
    theta <- runif(N) 
    genetics_multiple <-  sapply(rep(0.3,variant_number), function(x) generate_genetics(x, N))
    lia_ds <- rnorm(N) + snp2D * t2d * theta * genetics_multiple %*% rep(t2d, variant_number)
    disease <- 1 * (quantile(lia_ds, 0.8) < lia_ds) 
    # save all the SNP results to the 100 replicates
    idx_list <- pmin((variant_number*(rep-1) + 1), 100):pmin(variant_number*rep, 100)
    p_gxT[id, idx_list] <- 
      test_gxT(theta, disease, genetics_multiple)[1:length(idx_list)]
    p_gxT_T2[id, idx_list] <- 
      test_gxT_T2(theta, disease, genetics_multiple)[1:length(idx_list)]
  }
}
pw_thr <- 0.05 
gxT_mean <- rowMeans(p_gxT < pw_thr,na.rm = T)
gxT_se <- sqrt(gxT_mean*(1-gxT_mean)/dim(p_gxT)[2])
gT2_mean <- rowMeans(p_gxT_T2 < pw_thr,na.rm = T)
gT2_se <- sqrt(gT2_mean*(1-gT2_mean)/dim(p_gxT_T2)[2])
df_2 <- data.frame(variant_number = 2*(1:10),  gxT_mean= gxT_mean, gxT_se = gxT_se, gT2_mean=gT2_mean, gT2_se = gT2_se)
plt_gxT <- ggplot(df_2,aes(x = variant_number)) +
  geom_line(aes(y = gxT_mean), color = red, linetype = "dashed") + 
  geom_pointrange(aes(y = gxT_mean, ymin=gxT_mean-1.96*gxT_se, ymax = gxT_mean + 1.96*gxT_se, color = "G x Topic model"), size = .3) +   
  geom_line(aes(y = gT2_mean), color = blue, linetype = "dashed") + 
  geom_pointrange(aes(y = gT2_mean, ymin=gT2_mean-1.96*gT2_se, ymax = gT2_mean + 1.96*gT2_se, color = "Topic^2 + G x Topic model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "Topic^2 + G x Topic model" = blue)) + 
  labs(x ="Number of interacting variants", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","GxT_SNP_number.png"), plt_gxT, width = 4, height = 4)





##################################################
# additional simulation: multiple SNP testing using single or multiple model
##################################################
genetics_multiple <-  sapply(rep(0.3,20), function(x) generate_genetics(x, N))
lia_ds <- rnorm(N) + theta * genetics_multiple %*% rep(t2d, 20)
disease <- 1 * (quantile(lia_ds, 0.8) < lia_ds) 
p_gxT <- test_gxT(theta, disease, genetics_multiple)
gxT <- glm(disease ~ theta * genetics_multiple, family = binomial)
summary(gxT)

# 4. SNP -> T -> D & SNP -> D 
snp2T <- .2
snp2D <- .2
t2d <- 4
theta <- pnorm(rnorm(N) + genetics * snp2T) 
# should use a static threshold
# disease <- 1 * (runif(N) < pnorm(rnorm(N) + genetics * snp2D +  genetics^2 * snp2D +theta * t2d))
disease <- 1 * (0.8 < pnorm(rnorm(N) + genetics * snp2D +  genetics^2 * snp2D +theta * t2d)) 
p_gxT <- test_gxT(theta, disease, genetics)
p_gxD <- test_gxD(theta, disease, genetics)
gxT <- glm(disease ~ theta * genetics, family = binomial)
summary(gxT)
gxD <- lm(theta ~ disease * genetics)
summary(gxD)



############################################
# perform GxTopic analysis 
############################################
source("simulation_functions.R")

# 1. SNP -> T -> D
v_id <- 1
pw_thr <- 0.05 # threshold of power 
df_1 <- data.frame( causal_v2t = as.numeric(), gxT_mean= as.numeric(), gxT_se = as.numeric(),
                    gT2_mean= as.numeric(), gT2_se = as.numeric())
for(s_id in 1:10){
  print(s_id)
  load(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/genetic_simulation/", 
              "gen_simulation_v",v_id,"s", s_id,".RData" ))
  p_gxT <- NULL
  p_gT2 <- NULL
  p_gxD <- NULL
  for(rep_id in 1:length(genetic_simulation)){
    genetics_population <- genetic_simulation[[rep_id]][[1]]
    patient_loadings <- genetic_simulation[[rep_id]][[2]] 
    disease_matrix <- genetic_simulation[[rep_id]][[3]] 
    sim_para <- genetic_simulation[[rep_id]][[4]]
    # find which topic is genetically modulated
    tp_genetic <- colMeans(patient_loadings[which(disease_matrix[,1] == 1),]) %>% which.max
    loadings <- patient_loadings[,tp_genetic]
    for(ds_id in 1:((sim_para$disease_number-2)/2)){
    # for(ds_id in 1:1){
      p_gxT <- c(p_gxT, test_gxT(loadings, disease_matrix[,ds_id], genetics_population[,21:(20+sim_para$cont_v2t)]))
      p_gT2 <- c(p_gT2, test_gxT_T2(loadings, disease_matrix[,ds_id], genetics_population[,21:(20+sim_para$cont_v2t)]))
      # p_gxD <- c(p_gxD,test_gxD_simple(loadings, disease_matrix[,ds_id], genetics_population[,21:(20+sim_para$cont_v2t)]))
    }
  }
  gxT_mean <- mean(p_gxT < pw_thr,na.rm = T)
  gxT_se <- sqrt(gxT_mean*(1-gxT_mean)/length(p_gxT))
  gT2_mean <- mean(p_gT2 < pw_thr,na.rm = T)
  gT2_se <- sqrt(gT2_mean*(1-gT2_mean)/length(p_gT2))
  # gT2_mean <- mean(p_gxD < pw_thr,na.rm = T)
  # gT2_se <- sqrt(gT2_mean*(1-gxT_mean)/length(p_gxD))
  df_1 <- df_1 %>%
    add_row(causal_v2t = sim_para$cont_v2t,  gxT_mean= gxT_mean, gxT_se = gxT_se,
            gT2_mean= gT2_mean, gT2_se = gT2_se)
}
plt_g2t2d <- ggplot(df_1,aes(x = causal_v2t)) +
  geom_line(aes(y = gxT_mean), color = red, linetype = "dashed") + 
  geom_line(aes(y =  gT2_mean), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = gxT_mean, ymin=gxT_mean-1.96*gxT_se, ymax = gxT_mean + 1.96*gxT_se, color = "G x Topic model"), size = .3) + 
  geom_pointrange(aes(y = gT2_mean, ymin=gT2_mean-1.96*gT2_se, ymax = gT2_mean + 1.96*gT2_se, color = "G x Disease model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "G x Disease model" = green)) + 
  labs(x ="Number of causal SNP to topic", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
ggsave(paste0("../paper_writing/Production_figures/","interaction_g2t2d.png"), plt_g2t2d, width = 4, height = 4)


# 2. causal disease SNP->D->T
v_id <- 2
pw_thr <- 0.05 # threshold of power 
df_2 <- data.frame( disease2topic = as.numeric(), gxT_mean= as.numeric(), gxT_se = as.numeric(),
                    gT2_mean= as.numeric(), gT2_se = as.numeric())
for(s_id in 1:10){
  print(s_id)
  load(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/genetic_simulation/", 
              "gen_simulation_v",v_id,"s", s_id,".RData" ))
  p_gxT <- NULL
  p_gT2 <- NULL
  for(rep_id in 1:length(genetic_simulation)){
    genetics_population <- genetic_simulation[[rep_id]][[1]]
    patient_loadings <- genetic_simulation[[rep_id]][[2]] 
    disease_matrix <- genetic_simulation[[rep_id]][[3]] 
    sim_para <- genetic_simulation[[rep_id]][[4]]
    # find which topic is genetically modulated
    tp_genetic <- colMeans(patient_loadings[which(disease_matrix[,1] == 1),]) %>% which.max
    loadings <- patient_loadings[,tp_genetic]
    ds_id <- sim_para$disease_number + 1
    p_gxT <- c(p_gxT, test_gxT(loadings, disease_matrix[,ds_id], genetics_population[,1:20]))
    p_gT2 <- c(p_gT2,test_gxT_T2(loadings, disease_matrix[,ds_id], genetics_population[,1:20]))
  }
  gxT_mean <- mean(p_gxT < pw_thr,na.rm = T)
  gxT_se <- sqrt(gxT_mean*(1-gxT_mean)/length(p_gxT))
  gT2_mean <- mean(p_gT2 < pw_thr,na.rm = T)
  gT2_se <- sqrt(gT2_mean*(1-gxT_mean)/length(p_gT2))
  df_2 <- df_2 %>%
    add_row(disease2topic = sim_para$disease2topic,  gxT_mean= gxT_mean, gxT_se = gxT_se,
            gT2_mean= gT2_mean, gT2_se = gT2_se)
}
plt_g2d2t <- ggplot(df_2,aes(x = disease2topic)) +
  geom_line(aes(y = gxT_mean), color = red, linetype = "dashed") + 
  geom_line(aes(y =  gT2_mean), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = gxT_mean, ymin=gxT_mean-1.96*gxT_se, ymax = gxT_mean + 1.96*gxT_se, color = "G x Topic model"), size = .3) + 
  geom_pointrange(aes(y = gT2_mean, ymin=gT2_mean-1.96*gT2_se, ymax = gT2_mean + 1.96*gT2_se, color = "G x Disease model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "G x Disease model" = green)) + 
  labs(x ="Causal disease effect on topic", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
ggsave(paste0("../paper_writing/Production_figures/","interaction_g2d2t.png"), plt_g2d2t, width = 4, height = 4)

######################################################
# compute per-sd change in disease to topic

# 3. SNP*T->D
v_id <- 3
pw_thr <- 0.05 # threshold of power 
df_3 <- data.frame( itr_effect = as.numeric(), gxT_mean= as.numeric(), gxT_se = as.numeric(),
                    gT2_mean= as.numeric(), gT2_se = as.numeric())
for(s_id in 1:10){
  print(s_id)
  load(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/genetic_simulation/", 
              "gen_simulation_v",v_id,"s", s_id,".RData" ))
  p_gxT <- NULL
  p_gT2 <- NULL
  for(rep_id in 1:length(genetic_simulation)){
    genetics_population <- genetic_simulation[[rep_id]][[1]]
    patient_loadings <- genetic_simulation[[rep_id]][[2]] 
    disease_matrix <- genetic_simulation[[rep_id]][[3]] 
    sim_para <- genetic_simulation[[rep_id]][[4]]
    # find which topic is genetically modulated
    tp_genetic <- colMeans(patient_loadings[which(disease_matrix[,1] == 1),]) %>% which.max
    loadings <- patient_loadings[,tp_genetic]
    ds_id <- sim_para$disease_number + 2
    p_gxT <- c(p_gxT, test_gxT(loadings, disease_matrix[,ds_id], genetics_population[,41:60]))
    p_gT2 <- c(p_gT2,test_gxT_T2(loadings, disease_matrix[,ds_id], genetics_population[,41:60]))
    ############################
    # important observation
    # when simulating multiple interaction and testing for one -- the power is lost!
    ############################
    # summary(glm(disease_matrix[,ds_id] ~ genetics_population[,41:60]*loadings, family = binomial))
  }
  gxT_mean <- mean(p_gxT < pw_thr,na.rm = T)
  gxT_se <- sqrt(gxT_mean*(1-gxT_mean)/length(p_gxT))
  gT2_mean <- mean(p_gT2 < pw_thr,na.rm = T)
  gT2_se <- sqrt(gT2_mean*(1-gxT_mean)/length(p_gT2))
  df_3 <- df_3 %>%
    add_row(itr_effect = sim_para$itr_effect, gxT_mean= gxT_mean, gxT_se = gxT_se,
            gT2_mean= gT2_mean, gT2_se = gT2_se)
}
plt_gxt <- ggplot(df_3,aes(x = itr_effect)) +
  geom_line(aes(y = gxT_mean), color = red, linetype = "dashed") + 
  geom_line(aes(y =  gT2_mean), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = gxT_mean, ymin=gxT_mean-1.96*gxT_se, ymax = gxT_mean + 1.96*gxT_se, color = "G x Topic model"), size = .3) + 
  geom_pointrange(aes(y = gT2_mean, ymin=gT2_mean-1.96*gT2_se, ymax = gT2_mean + 1.96*gT2_se, color = "G x Disease model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "G x Disease model" = green)) + 
  labs(x ="Effect size of GxTopic interaction", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
ggsave(paste0("../paper_writing/Production_figures/","interaction_gxt.png"), plt_gxt, width = 4, height = 4)

# 4. SNP->D; SNP->T->D
v_id <- 4
pw_thr <- 0.05 # threshold of power 
df_4 <- data.frame( topic2disease = as.numeric(), gxT_mean= as.numeric(), gxT_se = as.numeric(),
                    gT2_mean= as.numeric(), gT2_se = as.numeric())
for(s_id in 1:10){
  print(s_id)
  load(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/genetic_simulation/", 
              "gen_simulation_v",v_id,"s", s_id,".RData" ))
  p_gxT <- NULL
  p_gT2 <- NULL
  for(rep_id in 1:length(genetic_simulation)){
    genetics_population <- genetic_simulation[[rep_id]][[1]]
    patient_loadings <- genetic_simulation[[rep_id]][[2]] 
    disease_matrix <- genetic_simulation[[rep_id]][[3]] 
    sim_para <- genetic_simulation[[rep_id]][[4]]
    # find which topic is genetically modulated
    tp_genetic <- colMeans(patient_loadings[which(disease_matrix[,1] == 1),]) %>% which.max
    loadings <- patient_loadings[,tp_genetic]
    ds_id <- sim_para$disease_number + 3
    p_gxT <- c(p_gxT, test_gxT(loadings, disease_matrix[,ds_id], genetics_population[,21:40]))
    p_gT2 <- c(p_gT2,test_gxT_T2(loadings, disease_matrix[,ds_id], genetics_population[,21:40]))
  }
  gxT_mean <- mean(p_gxT < pw_thr,na.rm = T)
  gxT_se <- sqrt(gxT_mean*(1-gxT_mean)/length(p_gxT))
  gT2_mean <- mean(p_gT2 < pw_thr,na.rm = T)
  gT2_se <- sqrt(gT2_mean*(1-gxT_mean)/length(p_gT2))
  df_4 <- df_4 %>%
    add_row(topic2disease = sim_para$topic2disease, gxT_mean= gxT_mean, gxT_se = gxT_se,
            gT2_mean= gT2_mean, gT2_se = gT2_se)
}
plt_pleiotropy <- ggplot(df_4,aes(x = topic2disease)) +
  geom_line(aes(y = gxT_mean), color = red, linetype = "dashed") + 
  geom_line(aes(y =  gT2_mean), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = gxT_mean, ymin=gxT_mean-1.96*gxT_se, ymax = gxT_mean + 1.96*gxT_se, color = "G x Topic model"), size = .3) + 
  geom_pointrange(aes(y = gT2_mean, ymin=gT2_mean-1.96*gT2_se, ymax = gT2_mean + 1.96*gT2_se, color = "G x Disease model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "G x Disease model" = green)) + 
  labs(x ="Effect of topic on disease", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
ggsave(paste0("../paper_writing/Production_figures/","pleitropy_gxt.png"), plt_pleiotropy, width = 4, height = 4)

# additional cases with non-linear genetic effect
v_id <- 4
pw_thr <- 0.05 # threshold of power 
df_5 <- data.frame( topic2disease = as.numeric(), gxT_mean= as.numeric(), gxT_se = as.numeric(),
                    gT2_mean= as.numeric(), gT2_se = as.numeric())
for(s_id in 1:10){
  print(s_id)
  load(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/genetic_simulation/", 
              "gen_simulation_v",v_id,"s", s_id,".RData" ))
  p_gxT <- NULL
  p_gT2 <- NULL
  for(rep_id in 1:length(genetic_simulation)){
    genetics_population <- genetic_simulation[[rep_id]][[1]]
    patient_loadings <- genetic_simulation[[rep_id]][[2]] 
    disease_matrix <- genetic_simulation[[rep_id]][[3]] 
    sim_para <- genetic_simulation[[rep_id]][[4]]
    # find which topic is genetically modulated
    tp_genetic <- colMeans(patient_loadings[which(disease_matrix[,1] == 1),]) %>% which.max
    loadings <- patient_loadings[,tp_genetic]
    ds_id <- sim_para$disease_number + 4
    p_gxT <- c(p_gxT, test_gxT(loadings, disease_matrix[,ds_id], genetics_population[,21:40]))
    p_gxD <- c(p_gxD,test_gxT_T2(loadings, disease_matrix[,ds_id], genetics_population[,21:40]))
  }
  gxT_mean <- mean(p_gxT < pw_thr,na.rm = T)
  gxT_se <- sqrt(gxT_mean*(1-gxT_mean)/length(p_gxT))
  gT2_mean <- mean(p_gxD < pw_thr,na.rm = T)
  gT2_se <- sqrt(gT2_mean*(1-gxT_mean)/length(p_gxD))
  df_5 <- df_5 %>%
    add_row(topic2disease = sim_para$topic2disease, gxT_mean= gxT_mean, gxT_se = gxT_se,
            gT2_mean= gT2_mean, gT2_se = gT2_se)
}
plt_nonlinrG <- ggplot(df_5,aes(x = topic2disease)) +
  geom_line(aes(y = gxT_mean), color = red, linetype = "dashed") + 
  geom_line(aes(y =  gT2_mean), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = gxT_mean, ymin=gxT_mean-1.96*gxT_se, ymax = gxT_mean + 1.96*gxT_se, color = "G x Topic model"), size = .3) + 
  geom_pointrange(aes(y = gT2_mean, ymin=gT2_mean-1.96*gT2_se, ymax = gT2_mean + 1.96*gT2_se, color = "G x Disease model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "G x Disease model" = green)) + 
  labs(x ="Effect of topic on disease", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
ggsave(paste0("../paper_writing/Production_figures/","interaction_nonlinear.png"), plt_nonlinrG, width = 4, height = 4)

# QQ plot
sort_logP <- sort(-log10(p_gxT))
sim_logP <- sort(-log10(runif(length(sort_logP))))
plot(sim_logP, sort_logP)

# D = T + T^2 non linear case
v_id <- 4
pw_thr <- 0.05 # threshold of power 
df_6 <- data.frame( topic2disease = as.numeric(), gxT_mean= as.numeric(), gxT_se = as.numeric(),
                    gT2_mean= as.numeric(), gT2_se = as.numeric())
for(s_id in 1:10){
  print(s_id)
  load(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/genetic_simulation/", 
              "gen_simulation_v",v_id,"s", s_id,".RData" ))
  p_gxT <- NULL
  p_gT2 <- NULL
  for(rep_id in 1:length(genetic_simulation)){
    genetics_population <- genetic_simulation[[rep_id]][[1]]
    patient_loadings <- genetic_simulation[[rep_id]][[2]] 
    disease_matrix <- genetic_simulation[[rep_id]][[3]] 
    sim_para <- genetic_simulation[[rep_id]][[4]]
    # find which topic is genetically modulated
    tp_genetic <- colMeans(patient_loadings[which(disease_matrix[,1] == 1),]) %>% which.max
    loadings <- patient_loadings[,tp_genetic]
    ds_id <- sim_para$disease_number + 5
    p_gxT <- c(p_gxT, test_gxT(loadings, disease_matrix[,ds_id], genetics_population[,21:40]))
    p_gT2 <- c(p_gT2,test_gxT_T2(loadings, disease_matrix[,ds_id], genetics_population[,21:40]))
  }
  gxT_mean <- mean(p_gxT < pw_thr,na.rm = T)
  gxT_se <- sqrt(gxT_mean*(1-gxT_mean)/length(p_gxT))
  gT2_mean <- mean(p_gT2 < pw_thr,na.rm = T)
  gT2_se <- sqrt(gT2_mean*(1-gxT_mean)/length(p_gT2))
  df_6 <- df_6 %>%
    add_row(topic2disease = sim_para$topic2disease, gxT_mean= gxT_mean, gxT_se = gxT_se,
            gT2_mean= gT2_mean, gT2_se = gT2_se)
}
plt_nonlinrG <- ggplot(df_6,aes(x = topic2disease)) +
  geom_line(aes(y = gxT_mean), color = red, linetype = "dashed") + 
  geom_line(aes(y =  gT2_mean), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = gxT_mean, ymin=gxT_mean-1.96*gxT_se, ymax = gxT_mean + 1.96*gxT_se, color = "G x Topic model"), size = .3) + 
  geom_pointrange(aes(y = gT2_mean, ymin=gT2_mean-1.96*gT2_se, ymax = gT2_mean + 1.96*gT2_se, color = "G x Disease model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "G x Disease model" = green)) + 
  labs(x ="Effect of topic on disease", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
ggsave(paste0("../paper_writing/Production_figures/","interaction_t2_nonlinear.png"), plt_nonlinrG, width = 4, height = 4)


############################################
# perform MR
############################################
source("simulation_functions.R")

# 1. SNP -> T -> D
v_id <- 1
pw_thr <- 0.05 # threshold of power 
df_1 <- data.frame( causal_v2t = as.numeric(), t2d_mean= as.numeric(), t2d_se = as.numeric(),
                    d2t_mean= as.numeric(), d2t_se = as.numeric())
for(s_id in 1:10){
  print(s_id)
  load(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/genetic_simulation/", 
              "gen_simulation_v",v_id,"s", s_id,".RData" ))
  p_t2d <- rep(1, length(genetic_simulation)*5)
  p_d2t <- rep(1, length(genetic_simulation)*5)
  for(rep_id in 1:length(genetic_simulation)){
    genetics_population <- genetic_simulation[[rep_id]][[1]]
    patient_loadings <- genetic_simulation[[rep_id]][[2]] 
    disease_matrix <- genetic_simulation[[rep_id]][[3]] 
    sim_para <- genetic_simulation[[rep_id]][[4]]
    # find which topic is genetically modulated
    tp_genetic <- colMeans(patient_loadings[which(disease_matrix[,1] == 1),]) %>% which.max
    loadings <- patient_loadings[,tp_genetic]
    for(ds_id in 1:5){
      MRresults <- eggerMR(loadings, disease_matrix[,ds_id], genetics_population, method = "ivw")
      try({
        p_t2d[(ds_id + (rep_id - 1)*5 )] <- MRresults[[1]]@Pvalue
        p_d2t[(ds_id + (rep_id - 1)*5 )] <- MRresults[[2]]@Pvalue
      })
    }
  }
  t2d_mean <- mean(p_t2d < pw_thr,na.rm = T)
  t2d_se <- sqrt(t2d_mean*(1-t2d_mean)/length(p_t2d))
  d2t_mean <- mean(p_d2t < pw_thr,na.rm = T)
  d2t_se <- sqrt(d2t_mean*(1-d2t_mean)/length(p_d2t))
  df_1 <- df_1 %>%
    add_row(causal_v2t = sim_para$cont_v2t,  t2d_mean= t2d_mean, t2d_se = t2d_se,
            d2t_mean= d2t_mean, d2t_se = d2t_se)
}
plt_g2t2d <- ggplot(df_1,aes(x = causal_v2t)) +
  geom_line(aes(y = t2d_mean), color = red, linetype = "dashed") + 
  geom_line(aes(y =  d2t_mean), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = t2d_mean, ymin=t2d_mean-1.96*t2d_se, ymax = t2d_mean + 1.96*t2d_se, color = "G x Topic model"), size = .3) + 
  geom_pointrange(aes(y = d2t_mean, ymin=d2t_mean-1.96*d2t_se, ymax = d2t_mean + 1.96*d2t_se, color = "G x Disease model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "G x Disease model" = green)) + 
  labs(x ="Number of causal SNP to topic", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","MR_g2t2d.png"), plt_g2t2d, width = 4, height = 4)


# 2. SNP -> T -> D
v_id <- 2
pw_thr <- 0.05 # threshold of power 
df_2 <- data.frame( disease2topic = as.numeric(), t2d_mean= as.numeric(), t2d_se = as.numeric(),
                    d2t_mean= as.numeric(), d2t_se = as.numeric())
for(s_id in 1:10){
  print(s_id)
  load(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/genetic_simulation/", 
              "gen_simulation_v",v_id,"s", s_id,".RData" ))
  p_t2d <- rep(1, length(genetic_simulation))
  p_d2t <- rep(1, length(genetic_simulation))
  for(rep_id in 1:length(genetic_simulation)){
    genetics_population <- genetic_simulation[[rep_id]][[1]]
    patient_loadings <- genetic_simulation[[rep_id]][[2]] 
    disease_matrix <- genetic_simulation[[rep_id]][[3]] 
    sim_para <- genetic_simulation[[rep_id]][[4]]
    # find which topic is genetically modulated
    tp_genetic <- colMeans(patient_loadings[which(disease_matrix[,1] == 1),]) %>% which.max
    loadings <- patient_loadings[,tp_genetic]
    MRresults <- eggerMR(loadings, disease_matrix[,sim_para$disease_number + 1], genetics_population, method = "ivw")
    try({
      p_t2d[rep_id] <- MRresults[[1]]@Pvalue
      p_d2t[rep_id] <- MRresults[[2]]@Pvalue
    })
  }
  t2d_mean <- mean(p_t2d < pw_thr,na.rm = T)
  t2d_se <- sqrt(t2d_mean*(1-t2d_mean)/length(p_t2d))
  d2t_mean <- mean(p_d2t < pw_thr,na.rm = T)
  d2t_se <- sqrt(d2t_mean*(1-d2t_mean)/length(p_d2t))
  df_2 <- df_2 %>%
    add_row(disease2topic = sim_para$disease2topic,  t2d_mean= t2d_mean, t2d_se = t2d_se,
            d2t_mean= d2t_mean, d2t_se = d2t_se)
}
plt_g2d2t <- ggplot(df_2,aes(x = disease2topic)) +
  geom_line(aes(y = t2d_mean), color = red, linetype = "dashed") + 
  geom_line(aes(y =  d2t_mean), color = green, linetype = "dashed") +
  geom_pointrange(aes(y = t2d_mean, ymin=t2d_mean-1.96*t2d_se, ymax = t2d_mean + 1.96*t2d_se, color = "G x Topic model"), size = .3) + 
  geom_pointrange(aes(y = d2t_mean, ymin=d2t_mean-1.96*d2t_se, ymax = d2t_mean + 1.96*d2t_se, color = "G x Disease model"), size = .3) + 
  scale_color_manual(values = c("G x Topic model" = red, "G x Disease model" = green)) + 
  labs(x ="Topic to disease effect size", y = "Power at pvalue=0.05") + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_hline(yintercept=pw_thr, linetype="dashed", color = grey, size = 1) +
  theme(legend.position = "None",panel.background=element_blank())
ggsave(paste0("../paper_writing/Production_figures/","MR_g2d2t.png"), plt_g2d2t, width = 4, height = 4)









##########################################
# previous task
# find a better statistics for identifying the number of subgroups in the posterior multinomials zn
##########################################
# step 1: simulations showing the recover of disease profiles
source("topic_functions.R")
sample_sz <- 10000 # 20000
topic_number <- 2 # 2
disease_number <- 20 # 30
degree_freedom <- 3
para <- simulate_age_topic_data(sample_sz, topic_number, disease_number, degree_freedom,ds_per_idv = 6.1)
para_org <- para
# make sure each disease is only kept for once
rec_data <- para$unlist_Ds_id %>% 
  rename(diag_icd10 = Ds_id) %>%
  select(eid, diag_icd10, age_diag) %>%
  group_by(eid, diag_icd10) %>% 
  arrange(age_diag, .by_group = T) %>% # keep only the first records of repeated diagnosis
  slice(1) %>%
  ungroup() %>%
  mutate(rowid = 1:n())
ds_list <- rec_data %>%
  group_by(diag_icd10) %>%
  summarise(occ = n())
para <- topic_init_age(rec_data, ds_list, topic_number, degree_freedom)
# infer topics on the simulated data set
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

################################
# save disease simulation examples
################################

# plot simulated beta
longData<-melt(para_org$true_beta[30:81,,1]) %>%
  mutate(Var1 = Var1+29, Var2 = factor(Var2, levels = 1:disease_number) )
plt <- ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Diseases", y="Age", title="") +
  theme_bw(base_size = 20) + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/example_topic_simulation.png"), plt, width = 6, height = 10)

longData<-melt(para$pi_beta_basis[30:81,,1]) %>%
  mutate(Var1 = Var1+29, Var2 = factor(Var2, levels = 1:disease_number) )
plt <- ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Diseases", y="Age", title="") +
  theme_bw(base_size = 20) + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
ggsave(paste0("~/Desktop/comorbidity/paper_writing/Production_figures/example_infertopic_simulation.png"), plt, width = 6, height = 10)

# testing two different distance:
# a <- c(1,0,0,0,0)
# b <- rep(1,5)/sum(rep(1,5))
# c <- c(0,0,1,0,0)
# eucl_dist <- c(sum((a - b)^2), sum((a - c)^2) ) 
# Dpq <- function(x, y ){
#   sum(x * log(2*x/(x + y + 10^(-256)) + 10^(-256)) + y * log(2*y/(x + y + 10^(-256)) + 10^(-256)))
# }
# Dpq_dist <- c(Dpq(a,b), Dpq(a,c))
#                

##########################################
# plotting the subtype simulation calibration results
##########################################
DIR <- "~/Desktop/comorbidity/Results/"
# pt <- paste0("CVB0_model_output_A2N_age_dependent_K", K,"_P",df_P, "_rep*")
pt <- paste0("^subtype_tests_simulation_rep*")
temp <- list.files(paste(DIR, sep=""), pattern=pt)
null_ps <- c()
for(rep_id in 1:length(temp)){
  load(paste0(DIR,temp[rep_id]))
  null_ps <- c(null_ps, as.vector(subtype_tests[[1]][1:8,1]))
}
expected_ps <- runif(length(null_ps))
df_null <- data.frame(expected_logP  = sort(expected_ps), Null_logP = sort(null_ps)) %>%
  mutate(expected_logP = -log(expected_logP), Null_logP = -log(Null_logP))
ggplot(data = df_null) +
  geom_point(aes(x = sort(expected_logP), y = sort(Null_logP))) +
  geom_abline(slope = 1, intercept = 0)

# get power of the two methods
DIR <- "~/Desktop/comorbidity/Results/"
# pt <- paste0("CVB0_model_output_A2N_age_dependent_K", K,"_P",df_P, "_rep*")
pt <- paste0("^subtype_tests_simulation_rep*")
ratio_list <- 0.02 * (0:25)
temp <- list.files(paste(DIR, sep=""), pattern=pt)
powers_agelda1 <- matrix(NA, nrow = length(temp), ncol = length(ratio_list))
powers_basiclda1 <- matrix(NA, nrow = length(temp), ncol = length(ratio_list))
powers_agelda2 <- matrix(NA, nrow = length(temp), ncol = length(ratio_list))
powers_basiclda2 <- matrix(NA, nrow = length(temp), ncol = length(ratio_list))
for(rep_id in 1:length(temp)){
  load(paste0(DIR,temp[rep_id]))
  powers_agelda1[rep_id,] <- subtype_tests[[1]][1,]
  powers_basiclda1[rep_id,] <- subtype_tests[[2]][1,]
  powers_agelda2[rep_id,] <- subtype_tests[[3]][1,]
  powers_basiclda2[rep_id,] <- subtype_tests[[4]][1,]
}
powers_basiclda1 <- apply(powers_basiclda1, 2, function(x) mean(x <= 0.05))
powers_agelda1 <- apply(powers_agelda1, 2, function(x) mean(x <= 0.05))
powers_basiclda1_var <- sqrt(powers_basiclda1*(1-powers_basiclda1)/100)
powers_agelda1_var <- sqrt(powers_agelda1*(1-powers_agelda1)/100)

powers_basiclda2 <- apply(powers_basiclda2, 2, function(x) mean(x <= 0.05))
powers_agelda2 <- apply(powers_agelda2, 2, function(x) mean(x <= 0.05))
powers_basiclda2_var <- sqrt(powers_basiclda2*(1-powers_basiclda2)/100)
powers_agelda2_var <- sqrt(powers_agelda2*(1-powers_agelda2)/100)
df_power <- data.frame(mixing_ratio = ratio_list, powers_basiclda1, powers_agelda1, powers_basiclda2, powers_agelda2, 
                       powers_basiclda1_var, powers_agelda1_var, powers_basiclda2_var, powers_agelda2_var)
ggplot(df_power) + 
  geom_line(aes(x=ratio_list, y = powers_basiclda1), color = blue) +
  geom_ribbon(aes(x=ratio_list, ymax = powers_basiclda1 + 1.96*powers_basiclda1_var, ymin = powers_basiclda1 - 1.96*powers_basiclda1_var ), alpha = 0.3) +
  geom_line(aes(x=ratio_list, y = powers_agelda1), color = red) +
  geom_ribbon(aes(x=ratio_list, ymax = powers_agelda1 + 1.96*powers_agelda1_var, ymin = powers_agelda1 - 1.96*powers_agelda1_var ), alpha = 0.3) +
  theme(legend.position = "none",panel.background=element_blank()) +
  ylab("Power")

ggplot(df_power) + 
  geom_line(aes(x=ratio_list, y = powers_basiclda2), color = blue) +
  geom_ribbon(aes(x=ratio_list, ymax = powers_basiclda2 + 1.96*powers_basiclda2_var, ymin = powers_basiclda2 - 1.96*powers_basiclda2_var ), alpha = 0.3) +
  geom_line(aes(x=ratio_list, y = powers_agelda2), color = red) +
  geom_ribbon(aes(x=ratio_list, ymax = powers_agelda2 + 1.96*powers_agelda2_var, ymin = powers_agelda2 - 1.96*powers_agelda2_var ), alpha = 0.3) +
  theme(legend.position = "none",panel.background=element_blank()) +
  ylab("Power")
