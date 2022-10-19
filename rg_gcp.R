source("topic_functions.R")
source("plotting_functions.R")


######################################################
# causal (genetic) analysis of diseases in the topics
######################################################
# step 1: filter the diseases / subtypes for each topics
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

load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
thre_pick <- 1000
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
write.table(common_disease_within_topics, "all_disease_topic_list.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)
write.table(common_disease_within_topics$disease, "causality_disease_list.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)

# save the disease list with age and sample size -- will be used to compare h2g across age and numbers 
common_disease_within_topics <- data.frame(disease = disease_idx, topic = topic_assign, mean_age = mean_age, number_case = number_case)
write.table(common_disease_within_topics, "disease_topic_list_size_meanage.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)

# write a table for only a couple of important diseasse subtypes 
common_disease_within_topics %>%
  filter(disease %in% c(250.2, 495, 401.1, 272.11)) %>%
  write.table("important_disease_topic_list.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)

trait_number <- data.frame(disease = disease_idx, topic = topic_assign, mean_age = mean_age, case_number = number_case)
write.csv(trait_number, "trait_number.csv", row.names = F)

# filtering the diseases based on z-score
h2g_disease_topics <- read.csv("causality_analysis/h2g_causality.csv")
filtered_heritable <- h2g_disease_topics %>% 
  filter(z_score > 4 | disease == 250.2) 
h2g_disease_topics %>%
  group_by(topic) %>%
  summarise(mean(h2g))
write.table(filtered_heritable, "trait1_LCV.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)
write.table(filtered_heritable, "trait2_LCV.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)

filtered_heritable %>%
  filter(disease %in% c(250.2, 495, 401.1, 272.11), topic != "all") %>%
  write.table("important_LCV1.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)
filtered_heritable %>%
  filter(disease %in% c(250.2, 495, 401.1, 272.11), topic != "all") %>%
  write.table("important_LCV2.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)


#######################################
# rg analysis of subtypes and diseases
#######################################
DIR <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/causality_analysis/LCV_results/"
temp <- list.files(paste(DIR, sep=""), pattern="^[^t].*RData")
LCV_data <- data.frame(trait1 = as.numeric(), trait2 = as.numeric(),topic1= as.character(), 
                       topic2 = as.character(), gcp= as.numeric(),gcp_se= as.numeric(),log10P= as.numeric(),
                       rho.zscore = as.numeric(), rho.est = as.numeric(), rho.err = as.numeric())
for(i in 1:length(temp)){
  load(paste0(DIR, temp[i]))
  row_data <- LCV_save[[1]]
  row_data$rho.zscore <- LCV_save[[2]]$rho.est/LCV_save[[2]]$rho.err
  row_data$rho.est <- LCV_save[[2]]$rho.est
  row_data$rho.err <- LCV_save[[2]]$rho.err
  # row_data <- row_data %>% 
  #   mutate(topic1 = as.character(topic1), topic2 = as.character(topic2))
  LCV_data <- LCV_data %>% 
    add_row(row_data)
}

########################################################
# compare the rg for topic specific and non-topic specific
#######################################################
heter_ds <- read.table("top_hererogeneous_disease.txt")
# output number of subtypes; the rest are diseases withouth subtypes
read.table("trait2_LCV.txt", header = F) %>%
  filter(V1 %in% heter_ds$V1, V2 != "all")

# below are the topic specific disease we are testing here
subtypes_ds <- read.csv("causality_analysis/h2g_causality.csv") %>%
  filter(topic != "all") %>%
  group_by(disease) %>%
  tally() %>%
  filter(n > 1)
# extract all the non-subtype disease
non_subtp <- LCV_data %>%
  filter(topic1 == "all" |  !(trait1 %in% heter_ds$V1), 
         topic2 == "all" |  !(trait2 %in% heter_ds$V1) ) %>%
  select(-topic1, -topic2) %>%
  rename(non_tp_gcp = gcp, non_tp_gcp_se = gcp_se, non_tp_log10P = log10P,
         non_tp_rho.zscore = rho.zscore, non_tp_rho.est = rho.est, non_tp_rho.err=rho.err)
subtp_rg <- LCV_data %>%
  filter(topic1 != "all", topic2 != "all") %>%
  left_join(non_subtp, by = c("trait1", "trait2")) %>%
  mutate(diff.zscore = (non_tp_rho.est - rho.est)/sqrt(rho.err^2 + non_tp_rho.err^2))

# plot the big matrix with all diseases
disease_order <- subtp_rg %>%
  group_by(topic1, trait1) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(ds_tp =  paste0(topic1,"_", trait1)) 

plot_rg <- subtp_rg %>%
  mutate(ds_tp1 =  factor(paste0(topic1,"_", trait1), levels = disease_order$ds_tp, ordered = TRUE), 
         ds_tp2 =  factor(paste0(topic2,"_", trait2), levels = disease_order$ds_tp, ordered = TRUE) ) %>%
  mutate(non_tp_rho.zscore = if_else(is.infinite(non_tp_rho.zscore),  # 
                                     max(subtp_rg$rho.zscore[!is.infinite(subtp_rg$rho.zscore)], na.rm = T),
                                     non_tp_rho.zscore)) %>%
  mutate(rg = if_else(ds_tp1 > ds_tp2, non_tp_rho.est, rho.est - non_tp_rho.est), 
         rg.zscore.abs = if_else(ds_tp1 > ds_tp2, abs(non_tp_rho.zscore), abs(diff.zscore)),
         shape_type = if_else(ds_tp1 > ds_tp2, "lower", "upper")) %>%
  mutate(P = if_else(ds_tp1 < ds_tp2, (1-pnorm(abs(diff.zscore)))*2, 9)) %>%
  filter(! ds_tp2 == ds_tp1, !is.na(rg), !is.na(P)) %>% 
        #  ds_tp1 > ds_tp2 | trait1 %in% subtypes_ds$disease | trait2 %in% subtypes_ds$disease ) %>%
  mutate(ds_tp2 =  factor(paste0(topic2,"_", trait2), levels = rev(disease_order$ds_tp), ordered = TRUE) )
######################
# p-values should only be computed for rho.err < 0.1 and one of the rho (either subtype or non subtype) has z-score > 4

# adjust p value
test_rg <- plot_rg %>%
  filter(rho.err < 0.1, non_tp_rho.err < 0.1, P < 1, 
         abs(rho.zscore) > 4 | abs(non_tp_rho.zscore) > 4, 
         trait1 %in% subtypes_ds$disease | trait2 %in% subtypes_ds$disease) %>%
  select(ds_tp1, ds_tp2, P) 
  
test_rg$Q <- p.adjust(test_rg$P, method = "fdr")
plot_rg <- plot_rg %>%
  left_join(select(test_rg, ds_tp1, ds_tp2, Q), by = c("ds_tp1", "ds_tp2"))
plot_rg$star_label = if_else(plot_rg$Q < 0.1, T, NA)

# only keep disease that has rho.zscore > 4 with at least one another disease.
disease_to_include <- plot_rg %>% 
  filter(abs(rho.zscore) > 4, trait1 != trait2) %>%
  group_by(trait1) %>%
  slice(1) %>%
  pull(trait1)

# cap the size of z-score in this figure to be 10
plot_rg_meaningful_ds <- plot_rg %>%
  filter(trait1 %in% disease_to_include, trait2 %in% disease_to_include) %>%
  mutate(rg.zscore.abs  = pmin(rg.zscore.abs , 10)) 
  
  
plt <- ggplot(plot_rg_meaningful_ds, aes(x = ds_tp1, y = ds_tp2)) +
  geom_tile(aes( width = 0.9, height = 0.9), fill = if_else(plot_rg_meaningful_ds$P<0.05, grey, "white")) + 
  geom_point(aes( fill = rg, size = rg.zscore.abs, shape = shape_type)) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0) + 
  scale_size_area(name = "|z-score|", max_size = 5) + 
  scale_shape_manual(values=c("lower"= 21, "upper" = 22)) + 
  geom_point(shape=8, size = if_else(plot_rg_meaningful_ds$star_label, 3, -0.1)) +
  labs(x = "", y = "") +
  theme_bw(base_size = 20)+
  theme(legend.position = "None", 
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_blank())

ggsave(paste0("../paper_writing/Production_figures/rg_all.png"), plt, width = 10, height = 10)

# save the topic pallete for visualization:
disease_order <- plot_rg_meaningful_ds %>%
  group_by(topic1, trait1) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(ds_tp =  paste0(topic1,"_", trait1)) 
color_data <-  matrix(rev(disease_order$topic1), ncol = 1) %>%
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
ggsave(paste0("../paper_writing/Production_figures/all_rg_topic_color.png"), plt_palette, width = 1.5, height = 10)



# save the data table
label_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CAD", "UGI", "LGI", "SRD")
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
plot_rg %>% 
  filter(P!=9, P!=1) %>%
  mutate(topic1 = label_name[as.numeric(topic1)], topic2 = label_name[as.numeric(topic2)]) %>%
  select(trait1,  trait2, topic1, topic2, rho.zscore, rho.est, rho.err, non_tp_rho.zscore, non_tp_rho.est, non_tp_rho.err, diff.zscore, rg, rg.zscore.abs,P, Q) %>%
  left_join(select(phe_phecode, phecode, phenotype), by = c("trait1" = "phecode")) %>%
  left_join(select(phe_phecode, phecode, phenotype), by = c("trait2" = "phecode")) %>%
  rename(phenotype.1 = phenotype.x, phenotype.2 = phenotype.y, 
         subtype.rho.zscore = rho.zscore, subtype.rho.est = rho.est, subtyp.rho.err = rho.err, 
         all.rho.zscore = non_tp_rho.zscore, all.rho.est = non_tp_rho.est, all.rho.err = non_tp_rho.err, 
         diff.rg = rg, diff.rg.zscore.abs = rg.zscore.abs) %>%
  write.csv("~/Desktop/comorbidity/paper_writing/supplementary_files/subtype_rg.csv", row.names = F)


############################
# Fig 5: rg for some topics
############################
label_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CAD", "UGI", "LGI", "SRD")
order_ds <- c(4, 5, 10, 7, 8, 9, 3, 6, 1, 2)
library(colorBlindness)
displayAvailablePalette(color="white")
colors_topic <- c(cbPalette, paletteMartin[c(8,11)])
names(colors_topic) <- as.character(1:10)
# plot_rg <- subtp_rg %>%
#   mutate(ds_tp1 =  factor(paste0(topic1,"_", trait1), levels = disease_order$ds_tp, ordered = TRUE), 
#          ds_tp2 =  factor(paste0(topic2,"_", trait2), levels = disease_order$ds_tp, ordered = TRUE) ) %>%
#   mutate(rg = if_else(ds_tp1 > ds_tp2, non_tp_rho.est, rho.est)) %>%
#   mutate(P = if_else(ds_tp1 < ds_tp2, (1-pnorm(abs(diff.zscore)))*2, 9)) %>%
#   mutate(ds_tp2 =  factor(paste0(topic2,"_", trait2), levels = rev(disease_order$ds_tp), ordered = TRUE))

disease_order <- plot_rg %>%
  filter(topic1 %in% c(5,6, 7), topic2 %in% c(5,6,7)) %>%
  filter(!(trait1 %in% c(550.10, 550.40, 594.10, 600.00, 728.71, 272.10, 427.20, 411.10, 428.20, 250.10, 318.00, 496.10, 496.21)),
         !(trait2 %in% c(550.10, 550.40, 594.10, 600.00, 728.71, 272.10, 427.20, 411.10, 428.20, 250.10, 318.00, 496.10, 496.21))) %>%
  group_by(topic1, trait1) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(ds_tp =  paste0(topic1,"_", trait1)) 

phe_phecode <- read.csv("info_phenotype_phecode.csv")
disease_list <- disease_order %>%
  left_join(select(phe_phecode, phecode, phenotype), by = c("trait1" = "phecode")) %>%
  select(phenotype)
write.table(disease_list, "../paper_writing/Production_figures/high_h2_disease_list.txt", row.names = F, col.names = F, quote = F)

plot_subrg <- plot_rg %>% filter(ds_tp1 %in% disease_order$ds_tp, ds_tp2 %in% disease_order$ds_tp)
Q_sig <- plot_subrg %>%
  filter(Q < 0.1)
P_nomi <- plot_subrg %>%
  filter(P < 0.01)

sub_rg_plt <- ggplot(plot_subrg, aes(x = ds_tp1, y = ds_tp2)) +
  # geom_tile(aes( width = 0.8, height = 0.8), fill = if_else(plot_subrg$P<0.05, grey, "white"), alpha = 0.4) + 
  geom_point(aes( fill = rg, size = rg.zscore.abs, shape = shape_type)) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0) + 
  scale_size_area(name = expression(-log[10](P)) , max_size = 10) + 
  scale_shape_manual(values=c("lower"= 21, "upper" = 22)) + 
  geom_point(shape=8, size = if_else(plot_subrg$star_label, 3, 0)) +
  labs(x = "", y = "") +
  theme_bw(base_size = 20)+
  theme(legend.position = "None", 
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_blank())

ggsave(paste0("../paper_writing/Production_figures/rg_main.png"), sub_rg_plt, width = 7, height = 7)

color_data <-  matrix(rev(disease_order$topic1), ncol = 1) %>%
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
ggsave(paste0("../paper_writing/Production_figures/rg_topic_color.png"), plt_palette, width = 1, height = 8)

# plot legend
plt_rg <- ggplot(plot_subrg, aes(x = ds_tp1, y = ds_tp2)) +
  geom_tile(aes( width = 0.9, height = 0.9), fill = if_else(plot_subrg$P<0.05, grey, "white")) + 
  geom_point(aes( fill = rg, size = rg.zscore.abs, shape = shape_type)) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0) + 
  scale_size_area(name = "|z-score|", max_size = 10) + 
  scale_shape_manual(values=c("lower"= 21, "upper" = 22)) + 
  geom_point(shape=8, size = if_else(plot_subrg$Q<0.1, 3, -0.1)) +
  labs(x = "", y = "") +
  theme_bw(base_size = 20)+
  theme( axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_blank())
legend_plt <- cowplot::get_legend(plt_rg)
grid.newpage()
legend_plt <- grid.draw(legend_plt)

#################################
# saving lower and upper triangle separately
#################################
# get the fdr z-score threshold:
p_max <- plot_subrg %>% # p = 0.008
  filter(Q < 0.1) %>% 
  summarise(max(P))
full_size <- abs(qnorm(p_max[[1]]/2))
plot_difference <- plot_subrg %>%
  filter(shape_type == "upper") %>%
  mutate(rg.zscore.abs = pmin(rg.zscore.abs, full_size)) %>%
  mutate(ds_tp1 =  factor(ds_tp1, levels = disease_order$ds_tp, ordered = TRUE),
         ds_tp2 =  factor(ds_tp2, levels = rev(disease_order$ds_tp), ordered = TRUE)) %>% 
  filter(rg != 0)

# create the box tiles in the plot
grouping <- expand.grid(ds_tp1 = disease_order$ds_tp, ds_tp2 = disease_order$ds_tp) %>%
  mutate(ds_tp1 =  factor(ds_tp1, levels = disease_order$ds_tp, ordered = TRUE),
         ds_tp2 =  factor(ds_tp2, levels = disease_order$ds_tp, ordered = TRUE)) %>%
  filter(ds_tp1 <= ds_tp2) %>%
  mutate(ds_tp2 =  factor(ds_tp2, levels = rev(disease_order$ds_tp), ordered = TRUE)) 

# add dark background to shallow 
black_block <- grouping %>%
  anti_join(plot_difference, by = c("ds_tp1", "ds_tp2"))

sub_rg_plt <- ggplot(plot_difference, aes(x = ds_tp1, y = ds_tp2)) +
  # geom_tile(aes( width = 0.8, height = 0.8), fill = if_else(plot_subrg$P<0.05, grey, "white"), alpha = 0.4) + 
  geom_tile(aes( fill = rg, width = sqrt(rg.zscore.abs/full_size), height =sqrt(rg.zscore.abs/full_size)) ) +
  geom_tile(data = grouping, aes(x = ds_tp1, y = ds_tp2), color = grey, alpha = 0) + 
  geom_tile(data = black_block, aes(x = ds_tp1, y = ds_tp2), color = grey, fill = "black", alpha = 0.5) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0) + 
  # scale_size_area(name = expression(-log[10](P)) , max_size = 10) + 
  # scale_shape_manual(values=c("lower"= 21, "upper" = 22)) + 
  geom_point(shape=8, size = if_else(plot_difference$star_label, 3, 0)) +
  labs(x = "", y = "") +
  theme_void() + #theme_bw(base_size = 20)+
  theme(legend.position = "None", 
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_blank()) +
  scale_y_discrete(limits=rev) # reverse the figure
ggsave(paste0("../paper_writing/Production_figures/rg_upper.png"), sub_rg_plt, width = 7, height = 7)

sub_rg_plt <- ggplot(plot_difference, aes(x = ds_tp1, y = ds_tp2)) +
  # geom_tile(aes( width = 0.8, height = 0.8), fill = if_else(plot_subrg$P<0.05, grey, "white"), alpha = 0.4) + 
  geom_tile(aes( fill = rg, width = rg.zscore.abs/full_size, height = rg.zscore.abs/full_size)) +
  geom_tile(data = grouping, aes(x = ds_tp1, y = ds_tp2), color = grey, alpha = 0) + 
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0) + 
  # scale_size_area(name = expression(-log[10](P)) , max_size = 10) + 
  # scale_shape_manual(values=c("lower"= 21, "upper" = 22)) + 
  labs(x = "", y = "") +
  theme_void() + #theme_bw(base_size = 20)+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_blank()) +
  scale_y_discrete(limits=rev) 
legend_plt <- cowplot::get_legend(sub_rg_plt)
grid.newpage()
legend_plt <- grid.draw(legend_plt)


# plot disease that have multiple subtypes
plot_difference %>% 
  filter(trait1 %in% heter_ds$V1) %>%
  group_by(trait1) %>%
  slice(1)

# plot the upper triangle of disease non-topic separately.
full_size <- 4
plot_lower <- plot_subrg %>%
  filter(shape_type == "lower") %>%
  filter(!(ds_tp1 %in% c("6_250.2", "7_250.2", "6_401.1", "7_401.1")),
         !(ds_tp2 %in% c("6_250.2", "7_250.2", "6_401.1", "7_401.1"))) %>%
  mutate(ds_tp1 =  factor(ds_tp1, levels = disease_order$ds_tp, ordered = TRUE),
         ds_tp2 =  factor(ds_tp2, levels = disease_order$ds_tp, ordered = TRUE))  %>%
  mutate(rg.zscore.abs = pmin(rg.zscore.abs, full_size)) 
list_ds <-  c(as.character(unique(plot_lower$ds_tp1)), "5_250.2")
grouping <- expand.grid(ds_tp1 = list_ds, ds_tp2 =list_ds) %>%
  mutate(ds_tp1 =  factor(ds_tp1, levels = disease_order$ds_tp, ordered = TRUE),
         ds_tp2 =  factor(ds_tp2, levels = disease_order$ds_tp, ordered = TRUE)) %>% 
  filter(ds_tp1 >= ds_tp2)

# add dark background to shallow 
black_block <- grouping %>%
  anti_join(plot_lower, by = c("ds_tp1", "ds_tp2"))

sub_rg_lower <- ggplot(plot_lower, aes(x = ds_tp1, y = ds_tp2)) +
  # geom_tile(aes( width = 0.8, height = 0.8), fill = if_else(plot_subrg$P<0.05, grey, "white"), alpha = 0.4) + 
  geom_point(aes( fill = rg, size = rg.zscore.abs/full_size), shape = 21) +
  geom_tile(data = grouping, aes(x = ds_tp1, y = ds_tp2), color = grey, alpha = 0) + 
  geom_tile(data = black_block, aes(x = ds_tp1, y = ds_tp2), color = grey, fill = "black", alpha = 0.5) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0) + 
  scale_size_area(name = "rg" , max_size = 20) +
  # scale_size_area(name = expression(-log[10](P)) , max_size = 10) + 
  # scale_shape_manual(values=c("lower"= 21, "upper" = 22)) + 
  labs(x = "", y = "") +
  theme_void() + #theme_bw(base_size = 20)+
  theme(legend.position = "None", 
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_blank()) +
  scale_y_discrete(limits=rev) # reverse the figure
ggsave(paste0("../paper_writing/Production_figures/rg_lower.png"), sub_rg_lower, width = 7, height = 7)
sub_rg_lower <- ggplot(plot_lower, aes(x = ds_tp1, y = ds_tp2)) +
  # geom_tile(aes( width = 0.8, height = 0.8), fill = if_else(plot_subrg$P<0.05, grey, "white"), alpha = 0.4) + 
  geom_point(aes( fill = rg, size = rg.zscore.abs/full_size), shape = 20) +
  geom_tile(data = grouping, aes(x = ds_tp1, y = ds_tp2), color = grey, alpha = 0) + 
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0) + 
  scale_size_area(name = "rg" , max_size = 13) +
  # scale_size_area(name = expression(-log[10](P)) , max_size = 10) + 
  # scale_shape_manual(values=c("lower"= 21, "upper" = 22)) + 
  labs(x = "", y = "") +
  theme_void() + #theme_bw(base_size = 20)+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_blank()) 
legend_plt <- cowplot::get_legend(sub_rg_lower)
grid.newpage()
legend_plt <- grid.draw(legend_plt)


##########################
# plot a small version of rg -- used in ppt presentation
##########################
  label_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CAD", "UGI", "LGI", "SRD")
  order_ds <- c(4, 5, 10, 7, 8, 9, 3, 6, 1, 2)
  library(colorBlindness)
  displayAvailablePalette(color="white")
  colors_topic <- c(cbPalette, paletteMartin[c(8,11)])
  names(colors_topic) <- as.character(1:10)
  # plot_rg <- subtp_rg %>%
  #   mutate(ds_tp1 =  factor(paste0(topic1,"_", trait1), levels = disease_order$ds_tp, ordered = TRUE), 
  #          ds_tp2 =  factor(paste0(topic2,"_", trait2), levels = disease_order$ds_tp, ordered = TRUE) ) %>%
  #   mutate(rg = if_else(ds_tp1 > ds_tp2, non_tp_rho.est, rho.est)) %>%
  #   mutate(P = if_else(ds_tp1 < ds_tp2, (1-pnorm(abs(diff.zscore)))*2, 9)) %>%
  #   mutate(ds_tp2 =  factor(paste0(topic2,"_", trait2), levels = rev(disease_order$ds_tp), ordered = TRUE))
  
  disease_order <- plot_rg %>%
    filter(topic1 %in% c(5,6, 7), topic2 %in% c(5,6,7)) %>%
    filter((trait1 %in% c(272.11, 495, 250.2, 401.1 )),
           (trait2 %in% c(272.11, 495, 250.2, 401.1 ))) %>%
    group_by(topic1, trait1) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(ds_tp =  paste0(topic1,"_", trait1)) 
  
  phe_phecode <- read.csv("info_phenotype_phecode.csv")
  disease_list <- disease_order %>%
    left_join(select(phe_phecode, phecode, phenotype), by = c("trait1" = "phecode")) %>%
    select(phenotype)
  write.table(disease_list, "../paper_writing/Production_figures/rg_ppt_disease_list.txt", row.names = F, col.names = F, quote = F)
  
  plot_subrg <- plot_rg %>% filter(ds_tp1 %in% disease_order$ds_tp, ds_tp2 %in% disease_order$ds_tp)
  Q_sig <- plot_subrg %>%
    filter(Q < 0.1)
  P_nomi <- plot_subrg %>%
    filter(P < 0.05)
  
  sub_rg_plt <- ggplot(plot_subrg, aes(x = ds_tp1, y = ds_tp2)) +
    geom_tile(aes( width = 0.8, height = 0.8), fill = if_else(plot_subrg$P<0.05, grey, "white"), alpha = 0.4) + 
    geom_point(aes( fill = rg, size = rg.zscore.abs, shape = shape_type)) +
    scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0) + 
    scale_size_area(name = expression(-log[10](P)), max_size = 10) + 
    scale_shape_manual(values=c("lower"= 21, "upper" = 22)) + 
    geom_point(shape=8, size = if_else(plot_subrg$star_label, 3, 0)) +
    labs(x = "", y = "") +
    theme_bw(base_size = 20)+
    theme(legend.position = "None", 
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          plot.title=element_blank())
  
  
  ggsave(paste0("../paper_writing/Production_figures/rg_ppt.png"), sub_rg_plt, width = 7, height = 7)

color_data <-  matrix(rev(disease_order$topic1), ncol = 1) %>%
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
ggsave(paste0("../paper_writing/Production_figures/rg_ppt_color.png"), plt_palette, width = 1, height = 8)

##############################################################################
# rg where topic loading are adjusted in the case-control matching procedure
##############################################################################
# loading the topic_matched data
DIR <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/causality_analysis/LCV_results/"
temp <- list.files(paste(DIR, sep=""), pattern="^topic_match")
LCV_topic_match <- data.frame(trait1 = as.numeric(), trait2 = as.numeric(),topic1= as.numeric(), 
                       topic2 = as.numeric(), gcp= as.numeric(),gcp_se= as.numeric(),log10P= as.numeric(),
                       rho.zscore = as.numeric(), rho.est = as.numeric(), rho.err = as.numeric())
for(i in 1:length(temp)){
  load(paste0(DIR, temp[i]))
  row_data <- LCV_save[[1]]
  row_data$rho.zscore <- LCV_save[[2]]$rho.est/LCV_save[[2]]$rho.err
  row_data$rho.est <- LCV_save[[2]]$rho.est
  row_data$rho.err <- LCV_save[[2]]$rho.err
  LCV_topic_match <- LCV_topic_match %>% 
    add_row(row_data)
}
# non_subtp need to be computed from previous section
subtp_rg <- LCV_topic_match %>%
  filter(!(trait1 == 401.1 & topic1 == 7), !(trait2 == 401.1 & topic2 == 7)) %>%
  left_join(non_subtp, by = c("trait1", "trait2")) %>%
  mutate(diff.zscore = (non_tp_rho.est - rho.est)/sqrt(rho.err^2 + non_tp_rho.err^2))

# plot the big matrix with all diseases
disease_order <- subtp_rg %>%
  group_by(topic1, trait1) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(ds_tp =  paste0(topic1,"_", trait1)) 

plot_topic_match_rg <- subtp_rg %>%
  mutate(ds_tp1 =  factor(paste0(topic1,"_", trait1), levels = disease_order$ds_tp, ordered = TRUE), 
         ds_tp2 =  factor(paste0(topic2,"_", trait2), levels = disease_order$ds_tp, ordered = TRUE) ) %>%
  mutate(non_tp_rho.zscore = if_else(is.infinite(non_tp_rho.zscore),  # 
                                     max(subtp_rg$rho.zscore[!is.infinite(subtp_rg$rho.zscore)], na.rm = T),
                                     non_tp_rho.zscore)) %>%
  mutate(rg = if_else(ds_tp1 > ds_tp2, non_tp_rho.est, rho.est - non_tp_rho.est), 
         rg.zscore.abs = if_else(ds_tp1 > ds_tp2, abs(non_tp_rho.zscore), 5 * abs(diff.zscore)),
         shape_type = if_else(ds_tp1 > ds_tp2, "lower", "upper")) %>%
  mutate(P = if_else(ds_tp1 < ds_tp2, (1-pnorm(abs(diff.zscore)))*2, 9)) %>%
  filter(! ds_tp2 == ds_tp1, !is.na(rg), !is.na(P)) %>% 
  #  ds_tp1 > ds_tp2 | trait1 %in% subtypes_ds$disease | trait2 %in% subtypes_ds$disease ) %>%
  mutate(ds_tp2 =  factor(paste0(topic2,"_", trait2), levels = rev(disease_order$ds_tp), ordered = TRUE) )




label_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CAD", "UGI", "LGI", "SRD")
order_ds <- c(4, 5, 10, 7, 8, 9, 3, 6, 1, 2)
library(colorBlindness)
displayAvailablePalette(color="white")
colors_topic <- c(cbPalette, paletteMartin[c(8,11)])
names(colors_topic) <- as.character(1:10)
# plot_rg <- subtp_rg %>%
#   mutate(ds_tp1 =  factor(paste0(topic1,"_", trait1), levels = disease_order$ds_tp, ordered = TRUE), 
#          ds_tp2 =  factor(paste0(topic2,"_", trait2), levels = disease_order$ds_tp, ordered = TRUE) ) %>%
#   mutate(rg = if_else(ds_tp1 > ds_tp2, non_tp_rho.est, rho.est)) %>%
#   mutate(P = if_else(ds_tp1 < ds_tp2, (1-pnorm(abs(diff.zscore)))*2, 9)) %>%
#   mutate(ds_tp2 =  factor(paste0(topic2,"_", trait2), levels = rev(disease_order$ds_tp), ordered = TRUE))


phe_phecode <- read.csv("info_phenotype_phecode.csv")
disease_list <- disease_order %>%
  left_join(select(phe_phecode, phecode, phenotype), by = c("trait1" = "phecode")) %>%
  select(phenotype)
write.table(disease_list, "../paper_writing/Production_figures/rg_topic_match_disease_list.txt", row.names = F, col.names = F, quote = F)

# plot_subrg <- plot_topic_match_rg %>% filter(ds_tp1 %in% disease_order$ds_tp, ds_tp2 %in% disease_order$ds_tp)
# Q_sig <- plot_subrg %>%
#   filter(Q < 0.1)
# P_nomi <- plot_subrg %>%
#   filter(P < 0.05)

sub_rg_plt <- ggplot(plot_topic_match_rg, aes(x = ds_tp1, y = ds_tp2)) +
  geom_tile(aes( width = 0.8, height = 0.8), fill = if_else(plot_topic_match_rg$P<0.05, grey, "white"), alpha = 0.4) + 
  geom_point(aes( fill = rg, size = rg.zscore.abs, shape = shape_type)) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0) + 
  scale_size_area(name = expression(-log[10](P)), max_size = 10) + 
  scale_shape_manual(values=c("lower"= 21, "upper" = 22)) + 
  # geom_point(shape=8, size = if_else(plot_subrg$star_label, 3, 0)) +
  labs(x = "", y = "") +
  theme_bw(base_size = 20)+
  theme(legend.position = "None", 
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_blank())


ggsave(paste0("../paper_writing/Production_figures/rg_topic_match.png"), sub_rg_plt, width = 7, height = 7)

color_data <-  matrix(rev(disease_order$topic1), ncol = 1) %>%
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
ggsave(paste0("../paper_writing/Production_figures/rg_topic_match_color.png"), plt_palette, width = 1, height = 8)

# comparing the values in two cases: 
# note the plot_rg are run from the big figure plotting

rg_difference_non_topic_match <- plot_rg %>% 
  filter(shape_type == "upper") %>%
  select(trait1, trait2, topic1, topic2, rg) %>%
  mutate(topic1 = as.numeric(topic1), topic2 = as.numeric(topic2)) %>%
  rename(non_match.rg = rg)

rg_match_nonmatch <- plot_topic_match_rg %>% 
  filter(shape_type == "upper") %>%
  select(trait1, trait2, topic1, topic2, rg) %>%
  rename(match.rg = rg) %>%
  left_join(rg_difference_non_topic_match, by = c("trait1", "trait2", "topic1", "topic2"))

# make the scatter plot
plt_matched_non_match <- ggplot(rg_match_nonmatch) +
  geom_point(aes(x = non_match.rg, y =  match.rg)) + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "rg after adjusting for topic loadings", y = "rg") +
  theme_bw(base_size = 20)+
  theme(legend.position = "None", 
        plot.title=element_blank())
ggsave(paste0("../paper_writing/Production_figures/rg_scatter_topic_match.png"), plt_matched_non_match, width = 6, height = 6)


#################
# GCP 
#################
plot_rg %>% 
  filter(abs(rho.zscore) > 1.96, rho.est < 0.9) %>% 
  filter(abs(gcp)/gcp_se > 3,  abs(non_tp_gcp)/non_tp_gcp_se < 1)
plot_rg %>% 
  filter(abs(rho.zscore) > 1.96, rho.est < 0.9) %>% 
  filter( abs(gcp - non_tp_gcp)/sqrt(non_tp_gcp_se^2 + gcp_se^2) >3 )

################################
# Fig 5 comparing GCP for subtypes
################################
ds_id <- "250.2"

df_gcp <- subtp_rg %>% 
  filter(abs(rho.zscore) > 1.96, rho.est < 0.9) %>%
  filter(trait1 == ds_id, -log10P > 3 | -non_tp_log10P > 3)

other_ds <- c(272.11, 495, 550.2, 574.1)
plt_gcp <- subtp_rg %>% 
  filter(trait1 == ds_id, trait2 %in% other_ds) %>%
  select(trait1, trait2, topic1, topic2, gcp, log10P, non_tp_gcp, non_tp_log10P) %>%
  left_join(select(phe_phecode, phecode, phenotype), by =c("trait2" = "phecode") ) %>%
  filter(!(trait2 == 495 & topic2 ==3))  %>%
  mutate(topic1 = label_name[as.numeric(topic1)], topic2 = label_name[as.numeric(topic2)])

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



# file: extract_ldsc_rg_intercept.R
# extract intercept for LDSC according to equation 16 in supplementary note of Bulik-Sullivan 2015 NG
# compare if the results are consistent and whether the se of rg is smaller
subtype_data <- read.csv("Association_analysis/subtypes_rg.csv") 
not_constraint_intercept_subtype <- read.csv("Association_analysis/rg_intercept_not_fix.csv")

compare_intercept <- data.frame(rg.1 = subtype_data$rg.1, serg.1 = subtype_data$serg.1, 
                                rg.1.NotIncptCnst = not_constraint_intercept_subtype$rg.1, 
                                serg.1.NotIncptCnst =not_constraint_intercept_subtype$serg.1) 

ggplot(data = compare_intercept) +
  geom_point(aes(x = rg.1, y = rg.1.NotIncptCnst), size = 1) +
  labs(x = "Constraint on intecept", y = "No constraints") +
  geom_abline(slope = 1) + 
  theme_bw() 
ggplot(data = compare_intercept) +
  geom_point(aes(x = serg.1, y = serg.1.NotIncptCnst), size = 1) +
  labs(x = "Constraint on intecept", y = "No constraints", title = "s.e. of rg estimates") +
  geom_abline(slope = 1) + 
  lims(x = c(0,1.5), y = c(0,1.5))+
  theme_bw() 


# filtering the diseases: each subtype has to have at least 200 cases
subtype_diseases_list <- read.csv("subtype_disesea_list.csv")
min_sz <- rep(NA, dim(subtype_diseases_list)[1])
for(idx in 1:dim(subtype_diseases_list)[1]){
  j <- subtype_diseases_list$disease_idx[idx]
  loadings <- para$unlist_zn[para$ds_list[[j]]$id,]
  km_fit <- kmeans(loadings, subtype_diseases_list$cluster_number[idx])
  min_sz[idx] <- min(sapply(1:subtype_diseases_list$cluster_number[idx], function(x) sum(km_fit$cluster == x)))
}
subtype_diseases_list$min_sz <- min_sz
subtype_diseases_list %>% 
  filter(min_sz > 500)

subtype_data <- read.csv("Association_analysis/subtypes_rg.csv") %>%
  left_join(select(subtype_diseases_list,diag_icd10, occ, phenotype, cluster_number, min_sz ) , by = c("disease" = "diag_icd10"))

# plot Fst P-value
Fst_p <- subtype_data %>%
  filter(!is.na(subtype_data$p_fst))
ggplot(Fst_p, aes(label=phenotype)) +
  geom_point(aes(x= sort(-log(runif(dim(Fst_p)[1]))), y = sort(-log(p_fst))) , size = 1) +
  geom_abline(slope = 1) + 
  theme_bw() 
sum(p.adjust(Fst_p$p_fst, method = "fdr") < 0.1)


######################
# topic versus trait 
######################
DIR <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/causality_analysis/LCV_topic/"
temp <- list.files(paste(DIR, sep=""), pattern=".RData")
LCV_topic <- list()
for(i in 1:length(temp)){
  load(paste0(DIR, temp[i]))
  LCV_topic[[i]] <- LCV_data
}
LCV_topic <- bind_rows(LCV_topic) 
filtered_topic <- LCV_topic  %>%
  filter(abs(rho.zscore) > 3, abs(rho.est) < 0.9)
filtered_topic$Q <- p.adjust(10^(filtered_topic$log10P), method = "fdr")



# make sure LCV assumption is not violated
# (1) check if there is model violation (2) constrain rg z-score to be 2.58 (P value < 0.01)
p_thre <- filtered_topic %>% filter(Q < 0.05) %>% summarise(max(log10P))
disease_list <- LCV_topic  %>%
  filter(abs(rho.zscore) > 2.58, log10P < p_thre[[1]],  abs(gcp/gcp_se) > 1.96,
         !(p.topicfullycasaul < 0.1 & p.traitfullycasaul < 0.1)) %>%
  pull(trait) %>% unique()


topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CAD", "UGI", "LGI", "SRD")
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
topic_list <- c(1,4,5,7,10)
# # data to save in a data file
df_gcp_data <- LCV_topic %>%
  filter( # topic %in% topic_list,
    abs(rho.zscore) > 2.58, 
    !(trait_topic != 7 & trait == 272.11),
    !(trait_topic != 7 & trait == 401.1),
    !(trait_topic != "all" & trait == 250.2),
    !(trait_topic != "all" & trait == 495)) %>%
  mutate(topic = topic_name[topic],
         log10P = ifelse(abs(gcp/gcp_se) > 1.96 & abs(rho.zscore) > 2.58 & log10P < p_thre[[1]],
                         log10P, 0)) 
df_gcp_data %>% dim
# data that are included
df_plot <- LCV_topic %>%
  filter( # topic %in% topic_list,
    abs(rho.zscore) > 1.96, # abs(rho.est) < 0.9,
         trait %in% disease_list, 
         !(trait_topic != 7 & trait == 272.11),
    !(trait_topic != 7 & trait == 401.1),
      !(trait_topic != "all" & trait == 250.2),
    !(trait_topic != "all" & trait == 495)) %>%
  mutate(topic = topic_name[topic],
         log10P = ifelse(abs(gcp/gcp_se) > 1.96 & abs(rho.zscore) > 2.58 & log10P < p_thre[[1]],
                         log10P, 0)) %>%
  filter(log10P < 0) %>%
  mutate(topic = factor(topic, levels = 
        c( "NRI", "CER", "SRD", "CAD", "UGI", "LGI", "FGND", "MGND","MDS", "ARP"))) %>%
  left_join(select(phe_phecode, phecode, phenotype), by = c("trait" = "phecode")) %>%
  mutate(phenotype = ifelse(phenotype == "Unstable angina (intermediate coronary syndrome)", "Unstable angina", phenotype)) %>%
  mutate(phenotype = ifelse(phenotype == "Essential hypertension", "Essential hypertension (CAD topic)", phenotype)) %>%
  mutate(phenotype = ifelse(phenotype == "Hypercholesterolemia", "Hypercholesterolemia (CAD topic)", phenotype)) 
# order disease for plotting
ds_order <- df_plot %>%
  arrange(trait) %>%
  group_by(trait) %>%
  pull(phenotype)
df_plot <- df_plot %>% 
  mutate(phenotype = factor(phenotype, levels = rev(unique(ds_order))))
# use color that are different
library(colorBlindness)
paletteMartin[[4]]
Brown2Blue12Steps[[10]]
plt <- ggplot(df_plot) +
  geom_point(aes(x = topic, y = phenotype, fill = -gcp, size = -log10P), shape = 21) +
  scale_fill_gradient2(low = "#32E3FF", mid = "white", high = "#ff6db6") + 
  scale_size_area(max_size = 10) + 
  labs(x = "", y = "") +
  theme_bw(base_size = 20)+
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/causal_disease.png", plt, width = 12, height = 5)
# plot legend
plt_legend <- ggplot(df_plot) +
  geom_point(aes(x = topic, y = phenotype, fill = -gcp, size = -log10P), shape = 21) +
  scale_fill_gradient2(low = "#32E3FF", mid = "white", high = "#ff6db6") + 
  scale_size_area(max_size = 10) + 
  labs(x = "", y = "") +
  theme_bw(base_size = 20)+
  theme( 
    legend.position="right",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
ggsave("~/Desktop/comorbidity/paper_writing/Production_figures/causal_disease_legend.png", plt_legend, width = 14, height = 5)


###############################################
# first: plot heritability differences
###############################################
phe_phecode <- read.csv("info_phenotype_phecode.csv")
h2g_disease_topics <- read.csv("causality_analysis/h2g_causality.csv")
subtype_list <- h2g_disease_topics %>% group_by(disease) %>% summarise(num = n()) %>% filter(num > 1)
h2g_disease_topics %>% filter(disease %in% subtype_list$disease)
df_subtype <- data.frame(disease = as.character(), age.1 =as.numeric(), age.2 =as.numeric(), age.3 =as.numeric(),
                         h2g.1=as.numeric(), h2g.2=as.numeric(), h2g.3=as.numeric(), h2g_se.1=as.numeric(), 
                         h2g_se.2=as.numeric(), h2g_se.3=as.numeric(), topic.1=as.numeric(), 
                         topic.2=as.numeric(), topic.3=as.numeric())
for(ds in subtype_list$disease){
  subtps <- h2g_disease_topics %>%
    filter(disease == ds) %>%
    arrange(age)
  df_subtype <- df_subtype %>%
    add_row(disease = as.character(ds), age.1 = subtps$age[1], age.2 =subtps$age[2], age.3 =subtps$age[3],
            h2g.1=subtps$h2g[1], h2g.2=subtps$h2g[2], h2g.3=subtps$h2g[3], h2g_se.1=subtps$seh2g[1], 
            h2g_se.2=subtps$seh2g[2], h2g_se.3=subtps$seh2g[3],  topic.1=subtps$topic[1], 
            topic.2=subtps$topic[2], topic.3=subtps$topic[3])
}



phe_phecode <- phe_phecode %>% 
  mutate(phecode = as.character(phecode))
df_plot <- df_subtype %>% 
  left_join(select(phe_phecode, phecode, phenotype), by = c("disease" = "phecode")) %>%
  mutate(z_h2g_diff = (h2g.1 - h2g.2)/sqrt(h2g_se.1^2 + h2g_se.2^2)) %>%
  mutate(p_h2g_diff = 2*pnorm(q=abs(z_h2g_diff), lower.tail=FALSE))
df_plot$q_fdr_h2g_diff <- p.adjust(df_plot$p_h2g_diff, method = "fdr")

ggplot(df_plot,aes(x= factor(phenotype), y = (h2g.1 - h2g.2), label=phenotype)) +
  geom_pointrange(aes(ymin=(h2g.1 - h2g.2) - sqrt(h2g_se.1^2 + h2g_se.2^2), ymax=(h2g.1 - h2g.2) + sqrt(h2g_se.1^2 + h2g_se.2^2), color = age.2  - age.1), size = 1) +
  scale_color_gradient2(low= "blue", mid="grey90",  high="red") +
  geom_label_repel(aes(label=ifelse(q_fdr_h2g_diff < 0.1,as.character(phenotype),'')), max.overlaps = 50) +
  theme_bw()

# plot h2g differences
subtype_h2g_diff <- subtype_data %>%
  filter(h2g.1/h2g_se.1 > 3 | h2g.2/h2g_se.2 > 3 | h2g.3/h2g_se.3 > 3) %>% 
  mutate(youngh2g = if_else(age.1 < age.2, h2g.1, h2g.2), 
         oldh2g= if_else(age.1 < age.2, h2g.2, h2g.1),
         youngh2g_se = if_else(age.1 < age.2, h2g_se.1 , h2g_se.2 ), 
         oldh2g_se = if_else(age.1 < age.2, h2g_se.2 , h2g_se.1 )) %>%
  mutate(h2g_diff = youngh2g - oldh2g, se_diff = sqrt(youngh2g_se^2 + oldh2g_se^2), disease = as.factor(disease)) %>%
  mutate(p_h2g_diff = (1-pnorm(abs(h2g_diff)/se_diff)) *2 )  %>% # need to multiply by 2 as it is a two side test
  arrange(h2g_diff)
subtype_h2g_diff$q_fdr_h2g_diff <- p.adjust(subtype_h2g_diff$p_h2g_diff, method = "fdr")
ggplot(subtype_h2g_diff,aes(x= disease, y = h2g_diff, label=phenotype)) +
  geom_pointrange(aes(ymin=h2g_diff - se_diff, ymax=h2g_diff + se_diff, color = rg.1), size = 1) +
  scale_color_gradient2(low= "blue", mid="grey90",  high="red") +
  geom_label_repel(aes(label=ifelse(q_fdr_h2g_diff < 0.1,as.character(phenotype),'')), max.overlaps = 50) +
  theme_bw() 

# second: plot genetic correlation 
subtype_rg <- subtype_data %>%
  filter(h2g.1/h2g_se.1 > 3 & (h2g.2/h2g_se.2 > 3 | h2g.3/h2g_se.3 > 3) ) %>% # filter rg > than target
  mutate(youngh2g = if_else(age.1 < age.2, h2g.1, h2g.2), 
         oldh2g= if_else(age.1 < age.2, h2g.2, h2g.1),
         youngh2g_se = if_else(age.1 < age.2, h2g_se.1 , h2g_se.2 ), 
         oldh2g_se = if_else(age.1 < age.2, h2g_se.2 , h2g_se.1 ))
# filter(h2g.1 - 1.96 * h2g_se.1 > 0 | h2g.2 - 1.96 * h2g_se.2 > 0) %>%


# plot rg v.s. se of rg
ggplot(subtype_rg,aes(x = serg.1, y = rg.1, label=phenotype)) +
  geom_point( size = 1) +
  geom_label_repel(aes(label=ifelse(abs(rg.1)< 0.3 & serg.1 < 0.3,as.character(phenotype),'')), max.overlaps = 50) +
  theme_bw() 

subtype_rg <- subtype_data %>%
  filter(cluster_number == 2, rg.1 > 0.5, serg.1 < 0.5) %>% # filter rg > than target
  mutate(youngh2g = if_else(age.1 < age.2, h2g.1, h2g.2), 
         oldh2g= if_else(age.1 < age.2, h2g.2, h2g.1),
         youngh2g_se = if_else(age.1 < age.2, h2g_se.1 , h2g_se.2 ), 
         oldh2g_se = if_else(age.1 < age.2, h2g_se.2 , h2g_se.1 ))
ggplot(subtype_rg,aes(x= oldh2g, y = youngh2g,  label=phenotype)) +
  geom_pointrange(aes(ymin=youngh2g - youngh2g_se, ymax=youngh2g + youngh2g_se, color = rg.1), size = 1) +
  geom_pointrange(aes(xmin=oldh2g - oldh2g_se, xmax=oldh2g + oldh2g_se, color = rg.1), size = 1) +
  scale_color_gradient2(low= "blue", mid="grey90",  high="red") +
  geom_label_repel(aes(label=ifelse(abs(youngh2g - oldh2g)> 0.1,as.character(phenotype),'')), max.overlaps = 50) +
  geom_abline(slope = 1) + 
  theme_bw() 

subtype_3tp <- subtype_data%>%
  filter(cluster_number == 3) %>%
  mutate(youngh2g = if_else(age.1 < age.2, h2g.1, h2g.2), 
         oldh2g= if_else(age.1 < age.2, h2g.2, h2g.1),
         youngh2g_se = if_else(age.1 < age.2, h2g_se.1 , h2g_se.2 ), 
         oldh2g_se = if_else(age.1 < age.2, h2g_se.2 , h2g_se.1 ))

ggplot(subtype_3tp,aes(x= oldh2g, y = youngh2g,  label=phenotype)) +
  geom_pointrange(aes(ymin=youngh2g - youngh2g_se, ymax=youngh2g + youngh2g_se, color = rg.1), size = 1) +
  geom_pointrange(aes(xmin=oldh2g - oldh2g_se, xmax=oldh2g + oldh2g_se, color = rg.1), size = 1) +
  scale_color_gradient2(low= "blue", mid="grey90",  high="red") +
  geom_label_repel(aes(label=ifelse(abs(youngh2g - oldh2g)> 0.1,as.character(phenotype),'')), max.overlaps = 50) +
  geom_abline(slope = 1) + 
  theme_bw() 




write.csv(select(subtype_rg, phenotype, phe1h2g, phe2h2g, corr_lmm), "Association_analysis/Disease_correlation_form.csv", row.names =   F)  

#############################################
# specific disease examples
#############################################
subtype_diseases_list <- read.csv(file = "subtype_disesea_list.csv")
# breast cancer v.s. malignant neoplasm
order_incidences <- order(apply(para$unlist_zn[para$ds_list[[j]]$id,], 1, which.max))
breast_cancer <- 17
malignant_neoplasm_breast <- 18
Sample_BC <- para$unlist_Ds_id[para$ds_list[[breast_cancer]]$id,]$eid
length(Sample_BC)
Sample_MNB <- para$unlist_Ds_id[para$ds_list[[malignant_neoplasm_breast]]$id,]$eid
length(Sample_MNB)
length(intersect(Sample_BC, Sample_MNB))



########################################################
# extract the GCP
########################################################
DIR <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/causality_analysis/LCV_results/"
temp <- list.files(paste(DIR, sep=""), pattern=".RData")
LCV_data <- data.frame(trait1 = as.numeric(), trait2 = as.numeric(),topic1= as.numeric(), 
                       topic2 = as.numeric(), gcp= as.numeric(),gcp_se= as.numeric(),log10P= as.numeric(),
                       rho.zscore = as.numeric(), rho.est = as.numeric())
for(i in 1:length(temp)){
  load(paste0(DIR, temp[i]))
  row_data <- LCV_save[[1]]
  row_data$rho.zscore <- LCV_save[[2]]$rho.est/LCV_save[[2]]$rho.err
  row_data$rho.est <- LCV_save[[2]]$rho.est
  LCV_data <- LCV_data %>% 
    add_row(row_data)
}
# load other variables
phe_phecode <- read.csv("info_phenotype_phecode.csv")
h2g_disease_topics <- read.csv("causality_analysis/h2g_causality.csv")
trait_number <- read.csv("trait_number.csv") 
filtered_disease <- h2g_disease_topics %>% 
  filter(z_score > 4)
filtered_disease <- filtered_disease %>%
  left_join(select(phe_phecode, phecode, phenotype), by = c("disease" = "phecode"))
LCV_data <- LCV_data %>%
  left_join(select(filtered_disease, disease, topic, age, phenotype), by = c("trait1" = "disease", "topic1" = "topic"))
LCV_data <- LCV_data %>%
  left_join(select(filtered_disease, disease, topic, age,phenotype), by = c("trait2" = "disease", "topic2" = "topic"))

plot_data <- LCV_data %>%
  filter(abs(rho.zscore) > 0,  gcp_se < 10) %>%
  mutate(trait1 = paste0("Topic:", topic1, " ",phenotype.x ) ,
         trait2 = paste0("Topic:", topic2, " ",phenotype.y ) )
ggplot(plot_data, aes(x = trait1, y = trait2)) + 
  geom_tile(aes(fill=rho.est)) + 
  # scale_fill_gradientn(colors =  c(blue, "white", "white", red), limits = c(-1, NA)) +
  scale_fill_gradientn(colors =  c(blue, "white", red), limits = c(-1, NA)) +
  labs(x=NULL, y=NULL, title="rg estimates") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=90),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

plot_data <- LCV_data %>%
  filter(abs(rho.zscore) > 1.96,  abs(rho.est) < .9) %>%
  left_join(select(trait_number, disease, topic, case_number), by = c("trait1" = "disease", "topic1" = "topic")) %>%
  mutate(trait1 = paste0("Topic:", topic1, " ",phenotype.x ) ,
         trait2 = paste0("Topic:", topic2, " ",phenotype.y ) )
fdr_gcp <- p.adjust(10^(plot_data$log10P), method = "fdr")
plot_data$fdr_gcp <- fdr_gcp
plot_data <- plot_data %>% 
  filter(fdr_gcp < 0.1)
ggplot(plot_data, aes(x = trait1, y = trait2)) + 
  geom_tile(aes(fill=gcp)) + 
  scale_fill_gradientn(colors =  c(blue, "white", "white", red), limits = c(-1, NA)) +
  labs(x=NULL, y=NULL, title="GCP estimates") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=90),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

# for a specific disease, plot all its causally related diseases
source("plotting_functions.R")
pal_tree <- colorRampPalette(c(red, grey))
pal_tree_vector <- rev(pal_tree(100))

######################
# plot a target topic 
######################
topic_id <- 5
# Keep a network that have relatively high GCP, and only one direction
edges <- plot_data  %>% 
  rename(source = trait1, target = trait2) %>%
  filter(gcp > 0.6) %>%
  mutate(age_gap <- (age.x - age.y)) %>%
  filter(topic1 == topic_id | topic2 == topic_id)

######################
# plot a target disease 
######################
# disease_id <- "Diaphragmatic hernia"  
disease_id <- "Type 2 diabetes"  
# disease_id <- "Essential hypertension" 
# Keep a network that have relatively high GCP, and only one direction
edges <- plot_data  %>% 
  rename(source = trait1, target = trait2) %>%
  filter(gcp > 0) %>%
  mutate(age_gap <- (age.x - age.y)) %>%
  filter(phenotype.x == disease_id | phenotype.y == disease_id)

filter_nodes <- data.frame(diag_icd10 = unique(c(as.character(edges$source), as.character(edges$target)) ) )  %>%
  left_join(select(plot_data, trait1, case_number), by = c("diag_icd10" = "trait1")) %>%
  group_by(diag_icd10) %>%
  slice(1)
comorbidity_graph <- graph_from_data_frame(d=edges, vertices=filter_nodes, directed=T) 
E(comorbidity_graph)$arrow.mode <- 2
E(comorbidity_graph)$width <- edges$gcp
E(comorbidity_graph)$color <- pal_tree_vector[floor(edges$gcp * 100 + 1)]
V(comorbidity_graph)$size <- sqrt(V(comorbidity_graph)$case_number)/20
V(comorbidity_graph)$color <- pal_tree_vector[pmin(floor(V(comorbidity_graph)$case_number/max(filter_nodes$case_number) * 100) + 1, 100)]
V(comorbidity_graph)$frame.color <- NA
# pdf(file = paste0("~/Desktop/comorbidity/figures/GCP_full_network.pdf"), width = 8, height = 8)
# layout: layout_as_tree, or others e.g. layout.fruchterman.reingold
plot(comorbidity_graph,  layout = layout.fruchterman.reingold, edge.arrow.size = .2, vertex.label.color="black",vertex.label.cex=.5, vertex.label.dist=0.2, vertex.label.font = 2, edge.curved = 0)
# dev.off()



# plot T2D and hypertension
plot_data <- LCV_data %>%
  filter(trait1 %in% c(401.10, 250.20), trait2 %in% c(401.10, 250.20)) %>%
  filter(trait1 != trait2 ) %>%
  mutate(trait1 = paste0("Topic:", topic1, " ",phenotype.x ) ,
         trait2 = paste0("Topic:", topic2, " ",phenotype.y ) )

plt <- ggplot(plot_data, aes(x = trait1, y = trait2)) + 
  geom_tile(aes(fill=gcp)) + 
  scale_fill_gradientn(colors =  c(blue, "white", "white", red), limits = c(-1, NA), na.value = grey) +
  labs(x=NULL, y=NULL, title="rg estimates") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=90),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

# compare the gcp from same topic v.s. different topic
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- h2g_disease_topics %>%
  left_join(phe_phecode, by = c("disease" = "phecode") )

cardio_disese <- LCV_data %>%
  filter(trait1 > 400, trait1 < 460, trait2 > 400, trait2 < 460) %>%
  mutate(trait1 = paste0("Topic:", topic1, " ",phenotype.x ) ,
         trait2 = paste0("Topic:", topic2, " ",phenotype.y ) )
ggplot(cardio_disese, aes(x = trait1, y = trait2)) + 
  geom_tile(aes(fill=abs(rho.est))) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x=NULL, y=NULL, title="GCP estimates") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=90),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11)) 

diff_topic <- LCV_data %>%
  filter(topic1 != topic2, trait1 != trait2, rho.zscore > 2, abs(rho.est) < .9,gcp_se < 0.3)
ggplot(diff_topic, aes(label=paste0("trait 1:",phenotype.x," trait 2:", phenotype.y))) +
  geom_point(aes(x= sort(-log10(runif(dim(diff_topic)[1]))), y = sort(-log10P)) , size = 1) +
  geom_abline(slope = 1) + 
  theme_bw() 

############################################
# extract diseaes set for prediction task
############################################
# prediction update: we should predict the risk between 50-60; 60-70 & > 70 years old

########################################################################
# 2022-10-10: compute rg using imputed UK BB genotypes 
########################################################################
DIR <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/imputed_SNP_results/"
temp <- list.files(paste(DIR, sep=""), pattern=".RData")
LCV_data <- data.frame(trait1 = as.numeric(), trait2 = as.numeric(),topic1= as.character(), 
                       topic2 = as.character(), gcp= as.numeric(),gcp_se= as.numeric(),log10P= as.numeric(),
                       rho.zscore = as.numeric(), rho.est = as.numeric(), rho.err = as.numeric())
for(i in 1:length(temp)){
  load(paste0(DIR, temp[i]))
  row_data <- LCV_save[[1]]
  row_data$rho.zscore <- LCV_save[[2]]$rho.est/LCV_save[[2]]$rho.err
  row_data$rho.est <- LCV_save[[2]]$rho.est
  row_data$rho.err <- LCV_save[[2]]$rho.err
  LCV_data <- LCV_data %>% 
    add_row(row_data)
}
# load other variables
plot_data <- LCV_data %>%
  filter(abs(rho.zscore) > 0,  gcp_se < 10) 


# extract all the non-subtype disease
non_subtp <- LCV_data %>%
  filter(topic1 == "all" , 
         topic2 == "all" ) %>%
  select(-topic1, -topic2) %>%
  rename(non_tp_gcp = gcp, non_tp_gcp_se = gcp_se, non_tp_log10P = log10P,
         non_tp_rho.zscore = rho.zscore, non_tp_rho.est = rho.est, non_tp_rho.err=rho.err)
subtp_rg <- LCV_data %>%
  filter(topic1 != "all", topic2 != "all") %>%
  left_join(non_subtp, by = c("trait1", "trait2")) %>%
  mutate(diff.zscore = (non_tp_rho.est - rho.est)/sqrt(rho.err^2 + non_tp_rho.err^2))

# plot the big matrix with all diseases
disease_order <- subtp_rg %>%
  group_by(topic1, trait1) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(ds_tp =  paste0(topic1,"_", trait1)) 

plot_rg <- subtp_rg %>%
  mutate(ds_tp1 =  factor(paste0(topic1,"_", trait1), levels = disease_order$ds_tp, ordered = TRUE), 
         ds_tp2 =  factor(paste0(topic2,"_", trait2), levels = disease_order$ds_tp, ordered = TRUE) ) %>%
  mutate(non_tp_rho.zscore = if_else(is.infinite(non_tp_rho.zscore),  # 
                                     max(subtp_rg$rho.zscore[!is.infinite(subtp_rg$rho.zscore)], na.rm = T),
                                     non_tp_rho.zscore)) %>%
  mutate(rg = if_else(ds_tp1 > ds_tp2, non_tp_rho.est, rho.est - non_tp_rho.est), 
         rg.zscore.abs = if_else(ds_tp1 > ds_tp2, abs(non_tp_rho.zscore), abs(diff.zscore)),
         shape_type = if_else(ds_tp1 > ds_tp2, "lower", "upper")) %>%
  mutate(P = if_else(ds_tp1 < ds_tp2, (1-pnorm(abs(diff.zscore)))*2, 9)) %>%
  filter(! ds_tp2 == ds_tp1, !is.na(rg), !is.na(P)) %>% 
  #  ds_tp1 > ds_tp2 | trait1 %in% subtypes_ds$disease | trait2 %in% subtypes_ds$disease ) %>%
  mutate(ds_tp2 =  factor(paste0(topic2,"_", trait2), levels = rev(disease_order$ds_tp), ordered = TRUE) )

within_disease <- plot_rg %>%
  filter(trait1 == trait2, topic1 != topic2) %>%
  filter(P < 1)

within_disease$FDR <-  p.adjust(within_disease$P, method = "fdr")

within_disease %>% filter(FDR < 0.1)
