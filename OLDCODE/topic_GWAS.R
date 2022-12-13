#####################################################
# create topic loading file as covariates for GWAS
#####################################################
# this file contains code for performing genetic analysis on the topic loadings 
# load(file = paste0("~/Desktop/comorbidity/Results/REP_3_Yidong_data_Run_treeLDA_2021-03-11.RData"))
# another version with alpha = .1
load(file = paste0("~/Desktop/comorbidity/Results/REP_3_Yidong_data_Run_treeLDA_2021-04-03.RData"))
posti_theta <- para$alpha_z/rowSums(para$alpha_z)
df_loadings <- data.frame(eid = para$eid, fid = para$eid, theta = posti_theta) 
df_loadings %>%
  write.table(paste0("~/Desktop/comorbidity/Results/alpha",para$alpha[1],"_loading_K",para$K,"_P",para$P,".txt"), sep="\t", col.names = FALSE, row.names = FALSE)

# save the third layer topics for Yidong
# load the betas and the lb
rep_id <- 3 # save the third VB run since it has the highest lower bound.
load(paste0("~/Desktop/comorbidity/Results/REP_",rep_id,"_Yidong_data_Run_treeLDA_2021-04-03.RData"))
filled_list <- para$list_above500occu$diag_icd10
third_layer_topics <- matrix(NA, para$K, 100)
for(topic_id in 1:para$K){
    nodes_set <- list()
    visual_matrix <- matrix(0, nrow = para$L, ncol = para$D)
    for(l in 1:para$L){
      # put a value at each beta
      visual_matrix[l,] <- apply( as.matrix( sapply(1:l, function(l) para$treeBeta[[l]][para$NodeDS[l,],topic_id,drop = F]) ),
                                  1, prod)
      
    }
    nodes_l <- unlist(lapply(filled_list, function(x) (substring(x, 1,3)))) 
    node_set <- unique(nodes_l)
    
    third_layer_topics[topic_id,] <- visual_matrix[3,match(1:para$Cl[[3]], para$NodeDS[3,])]
}
names(third_layer_topics) <- node_set
third_layer_topics %>%
  write.table(paste0("~/Desktop/comorbidity/Results/alpha",para$alpha[1],"topic_for_Yidong.txt"), sep="\t", col.names = T, row.names = T)
################################################################################
# compute the age correlation of the loadings (using the posterior of z)
################################################################################
# first perform this over top 10 popoular diseases 

list_above500occu <- read.csv(paste0("listAbove500.csv")) 

lst_ds_high <- list_above500occu %>% # get the most prevelant disease
  filter(occ > 10000)


load("~/Desktop/comorbidity/Results/Run_treeLDA_2021-02-16.RData")
code2id <- function(x){
  return( match(x, para$list_above500occu$diag_icd10))
}


for(id in 1:length(lst_ds_high$diag_icd10)){
  ds_id <- code2id(lst_ds_high$diag_icd10[id]) 
  
  # here I am rounding the disease time to year for computation efficiency
  para$unlist_Ds_id <- first_icidence_age %>%
    mutate(Ds_id = code2id(diag_icd10)) %>%
    select(-diag_icd10) %>%
    mutate(age_diag = round(age_diag)) 
  
  df_Ez_n_zn <- data.frame(age = para$ds_list[[ds_id]]$age_diag, poste_zn = para$unlist_zn[para$ds_list[[ds_id]]$id, ]) %>%
    group_by(age) %>%
    summarise_all(mean) %>%
    arrange(age)
  
  library(reshape2)
  library(ggplot2)
  
  longData<-melt(as.matrix(df_Ez_n_zn)[,2:11])
  longData<-longData[longData$value!=0,]
  
  plt <- ggplot(longData, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill=value)) + 
    scale_fill_gradient(low="grey90", high="red") +
    labs(x="topics", y="age", title=as.character(lst_ds_high$diag_icd10[id])) +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))
  
  ggsave(paste0("~/Desktop/comorbidity/figures/",as.character(lst_ds_high$diag_icd10[id]),"age_topic_distribution.png"), 
         plt, width = 4, height = 4)
  
}

# for all disease, 
df_Ez_n_zn <- data.frame(age = para$unlist_Ds_id$age_diag, poste_zn = para$unlist_zn) %>%
  group_by(age) %>%
  summarise_all(mean) %>%
  arrange(age)
library(reshape2)
library(ggplot2)

longData<-melt(as.matrix(df_Ez_n_zn)[,2:11])
longData<-longData[longData$value!=0,]

plt <- ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="topics", y="age", title="posterior of all topics") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
ggsave(paste0("~/Desktop/comorbidity/figures/all_disease_age_topic_distribution.png"), 
       plt, width = 4, height = 4)
