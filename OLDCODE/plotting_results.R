library(topicmodels)
library(dplyr)
library(ggplot2)
require(survival)
library(stringr)
library(tidyverse)
library(gtools)
library(maxLik)
library(reshape2)
require(umap)

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

################################################################################
# below we will try to load the tree into a igraph network for visualization
################################################################################
load("~/Desktop/comorbidity/Results/REP_3_Yidong_data_Run_treeLDA_2021-03-11.RData")
library(igraph)
# there are three metrics to compute: 1. NodeDS index matrix of (L)-by-D: elements refereing parent for disease j at layer l
# 2. a disease list: treeDS[[l]][[cl]] the records associated with each node, tree equivalent of para$ds_list (using the list of diseases code from 1:para$D that are under that node)
# 3. the beta list: treeBeta[[l]] is the vector (matrix to include i) of betas at this layer 
# 4. para$Cl[[l]] save the number of nodes at each layer
# using ZZZZ for the death events
filled_list <- para$list_above500occu$diag_icd10
# make sure death is separate from the tree
if("Death" %in% filled_list){
  filled_list[which(filled_list == "Death")] <- "ZZZZ"}
for(topic_id in 1:10){
  sources <- list()
  ends <- list()
  edges <- list()
  nodes_set <- list()
  
  
  visual_matrix <- matrix(0, nrow = para$L, ncol = para$D)
  for(l in 1:para$L){
    # put a value at each beta
    visual_matrix[l,] <- apply( as.matrix( sapply(1:l, function(l) para$treeBeta[[l]][para$NodeDS[l,],topic_id,drop = F]) ),
                                1, prod)
    
  }
  
  for(l in 1:para$L){
    if(l == para$L){
      nodes_l <- unlist(lapply(filled_list, function(x) (substring(x, 1,l+4)))) # longest ICD10 is 7 character
    }else{
      nodes_l <- unlist(lapply(filled_list, function(x) (substring(x, 1,l)))) 
    }
    node_set <- unique(nodes_l)
    if(l == 1){
      sources[[1]] <- rep("root", length(node_set))
      ends[[1]] <- node_set
    }else{
      sources[[l]] <- unlist(lapply(node_set, function(x) (substring(x, 1,l-1))))
      ends[[l]] <- node_set
    }
    # edges[[l]] <- data.frame(source = sources[[l]], target = ends[[l]], beta = para$treeBeta[[l]][,topic_id]) # directly using beta
    edges[[l]] <- data.frame(source = sources[[l]], target = ends[[l]], beta = visual_matrix[l,match(1:para$Cl[[l]], para$NodeDS[l,])]) # using the product of multinomials
    nodes_set[[l]] <- node_set
  }
  
  pal_tree <- colorRampPalette(c(red, grey))
  pal_tree_vector <- rev(pal_tree(100))
  
  edges <- bind_rows(edges) %>% # remove edges that have same source and target
    filter(! as.character(source) == as.character(target) )
  filter_edges <- edges %>%
    filter(beta > 10/para$D)
  # filter_nodes <- unique(c(as.character(filter_edges$source), as.character(filter_edges$target)) )
  filter_nodes <- filter_edges %>%
    select(target, beta) %>%
    add_row(target = "root", beta = 1)
  
  significant_icdtree <- graph_from_data_frame(d=filter_edges, vertices=filter_nodes, directed=T) 
  E(significant_icdtree)$arrow.mode <- 0
  E(significant_icdtree)$width <- filter_edges$beta*100
  E(significant_icdtree)$color <- pal_tree_vector[floor(filter_edges$beta/max(filter_edges$beta) * 100)]
  V(significant_icdtree)$size <- sqrt(V(significant_icdtree)$beta) * 10
  V(significant_icdtree)$color <- pal_tree_vector[pmin(floor(V(significant_icdtree)$beta/max(filter_edges$beta) * 100), 100)]
  V(significant_icdtree)$frame.color <- NA
  png(filename = paste0("~/Desktop/comorbidity/figures/tree_topics",topic_id,".png"), width = 3200, height = 2400)
  plot(significant_icdtree, vertex.label.color="black",vertex.label.cex=2.5, vertex.label.dist=.8, vertex.label.font = 2, edge.curved = 0, layout = layout_as_tree)
  dev.off()
}
################################################################################
# plot the convergence analysis for Yidong's Gibbs and VB
################################################################################
# saving the third layer separately
topics_yidong_gibbs <- read.table("UKB_data_yidong/11.phi.ave.txt", header = T, sep = ",")

# normalise Yidong's topics
# topics_yidong_gibbs <- apply(topics_yidong_gibbs,1, function(x) x/sum(x)) %>% t

library(Rtsne)
source("plotting_functions.R")
# load the betas and the lb
lb_lst <- list()
beta_lst <- list()
for(rep_id in 1:5){
  load(paste0("~/Desktop/comorbidity/Results/REP_",rep_id,"_Yidong_data_Run_treeLDA_2021-03-11.RData"))
  lb_lst[[rep_id]] <- para$lb
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
  beta_lst[[rep_id]] <- third_layer_topics
}
# compare the convergence
covrg <- matrix(NA, nrow = 5, ncol = 100)
# order beta_lst by lower bound 
lb_cvg <- rep(NA, 5)
for(rep_id in 1:5){
  lb_cvg[rep_id] <- lb_lst[[rep_id]] %>% 
    slice_tail %>% 
    pull(2)
  covrg[rep_id, ] <- lb_lst[[rep_id]] %>% 
    slice_tail(n=100) %>% 
    pull(2)
}
df_lb <- data_frame(itr = 1:100, lower_bound = t(covrg))
ggplot(df_lb) + 
  geom_line(aes(x = itr, y = lower_bound[,1])) + 
  geom_line(aes(x = itr, y = lower_bound[,2])) + 
  geom_line(aes(x = itr, y = lower_bound[,3])) + 
  geom_line(aes(x = itr, y = lower_bound[,4])) + 
  geom_line(aes(x = itr, y = lower_bound[,5])) 
# order the topics according to hierarchical clustering results
ordered_lst <- beta_lst[order(lb_cvg)[5:3]]
ordered_lst[[4]]  <- unname(data.matrix(topics_yidong_gibbs[2:11,, drop = F])) # exclude the first topic which is a healthy topic
beta_for_clustering <- do.call(rbind, ordered_lst)
# create labels for each run
runs_vb_gibs <- floor((1:dim(beta_for_clustering)[1]-1)/10) + 1
label_for_each_row <- sapply(runs_vb_gibs, function(x) ifelse(x == 4, "Gibbs", paste0("VB",x)))

ddrg <- dendro_plot(beta_for_clustering, label_for_each_row)
dg_topics <- ddrg[[2]]
plt <- ddrg[[2]]
ggsave("~/Desktop/comorbidity/figures/consine_dist_Yidong.png",plt,width = 6, height =6) 

# Extract the order of the tips in the dendrogram
beta.order <- order.dendrogram(dg_topics)
beta.long <- melt(beta_for_clustering[beta.order, ]) %>%
  rename("Disease" = Var2 , "Topics" = Var1, "weight" = value)

# Create heatmap plot
heatmap.plot <- ggplot(data = beta.long, aes(x = Disease, y = Topics)) +
  geom_tile(aes(fill = weight)) +
  scale_fill_gradient2() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top")

# All together
grid.newpage()
par(cex=5)
print(heatmap.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.43, width = 0.2, height = 0.92))


# below provide a handle for the cluster membership of each topic
beta_for_clustering %>%
  dist() %>%
  clust() %>% 
  cutree(k = 10)
################################################################################
# plot the lower bound distribution for curve degree of freedom 
################################################################################
rep_number <- 10
maxP <- 6
df_lb_P_K <- data_frame(df_P = as.integer(),df_K = as.integer(), lower_bound = as.numeric())
for(K in c(5:20, 25, 30, 35, 40, 45, 50) ){
  lb_lst <- list()
  for(df_P in 2:maxP){
    for(rep_id in 1:rep_number){
      try({# load(paste0("~/Desktop/comorbidity/Results/strongRegularization_model_output_PheCode_age_dependent_K", K,"_P",df_P, "_rep",rep_id,".RData"))
         load(paste0("~/Desktop/comorbidity/Results/rec2CVB0_model_output_PheCode_age_dependent_K", K,"_P",df_P, "_rep",rep_id,".RData"))
      cvrg_lb <-  model_output[[2]] %>% 
        filter(!is.infinite(Lower_bound)) %>%
        slice_tail %>% 
        pull(2)
      df_lb_P_K  <- df_lb_P_K %>% 
        add_row(df_P = df_P, df_K = K, lower_bound = cvrg_lb)
      # if(cvrg_lb < - 9 * 10^(6)){ # avoid those old files which computed the wrong lower bound
      #   df_lb_P_K  <- df_lb_P_K %>% 
      #     add_row(df_P = df_P, df_K = K, lower_bound = cvrg_lb)}
      })
    }
  }
}

# plot a box plot for selecting best topic number and degree of freedom
df_boxplot <- df_lb_P_K %>%
  # filter(df_K  < 21, df_K > 4) %>% # help visualise one or two Ks
  mutate(df_P = as.character(df_P), df_K = as.factor(df_K)) 
  
plt <- ggplot(data=df_boxplot,aes(x=df_K, y=lower_bound, fill = df_P)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width=.8,position=position_dodge(width=0.85)) +
  geom_point(position=position_jitterdodge(jitter.width = .02, dodge.width=0.85), size = 0.3, alpha=0.6) +
  # geom_jitter(width=0.1, size = 0.2, alpha=0.4) + 
  labs(x = "Number of topics", y = "Lower bound") + 
  scale_fill_manual(values=cbPalette[2:7]) + 
  theme(panel.background=element_blank()) 
ggsave(paste0("~/Desktop/comorbidity/figures/CVB0_lb_change_with_topic_number_degree_freedom.png"), plt, width = 10, height = 4)

##################################
# plot some age topics 
##################################
source("plotting_functions.R")
# first find the best rep
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

topic_id <- 7
para$D <- model_output[[3]]
para$list_above500occu <- read.csv("listAbove1000_PheCode.csv") 
# focus on 30-81 years old
dominant_ds_id <- which(sapply(1:para$D, function(j) mean(model_output[[1]][30:81,j,topic_id]) > 10/para$D))
para$list_above500occu$diag_icd10[dominant_ds_id]
filtered_topics <- model_output[[1]][30:81,dominant_ds_id,topic_id]

longData<-melt(filtered_topics)
longData<-longData[longData$value!=0,, drop = F]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

para$D <- model_output[[3]]
pal_age <- colorRampPalette(c(blue, red, orange, purple, green))
pal_age_vector <- pal_age(para$D)
thre_pick <- 10/para$D
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- para$list_above500occu %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") )
icd10 <- ds_list$phenotype
# using the larger results file: ordering are changed in "model_output"
betas <- para$pi_beta_basis
for(topic_id in 1:K){
  trajs <- betas[30:80,,topic_id] # trajectories
  plot_title <- paste0("Topic: ", topic_id)
  plt <- plot_age_topics(icd10, trajs, thre_pick, pal_age_vector, plot_title, start_age = 30)
  ggsave(paste0("~/Desktop/comorbidity/figures/topics/rep", rep_id, "K",K,"P",df_P,"age_topics",topic_id,".png"), 
         plt, width = 8, height = 8)
}


################################################################################
# clustering of topics of age profiles for different topic numbers
################################################################################
## file order_topic_output.R
# order the topics by posterior
args <- commandArgs(trailingOnly = TRUE) # args[1] defines the number of topics; args[2] defines degree of freedom for the curves; args[3] is the repitation id
K <- as.numeric(args[1]) 
P <- as.numeric(args[2]) # degrees of freedom
rep_id <- args[3]

var <- load("Results/Run_A2N_age_dependent_K",K,"_P",P,"_rep",rep_id, ".RData")
order_by_post_z <- order(colMeans(para$alpha_z), decreasing = T)
ordered_pi_beta <- para$pi_beta_basis[,,order_by_post_z]
model_output <- list(ordered_pi_beta, para$lb, para$D, para$M, para$K, para$P)
save(model_output, file = paste0("Results/","model_output_A2N_age_dependent_K",K,"_P",P,"_rep",rep_id, ".RData"))

# Step1 (on the cluster): order the topics by posterior distribution in the population 
# Step2: perform clustering and identify those topics that are robust
# repeat this procedure for all Ks and Ps, look at the pattern, is there a best number of cluster? 
# repeat the procedure with normal LDA package, see if normal LDA in the package is consistent
beta_lst <- list()
K <- 10
P <- 5 # quadratic polynomial
DIR <- "~/Desktop/comorbidity/Results/"
# pt <- paste0("^strongRegularization_model_output_PheCode_age_dependent_K", K,"_P", P, "_rep*")
# pt <- paste0("^rec2CVB0_model_output_A2N_age_dependent_K", K,"_P", P, "_rep*")
pt <- paste0("^rec2CVB0_model_output_PheCode_age_dependent_K", K,"_P", P, "_rep*")
#pt <- paste0("^model_output_A2N_age_dependent_K", K,"_P", P, "_rep*")
temp <- list.files(paste(DIR, sep=""), pattern=pt)
covrg <- matrix(NA, nrow = length(temp), ncol = 1)
for(rep_id in 1:length(temp)){
  load(paste0(DIR, temp[rep_id]))
  covrg[rep_id, ] <- model_output[[2]] %>% 
    filter(!is.infinite(Lower_bound)) %>%
    slice_tail %>% 
    pull(2)
  # collapse the topic matrices into vectors for comparison
  beta_lst[[rep_id]] <- matrix(model_output[[1]][30:80,,], prod(dim(model_output[[1]][30:80,,])[1:2]), dim(model_output[[1]])[3]) %>% t
}
# order the runs by final lower bound; the topics are ordered by average posterior distribution in the population
ordered_lst <- beta_lst[order(covrg, decreasing = T)[1:5]] #[1:length(temp)]]
beta_for_clustering <- do.call(rbind, ordered_lst)
plt <- cosine_dist_plot(beta_for_clustering)

ggsave(paste0("~/Desktop/comorbidity/figures/cosine_similarity_K",K,"P",P,"age_topics",topic_id,".png"), 
       plt[[1]], width = 4, height = 4)

# perform hierarchical clustering
runs_vb_gibs <- floor((1:dim(beta_for_clustering)[1]-1)/K) + 1
label_for_each_row <- sapply(runs_vb_gibs, function(x) paste0("VB",x))

ddrg <- dendro_plot(beta_for_clustering, label_for_each_row)
dg_topics <- ddrg[[2]]
plt <- ddrg[[2]]

# perform kmeans
sum_variance <- matrix(nrow = 20)
sum_variance[1] <- (nrow(beta_for_clustering)-1)*sum(apply(beta_for_clustering,2,var))
for (i in 2:20){
  sum_variance[i] <- sapply(1:3, function(x) sum(kmeans(beta_for_clustering, centers=i)$withinss)) %>%
    mean()
}
# plot the kmeans changes for choice of clusters 
plot(1:19, (sum_variance[1:19] - sum_variance[2:20]), type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

K_best <- 15 # or 15
fit <- kmeans(beta_for_clustering, K_best)

# visualise the kmeans

# Dimension reduction using PCA
res.pca <- prcomp(beta_for_clustering,  scale = TRUE)
summary(res.pca)
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(res.km$cluster)
# Add Species groups from the original data sett
ind.coord$Species <- df$Species



# first checking if matrix and array will turn not be of same size
a <- matrix(model_output[[1]], prod(dim(model_output[[1]])[1:2]), dim(model_output[[1]])[3]) %>% t
stopifnot(model_output[[1]] == array(t(a), dim = dim(model_output[[1]])))

# get cluster means 
centers_ds <- array(t(fit$centers), dim = dim(model_output[[1]])) 
# plot these topics
num_D <- dim(centers_ds)[2] # number of diseases
pal_age <- colorRampPalette(c(blue, red, orange, purple, green))
pal_age_vector <- pal_age(num_D)
thre_pick <- 20/num_D
icd10 <- read.csv("listAbove500include_deaths_ICDA2N.csv") %>% pull(1)
for(topic_id in 1:K){
  trajs <- centers_ds[30:80,,topic_id] # trajectories
  plt <- plot_age_topics(icd10, trajs, thre_pick, pal_age_vector, start_age = 30)
  ggsave(paste0("~/Desktop/comorbidity/figures/topics/kmean_centers_K",K,"P",df_P,"age_topics",topic_id,".png"), 
         plt, width = 4, height = 4)
}




#######################################################################
# 2021-07-02: what are those non-consistent topics? 
# Are they formed by the individuals who only have one records?
###########################################################
# first step: find all those topics that are not consistent. 
beta_lst <- list()
K <- 11
P <- 5 # quadratic polynomial
DIR <- "~/Desktop/comorbidity/Results/"
pt <- paste0("^rec2CVB0_model_output_A2N_age_dependent_K", K,"_P", P, "_rep*")
temp <- list.files(paste(DIR, sep=""), pattern=pt)
covrg <- matrix(NA, nrow = length(temp), ncol = 1)
for(rep_id in 1:length(temp)){
  load(paste0(DIR, temp[rep_id]))
  covrg[rep_id, ] <- model_output[[2]] %>% 
    filter(!is.infinite(Lower_bound)) %>%
    slice_tail %>% 
    pull(2)
  # collapse the topic matrices into vectors for comparison
  beta_lst[[rep_id]] <- matrix(model_output[[1]][30:80,,], prod(dim(model_output[[1]][30:80,,])[1:2]), dim(model_output[[1]])[3]) %>% t
}
# order the runs by final lower bound; the topics are ordered by average posterior distribution in the population
ordered_lst <- beta_lst[order(covrg, decreasing = T)[1:2]] #[1:length(temp)]]
beta_for_clustering <- do.call(rbind, ordered_lst)
rslt <- cosine_dist_plot(beta_for_clustering)
plt <- rslt[[1]]
C <- rslt[[2]]

#####################################
# get the prediction likelihoood
#####################################

rep_number <- 10
maxP <- 7
df_predict_lik_P_K <- data_frame(df_P = as.integer(),df_K = as.integer(), likelihood = as.numeric())
for(K in 5:20){
  lb_lst <- list()
  for(df_P in 2:maxP){
    for(rep_id in 1:rep_number){
      try({load(paste0("../Results/","prediction_age_K",K,"_P",df_P,"_rep",rep_id, ".RData"))
        df_predict_lik_P_K  <- df_predict_lik_P_K %>% 
            add_row(df_P = df_P, df_K = K, likelihood = predictions[[1]]$test_log_lik)
        })
    }
  }
}

# plot a box plot for selecting best topic number and degree of freedom
df_boxplot <- df_predict_lik_P_K %>%
  filter(df_K  < 21, df_K > 4) %>% # help visualise one or two Ks
  mutate(df_P = as.character(df_P), df_K = as.factor(df_K)) 

plt <- ggplot(data=df_boxplot,aes(x=df_K, y=likelihood, fill = df_P)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width=.8,position=position_dodge(width=0.85)) +
  geom_point(position=position_jitterdodge(jitter.width = .02, dodge.width=0.85), size = 0.3, alpha=0.6) +
  # geom_jitter(width=0.1, size = 0.2, alpha=0.4) + 
  labs(x = "Number of topics", y = "Lower bound") + 
  scale_fill_manual(values=cbPalette[2:7]) + 
  theme(panel.background=element_blank()) 
ggsave(paste0("~/Desktop/comorbidity/figures/prediction_likelihood_with_topic_number_degree_freedom.png"), plt, width = 10, height = 4)

# check if the algorithms simply predict high incidence diseases
incidence <- rec_data %>% 
  group_by(diag_icd10) %>%
  summarise(n()) %>% 
  arrange(desc(`n()`))
sum(incidence$`n()`[1:floor(para$D/10)])/dim(rec_data)[1]
  
########################################################
# analysis: compute the age assignments of each topic
########################################################
load("../Results/Run_2rec_A2N_age_dependent_K11_P5_rep2.RData")
# compute the average age distribution of diseases
age_distribution <- sapply(1:max(para$unlist_Ds_id$age_diag), function(x) colSums(para$unlist_zn[para$unlist_Ds_id$age_diag == x, ,drop = F]) )
# normalise the age distribution
age_distribution <- sapply(1:para$K, function(x) age_distribution[x,]/sum(age_distribution[x,]))
mean_age <- colSums(sapply(1:para$K, function(x) age_distribution[,x] * (1:dim(age_distribution)[1])))

########################################################
# estimate ontology heterogeneity in topics over age
########################################################
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
ds_list <- read.csv("listAbove1000include_deaths_PheCode.csv")
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode") )
assignment_per_ds <- sapply(1:para$D, function(j) colMeans(para$unlist_zn[para$ds_list[[j]]$id,]) ) %>% t
# then compute the average age for each loading
age_topic_assoc_per_ds <- matrix(NA, nrow = para$D, ncol = para$K)
for(j in 1:para$D){
  loadings <- para$unlist_zn[para$ds_list[[j]]$id,]
  sum_loadings <- colSums(loadings)
  loadings <- sweep(loadings, 2, sum_loadings, FUN = "/")
  age <- para$unlist_Ds_id[para$ds_list[[j]]$id,]$age_diag
  age_topic_assoc_per_ds[j,] <- sapply(1:para$K, function(i) age %*% loadings[,i])
}
j <- 4
plot(age_topic_assoc_per_ds[j,],assignment_per_ds[j,])

# plot the assignments over age
ds_id <- 174.11  # 39, 30
j <- match(ds_id, para$list_above500occu$diag_icd10)
age_topic_matrix <- matrix(NA, nrow = 81, ncol = para$K)
for(ag in 31:81){
  ids <- para$ds_list[[j]] %>% 
    filter(age_diag == ag) %>%
    pull(3)
  if(length(ids)){
    age_topic_matrix[ag,] <- colMeans(para$unlist_zn[ids,,drop = F])
  }
}
longData<-melt(age_topic_matrix)
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Topics", y="Age", title=paste0("ICD10: ", para$list_above500occu$diag_icd10[j])) +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

# testing for heterogeneity: null multinomial: the probability of all ages
# test for the first half of age and second half of age

##########################################
# using k-means and do permutation test to identify the set of disease with heterogeneity
##########################################
rs_kmean <- function(loadings, center_num){
  kfit <- kmeans(loadings, centers=center_num)
  return(kfit$betweenss)
}

heterogeneity_age_test <- matrix(NA, nrow = para$D, ncol = 4)
for(j in 1:para$D){
  print(paste0("disease id: ", j))
  loadings <- para$unlist_zn[para$ds_list[[j]]$id,]
  tot_sz <- dim(para$unlist_zn)[1]
  target_sz <- dim(loadings)[1]
  
  perm_sz <- 1000 # how many permutation samples to collect
  stats_target <- matrix(0, nrow = 5, ncol = 1) 
  permutate_stats <- matrix(0, nrow = 5, ncol = perm_sz) 
  for(k_num in 2:5){
    stats_target[k_num] <- rs_kmean(loadings, k_num) # get the target test statistic
    
    permutate_stats[k_num, ] <- sapply(1:perm_sz, function(x) 
      rs_kmean(para$unlist_zn[sample(1:tot_sz, target_sz), ], k_num) )
    heterogeneity_age_test[j,(k_num - 1)] <- (1 + sum((stats_target[k_num] - stats_target[k_num-1]) < 
                                                        (permutate_stats[k_num, ] - permutate_stats[k_num-1, ]) ))/perm_sz
  }

}

plot(-log(p.adjust(heterogeneity_age_test[,1], method = "fdr")))

qqdf <- data.frame(p_true = sort(-log(heterogeneity_age_test[,1])), p_sim = sort(- log( runif(length(heterogeneity_age_test))) ))  
ggplot(qqdf) + 
  geom_point(aes(x = p_sim, y = p_true))

df_pvalues <- data.frame(para$list_above500occu$diag_icd10, heterogeneity_age_test)
write.csv(df_pvalues, "test_disease_heterogeneity_pvalues.csv")

top_ds <- which(rank(log(p.adjust(heterogeneity_age_test, method = "fdr"))) <= 50)
write.table(para$list_above500occu$diag_icd10[top_ds], "top50_hererogeneous_disease.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)

# # the test below is replaced by a better one above. 
# library(stats)
# heterogeneity_age_test <- matrix(NA, nrow = para$D, ncol = 2)
# for(j in 1:para$D){
#   loadings <- para$unlist_zn[para$ds_list[[j]]$id,]
#   mean_loadings <- colMeans(loadings)
#   ids_first_half <- para$ds_list[[j]] %>% 
#     filter(age_diag <= 55) %>%
#     pull(3)
#   ids_second_half <- para$ds_list[[j]] %>% 
#     filter(age_diag > 55) %>%
#     pull(3)
#   first_half <- as.integer(colSums(para$unlist_zn[ids_first_half,,drop = F]))
#   second_half <- as.integer(colSums(para$unlist_zn[ids_second_half,,drop = F]))
#   if(length(ids_first_half) & length(ids_second_half)){
#     p_first_half <- chisq.test(x = first_half, p = mean_loadings)$p.value
#     p_second_half <- chisq.test(x = second_half, p = mean_loadings)$p.value
#   }else{
#     p_first_half <- 1
#     p_second_half <- 1
#   }
#   
#   heterogeneity_age_test[j,] <- c(p_first_half, p_second_half)
# }
# 
# plot(-log(p.adjust(heterogeneity_age_test[,2], method = "fdr")))
# 
# # plot a qq plot
# qqdf <- data.frame(p_true = sort(-log(heterogeneity_age_test[,2])), p_sim = sort(- log( runif(length(heterogeneity_age_test[,2]))) ))  
# ggplot(qqdf) + 
#   geom_point(aes(x = p_sim, y = p_true))
# 
# # plot each disease assignments over age
# top_ds <- which(rank(log(p.adjust(heterogeneity_age_test[,2], method = "fdr"))) <= 50)
# write.table(para$list_above500occu$diag_icd10[top_ds], "top50_hererogeneous_disease.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)

##########################################
# save a keep for 2 records individuals to reduce file size (for plink)
##########################################
rec_data <- read.csv("rec2subjectAbove500occur_include_death_ICDA2N.csv")
individuals <- rec_data %>%
  group_by(eid) %>%
  summarise() %>%
  select(eid) %>%
  mutate(fid = eid)
write.table(individuals, paste0("keep2DiseaseIndividuals.txt"), sep="\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#####################################
# compute the comorbidity profiles
#####################################
topic_num <- 10
degree_free <- 5
set.seed(19940110)
ds_list <- read.table("top50_hererogeneous_disease.txt")
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("V1" = "phecode"))
phecodeList <- para$list_above500occu %>%
  left_join(phe_phecode, by = c("diag_icd10" = "phecode"))

# find the best convergence lower bound
DIR <- "~/Desktop/comorbidity/Results/"
pt <- paste0("^rec2CVB0_model_output_PheCode_age_dependent_K", topic_num,"_P", degree_free, "_rep*")
temp <- list.files(paste(DIR, sep=""), pattern=pt)
covrg <- matrix(NA, nrow = length(temp), ncol = 1)
for(rep_id in 1:length(temp)){
  load(paste0(DIR, temp[rep_id]))
  covrg[rep_id, ] <- model_output[[2]] %>% 
    filter(!is.infinite(Lower_bound)) %>%
    slice_tail %>% 
    pull(2)
}
# order the runs by final lower bound; the topics are ordered by average posterior distribution in the population
repid <- order(covrg, decreasing = T)[1]
load(paste0("../Results/Run_2rec_PheCode_age_dependent_K",topic_num,"_P",degree_free,"_rep",repid, ".RData"))

k_best_age_test <- matrix(NA, nrow = length(ds_list$V1), ncol = 4)
for(id in 1:length(ds_list$V1)){
  ds_id <- as.character(ds_list$V1[id])
  j <- match(ds_id, para$list_above500occu$diag_icd10)
  loadings <- para$unlist_zn[para$ds_list[[j]]$id,]
  longData<-melt(loadings)
  longData<-longData[longData$value!=0,]
  
  ggplot(longData, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill=value)) + 
    scale_fill_gradient(low="grey90", high="red") +
    labs(x="Topics", y="incidence", title=paste0("Topic assignments: ", ds_list$phenotype[id])) +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))
  # using k-means cluster to identify cluster numbers
  # sum_variance <- matrix(nrow = 10)
  # sum_variance[1] <- (nrow(loadings)-1)*sum(apply(loadings,2,var))
  # for(i in 2:10){
  #   sum_variance[i] <- sapply(1:5, function(x) kmeans(loadings, centers=i)$betweenss) %>%
  #     mean()
  # }
  # # plot the kmeans changes for choice of clusters 
  # plot(1:10, c(sum_variance[1], (sum_variance[1] - sum_variance[2:10])), type="b", xlab="Number of Clusters",
  #      ylab="Within groups sum of squares")
    # using a permutation strategy for k_best
  tot_sz <- dim(para$unlist_zn)[1]
  target_sz <- dim(loadings)[1]
  perm_sz <- 1000 # how many permutation samples to collect
  permutate_samples <- sapply(1:perm_sz, function(x)
    sample(1:tot_sz, target_sz))
  stats_target <- sapply(2:5, function(k_num) rs_kmean(loadings, k_num)) # get the target test statistic
  permutate_stats <- sapply(1:perm_sz, function(sp) 
    sapply(2:5, function(k_num) rs_kmean(para$unlist_zn[permutate_samples[,sp], ], k_num) ) )
  
  stats_target <- stats_target - c(0, stats_target[1:3])
  permutate_stats <- permutate_stats - rbind(rep(0, perm_sz), permutate_stats[1:3,])
  for(k in 1:4){ # four tests here
    k_best_age_test[id,k] <- (1 + sum(stats_target[k] < permutate_stats[k,]))/perm_sz
  }
  K_best <- length(which(k_best_age_test[id,] < 0.005)) + 1
  # get the elbow point: use the ratio of newly explained variance, with 1/K as threshold (K is number of topics available)
  # K_best <- which((sum_variance[2:10] - c(0,sum_variance[2:9]))/(sum_variance[1:9] < 1/para$K)[1] 
  print(paste0(ds_list$V1[id]," ",ds_list$phenotype[id], " number of clusters ", K_best ))
  kmfit <- kmeans(loadings, K_best)
  longData<-melt(kmfit$centers)
  longData<-longData[longData$value!=0,]
  
  ggplot(longData, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    scale_fill_gradient(low="grey90", high="red") +
    labs(x="Topics", y="cluster", title=paste0("ICD10: ", ds_id)) +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))
  
  # using the centers to get the disease profile
  cases_eid <- list()
  for(k in 1:K_best){
    cases_id_k <- para$ds_list[[j]]$id[which(kmfit$cluster == k)]
    cases_eid[[k]] <- para$unlist_Ds_id[cases_id_k,] %>%
      mutate(disease_group = k) %>%
      mutate(Ds_id = para$list_above500occu[Ds_id, 1])
  }
  cases_eid <- bind_rows(cases_eid)
  save(cases_eid, file = paste0("../association_results/","ds",ds_id,"_group_K",para$K,"_P",para$P,"_rep",para$rep_ID, ".RData"))
}

#########################################################
####### now compute the disease comorbidity for each group 
#########################################################
# within each individual that has this disease (assigned to clusters, we compute the weight as the cosine similarity between the weight of specific disease and the target disease)
for(id in 1:length(ds_list$V1)){
  ds_id <- as.character(ds_list$V1[id])
  j <- match(ds_id, para$list_above500occu$diag_icd10)
  topic_num <- 10
  degree_free <- 5
  load(paste0("../association_results/","ds",ds_id,"_group_K",topic_num,"_P",degree_free,"_rep",repid, ".RData"))
  K_best <- max(cases_eid$disease_group)
  loadings <- para$unlist_zn[para$ds_list[[j]]$id,]
  kmfit <- kmeans(loadings, K_best)
  
  # plot histogram for each group
  plt <- ggplot(cases_eid,aes(x = age_diag, fill = as.factor(disease_group), colour = as.factor(disease_group))) +
    geom_histogram(alpha = 0.5, position = "identity", binwidth = 1) + 
    scale_fill_manual(values = c(red, blue)) + 
    scale_color_manual(values = c(red, blue)) +
    theme(legend.position = "None") + 
    labs(x="Age (years)", y="incidence count", title=paste0("Disease: ", phecodeList$phenotype[j])) 
  ggsave(paste0("../figures/","age_distribution_",ds_id,"_group_K",topic_num,"_P",degree_free,"_rep",repid, ".png"), plt, width = 6, height = 4)
  
  
  profile_k <- list()
  for(k in 1:K_best){
    # find all individuals of this group
    cases_k <- cases_eid %>%
      filter(disease_group == k)
    ds_same_individual <- sapply(match(cases_k$eid, para$eid), 
                                 function(x) para$unlist_Ds_id[para$patient_lst[[x]],], 
                                 simplify = F) %>% 
      bind_rows() %>%
      mutate(id = row_number())
    zn_same_individual <- sapply(match(cases_k$eid, para$eid), 
                                 function(x) para$unlist_zn[para$patient_lst[[x]],], 
                                 simplify = F) 
    zn_same_individual <- do.call(rbind, zn_same_individual)
    # if there is no incidence, use 0
    profile_k[[k]] <- matrix(0, nrow = para$D, ncol = 1)
    target_z <- kmfit$centers[k,]
    for(ds in 1:para$D){
      ds_row_id <- ds_same_individual %>% 
        filter(Ds_id == ds)
      if(length(ds_row_id$id)){ # some disease might not have comorbidity
        cos_sim <- (zn_same_individual[ds_row_id$id,,drop = F] %*% target_z)/
          sqrt(rowSums(zn_same_individual[ds_row_id$id,,drop = F]^2) * sum(target_z^2))
        profile_k[[k]][ds] <- sum(cos_sim)
      }
    }
    # compute relative risk (divided by prevalence)
    profile_k[[k]] <- profile_k[[k]]/para$list_above500occu$occ
    # normalise within each profile
    profile_k[[k]] <- profile_k[[k]]/mean(profile_k[[k]], na.rm = T)
  }
  profile_k <- do.call(cbind, profile_k)
  longData <- melt(profile_k) %>% 
    filter(Var1 != j) %>%
    mutate(Var1 = as.character(phecodeList$phenotype[Var1])) 
  # threshold, number of records associated with the disease divided by number of distinct disease
  longData<-longData[longData$value > 5,]
  
  plt <- ggplot(longData, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill=value)) + 
    scale_fill_gradient(low="grey90", high="red") +
    labs(x="clusters", y="comorbidity", title=paste0("Phecode: ", para$list_above500occu$diag_icd10[j], " ", phecodeList$phenotype[j])) +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))
  ggsave(paste0("../figures/","ds",ds_id,"_group_K",topic_num,"_P",degree_free,"_rep",repid, ".png"), plt, width = 8, height = 6)
}
# compute the age profiles of each profile: it should contain the age profile of the disease itself as well.
age <- para$unlist_Ds_id[para$ds_list[[j]]$id,]$age_diag
mean(age[which(kmfit$cluster == 1)])
mean(age[which(kmfit$cluster == 2)])

####################################################################
# check if PRS score is significantly different within two subgroups
####################################################################
DIR <- "~/Desktop/comorbidity/association_results/"
ds_list <- read.table("top50_hererogeneous_disease.txt")
# playing with the code
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("V1" = "phecode"))
repid <- 2
topic_num <- 10
degree_free <- 5
regression_md <- list()
p_values <- matrix(nrow = length(ds_list$V1), ncol = 1)
data_percentile <- list()
for(id in 1:length(ds_list$V1)){
  ds_id <- as.character(ds_list$V1[id])
  # check the PRS data are indeed correct
  ds_eid <- rec_data %>% 
    filter(diag_icd10  == ds_id) %>%
    select(eid)
  PRS_profile <- read.table(paste0(DIR, ds_id, ".profile"), header =T)
  
  # get the disease of different groups
  load(file = paste0("../association_results/","ds",ds_id,"_group_K", topic_num,"_P",degree_free,"_rep",repid, ".RData"))
  K_best <- max(cases_eid$disease_group)
  eid_subgroups <- list()
  for(k in 1:K_best){
    eid_subgroups[[k]] <- cases_eid %>%
      filter(disease_group == k) %>%
      select(eid)
  }
  
  # plot percentile PRS with incidence rate
  percentile <- PRS_profile %>% mutate(percentiles = ntile(SCORE,100))
  prevalence <- rep(NA, 100)
  prevalence_subgroup <- list()
  prevalence_subgroup <-  matrix(NA, nrow = 100, ncol = K_best)
  for(pct in 1:100){
    pct_eid <- percentile %>%
      filter(percentiles == pct) %>%
      pull(1)
    prevalence[pct] <- mean(pct_eid %in% ds_eid$eid )
    for(k in 1:K_best){
      prevalence_subgroup[pct, k] <- mean(pct_eid %in% eid_subgroups[[k]]$eid ) 
    }
  }
  df_plt <- data.frame(percentile = 1:100, prevalence, prevalence_subgroup = prevalence_subgroup)
  data_percentile[[id]] <- df_plt
  # ggplot(df_plt) +
  #   geom_point(aes(x = percentile, y = prevalence), colour = grey)
  regression_df <- cases_eid %>% 
    left_join(select(PRS_profile, FID, SCORE), by = c("eid" = "FID"))
  try({
    regression_md[[id]] <- glm(SCORE ~ disease_group, data = regression_df)
    p_values[id] <- summary(regression_md[[id]])$coefficients[2,4] 
    print(paste(ds_id, " has a p-value ", p_values[id], " and ", K_best, " clusters"))
  })
}
ds_list$p_true <- -log(p_values)
ds_list <- ds_list %>%
  arrange(p_true)
ds_list$p_sim = sort(- log( runif(length(p_values))) )
ggplot(data = ds_list) + 
  geom_point(aes(x = p_sim, y = p_true), color = grey, size = 2) + 
  geom_abline(slope = 1) + 
  geom_label_repel(aes(x = p_sim, y = p_true, label=ifelse(p_true> 3,as.character(phenotype),''))) 
  
ds_list$phenotype[which(p_values < 0.05)] 
# do two plotting for each disease: one age of clustering, one cluster centers
  


#################################################################################
# play with data: include more people and then plot the number of unique types
#################################################################################
all_eid <- rec_data %>%
  group_by(eid) %>%
  summarise() 
num_uniq <- c()
for(num_eid in (1:100)*10){
  rec_unique <- all_eid %>% 
    sample_n(size  = num_eid) %>%  
    left_join(rec_data, by = "eid")
  num_uniq <- c(num_uniq,length(unique(rec_unique$diag_icd10)))
}
plot_df <- data_frame(sample_size =  (1:100)*10, unique_rec = num_uniq)
ggplot(plot_df) + 
  geom_point(aes(x = sample_size, y = unique_rec))









#########################################
# plot Yidong's Phewas
#########################################
p_values <- read.csv("result-all.new.anno.txt", header = T)
ggplot(p_values,aes(x = BETA, y = -log(P))) + 
  geom_point()+
  geom_label_repel(aes(label=ifelse( (-log(P) > 100 & abs(BETA) > 0.01) ,as.character(Field),''))) +
  theme_bw() 

##########################################
# EM for mixture of dirichlet distribution
##########################################
source("topic_functions.R")
set.seed(19940110)
# test power calibration under a simple dirichlet sample 
LRT <- rep(NA, 100)
BIC1 <- rep(NA, 100)
BIC2 <- rep(NA, 100)
for(i in 1:100){
  loadings1 <-  rdirichlet(1000, c(0.1,0.1,3,0.1,0.1)) 
  parMixDir1 <- fit_MixDir(loadings1, 1)
  parMixDir2 <- fit_MixDir(loadings1, 2)
  LRT[i] <- -2 * ( parMixDir1$lb$Lower_bound[dim(parMixDir1$lb)[1]] -
                     parMixDir2$lb$Lower_bound[dim(parMixDir2$lb)[1]] )
  BIC1[i] <- (parMixDir1$M - 1)*log(dim(loadings1)[1]) -2 *  parMixDir1$lb$Lower_bound[dim(parMixDir1$lb)[1]]
  BIC2[i] <- (2*parMixDir2$M - 1)*log(dim(loadings1)[1]) -2 *  parMixDir2$lb$Lower_bound[dim(parMixDir2$lb)[1]]
}
degree_free <- parMixDir1$M # the degree of freedom is only M-1 + pi which is 1
df_calibration <- data.frame(p_sim = sort(runif(100)) , p_true = sort(1 - pchisq(LRT, df = degree_freedom)))
ggplot(df_calibration) + 
  geom_point(aes(x = p_sim, y = p_true)) + 
  geom_abline(slope = 1)

x <- runif(1000)
sp2 <- rdirichlet(1000, c(0.1,0.1,0.1,0.1,3))
loadings2 <- rdirichlet(1000, c(0.1,0.1,3,0.1,0.1))
loadings2[x < 0.2, ] <- sp2[x < 0.2, ]
parMixDir <- fit_MixDir(loadings2, 2)

# test the calibration using the simulated topics
set.seed(19940110)
LRT <- rep(NA, para$D)
BIC1 <- rep(NA, para$D)
BIC2 <- rep(NA, para$D)
for(i in 1:para$D){
  loadings1 <-  para$unlist_zn[para$ds_list[[j]]$id,]
  parMixDir1 <- fit_MixDir(loadings1, 1)
  parMixDir2 <- fit_MixDir(loadings1, 2)
  LRT[i] <- -2 * ( parMixDir1$lb$Lower_bound[dim(parMixDir1$lb)[1]] -
                     parMixDir2$lb$Lower_bound[dim(parMixDir2$lb)[1]] )
  BIC1[i] <- (parMixDir1$M - 1)*log(dim(loadings1)[1]) -2 *  parMixDir1$lb$Lower_bound[dim(parMixDir1$lb)[1]]
  BIC2[i] <- (2*parMixDir2$M - 1)*log(dim(loadings1)[1]) -2 *  parMixDir2$lb$Lower_bound[dim(parMixDir2$lb)[1]]
}

j1 <- 4
j2 <- 13
LRT <- rep(NA, 1)
for(i in (1:10)){
  loadings1 <-  para$unlist_zn[para$ds_list[[j1]]$id,]
  loadings2 <-  para$unlist_zn[para$ds_list[[j2]]$id,][1:(i*20),]
  loadings <- rbind(loadings1, loadings2)
  parMixDir1 <- fit_MixDir(loadings, 1)
  parMixDir2 <- fit_MixDir(loadings, 2)
  
  LRT[i] <- -2 * ( parMixDir1$lb$Lower_bound[dim(parMixDir1$lb)[1]] -
                     parMixDir2$lb$Lower_bound[dim(parMixDir2$lb)[1]] )
}
degree_free <- parMixDir1$M # the degree of freedom is only M-1 + pi which is 1
df_power <- data.frame(proportion_heterogeneity = 20*(1:100)/(dim(loadings1)[1] + 20*(1:100)) , p_true = (1 - pchisq(LRT, df = degree_freedom)))
ggplot(df_power) + 
  geom_line(aes(x = proportion_heterogeneity, y = p_true)) 

degree_free <- parMixDir1$M # the degree of freedom is only M-1 + pi which is 1
df_calibration <- data.frame(p_sim = sort(runif(length(LRT))) , p_true = sort(1 - pchisq(LRT, df = degree_freedom)))
ggplot(df_calibration) + 
  geom_point(aes(x = p_sim, y = p_true)) + 
  geom_abline(slope = 1)

