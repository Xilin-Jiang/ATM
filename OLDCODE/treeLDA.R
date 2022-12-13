# tree LDA
library(topicmodels)
library(dplyr)
library(ggplot2)
require(survival)
library(stringr)
library(tidyverse)
library(gtools)
library(maxLik)

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

##########################
# load the packages
##########################
source("topic_functions.R")
rec_data <- read.csv("DiseaseAbove500occur_include_death_ICDA2N.csv")
ds_list <- read.csv("listAbove500include_deaths_ICDA2N.csv")
topic_num <- 10
tree_leveles <- 4
para <- topic_init_tree(rec_data, ds_list, topic_num, tree_leveles)


############################
# start optimization
############################
# set the number of update 
para$max_itr <- 500
para$lb <- data.frame("Iteration" = 0,"Lower_bound" = comp_lda_lb(para))
para$alpha <- rep(.1, para$K)
para$tol <- 10^(-6) # tolerance of lower bound step 
para$itr_check_lb <- 5
for(itr in 1:para$max_itr){
  print(paste0("Interation: ",itr))
  para <- CVB0_E_zn(para)
  para <- update_beta_treeLDA(para)
  # print(paste0("Current Lower bound at beta ", comp_lda_lb(para)))
  # para <- update_alpha(para) # we do not need to update alpha
  # print(paste0("Current Lower bound at alpha ", comp_lda_lb(para)))
  if(mean(para$alpha) < 1){ # in this case, abandone the approximated lower bound, using the Mean-field lower bound.
    para <- comp_E_lntheta(para)
    para$lb[nrow(para$lb) + 1,] <- c(itr, comp_lda_lb(para))
  }else{
    para$lb[nrow(para$lb) + 1,] <- c(itr, CVB_lb(para))
  }
  if(itr %% para$itr_check_lb  ==0){
    print(paste0("Current Lower bound ", pull(filter(para$lb, Iteration == itr), Lower_bound), " at iteration: ",itr))
    curr_lb <- pull(filter(para$lb, Iteration == itr), Lower_bound) 
    prev_lb <- pull(filter(para$lb, Iteration == (itr -para$itr_check_lb )), Lower_bound) 
    if(is.finite((curr_lb - prev_lb)) & abs(curr_lb - prev_lb)/abs(prev_lb) < para$tol ){
      print(paste0("Optimization converged at step ", itr))
      break
    }
  }
}

# save the parameter
save(para, file = paste0("~/Desktop/comorbidity/Results/","Run_treeLDA_",Sys.Date(),".RData"))

# compare the lower bound 
load("~/Desktop/comorbidity/Results/Run_treeLDA_2021-02-12.RData")
lb_tree <- para$lb
load("~/Desktop/comorbidity/Results/Run_2021-02-10.RData")
lb_flat <- para$lb
df_log_lb <- data_frame(number_iteration = seq(0,300,5), lb_tree = lb_tree$Lower_bound[seq(1,301,5)], lb_flat = lb_flat$Lower_bound[1:61])
plt <- ggplot(df_log_lb) + 
  geom_line(aes(x =number_iteration, y = lb_tree, color = "TreeLDA lowerbound")) + 
  geom_line(aes(x = number_iteration, y = lb_flat, color = "FlatLDA lowerbound")) + 
  # geom_point(aes(x = number_topics, y = constant_effect_at_50, color = "lowerbound without age effect"), size = 1) + 
  # geom_point(aes(x = number_topics, y = age_effect_at_30, color = "lower bound fitted using age effect"), size = 1) +
  scale_color_manual(name="Model Type", values=c("TreeLDA lowerbound" = red, "FlatLDA lowerbound"  = green)) +
  labs(x = "iteration time", y = "Lower bound of log likelihood")
ggsave("~/Desktop/comorbidity/figures/lower_bound_tree_flat.png",plt,width = 6, height =4) 

# visualization:
for(topic_id in 1:para$K){
  vsual_matrix <- matrix(0, nrow = para$L, ncol = para$D)
  for(l in 1:para$L){
    for(j in 1:para$D){
      vsual_matrix[l,j] <- para$treeBeta[[l]][para$NodeDS[l,j],topic_id]
    }
  }
  library(reshape2)
  library(ggplot2)
  longData<-melt(vsual_matrix)
  longData<-longData[longData$value!=0,]
  
  plt <- ggplot(longData, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill=value)) + 
    scale_fill_gradient(low="grey90", high="red") +
    labs(x="levels", y="disease", title=paste0("Tree topics ", topic_id) ) +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))
  ggsave(paste0("~/Desktop/comorbidity/figures/tree_topics",topic_id,".png"), 
         plt, width = 6, height = 4)
}

library(reshape2)
library(ggplot2)
topic_id <- 1
plot_tree <- function(para, topic_id, keep_tree_threshold){
  visual_matrix <- matrix(0, nrow = para$L, ncol = para$D)
  for(l in 1:para$L){
    # put a value at each beta
    visual_matrix[l,] <- apply( as.matrix( sapply(1:l, function(l) para$treeBeta[[l]][para$NodeDS[l,],topic_id,drop = F]) ),
                               1, prod)

  }
  # filter visual_matrix over the threshold
  normalised_tree <- apply(visual_matrix, 1, function(x) x/max(x)) %>% t
  library(reshape2)
  library(ggplot2)
  longData<-melt(normalised_tree)
  longData<-longData[longData$value!=0,]
  
  plt <- ggplot(longData, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill=value)) + 
    scale_fill_gradient(low="grey90", high="red") +
    labs(x="levels", y="disease", title=paste0("Tree topics ", topic_id) ) +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))
  ggsave(paste0("~/Desktop/comorbidity/figures/tree_topics",topic_id,".png"), 
         plt, width = 6, height = 4)
} 

################################################################################
# below we will try to load the tree into a igraph network for visualization
################################################################################
load("~/Desktop/comorbidity/Results/Run_treeLDA_2021-02-16.RData")
library(igraph)
# there are three metrics to compute: 1. NodeDS index matrix of (L)-by-D: elements refereing parent for disease j at layer l
# 2. a disease list: treeDS[[l]][[cl]] the records associated with each node, tree equivalent of para$ds_list (using the list of diseases code from 1:para$D that are under that node)
# 3. the beta list: treeBeta[[l]] is the vector (matrix to include i) of betas at this layer 
# 4. para$Cl[[l]] save the number of nodes at each layer
# using ZZZZ for the death events
filled_list <- as.character(para$list_above500occu$diag_icd10)
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
  png(filename = paste0("~/Desktop/comorbidity/figures/Topics/tree_topics",topic_id,"_",Sys.Date(),".png"), width = 3200, height = 2400)
  plot(significant_icdtree, vertex.label.color="black",vertex.label.cex=2.5, vertex.label.dist=.8, vertex.label.font = 2, edge.curved = 0, layout = layout_as_tree)
  dev.off()
}


