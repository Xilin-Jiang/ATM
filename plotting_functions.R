library(dplyr)
library(ggplot2)
library(stringr)
library(tidyverse)
library("ggrepel")
library("grid")
library("dendextend")
library("ggdendro")
library(reshape2)
library(igraph)
library(colorBlindness)
library(gridExtra) 

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
# set the commonly used color
blue <- cbPalette[6]
red <- cbPalette[7]
green <- cbPalette[4]
orange <- cbPalette[2]
grey <- cbPalette[1]
yellow <- cbPalette[5]
purple <- cbPalette[8]
skyblue <- cbPalette[3]
##############################
# functions related to treeLDA
##############################

#########################################
# functions for clustering of topics
#########################################
cos.sim <- function(idx, X){
  A = X[idx[1],]
  B = X[idx[2],]
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}
cosine_dist_plot <- function(beta_for_clustering){
  # compute cosine similarity of the two matrices
  n <- nrow(beta_for_clustering) 
  cmb <- expand.grid(i=1:n, j=1:n) 
  
  C <- matrix(apply(cmb,1,function(x) cos.sim(x, beta_for_clustering)),n,n)
  
  longData<-melt(C)
  longData<-longData[longData$value!=0,]
  plt <- ggplot(longData, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill=value)) + 
    scale_fill_gradient(low="grey90", high="red") +
    labs(x="topics", y="topics", title="cosine similarity") +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))
  return(list(plt, C))
}

dendro_plot <- function(beta_for_clustering, label_for_each_row){
  beta_for_clustering <- apply(beta_for_clustering,1, function(x) x/sum(x)) %>% t
  dg_topics <- beta_for_clustering %>%
    dist() %>%
    hclust() %>%
    as.dendrogram()
  # change the labels to runs
  labels(dg_topics) <- label_for_each_row[order.dendrogram(dg_topics)]
  
  dendro.plot <- ggdendrogram(data = dg_topics, rotate = TRUE) + 
    theme(axis.text.y = element_text(size = 6))
  return(list(dg_topics, dendro.plot))
}
  


##############################
# functions related to age LDA
##############################
# first function plot
plot_age_topics <- function(icd10, trajs,  pal_age_vector, plot_title, start_age = 30, top_ds = 10){
  # dominant_ds_id <- which(sapply(1:dim(trajs)[2], function(j) mean(trajs[,j], na.rm = T) > thre_pick))
  dominant_ds_id <- order(sapply(1:dim(trajs)[2], function(j) mean(trajs[,j], na.rm = T)), decreasing = T)[1:top_ds]
  # if nothing pass the threshold
  if(!(length(dominant_ds_id) > 0 ))  return(ggplot())
  filtered_topics <- trajs[,dominant_ds_id, drop = F]
  print(paste0("Number of diseases selected: ",length(dominant_ds_id)))
  df_topics <- data.frame(age=start_age:(start_age + dim(trajs)[1] - 1), inferred_topics =  filtered_topics)
  legend_xpos <- apply(filtered_topics, 2, which.max)
  legend_ypos <- sapply(1:length(dominant_ds_id), function(x) filtered_topics[legend_xpos[x],x])
  legend <- data.frame( ds_label = as.character(icd10[dominant_ds_id]),
                        x_pos = legend_xpos + start_age ,
                        y_pos = legend_ypos )
  col_ds <- pal_age_vector[dominant_ds_id]
  if(length(dominant_ds_id) == 1 ){
    plt <- ggplot(data = df_topics, aes(x = age)) +
      geom_label_repel(data = legend, aes(x = x_pos, y = y_pos, label = ds_label,fontface = "bold"),color = col_ds ) + # , vjust = "inward", hjust = "inward") +
      theme_bw(base_size = 20) +
      geom_line(aes(x = age, y = inferred_topics), color = col_ds[1], size = 1.5) +
      labs(x="Age", y="Multinomial probability", title=plot_title)
  }else{
    plt <- ggplot(data = df_topics, aes(x = age)) +
      theme_bw(base_size = 20) +
      labs(x="Age", y="Multinomial probability", title=plot_title)
    for(line_id in 1:length(dominant_ds_id)){
      plt <- plt + 
        geom_line(aes_string(y = paste0("inferred_topics.",line_id) ), size = 1.5, color = col_ds[line_id])
    }
    plt <- plt + 
      geom_label_repel(data = legend, aes(x = x_pos, y = y_pos, label = ds_label,fontface = "bold", fontsize = 3),size = 5 ,color = col_ds)  # , vjust = "inward", hjust = "inward") +
      
  }
  return(plt)
}

# first function plot
plot_age_topics_specify_id <- function(icd10, trajs, dominant_ds_id, pal_age_vector, plot_title, start_age = 30){
  # dominant_ds_id <- which(sapply(1:dim(trajs)[2], function(j) mean(trajs[,j], na.rm = T) > thre_pick))
  # if nothing pass the threshold
  if(!(length(dominant_ds_id) > 0 ))  return(ggplot())
  filtered_topics <- trajs[,dominant_ds_id, drop = F]
  print(paste0("Number of diseases selected: ",length(dominant_ds_id)))
  df_topics <- data.frame(age=start_age:(start_age + dim(trajs)[1] - 1), inferred_topics =  filtered_topics)
  legend_xpos <- apply(filtered_topics, 2, which.max)
  legend_ypos <- sapply(1:length(dominant_ds_id), function(x) filtered_topics[legend_xpos[x],x])
  legend <- data.frame( ds_label = as.character(icd10[dominant_ds_id]),
                        x_pos = legend_xpos + start_age ,
                        y_pos = legend_ypos )
  col_ds <- pal_age_vector[dominant_ds_id]
  if(length(dominant_ds_id) == 1 ){
    plt <- ggplot(data = df_topics, aes(x = age)) +
      geom_label_repel(data = legend, aes(x = x_pos, y = y_pos, label = ds_label,fontface = "bold"),color = col_ds ) + # , vjust = "inward", hjust = "inward") +
      theme_bw(base_size = 20) +
      geom_line(aes(x = age, y = inferred_topics), color = col_ds[1], size = 1.5) +
      labs(x="Age", y="Multinomial probability", title=plot_title)
  }else{
    plt <- ggplot(data = df_topics, aes(x = age)) +
      theme_bw(base_size = 20) +
      labs(x="Age", y="Multinomial probability", title=plot_title)
    for(line_id in 1:length(dominant_ds_id)){
      plt <- plt + 
        geom_line(aes_string(y = paste0("inferred_topics.",line_id) ), size = 1.5, color = col_ds[line_id])
    }
    plt <- plt + 
      geom_label_repel(data = legend, aes(x = x_pos, y = y_pos, label = ds_label,fontface = "bold", fontsize = 3),size = 5 ,color = col_ds)  # , vjust = "inward", hjust = "inward") +
    
  }
  return(plt)
}



# function below recursively find all childrens
recurs_children <- function(node, l_x, l_y, edges,filter_nodes){
  x_pos <- l_x[which(filter_nodes$diag_icd10 == node)]
  y_pos <- l_y[which(filter_nodes$diag_icd10 == node)]
  child <- edges %>% 
    filter(source == node)
  child <- child$target
  if(length(child) > 0){
    cnt <- 0
    for(childnd in child){
      if(is.na(l_x[which(filter_nodes$diag_icd10 == childnd)])){
        child_pos <- which(filter_nodes$diag_icd10 == childnd)
        l_x[child_pos] <- x_pos + filter(edges, source == node, target == childnd)$age_gap
        # l_y[child_pos] <- y_pos + (floor(length(child)/2) - cnt)/pmin(length(l_y),pmax(1, max(l_y,na.rm = T) - min(l_y,na.rm = T)))
        l_y[child_pos] <- sum(!is.na(l_x))
        # also recurse through parents of children
        rslt <- recurs_parent(childnd, l_x, l_y, edges,filter_nodes)
        l_x <- rslt[[1]]
        l_y <- rslt[[2]]
      }
      rslt <- recurs_children(childnd, l_x, l_y, edges,filter_nodes)
      l_x <- rslt[[1]]
      l_y <- rslt[[2]]
      cnt <- cnt + 1
    }
    return(list(l_x, l_y))
  }else{
    return(list(l_x, l_y))
  }
}

# function below recursively find all parents
recurs_parent <- function(node, l_x, l_y, edges,filter_nodes){
  x_pos <- l_x[which(filter_nodes$diag_icd10 == node)]
  y_pos <- l_y[which(filter_nodes$diag_icd10 == node)]
  parent <- edges %>% 
    filter(target == node)
  parent <- parent$source
  if(length(parent) > 0){
    cnt <- 0
    for(parentnd in parent){
      if(is.na(l_x[which(filter_nodes$diag_icd10 == parentnd)])){
        parent_pos <- which(filter_nodes$diag_icd10 == parentnd)
        l_x[parent_pos] <- x_pos - filter(edges, source == parentnd, target == node)$age_gap
        # l_y[parent_pos] <- y_pos + (floor(length(child)/2) - cnt)
        l_y[parent_pos] <- sum(!is.na(l_x))
        # also recurse through children of parent
        rslt <- recurs_children(parentnd, l_x, l_y, edges,filter_nodes)
        l_x <- rslt[[1]]
        l_y <- rslt[[2]]
      }
      rslt <- recurs_parent(parentnd, l_x, l_y, edges,filter_nodes)
      l_x <- rslt[[1]]
      l_y <- rslt[[2]]
      cnt <- cnt + 1
    }
    return(list(l_x, l_y))
  }else{
    return(list(l_x, l_y))
  }
}
# use this to create a network object for comorbidity pairs
create_net <- function(matrix_year_diff, matrix_coocurre, age_gap_min,comorbidity_rate_thre){
  disease_pair <- melt(matrix_year_diff) %>%
    rename(age_gap = value) %>%
    left_join(melt(matrix_coocurre), by = c("Var1", "Var2")) %>%
    rename(comorbidity_rate = value) %>%
    filter(age_gap > age_gap_min, comorbidity_rate > comorbidity_rate_thre) %>%
    mutate(Disease1 = para$list_above500occu$diag_icd10[Var1], Disease2 = para$list_above500occu$diag_icd10[Var2]) %>%
    select(-Var1,-Var2)
  
  pal_tree <- colorRampPalette(c(red, grey))
  pal_tree_vector <- rev(pal_tree(100))
  
  edges <- disease_pair  %>% # remove edges that have same source and target
    rename(source = Disease2, target = Disease1) %>%
    select(source, target, age_gap, comorbidity_rate)
  
  filter_nodes <- data.frame(diag_icd10 = unique(c(as.character(edges$source), as.character(edges$target)) ) )  %>%
    left_join(para$list_above500occu, by = "diag_icd10")
  
  comorbidity_graph <- graph_from_data_frame(d=edges, vertices=filter_nodes, directed=T) 
  E(comorbidity_graph)$arrow.mode <- 2
  E(comorbidity_graph)$width <- edges$comorbidity_rate*10
  E(comorbidity_graph)$color <- pal_tree_vector[floor(edges$comorbidity_rate/max(edges$comorbidity_rate) * 100)]
  V(comorbidity_graph)$size <- sqrt(V(comorbidity_graph)$occ)/20
  V(comorbidity_graph)$color <- pal_tree_vector[pmin(floor(V(comorbidity_graph)$occ/max(filter_nodes$occ) * 100) + 1, 100)]
  V(comorbidity_graph)$frame.color <- NA
  return(comorbidity_graph)
}


################################
# GWAS results plotting
################################
# manhattan plot
plt_manhattan <- function(assoc.logistic, title){
  don <- assoc.logistic %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(assoc.logistic, ., by=c("CHR"="CHR")) %>%
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
    geom_label_repel(aes( label=ifelse(-log10(P) > 7.3,as.character(SNP),'')), max.overlaps = 50) +
    geom_hline(aes(yintercept = 7.3), color = "red", linetype = "dashed") +
    labs(x = "Chromosome number", y = "-log10(P)", title = title) +
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  return(plt)
}

