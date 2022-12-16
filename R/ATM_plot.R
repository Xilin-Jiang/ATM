# plot topic loadings
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
