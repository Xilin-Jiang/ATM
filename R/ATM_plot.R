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
# plot topic loadings
#' Title plot the topic loadings across age.
#'
#' @param disease_names the list of disease names, ordered as the topic.
#' @param trajs one disease topic, which should be a matrix of age-by-disease.
#' @param plot_title the title of the figure.
#' @param start_age starting age of the matrix, default 30.
#' @param top_ds How many disease to show, default is 10. This will filter the disease by
#' the average topic laodings across age and pick the top.
#'
#' @return a ggplot object of the topic loading.
#' @export
#'
#' @examples disease_list <- UKB_349_disease %>%
#' left_join(disease_info_phecode_icd10, by = c("diag_icd10"="phecode" )) %>%
#' pull(phenotype)
#' topic_id <- 1 # plot the first topic
#' plot_age_topics(disease_names = disease_list,
#'         trajs = UKB_HES_10topics[30:80,,topic_id],
#'         plot_title = paste0("topic ", topic_id),
#'         top_ds = 7)
plot_age_topics <- function(disease_names, trajs,  plot_title = "", start_age = 30, top_ds = 10){
  pal_age <- colorRampPalette(c(blue, red, orange, purple, green))
  pal_age_vector <- pal_age(length(disease_names))
  # dominant_ds_id <- which(sapply(1:dim(trajs)[2], function(j) mean(trajs[,j], na.rm = T) > thre_pick))
  dominant_ds_id <- order(sapply(1:dim(trajs)[2], function(j) mean(trajs[,j], na.rm = T)), decreasing = T)[1:top_ds]
  # if nothing pass the threshold
  if(!(length(dominant_ds_id) > 0 ))  return(ggplot())
  filtered_topics <- trajs[,dominant_ds_id, drop = F]
  print(paste0("Number of diseases selected: ",length(dominant_ds_id)))
  df_topics <- data.frame(age=start_age:(start_age + dim(trajs)[1] - 1), inferred_topics =  filtered_topics)
  legend_xpos <- apply(filtered_topics, 2, which.max)
  legend_ypos <- sapply(1:length(dominant_ds_id), function(x) filtered_topics[legend_xpos[x],x])
  legend <- data.frame( ds_label = as.character(disease_names[dominant_ds_id]),
                        x_pos = legend_xpos + start_age ,
                        y_pos = legend_ypos )
  col_ds <- pal_age_vector[dominant_ds_id]
  if(length(dominant_ds_id) == 1 ){
    plt <- ggplot(data = df_topics, aes(x = age)) +
      ggrepel::geom_label_repel(data = legend, aes(x = x_pos, y = y_pos, label = ds_label,fontface = "bold"),color = col_ds ) + # , vjust = "inward", hjust = "inward") +
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
      ggrepel::geom_label_repel(data = legend, aes(x = x_pos, y = y_pos, label = ds_label,fontface = "bold", fontsize = 3),size = 5 ,color = col_ds)  # , vjust = "inward", hjust = "inward") +

  }
  return(plt)
}
