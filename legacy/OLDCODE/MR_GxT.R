# filter all diseases that are heritable 
h2g_disease_topics <- read.csv("causality_analysis/h2g_causality.csv")
h2g_disease_topics %>% 
  filter(z_score > 4) %>%
  group_by(disease) %>%
  slice(1) %>%
  pull(disease) %>%
  write.table("GxTopic_disease_list.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = F)

#######################################
# making the qq plot and manhattan plot
#######################################
# compute qq plot for D = G + T + G*T model
cc_dir <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/GxTopic/"
ds_list <- read.table("GxTopic_disease_list.txt")
# playing with the code
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("V1" = "phecode"))
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
for(idx in 1:length(ds_list$V1)){
  ds_id <- ds_list$V1[idx]
  print(ds_id)
  ds_name <- ds_list$phenotype[idx]
  for(ccgwas in 1:10){
    # only load non-sex chromosome
      gxtopic <- read.table(paste0(cc_dir, ds_id, "GxTopic_topic", ccgwas, ".gxtopic"), header = T) %>%
        mutate(logP = -log10(P)) %>%
        filter(!is.na(logP), CHR <= 22)  # only load non-sex chromosome
      df_qq.logistic <- gxtopic %>%
        arrange(logP)
      df_qq.logistic$sim_logP <- sort( - log10( runif(dim(df_qq.logistic)[1])))
      plt_qq <- ggplot(df_qq.logistic) + 
        geom_point(aes(x = sim_logP, y = logP), size = 0.5, alpha = 0.5) + 
        geom_abline(slope = 1, linetype = "dashed", color = red) +
        theme(legend.position = "none",panel.background=element_blank()) + 
        xlab(expression("Expected" * -log[10](P)) ) + ylab(expression("Observed" * -log[10](P))) +
        ggtitle(paste0(ds_name, " x topic: ", topic_name[ccgwas]))
      
      ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/GxT_simple_model/", ds_id, "_topic", 
                    topic_name[ccgwas],"qq_GxTopic.png"), plt_qq, width = 4, height = 4)
      
      # plot manhatten plot
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
      plt_mht <- ggplot(don, aes(x=BPcum, y=logP)) +
        # Show all points
        geom_point( aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
        scale_color_manual(values = rep(c(blue, grey), 22 )) +
        # custom X axis:
        scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
        scale_y_continuous(expand = c(0.1, 0) ) +     # remove space between plot area and x axis
        geom_hline(aes(yintercept = 7.3), color = "red", linetype = "dashed") +
        labs(x = "Chromosome number", y = expression(-log[10](P)), title = paste0(ds_name, " x ", topic_name[ccgwas])) +
        # Custom the theme:
        theme_bw() +
        theme( 
          legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
        )
      
      ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/GxT_simple_model/manhatten_",ds_id,"_topic", topic_name[ccgwas],".png"), plt_mht, width = 10, height = 4)
  }
}

# GxTopic for BMI
subtype <- read.table(paste0(cc_dir,"278.1topic_id_subtp.txt"))
for(ccgwas in 1:3){
  gxtopic <- read.table(paste0(cc_dir, "BMI_allcov_GxTopic_", ccgwas, ".assoc.linear"), header = T) %>%
    filter(TEST == "ADDxCOV1") %>%
    mutate(logP = -log10(P)) %>%
    filter(!is.na(logP), CHR <= 22) 
  df_qq.logistic <- gxtopic %>%
    arrange(logP)
  df_qq.logistic$sim_logP <- sort( - log10( runif(dim(df_qq.logistic)[1])))
  plt <- ggplot(df_qq.logistic) + 
    geom_point(aes(x = sim_logP, y = logP), size = 0.5, alpha = 0.5) + 
    geom_abline(slope = 1, linetype = "dashed", color = red) +
    theme(legend.position = "none",panel.background=element_blank()) + 
    xlab(expression("Expected" * -log[10](P)) ) + ylab(expression("Observed" * -log[10](P))) +
    ggtitle(paste0("BMI x topic: ", topic_name[subtype$V3[ccgwas]]))
  
  ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/GxT_simple_model/BMI_topic", 
                topic_name[subtype$V3[ccgwas]],"qq_GxTopic.png"), width = 4, height = 4,plt)
  
  # plot manhatten plot
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
  plt_mht <- ggplot(don, aes(x=BPcum, y=logP)) +
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
    scale_color_manual(values = rep(c(blue, grey), 22 )) +
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0.1, 0) ) +     # remove space between plot area and x axis
    geom_hline(aes(yintercept = 7.3), color = "red", linetype = "dashed") +
    labs(x = "Chromosome number", y = expression(-log[10](P)), title = paste0("BMI x ", subtype$V3[ccgwas])) +
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/GxT_simple_model/manhatten_BMI_topic", topic_name[subtype$V3[ccgwas]],".png"), plt_mht, width = 10, height = 4)
}
# also save the effect on BMI controlled the topic loading
subtype <- read.table(paste0(cc_dir,"278.1topic_id_subtp.txt"))
for(ccgwas in 1:3){
  gxtopic <- read.table(paste0(cc_dir, "BMI_allcov_GxTopic_", ccgwas, ".assoc.linear"), header = T) %>%
    filter(TEST == "ADD") %>%
    mutate(logP = -log10(P), CHR <= 22) %>%
    filter(!is.na(logP)) 
  df_qq.logistic <- gxtopic %>%
    arrange(logP)
  df_qq.logistic$sim_logP <- sort( - log10( runif(dim(df_qq.logistic)[1])))
  plt <- ggplot(df_qq.logistic) + 
    geom_point(aes(x = sim_logP, y = logP), size = 0.5, alpha = 0.5) + 
    geom_abline(slope = 1, linetype = "dashed", color = red) +
    theme(legend.position = "none",panel.background=element_blank()) + 
    xlab(expression("Observed" * -log[10](P)) ) + ylab(expression("Expected" * -log[10](P))) +
    ggtitle(paste0("BMI controled by topic: ", topic_name[subtype$V3[ccgwas]]))
  
  ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/GxT_simple_model/BMI_topic", 
                topic_name[subtype$V3[ccgwas]],"qq_G.png"), width = 4, height = 4,plt)
}

########################################
# sex specific GxTopic analysis
#######################################
# male
cc_dir <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/GxTopic/"
ds_list <- read.table("GxTopic_disease_list.txt")
# playing with the code
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("V1" = "phecode"))
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
for(idx in 1:length(ds_list$V1)){
  ds_id <- ds_list$V1[idx]
  print(ds_id)
  ds_name <- ds_list$phenotype[idx]
  for(ccgwas in 1:10){
    # only load non-sex chromosome
    gxtopic <- read.table(paste0(cc_dir, ds_id, "GxTopic_men_topic", ccgwas, ".gxtopic"), header = T) %>%
      mutate(logP = -log10(P)) %>%
      filter(!is.na(logP))  # only load non-sex chromosome
    df_qq.logistic <- gxtopic %>%
      arrange(logP)
    df_qq.logistic$sim_logP <- sort( - log10( runif(dim(df_qq.logistic)[1])))
    plt_qq <- ggplot(df_qq.logistic) + 
      geom_point(aes(x = sim_logP, y = logP), size = 0.5, alpha = 0.5) + 
      geom_abline(slope = 1, linetype = "dashed", color = red) +
      theme(legend.position = "none",panel.background=element_blank()) + 
      xlab(expression("Expected" * -log[10](P)) ) + ylab(expression("Observed" * -log[10](P))) +
      ggtitle(paste0(ds_name, " x topic: ", topic_name[ccgwas]))
    
    ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/GxT_men/", ds_id, "_men_topic", 
                  topic_name[ccgwas],"qq_GxTopic.png"), plt_qq, width = 4, height = 4)
    
    # plot manhatten plot
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
    plt_mht <- ggplot(don, aes(x=BPcum, y=logP)) +
      # Show all points
      geom_point( aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
      scale_color_manual(values = rep(c(blue, grey), 22 )) +
      # custom X axis:
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
      scale_y_continuous(expand = c(0.1, 0) ) +     # remove space between plot area and x axis
      geom_hline(aes(yintercept = 7.3), color = "red", linetype = "dashed") +
      labs(x = "Chromosome number", y = expression(-log[10](P)), title = paste0(ds_name, " x ", topic_name[ccgwas])) +
      # Custom the theme:
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
    
    ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/GxT_men/manhatten_",ds_id,"_men_topic", topic_name[ccgwas],".png"), plt_mht, width = 10, height = 4)
  }
}

# female
cc_dir <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/GxTopic/"
ds_list <- read.table("GxTopic_disease_list.txt")
# playing with the code
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("V1" = "phecode"))
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
for(idx in 1:length(ds_list$V1)){
  ds_id <- ds_list$V1[idx]
  print(ds_id)
  ds_name <- ds_list$phenotype[idx]
  for(ccgwas in 1:10){
    # only load non-sex chromosome
    gxtopic <- read.table(paste0(cc_dir, ds_id, "GxTopic_women_topic", ccgwas, ".gxtopic"), header = T) %>%
      mutate(logP = -log10(P)) %>%
      filter(!is.na(logP), CHR <= 22)  # only load non-sex chromosome
    df_qq.logistic <- gxtopic %>%
      arrange(logP)
    df_qq.logistic$sim_logP <- sort( - log10( runif(dim(df_qq.logistic)[1])))
    plt_qq <- ggplot(df_qq.logistic) + 
      geom_point(aes(x = sim_logP, y = logP), size = 0.5, alpha = 0.5) + 
      geom_abline(slope = 1, linetype = "dashed", color = red) +
      theme(legend.position = "none",panel.background=element_blank()) + 
      xlab(expression("Expected" * -log[10](P)) ) + ylab(expression("Observed" * -log[10](P))) +
      ggtitle(paste0(ds_name, " x topic: ", topic_name[ccgwas]))
    
    ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/GxT_women/", ds_id, "_women_topic", 
                  topic_name[ccgwas],"qq_GxTopic.png"), plt_qq, width = 4, height = 4)
    
    # plot manhatten plot
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
    plt_mht <- ggplot(don, aes(x=BPcum, y=logP)) +
      # Show all points
      geom_point( aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
      scale_color_manual(values = rep(c(blue, grey), 22 )) +
      # custom X axis:
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
      scale_y_continuous(expand = c(0.1, 0) ) +     # remove space between plot area and x axis
      geom_hline(aes(yintercept = 7.3), color = "red", linetype = "dashed") +
      labs(x = "Chromosome number", y = expression(-log[10](P)), title = paste0(ds_name, " x ", topic_name[ccgwas])) +
      # Custom the theme:
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
    
    ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/GxT_women/manhatten_",ds_id,"_women_topic", topic_name[ccgwas],".png"), plt_mht, width = 10, height = 4)
  }
}







#######################################
# making nonlinear qq plot
#######################################
# compute qq plot for D = G + T + G*T model
cc_dir <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/GxTopic/"
ds_list <- read.table("GxTopic_disease_list.txt")
# playing with the code
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  left_join(phe_phecode, by = c("V1" = "phecode"))
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
for(idx in 1:length(ds_list$V1)){
  ds_id <- ds_list$V1[idx]
  print(ds_id)
  ds_name <- ds_list$phenotype[idx]
  for(ccgwas in 1:10){
    # only load non-sex chromosome
    gxtopic <- read.table(paste0(cc_dir, ds_id, "GxT_nonlinear_topic", ccgwas, ".gxtopic"), header = T) %>%
      mutate(logP = -log10(P)) %>%
      filter(!is.na(logP), CHR <= 22)  # only load non-sex chromosome
    df_qq.logistic <- gxtopic %>%
      arrange(logP)
    df_qq.logistic$sim_logP <- sort( - log10( runif(dim(df_qq.logistic)[1])))
    plt_qq <- ggplot(df_qq.logistic) + 
      geom_point(aes(x = sim_logP, y = logP), size = 0.5, alpha = 0.5) + 
      geom_abline(slope = 1, linetype = "dashed", color = red) +
      theme(legend.position = "none",panel.background=element_blank()) + 
      xlab(expression("Expected" * -log[10](P)) ) + ylab(expression("Observed" * -log[10](P))) +
      ggtitle(paste0(ds_name, " x topic: ", topic_name[ccgwas]))
    
    ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/GxT_nonlinear/", ds_id, "_topic", 
                  topic_name[ccgwas],"qq_GxT_nonlinear.png"), plt_qq, width = 4, height = 4)
    
    # # plot manhatten plot
    # don <- gxtopic %>% 
    #   # Compute chromosome size
    #   group_by(CHR) %>% 
    #   summarise(chr_len=max(BP)) %>% 
    #   # Calculate cumulative position of each chromosome
    #   mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    #   select(-chr_len) %>%
    #   # Add this info to the initial dataset
    #   left_join(gxtopic, ., by=c("CHR"="CHR")) %>%
    #   # Add a cumulative position of each SNP
    #   arrange(CHR, BP) %>%
    #   mutate( BPcum=BP+tot)
    # axisdf <- don %>% 
    #   group_by(CHR) %>% 
    #   summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
    # plt_mht <- ggplot(don, aes(x=BPcum, y=logP)) +
    #   # Show all points
    #   geom_point( aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
    #   scale_color_manual(values = rep(c(blue, grey), 22 )) +
    #   # custom X axis:
    #   scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    #   scale_y_continuous(expand = c(0.1, 0) ) +     # remove space between plot area and x axis
    #   geom_hline(aes(yintercept = 7.3), color = "red", linetype = "dashed") +
    #   labs(x = "Chromosome number", y = expression(-log[10](P)), title = paste0(ds_name, " x ", topic_name[ccgwas])) +
    #   # Custom the theme:
    #   theme_bw() +
    #   theme( 
    #     legend.position="none",
    #     panel.border = element_blank(),
    #     panel.grid.major.x = element_blank(),
    #     panel.grid.minor.x = element_blank()
    #   )
    # 
    # ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/GxT_simple_model/manhatten_",ds_id,"_topic", topic_name[ccgwas],".png"), plt_mht, width = 10, height = 4)
  }
}



########################################
# focusing on the SNPs that are significant GWAS loci
########################################
# example
ds_id <- 495
thre_gwas <- 5 * 10^(-8)
idx <- which(ds_list$V1 == ds_id)
ds_name <- ds_list$phenotype[idx]
gwas_hits <- read.table(paste0("CC_gwas/",ds_id, ".assoc.logistic"), header = T) %>%
  filter(P <= thre_gwas)
gxt_hits_filter <- list()
for(topic_id in 1:10){
  gxt_hits_filter[[topic_id]] <- read.table(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/GxTopic/",ds_id, "GxT_nonlinear_topic", topic_id,".gxtopic"), header = T) %>%
    filter(TEST == "ADDxCOV1", SNP %in% gwas_hits$SNP) %>%
    mutate(topic = topic_id)
}
gxt_hits_filter <- bind_rows(gxt_hits_filter) 
sum(p.adjust(gxt_hits_filter$P, method = "fdr") < 0.1)
plt_df <- data.frame(y = sort(-log10(gxt_hits_filter$P)), x = sort(-log10(runif(length(gxt_hits_filter$P)))))
ggplot(plt_df) +
  geom_point(aes(x = x, y = y)) +
  geom_abline( slope = 1, intercept = 0)

#############################################
# count the number of GxTopic hits
#############################################
cc_dir <- "/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/GxTopic/"
ds_list <- read.table("GxTopic_disease_list.txt")
subtype_ds <- read.table("subtype_disease_topic_list.txt", sep = "\t") %>%
  pull(V1) %>% unique
# playing with the code
phe_phecode <- read.csv("info_phenotype_phecode.csv") 
ds_list <- ds_list %>%
  filter(V1 %in% subtype_ds) %>%
  left_join(phe_phecode, by = c("V1" = "phecode"))
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
df_subtype <- data.frame(disease = as.numeric(), 
                         phenotypes = as.character(),
                         number_hits_gwas = as.numeric(),
                         number_hits_gxt = as.numeric())
GxTopic_results <- list()
for(idx in 1:length(ds_list$V1)){
  try({
    ds_id <- ds_list$V1[idx]
    description <- ds_list$phenotype[idx]
    print(ds_id)
    # load genetic data: note the genotype data are LD-clumped for r2>0.6 (only exclude those that are in strong LD)
    genotype <- read.table(paste0("GxTopic/",ds_id, "_GWAS_SNP.raw"), header = T )
    snps_independent <- sapply(names(genotype), function(x) str_split(x, "_")[[1]][1])
    
    gxt_hits_filter <-  read.table(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/GxTopic/",ds_id, "_filter_nonlinear.gxtopic"), header = T)%>%
      filter(SNP %in% snps_independent) 
    gxt_hits_filter$GxT_FDR <- p.adjust(gxt_hits_filter$GxT_P, method = "fdr") 
    sig_hits <- sum(gxt_hits_filter$FDR < 0.1)
    number_gwas <- dim(gxt_hits_filter)[1]/10
    df_subtype <- df_subtype %>%
      add_row(disease = ds_id, 
              phenotypes = ds_list$phenotype[idx],
              number_hits_gwas = number_gwas,
              number_hits_gxt = sig_hits)
    GxTopic_results[[idx]] <- gxt_hits_filter %>% 
      mutate(Phecode = ds_id, disease = description) %>%
      select(-TEST, -NMISS, -L95, -U95, -A1)
  })
}
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")
GxTopic_results <- bind_rows(GxTopic_results) %>% 
  filter(!is.na(GxT_P)) %>%
  mutate(topic = topic_name[topic])

subtype_GxToic_results <- read.table("subtype_disease_topic_list.txt", sep = "\t") %>%
  select(V1, V2) %>%
  rename(Phecode = V1, topic = V2) %>%
  filter(Phecode %in% ds_list$V1) %>% 
  left_join(GxTopic_results, by = c("Phecode", "topic")) %>%
  mutate(topic = topic_name[topic]) %>%
  filter(!is.na(GxT_OR))
  
subtype_GxToic_results$studywise_FDR <- p.adjust(subtype_GxToic_results$GxT_P, "fdr")
subtype_GxToic_results %>% filter(studywise_FDR < 0.1) %>% arrange(Phecode, topic) %>%
  write.csv("~/Desktop/comorbidity/paper_writing/supplementary_files/significant_GxT.csv", row.names = F)
subtype_GxToic_results %>%
  arrange(Phecode, topic) %>%
  write.csv("~/Desktop/comorbidity/paper_writing/supplementary_files/all_GWAS_GxT.csv", row.names = F)

plt_df <- data.frame(y = sort(-log10(subtype_GxToic_results$GxT_P)), x = sort(-log10(runif(length(subtype_GxToic_results$GxT_P)))))
plt_qq <- ggplot(plt_df) + 
  geom_point(aes(x = x, y = y), size = 0.5, alpha = 0.5) + 
  geom_abline(slope = 1, linetype = "dashed", color = red) +
  theme(legend.position = "none",panel.background=element_blank()) + 
  xlab(expression("Expected" * -log[10](P)) ) + ylab(expression("Observed" * -log[10](P))) +
  ggtitle(paste0("Diseasex subtype-topic interaction"))

ggsave(paste0("~/Desktop/comorbidity/paper_writing/supplementary_files/QQ_topicall_GxTopic.png"), plt_qq, width = 4, height = 4)

# no subtype diseases
no_subtypes <- GxTopic_results %>%
  anti_join(subtype_GxToic_results, by = c("Phecode", "topic", "SNP"))
plt_df <- data.frame(y = sort(-log10(no_subtypes$GxT_P)), x = sort(-log10(runif(length(no_subtypes$GxT_P)))))
plt_qq <- ggplot(plt_df) + 
  geom_point(aes(x = x, y = y), size = 0.5, alpha = 0.5) + 
  geom_abline(slope = 1, linetype = "dashed", color = red) +
  theme(legend.position = "none",panel.background=element_blank()) + 
  xlab(expression("Expected" * -log[10](P)) ) + ylab(expression("Observed" * -log[10](P))) +
  ggtitle(paste0("Disease x nonsubtype-topic interaction "))
ggsave(paste0("~/Desktop/comorbidity/paper_writing/supplementary_files/QQ_nonsubtype_GxTopic.png"), plt_qq, width = 4, height = 4)
# exploration code below
ds_id <- 174.11 # 250.2
idx <- which(ds_list$V1 == ds_id)
ds_name <- ds_list$phenotype[idx]
# topic_interest <- c(5,7)
gxt_hits_filter <-  read.table(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/GxTopic/",ds_id, "_filter_nonlinear.gxtopic"), header = T) 
gxt_hits_filter$FDR <- p.adjust(gxt_hits_filter$GxT_P, method = "fdr")
plt_df <- data.frame(y = sort(-log10(gxt_hits_filter$GxT_P)), x = sort(-log10(runif(length(gxt_hits_filter$GxT_P)))))
ggplot(plt_df) +
  geom_point(aes(x = x, y = y)) +
  geom_abline( slope = 1, intercept = 0)

topic_id <- 7
gxt_topic <- gxt_hits_filter %>%
  filter(topic == topic_id)

gxt_topic %>% arrange(FDR)

plt_df <- data.frame(y = sort(-log10(gxt_topic$GxT_P)), x = sort(-log10(runif(length(gxt_topic$GxT_P)))))
ggplot(plt_df) +
  geom_point(aes(x = x, y = y)) +
  geom_abline( slope = 1, intercept = 0)
plot(-log10(gxt_topic$GxT_P))


plt_qq <- ggplot(plt_df) + 
  geom_point(aes(x = x, y = y), size = 0.5, alpha = 0.5) + 
  geom_abline(slope = 1, linetype = "dashed", color = red) +
  theme(legend.position = "none",panel.background=element_blank()) + 
  xlab(expression("Expected" * -log[10](P)) ) + ylab(expression("Observed" * -log[10](P))) +
  ggtitle(paste0(ds_name, " x topic: "))

ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/GWAS_hits_GxT/", ds_id, "_topicall","qq_GxTopic.png"), plt_qq, width = 4, height = 4)

# plotting example: manhatten plot with different color refering to different topic;
# all other topics are treated as background 
# plot manhhatten plot for each disease

for(idx in 1:length(ds_list$V1)){
  try({
  ds_id <- ds_list$V1[idx]
  # idx <- which(ds_list$V1 == ds_id)
  ds_name <- ds_list$phenotype[idx]
  # topic_interest <- c(5,9)
  colors_topic <- c(cbPalette, paletteMartin[c(8,11)])
  names(colors_topic) <- as.character(1:10)
  
  gxt_hits_filter <-  read.table(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/GxTopic/",ds_id, "_filter_nonlinear.gxtopic"), header = T) 
  gxt_hits_filter$FDR <- p.adjust(gxt_hits_filter$GxT_P, method = "fdr")
  hits <- gxt_hits_filter %>% filter(FDR < 0.1) %>%
    group_by(topic) %>% tally()
  
  
  min_fdr_p <- gxt_hits_filter %>%
    filter(FDR < 0.1) %>%
    summarise(max(GxT_P)) %>%
    pull(1)
  for(topic_id in hits$topic){
    gxt_topic <- gxt_hits_filter %>%
      filter(topic == topic_id)
    
    # input_gxt <- gxt_hits_filter
    input_gxt <- gxt_topic
    
    
    don <- input_gxt %>% 
      # Compute chromosome size
      group_by(CHR) %>% 
      summarise(chr_len=max(BP)) %>% 
      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
      select(-chr_len) %>%
      # Add this info to the initial dataset
      left_join(input_gxt, ., by=c("CHR"="CHR")) %>%
      # Add a cumulative position of each SNP
      arrange(CHR, BP) %>%
      mutate( BPcum=BP+tot)
    axisdf <- don %>% 
      group_by(CHR) %>% 
      summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
    plt_mht <- ggplot(don, aes(x=BPcum, y= -log10(GxT_P) )) +
      # Show all points
      geom_point( aes(color=factor(topic, levels = 1:10)), alpha=1, size=1.3) +
      scale_color_manual(values =colors_topic) +
      # custom X axis:
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
      scale_y_continuous(expand = c(0.1, 0) ) +     # remove space between plot area and x axis
      geom_hline(aes(yintercept = -log10(min_fdr_p)), color = "red", linetype = "dashed") +
      labs(x = "Chromosome number", y = expression(-log[10](P)), title = paste0(ds_name, " x ", topic_name[topic_id])) +
      # Custom the theme:
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
    ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/GWAS_hits_GxT/gwas_hits_manhatten_",ds_id,"_topic", topic_name[topic_id],".png"), plt_mht, width = 10, height = 4) 
    
  }
  # a piece of code ajusting the manhattan plot location
  # gxt_topic <- gxt_hits_filter %>%
  #   filter(topic == topic_id) %>%
  #   # Compute gap between hits
  #   group_by(CHR) %>% 
  #   mutate(gap_hits = pmin(BP - lag(BP), 2*10^5) ) %>%
  #   mutate(gap_hits = if_else(is.na(gap_hits), 10^6, gap_hits) ) %>%
  #   # Calculate cumulative position of each chromosome
  #   mutate(BP=cumsum(gap_hits)  ) %>%
  #   ungroup()
  })
}


##############################################################################
# compute effect sizes over different percentile of topic weights
##############################################################################
load("../Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
# compute the patient loadings 
patient_loadings <- sweep((para$alpha_z-1), 1, rowSums(para$alpha_z - 1), FUN="/")
patient_loadings <- data.frame(eid = para$eid, loadings = patient_loadings)
topic_name <- c("MDS", "ARP", "FGND", "NRI", "CER", "MGND", "CVD", "UGI", "LGI", "SRD")

ds_id <- "250.2"
topic_id <- 7
# variant_name <- "rs1412829" # rs343092
variant_name <- "rs1063192" 
variant_name <- "rs1042725" 

# for asthma
ds_id <- "495"
topic_id <- 10
read.table(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/GxTopic/",ds_id, "_filter_nonlinear.gxtopic"), header = T) %>%
  arrange(GxT_P)
variant_name <- "rs72836344" 
variant_name <- "rs1837253" 


# systematically plot the SNP disease
subtype_GxToic_results <- read.csv("~/Desktop/comorbidity/paper_writing/supplementary_files/significant_GxT.csv") %>%
  mutate(topic = match(topic, topic_name))

ds_list_gxt <- subtype_GxToic_results %>%
  group_by(Phecode) %>%
  slice(1) %>%
  pull(Phecode)

# read the gwas_catalog data 
gwas_catalog <- read.table("../UKBB_interim_files/GWAS_catalog_data", header = T, sep = "\t",fill = T)
v2gene <- gwas_catalog %>%
  select(SNPS, REGION, CHR_ID, MAPPED_GENE) %>%
  group_by(SNPS) %>%
  slice(1) %>%
  ungroup()

rec_data <- read.csv("rec2subjectAbove1000occur_include_death_PheCode.csv")
all_eid <- rec_data %>%
  group_by(eid) %>%
  summarise() 

# set the number of tiles
number_tiles <- 4
tile_effect <- list()
cnt <- 1
for(ds_id in ds_list_gxt){
  print(paste0("Disease: ", ds_id))
  idx <- which(ds_list$V1 == ds_id)
  ds_name <- ds_list$phenotype[idx]

  # load genetic data: note the genotype data are LD-clumped for r2>0.6 (only exclude those that are in strong LD)
  genotype <- read.table(paste0("GxTopic/",ds_id, "_GWAS_SNP.raw"), header = T )
  snps_independent <- sapply(names(genotype), function(x) str_split(x, "_")[[1]][1])
  
  
  # identify the list of SNPs
  SNP_topic_list <- read.table(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/GxTopic/",ds_id, "_filter_nonlinear.gxtopic"), header = T) %>%
    filter(SNP %in% snps_independent)
  SNP_topic_list$FDR <- p.adjust(SNP_topic_list$GxT_P, method = "fdr")
  SNP_topic_list <- subtype_GxToic_results  %>% 
    filter(Phecode ==  ds_id) %>%
    left_join(v2gene, by = c("SNP" = "SNPS"))
  
  topic_list <- SNP_topic_list %>%
    pull(topic) %>%
    unique()
  
  for(topic_id in topic_list){
    print(paste0("topic: ", topic_id))
    snp_list <- SNP_topic_list %>%
      filter(topic == topic_id)
    
    # casectrl <- read.table(paste0("GxTopic/",ds_id, "_matching_topic_casectrl.txt"), header = T)
    # topic_casectrl <-  casectrl %>%
    #   left_join(select(patient_loadings, eid, paste0("loadings.", topic_id)), by = c("case" = "eid") )
    # topic_casectrl <- topic_casectrl %>%
    #   mutate(loading_percentile = ntile(topic_casectrl[[paste0("loadings.",topic_id)]], number_tiles))
    
    casectrl <- read.table(paste0("GxTopic/",ds_id, "_matching_topic_casectrl.txt"), header = T)
    cases <- casectrl %>% 
      select(case) %>%
      group_by(case) %>%
      slice(1) %>%
      ungroup() %>%
      rename(eid = case) %>%
      mutate(phenotype = 1) %>%
      left_join(select(patient_loadings, eid, paste0("loadings.", topic_id)), by = c("eid") ) 
    cases <- cases %>%
      mutate(loading_percentile = ntile(cases[[paste0("loadings.",topic_id)]], number_tiles))  %>%
      rename(loading = paste0("loadings.", topic_id))
    controls <- all_eid %>%
      anti_join(cases, by = c("eid")) %>%
      left_join(select(patient_loadings, eid, paste0("loadings.", topic_id)), by = c("eid") ) %>%
      rename(loading = paste0("loadings.", topic_id))
    # 2022-06-30: sampling controls that have same topic loadings as cases. To avoid signal on topic loadings
    control_list <- list()
    for(tls in 1:number_tiles){
      tls_case <- cases %>% 
        filter(loading_percentile == tls)
      min_tls <- tls_case$loading %>%
        min(na.rm = T)
      max_tls <- tls_case$loading %>%
        max(na.rm = T)
      control_sz <- controls %>%
        filter(loading < max_tls, loading > min_tls) %>%
        dim
      if(control_sz[1] < 2*(dim(tls_case)[1])){
        control_list[[tls]] <- controls %>%
          filter(loading < max_tls, loading > min_tls) %>%
          sample_n(2*(dim(tls_case)[1]), replace = T) %>%
          mutate(loading_percentile = tls, phenotype = 0)
      }else{
          control_list[[tls]] <- controls %>%
            filter(loading < max_tls, loading > min_tls) %>%
            sample_n(2*(dim(tls_case)[1])) %>%
            mutate(loading_percentile = tls, phenotype = 0)
      }

    }
    controls <- bind_rows(control_list)
    topic_casectrl <- bind_rows(cases, controls)
    # controls <- casectrl %>%
    #   select(control) %>%
    #   rename(eid = control) %>%
    #   mutate(phenotype = 0)
    # topic_casectrl <-  bind_rows(cases, controls) %>%
    #   left_join(select(patient_loadings, eid, paste0("loadings.", topic_id)), by = c("eid") )
    # topic_casectrl <- topic_casectrl %>%
    #   mutate(loading_percentile = ntile(topic_casectrl[[paste0("loadings.",topic_id)]], number_tiles))
    # topic_casectrl %>%
    #   group_by(loading_percentile) %>%
    #   summarise(mean(phenotype))
    #  # also separate the controls
    # topic_ctrl <-  casectrl %>%
    #   left_join(select(patient_loadings, eid, paste0("loadings.", topic_id)), by = c("control" = "eid") )
    # topic_ctrl <- topic_ctrl %>%
    #   mutate(loading_percentile = ntile(topic_ctrl[[paste0("loadings.",topic_id)]], number_tiles))
    
    ############################################
    # compute non-interacting SNP effect changes 
    ############################################
    gxt_hits_filter <-  read.table(paste0("/Users/xilin/Desktop/comorbidity/Multi-morbidity_biobank/GxTopic/",ds_id, "_filter_nonlinear.gxtopic"), header = T) %>%
      filter(topic == topic_id, GxT_P > 0.05)
    SNP_list_background <- gxt_hits_filter$SNP
    
    idx_background <- sapply(names(genotype), function(x) any(str_detect(x, SNP_list_background)), simplify = T ) 
    variant_background <- names(genotype)[idx_background]
    genotype_background <- genotype %>%
      select( c("FID", variant_background) )
    
    background_effect_estimate <- list()
    for(percentile_of_interest in 1:number_tiles){
      print(paste0("tile:", percentile_of_interest) ) 
      # cases <- topic_casectrl %>%
      #   filter(loading_percentile == percentile_of_interest) %>%
      #   select(case) %>%
      #   group_by(case) %>%
      #   slice(1) %>%
      #   ungroup() %>%
      #   rename(eid = case) %>%
      #   mutate(phenotype = 1)
      # controls <- topic_ctrl %>%
      #   filter(loading_percentile == percentile_of_interest) %>%
      #   select(control) %>%
      #   rename(eid = control) %>%
      #   mutate(phenotype = 0)
      # 
      gwas_data <- topic_casectrl %>% 
        filter(loading_percentile == percentile_of_interest) %>%
        left_join(genotype_background, by = c("eid" = "FID") ) %>%
        left_join(select(patient_loadings, eid, paste0("loadings.", topic_id)), by = c("eid") )
      # print(mean(gwas_data$rs1412829_G, na.rm  = T))
      effect_sizes <- c()
      for(vrt in variant_background){
        lm_model <- glm(data = gwas_data, formula = as.formula(paste0("phenotype ~", vrt) ), family = binomial)
        effect_sizes <- c(effect_sizes, summary(lm_model)$coefficients[2,1])
      }
      background_effect_estimate[[percentile_of_interest]] <- effect_sizes
    }
    background_effect_estimate <- do.call(rbind, background_effect_estimate)
    # names(background_effect_estimate) <- 1:number_tiles
    for(rw in 1:dim(background_effect_estimate)[2]){ 
      if(mean(background_effect_estimate[,rw]) < 0){
        background_effect_estimate[,rw] <- -background_effect_estimate[,rw]
      }
    }
    
    # perform a regression on background SNPs
    background_data <- data.frame(tiles = 1:number_tiles, effect = background_effect_estimate) %>%
      pivot_longer(cols = !tiles, names_to = "SNP", values_to = "effect")
    
    rg_model <- lm(data = background_data, formula = effect ~ tiles)
    summary(rg_model)
    p_background <- summary(rg_model)$coefficients[2,4]
    
    for(variant_id in 1:dim(snp_list)[1]){
      variant_name <- snp_list$SNP[variant_id]
      # break if the SNP is no in the genotype data (in ld > 0.99 with another snp)
      # if(!sum(sapply(names(genotype), function(x) str_detect(x, variant_name) ))){
      #   next
      # }
      idx <- which( sapply(names(genotype), function(x) str_detect(x, variant_name) ) )
      variant <- names(genotype)[idx]
      genotype_data <- genotype %>%
        select(FID, variant)
      effect_sizes <- c()
      se_effect <- c()
      for(percentile_of_interest in 1:number_tiles){
        # cases <- topic_casectrl %>%
        #   filter(loading_percentile == percentile_of_interest) %>%
        #   select(case) %>%
        #   group_by(case) %>%
        #   slice(1) %>%
        #   ungroup() %>%
        #   rename(eid = case) %>%
        #   mutate(phenotype = 1)
        # controls <- topic_ctrl %>%
        #   filter(loading_percentile == percentile_of_interest) %>%
        #   select(control) %>%
        #   rename(eid = control) %>%
        #   mutate(phenotype = 0)
        
        gwas_data <- topic_casectrl %>% 
          filter(loading_percentile == percentile_of_interest) %>% 
          left_join(genotype_data, by = c("eid" = "FID") ) %>%
          left_join(select(patient_loadings, eid, paste0("loadings.", topic_id)), by = c("eid") )
        # print(mean(gwas_data$rs1412829_G, na.rm  = T)) # check if the MAF are the same in different quartiles
        lm_model <- glm(data = gwas_data, formula = as.formula(paste0("phenotype ~", variant) ), family = binomial)
        effect_sizes <- c(effect_sizes, summary(lm_model)$coefficients[2,1])
        se_effect <- c(se_effect, summary(lm_model)$coefficients[2,2])
      }

      
      # computing the population level effect sizes
      # cases_all <- topic_casectrl %>%
      #   select(case) %>%
      #   group_by(case) %>%
      #   slice(1) %>%
      #   ungroup() %>%
      #   rename(eid = case) %>%
      #   mutate(phenotype = 1)
      # controls_all <- topic_ctrl %>%
      #   select(control) %>%
      #   rename(eid = control) %>%
      #   mutate(phenotype = 0)
      gwas_data_all <-  topic_casectrl %>% 
        left_join(genotype_data, by = c("eid" = "FID") ) %>%
        left_join(select(patient_loadings, eid, paste0("loadings.", topic_id)), by = c("eid") )
      lm_model <- glm(data = gwas_data_all, formula = as.formula(paste0("phenotype ~", variant) ), family = binomial)
      effect_population <- summary(lm_model)$coefficients[2,1]
      
      # flip the effect sizes
      if(mean(effect_sizes) < 0){
        effect_sizes <- - effect_sizes
        effect_population <- -effect_population
      }
      plot_df <- data.frame(quantile = 1:number_tiles, mean_effect = effect_sizes, se_effect = se_effect)
      p_diff <- (1-pnorm(abs(plot_df$mean_effect[1] - plot_df$mean_effect[number_tiles])/
                           sqrt(plot_df$se_effect[1]^2 + plot_df$se_effect[number_tiles]^2) ))*2
      
      # save the data
      tile_effect[[cnt]] <- plot_df %>%
        mutate(disease = ds_name , topic = topic_name[topic_id], SNP = variant_name, GWAS_catalog_gene = snp_list$MAPPED_GENE[variant_id], p_top_bottom = p_diff)
      cnt <- cnt+1
      
      #################### plot with both dots and the regression line as a background 
      # pretty complicated to determing the upper bound and lower bound of the figure
      plt_SNP <- ggplot() + 
        geom_jitter(data = background_data, 
                    mapping=aes(x=tiles, y=effect),width=0.2, size = 0.1, alpha=0.2) + 
        lims(y = c( min(c(effect_sizes - se_effect, effect_population, quantile(background_data$effect, 0.05))) - 0.05,
                    max(c(effect_sizes + se_effect, effect_population, quantile(background_data$effect, 0.95))) + 0.05)) + 
        geom_smooth(data = background_data, mapping=aes(x=tiles, y=effect), color = grey, method='lm') + 
        geom_pointrange(data = plot_df, aes(x=quantile, y=mean_effect, ymin=mean_effect - se_effect, ymax=mean_effect + se_effect), color = blue) + 
        theme(legend.position = "none",panel.background=element_blank(), plot.title = element_text(size = 12)) + 
        xlab( paste0("Quartile of ", topic_name[topic_id], " Topic") ) + ylab("Log odds ratio") +
        ggtitle(paste0(variant_name, " x ", topic_name[topic_id], " for " ,ds_name, "\n",
                       "chr", snp_list$CHR[variant_id],": ", snp_list$REGION[variant_id], " (nearest gene:", str_split(snp_list$MAPPED_GENE[variant_id], ",")[[1]][1], ")\n",
                       "P = ", snp_list$GxT_P[variant_id], " for interaction \n")) +
        # ggtitle(paste0(ds_name, ": ", variant_name, "\n",
        #                "Mapped gene:", snp_list$MAPPED_GENE[variant_id], "\n",
        #                "GWAS P-value: ", snp_list$P[variant_id], "\n",
        #                "GxTopic P-value: ", snp_list$GxT_P[variant_id], "\n",
        #                "Topic: ", topic_name[topic_id])) +
        # annotate(geom = 'text', label = paste0("Top vs. bottom quartile P=",round(p_diff, digits = 5)), x = -Inf, y = Inf, hjust = 0, vjust = 1) +
        geom_hline(aes(yintercept = effect_population), color = red, linetype = "dashed") 
      ggsave(paste0("~/Desktop/comorbidity/figures/GxTopic/GxT_effect_changes/single_SNP_effect_",ds_id,"_",topic_name[topic_id], "_",variant_name, ".png"), plt_SNP, width = 4, height = 4)
      
    }
    
  }
  
}
tile_effect <- bind_rows(tile_effect)
tile_effect %>% 
  write.csv("~/Desktop/comorbidity/paper_writing/supplementary_files/effect_sizes_topic_tile.csv", row.names = F)

##############################################################################
# 2022-06-30: checking if the case control topic loading are indeed balanced
##############################################################################
ds_id <- "250.2"
topic_id <- 7


casectrl <- read.table(paste0("GxTopic/",ds_id, "_matching_topic_casectrl.txt"), header = T)
topic_casectrl <-  casectrl %>%
  left_join(select(patient_loadings, eid, paste0("loadings.", topic_id)), by = c("case" = "eid") )
topic_casectrl <- topic_casectrl%>%
  mutate(loading_percentile = ntile(topic_casectrl[[paste0("loadings.",topic_id)]], number_tiles))

# also separate the controls
topic_ctrl <-  casectrl %>%
  left_join(select(patient_loadings, eid, paste0("loadings.", topic_id)), by = c("control" = "eid") )
topic_ctrl <- topic_ctrl %>%
  mutate(loading_percentile = ntile(topic_ctrl[[paste0("loadings.",topic_id)]], number_tiles))

hist(topic_casectrl$loadings.7)
# looking at the highest percentile, since the distribution of topic loading has a heavy tail
percentile_of_interest <- 4
cases <- topic_casectrl %>%
  filter(loading_percentile == percentile_of_interest) %>%
  select(case) %>%
  group_by(case) %>%
  slice(1) %>%
  ungroup() %>%
  rename(eid = case) %>%
  left_join(select(patient_loadings, eid, paste0("loadings.", topic_id)), by = c("eid") ) %>%
  mutate(phenotype = 1)
controls <- topic_ctrl %>%
  filter(loading_percentile == percentile_of_interest) %>%
  select(control) %>%
  rename(eid = control) %>%
  left_join(select(patient_loadings, eid, paste0("loadings.", topic_id)), by = c("eid") ) %>%
  mutate(phenotype = 0)
hist_data <- bind_rows(cases, controls) %>%
  rename(loading =  paste0("loadings.", topic_id))
ggplot(data = hist_data) + 
  geom_histogram(aes(x = loading, fill = phenotype))


genotype <- read.table(paste0("GxTopic/",ds_id, "_GWAS_SNP.raw"), header = T )
variant_name <- "rs1063192"
idx <- which( sapply(names(genotype), function(x) str_detect(x, variant_name) ) )
variant <- names(genotype)[idx]
genotype_data <- genotype %>%
  select(FID, variant)
effect_sizes <- c()
se_effect <- c()
for(percentile_of_interest in 1:number_tiles){
  cases <- topic_casectrl %>%
    filter(loading_percentile == percentile_of_interest) %>%
    select(case) %>%
    group_by(case) %>%
    slice(1) %>%
    ungroup() %>%
    rename(eid = case) %>%
    mutate(phenotype = 1)
  controls <- topic_ctrl %>%
    filter(loading_percentile == percentile_of_interest) %>%
    select(control) %>%
    rename(eid = control) %>%
    mutate(phenotype = 0)
  
  gwas_data <- bind_rows(cases, controls) %>% 
    left_join(genotype_data, by = c("eid" = "FID") ) %>%
    left_join(select(patient_loadings, eid, paste0("loadings.", topic_id)), by = c("eid") )
  # print(mean(gwas_data$rs1412829_G, na.rm  = T)) # check if the MAF are the same in different quartiles
  # lm_model <- glm(data = gwas_data, formula = as.formula(paste0("phenotype ~", variant, "+", "loadings.", topic_id) ), family = binomial)
  lm_model <- glm(data = gwas_data, formula = as.formula(paste0("phenotype ~", variant)), family = binomial)
  effect_sizes <- c(effect_sizes, summary(lm_model)$coefficients[2,1])
  se_effect <- c(se_effect, summary(lm_model)$coefficients[2,2])
}

