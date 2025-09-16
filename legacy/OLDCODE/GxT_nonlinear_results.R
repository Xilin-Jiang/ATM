# extract topic loading from the target results
# extract age distribution
setwd("/well/mcvean/xilin/Multimorbidity-biobank/")
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)                                  
ds_id <- args[1]
LOC <- paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/GxTopic/", ds_id, "/")

# filter only the significant SNPs 
thre_gwas <- 5*10^(-8)
GWASloc <- paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Association_analysis/", ds_id, "/")
# gwas_hits <- read.table(paste0(GWASloc, ds_id, ".assoc.logistic"), header = T) %>%
#   filter(P <= thre_gwas)
# rerun the analysis using the set of heritable diseases 
gwas_hits <- read.table(paste0(LOC, ds_id, ".assoc.logistic"), header = T) %>%
  filter(P <= thre_gwas)

# save hits from two model
StratifyDisease_hit <- c()
GxTopic_hit <- c()
gxt_hits_filter <- list()
for(topic_id in 1:10){
  gxt_hits <- read.table(paste0(LOC, ds_id, "GxT_nonlinear_topic", topic_id,".assoc.logistic"), header = T) %>%
    filter(TEST == "ADDxCOV1")
  write.table(gxt_hits, paste0(LOC,ds_id, "GxT_nonlinear_topic", topic_id,".gxtopic"), sep="\t", col.names = T, row.names = FALSE, quote = F)
  GxTopic_hit[topic_id] <- gxt_hits %>% 
    filter(P < 5*10^(-8)) %>%
    tally() %>%
    pull
    
  # save stratified method
  strf_hits <- read.table(paste0(LOC, ds_id, "GxT_DiseaseStratify_topic", topic_id,".qassoc.gxe"), header = T) 
  StratifyDisease_hit[topic_id] <- strf_hits %>% 
    filter(P_GXE < 5*10^(-8)) %>%
    tally() %>%
    pull
  
  # filter only genome-wide significant loci
  gxt_hits_filter[[topic_id]] <- read.table(paste0(LOC, ds_id, "GxT_nonlinear_topic", topic_id,".assoc.logistic"), header = T) %>%
    filter(TEST == "ADDxCOV1", SNP %in% gwas_hits$SNP) %>%
    mutate(topic = topic_id)
}
data.frame(topic = 1:10, StratifyDisease_hit = StratifyDisease_hit, 
           GxTopic_hit = GxTopic_hit) %>% mutate(disease = ds_id) %>%
  write.table(paste0(LOC,ds_id, "_nonlinear_hit.txt"), sep="\t", col.names = T, row.names = FALSE, quote = F)

gxt_hits_filter <- bind_rows(gxt_hits_filter) %>%
  rename(GxT_OR = OR, GxT_STAT = STAT, GxT_P = P) %>%
  select(-TEST, -NMISS)
gxt_hits_filter <- gwas_hits %>%
  left_join(gxt_hits_filter, by = c("CHR", "SNP", "BP", "A1"))

write.table(gxt_hits_filter, paste0(LOC,ds_id, "_filter_nonlinear.gxtopic"), sep="\t", col.names = T, row.names = FALSE, quote = F)
