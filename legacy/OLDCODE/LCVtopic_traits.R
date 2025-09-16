require(dplyr)
require(stringr)
# LCV example script written by Katie Siewert
args <- commandArgs(trailingOnly = TRUE)                                  
# topic <- as.numeric(args[1])
trait <- as.numeric(args[1])

old_group <- read.table("/users/mcvean/xilin/xilin/Multimorbidity-biobank/trait2_LCV.txt",quote = '', sep = '\t', header = F)
# using the function of package pdist
rec_data <- read.csv("/users/mcvean/xilin/xilin/Multimorbidity-biobank/rec2subjectAbove1000occur_include_death_PheCode.csv")

# create subtype files
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/causality_analysis/"
common_disease_within_topics <- read.table("/users/mcvean/xilin/xilin/Multimorbidity-biobank/all_disease_topic_list.txt", header = F)
load("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")

# ######################################
# # testing: 250.2_topic5; 401_topic5
# trait1 <- 13
# trait <- 23

ds_id <- old_group$V1[trait]
topic_id <- old_group$V2[trait]
#Start with data munged using the ldsc package
traitFile=paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/causality_analysis/", ds_id,"_topic", topic_id,"/", ds_id, "_topic", topic_id, ".sumstats.gz",sep="")
#Load trait 2 data and calculate Zs
d2 = na.omit(read.table(gzfile(traitFile),header=TRUE,sep="\t",stringsAsFactors = FALSE))

#Load LD scores
DIR_ldscore <- "/users/mcvean/xilin/xilin/UK_biobank/eur_w_ld_chr/"
temp <- list.files(paste(DIR_ldscore, sep=""), pattern=".ldscore.gz")
ldscores_data <- list()
for(i in 1:length(temp)){
  ldscores_data[[i]] <- read.table(gzfile(paste0(DIR_ldscore, temp[i])),header=TRUE,sep='\t',stringsAsFactors=FALSE)
}
ldscores_data <- bind_rows(ldscores_data)

#Run LCV-need to setwd to directory containing LCV package
setwd("/users/mcvean/xilin/xilin/Multimorbidity-biobank/")
source("RunLCV.R")

LCV_data <- data.frame(topic = as.numeric(), trait = as.numeric(),trait_topic= as.character(), 
                       gcp= as.numeric(),gcp_se= as.numeric(),log10P= as.numeric(),
                       rho.zscore = as.numeric(), rho.est = as.numeric(), rho.err = as.numeric(),
                       p.topicfullycasaul = as.numeric(), p.traitfullycasaul = as.numeric())
for(topic in 1:para$K){
  print(paste0("topic id ", topic))
  #Load trait 1 data and calculate Zs
  topicFile=paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/causality_analysis/topic", topic,"/topic", topic, ".sumstats.gz",sep="")
  d1 = na.omit(read.table(gzfile(topicFile),header=TRUE,sep="\t",stringsAsFactors = FALSE))
    
  #Merge
  m = merge(ldscores_data,d1,by="SNP")
  data = merge(m,d2,by="SNP")
  
  #Sort by position 
  data = data[order(data[,"CHR"],data[,"BP"]),]
  
  #Flip sign of one z-score if opposite alleles-shouldn't occur with UKB data
  #If not using munged data, will have to check that alleles match-not just whether they're opposite A1/A2
  mismatch = which(data$A1.x!=data$A1.y,arr.ind=TRUE)
  data[mismatch,]$Z.y = data[mismatch,]$Z.y*-1
  data[mismatch,]$A1.y = data[mismatch,]$A1.x
  data[mismatch,]$A2.y = data[mismatch,]$A2.x
  
  LCV = RunLCV(data$L2,data$Z.x,data$Z.y)
  LCV_data <- LCV_data %>% 
    tibble::add_row(topic = topic, trait = ds_id, trait_topic = topic_id, gcp = LCV$gcp.pm, 
                    gcp_se = LCV$gcp.pse, log10P = log(LCV$pval.gcpzero.2tailed)/log(10),
                    rho.zscore = LCV$rho.est/LCV$rho.err, rho.est = LCV$rho.est, rho.err = LCV$rho.err,
                    p.topicfullycasaul = LCV$pval.fullycausal[1], p.traitfullycasaul = LCV$pval.fullycausal[2])
}

save(LCV_data, 
     file = paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/causality_analysis/LCV_topic/", 
                   ds_id, "_topic", topic_id, "_topicLCV.RData"))

