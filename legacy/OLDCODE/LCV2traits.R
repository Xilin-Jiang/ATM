require(dplyr)
require(stringr)
# LCV example script written by Katie Siewert
args <- commandArgs(trailingOnly = TRUE)                                  
trait1 <- as.numeric(args[1])
trait2 <- as.numeric(args[2])

young_group <- read.table("/users/mcvean/xilin/xilin/Multimorbidity-biobank/trait1_LCV.txt",quote = '', sep = '\t', header = F)
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
# trait2 <- 8
ds_id1 <- young_group$V1[trait1]
topic_id1 <- young_group$V2[trait1]
ds_id2 <- old_group$V1[trait2]
topic_id2 <- old_group$V2[trait2]
#Start with data munged using the ldsc package
trait1File=paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/causality_analysis/", ds_id1,"_topic", topic_id1,"/", ds_id1, "_topic", topic_id1, ".sumstats.gz",sep="")

trait2File=paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/causality_analysis/", ds_id2,"_topic", topic_id2,"/", ds_id2, "_topic", topic_id2, ".sumstats.gz",sep="")

#Load trait 1 data and calculate Zs
d1 = na.omit(read.table(gzfile(trait1File),header=TRUE,sep="\t",stringsAsFactors = FALSE))

#Load trait 2 data and calculate Zs
d2 = na.omit(read.table(gzfile(trait2File),header=TRUE,sep="\t",stringsAsFactors = FALSE))

#Load LD scores
DIR_ldscore <- "/users/mcvean/xilin/xilin/UK_biobank/eur_w_ld_chr/"
temp <- list.files(paste(DIR_ldscore, sep=""), pattern=".ldscore.gz")
ldscores_data <- list()
for(i in 1:length(temp)){
  ldscores_data[[i]] <- read.table(gzfile(paste0(DIR_ldscore, temp[i])),header=TRUE,sep='\t',stringsAsFactors=FALSE)
}
ldscores_data <- bind_rows(ldscores_data)

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


#Run LCV-need to setwd to directory containing LCV package
setwd("/users/mcvean/xilin/xilin/Multimorbidity-biobank/")
source("RunLCV.R")

LCV = RunLCV(data$L2,data$Z.x,data$Z.y)
LCV_result <- data.frame(trait1 = ds_id1, trait2 = ds_id2, topic1 = topic_id1, topic2 = topic_id2, gcp = LCV$gcp.pm, gcp_se = LCV$gcp.pse, log10P = log(LCV$pval.gcpzero.2tailed)/log(10) ) 
LCV_save <- list(LCV_result, LCV)
save(LCV_save, 
     file = paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/causality_analysis/LCV_results/", 
                   ds_id1, "_topic", topic_id1, "_", ds_id2, "_topic", topic_id2, ".RData"))
sprintf("Estimated posterior gcp=%.2f(%.2f), log10(p)=%.1f; estimated rho=%.2f(%.2f)",LCV$gcp.pm, LCV$gcp.pse, log(LCV$pval.gcpzero.2tailed)/log(10), LCV$rho.est, LCV$rho.err)

