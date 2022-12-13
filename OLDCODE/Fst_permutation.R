library(dplyr)
library(stringr)
##########################################
# extract LDSC results
##########################################
args <- commandArgs(trailingOnly = TRUE)  
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/Association_analysis/"
ds_id <- args[1]
rep_id <- as.numeric(args[2])
set.seed(19940110 + rep_id)
if(rep_id == 1){
  permutation_Fst <- data.frame(weighted_Fst = as.numeric())
  write.csv(permutation_Fst, paste(DIR, ds_id,"/",ds_id,"_permutation_Fst.csv",sep=""), row.names = F)
}else{
  permutation_Fst <- read.csv(paste(DIR, ds_id,"/",ds_id,"_permutation_Fst.csv",sep=""))
  fst <- readLines(paste(DIR, ds_id,"/",ds_id,"_permutation_Fst.log",sep=""))
  permutation_Fst <- permutation_Fst %>%
    add_row(weighted_Fst = as.numeric(str_split(fst[length(fst) - 2], "\\s+")[[1]][4]))
  write.csv(permutation_Fst, paste(DIR, ds_id,"/",ds_id,"_permutation_Fst.csv",sep=""), row.names = F)
}
subtypes <- read.table(paste(DIR, ds_id,"/",ds_id,"subtypes.txt",sep=""), header = F)
permutation_subtypes <- subtypes
permutation_subtypes$V3 <- sample(subtypes$V3)
write.table(permutation_subtypes, paste(DIR, ds_id,"/",ds_id,"_permutation_subtypes.txt",sep=""), sep="\t", col.names = F, row.names = FALSE ,quote = F)
