# extract intercept for LDSC according to equation 16 in supplementary note of Bulik-Sullivan 2015 NG
require(dplyr)
ds_list <- read.table("/users/mcvean/xilin/xilin/Multimorbidity-biobank/top_hererogeneous_disease.txt") %>%
  mutate(V1 = as.character(V1))
intercept_rg <- matrix(NA, length(ds_list[[1]]), 3)
for(i in 1:length(ds_list[[1]])){
  ds_id <- ds_list[[1]][i]
  DIR <- paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Association_analysis/", ds_id, "/")
  pt <- paste( ds_id, "subtp_[0123456789]_pheno_over_age.txt",sep="")
  temp <- list.files(paste(DIR, sep=""), pattern=pt)
  
  N <- rep(NA, length(temp))
  P <- rep(NA, length(temp))
  Ncontrol <- rep(NA, length(temp))
  controls <- list()
  Ncc <- rep(NA, length(temp))
  for(subtp in 1:length(temp)){
    DataGwas <- read.table(paste0(DIR, temp[subtp]), header = F)
    N[subtp] <- dim(DataGwas)[1]
    controls[[subtp]] <- DataGwas %>%
      filter(V3 ==1) %>%
      pull(V1)
    Ncontrol[subtp] <- length(controls[[subtp]])
    P[subtp] <- (1 - Ncontrol[subtp]/N[subtp])
    Ncc[subtp] <- length(intersect(controls[[subtp]], controls[[1]])) # overlapped controls
  }
  intercept_rg[i,1:length(temp)] <- sqrt(as.numeric(N[1]) * as.numeric(N) * P[1] * (1- P[1]) * P * (1- P)) * Ncc/(as.numeric(Ncontrol[1]) * as.numeric(Ncontrol))
  ds_list$intercept_rg <- intercept_rg
}

write.table(ds_list, paste0("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Association_analysis/ldsc_intercept.txt"), sep="\t", col.names = FALSE, row.names = FALSE, quote = F)