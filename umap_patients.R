###########################################################
# need to add small amount of noise to avoid numeric error
# umap_patients.R
###########################################################
require(umap)
load("/users/mcvean/xilin/xilin/Multimorbidity-biobank/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
set.seed(19940110)
perturbed_data <- rbind(patient_loadings, onehot_anchors) + matrix(0.002 * rnorm(para$K * (para$M + 10)), nrow = para$M + 10)
umap.fit <- umap(perturbed_data)
save(umap.fit, file = "umap.fit.RData")