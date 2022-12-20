test_that("ATM object initialization", {
  #
})

test_that("each step of the inference", {
  #
})

test_that("inference wrapper", {
  inference_results <- wrapper_ATM(para_sim$rec_data, topic_num = 3, CVB_num = 1)
})

test_that("estimate topic weights from fixed topic loadings", {
  para_loadings <- loading2weights(HES_age_example)
  new_weights <- sweep((para_loadings$alpha_z - 1), 1, rowSums(para_loadings$alpha_z -1), FUN="/")
  new_weights <- data.frame(eid = para_loadings$eid, loading = new_weights)
  load("~/Desktop/comorbidity/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
  patient_weights <- sweep((para$alpha_z - 1), 1, rowSums(para$alpha_z -1), FUN="/")
  originial_weights <- data.frame(eid = para$eid, loading = patient_weights)
  corr_topics <- c()
  for(i in 1:para$K){
    combined_weights <- new_weights %>%
      select(eid, paste0("loading.", i)) %>%
      left_join(select(originial_weights, eid, paste0("loading.", i)), by = c("eid"))
    corr_topics[i] <- cor(combined_weights[,2], combined_weights[,3])
  }
  expect_gt(mean(corr_topics), 0.7)
})

