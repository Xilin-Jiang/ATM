test_that("Simulate non-genetic data", {
  para_sim <- simulate_age_topic_data(sample_sz = 1000, topic_number=3,
                                      disease_number=30, degree_freedom = 3, overlap = 2,
                                      ds_per_idv = 6.1)
  expect_equal(para_sim$M, 1000)
  expect_equal(para_sim$K, 3)
  expect_equal(para_sim$D, 30)
  expect_gt(mean(para_sim$Ns), 5)
  expect_lt(mean(para_sim$Ns), 7)
  expect_gt(dim(para_sim$rec_data)[1]/para_sim$M, 4)
  expect_lt(dim(para_sim$rec_data)[1]/para_sim$M, 6)
  invisible(capture.output(inference_results <- wrapper_ATM(para_sim$rec_data, topic_num = 3, CVB_num = 1)))
})

test_that("Simulate genetic data", {

  cont_v2t <- 20 # 1. number of variants contributing to topic
  disease2topic <- 0.1 # 2. causal disease SNP->D->T
  itr_effect <- 2 # 3. SNP*T->D
  topic2disease <- 2 # 4. SNP->D; SNP->T->D; topic to disease effect
  disease_number <- 20
  rslts <- simulate_topics(topic_number = 2, disease2topic  = disease2topic, v2t = cont_v2t, pop_sz = 1000 )
  para_sim <- rslts[[1]]
  genetics_population <- rslts[[2]]
  causal_disease <- rslts[[3]]
  reslt_ds <- simulate_genetic_disease_from_topic(para_sim, genetics_population,causal_disease,
                                                  disease_number = disease_number, itr_effect = itr_effect,
                                                  topic2disease = topic2disease, v2t = cont_v2t)
  rec_data <- reslt_ds[[1]]
  disease_data <- longdata2diseasematrix(rec_data)
  # testing the SNP -> Disease
  ds1 <- disease_data %>%
    select(eid, disease_number + 2) %>%
    arrange(eid)
  md1 <- lm(ds1[2][[1]] ~ genetics_population[,1:20])
  # test if we could discover 25% of the disease causing variant at P=0.05
  expect_gt(sum(summary(md1)$coefficients[2:21,4] < 0.05), 5)
})
