test_that("test prediction odds ratio", {
  topic_loadings <- ATM::UKB_HES_10topics
  disease_list <- ATM::UKB_349_disease
  testing_data <- ATM::HES_age_example %>%
    slice(1:10000)
  load("../../comorbidity/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
  para_training <- para
  output_old <- prediction_onebyone(testing_data, ds_list = disease_list, para_training, max_predict = 5)
  new_output <- prediction_OR(testing_data, ds_list = disease_list, topic_loadings =  topic_loadings, max_predict = 5)
  expect_equal(topic_loadings, para_training$pi_beta_basis)
  expect_equal(output_old, new_output$prediction_precision, tolerance = 0.01)

  ##### need to test with non_age LDA as well
  para_training$P <- NULL
  para_training$beta <- topic_loadings[50,,]
  output_old <- prediction_onebyone(testing_data, ds_list = disease_list, para_training, max_predict = 5)
  new_output <- prediction_OR(testing_data, ds_list = disease_list, topic_loadings =  para_training$beta, max_predict = 5)
  expect_equal(output_old, new_output, tolerance = 0.01)
})

