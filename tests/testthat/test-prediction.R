test_that("test prediction odds ratio", {
  set.seed(19940110)
  topic_loadings <- AgeTopicModels::UKB_HES_10topics
  disease_list <- AgeTopicModels::UKB_349_disease
  testing_data <- AgeTopicModels::HES_age_example %>%
    slice(1:10000)
  para_training <- topic_init_age(AgeTopicModels::HES_age_example, ds_list = AgeTopicModels::UKB_349_disease,topic_num = 10, degree_free_num = 5)
  para_training$pi_beta_basis <- AgeTopicModels::UKB_HES_10topics
  output_old <- prediction_onebyone(testing_data, ds_list = disease_list, para_training, max_predict = 5)
  new_output <- prediction_OR(testing_data, ds_list = disease_list, topic_loadings =  topic_loadings, max_predict = 5)
  expect_equal(topic_loadings, para_training$pi_beta_basis)
  expect_equal(output_old, new_output$prediction_precision, tolerance = 0.01)

  ##### need to test with non_age LDA as well
  para_training$P <- NULL
  para_training$beta <- topic_loadings[50,,]
  output_old <- prediction_onebyone(testing_data, ds_list = disease_list, para_training, max_predict = 5)
  new_output <- prediction_OR(testing_data, ds_list = disease_list, topic_loadings =  para_training$beta, max_predict = 5)
  expect_equal(output_old, new_output$prediction_precision, tolerance = 0.01)
})

