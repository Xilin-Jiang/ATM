test_that("test prediction odds ratio", {
  topic_loadings <- ATM::UKB_HES_10topics
  disease_list <- ATM::UKB_349_disease
  testing_data <- ATM::HES_age_example %>%
    slice(1:10000)
  load("../../comorbidity/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
  para_training <- para
  output_old <- prediction_OR_onebyone(testing_data, ds_list = disease_list, para_training, max_predict = 5)

})
