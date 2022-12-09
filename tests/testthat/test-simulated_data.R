test_that("Simulate data", {
  para_sim <- simulate_age_topic_data(sample_sz = 10000, topic_number=3,
                                      disease_number=30, degree_freedom = 3, overlap = 2,
                                      ds_per_idv = 6.1)
})
