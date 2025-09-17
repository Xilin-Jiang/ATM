test_that("each step of the inference", {
  # ATM object initialization
  set.seed(19940110)
  para <- topic_init_age(rec_data = HES_age_example, ds_list=UKB_349_disease, topic_num=10, degree_free_num= 5)
  # update the zn; each step should increase the lower bound
  lb1 <- CVB_lb(para)
  para <- CVB0_E_zn(para)
  lb2 <- CVB_lb(para)
  expect_gt(lb2, lb1)
  # update the beta
  para <- fast_update_age_depend_lda(para)
  lb3 <- CVB_lb(para)
  expect_gt(lb3, lb2)
})

test_that("inference wrapper", {
  set.seed(19940110)
  HES_age_small_sample <- HES_age_example %>%
    dplyr::slice_sample(prop = 0.1)
  inference_results <- wrapper_ATM(HES_age_small_sample, topic_num = 10, CVB_num = 1)
  expect_gt(inference_results$ELBO_convergence$Lower_bound[2], inference_results$ELBO_convergence$Lower_bound[1])
  expect_gt(inference_results$ELBO_convergence$Lower_bound[5], inference_results$ELBO_convergence$Lower_bound[4])
  disease_list <- inference_results$ds_list %>%
    left_join(disease_info_phecode_icd10, by = c("diag_icd10"="phecode" )) %>%
    dplyr::pull(phenotype)
  topic_id <- 5 # plot the first topic
  plt <- plot_age_topics(disease_names = disease_list,
                  trajs = inference_results$topic_loadings[35:75,,topic_id],
                  plot_title = paste0("topic ", topic_id),
                  top_ds = 7)
  expect_equal(class(plt)[2], class(ggplot())[2])
})

# test_that("estimate topic weights from fixed topic loadings", {
#   topic_weights_results <- loading2weights(HES_age_example)
#   new_weights <- topic_weights_results$topic_weights
#   load("~/Desktop/PROJECTS/Jiang_2023_NG_ATM/Results/Run_2rec_PheCode_age_dependent_K10_P5_rep10.RData")
#   patient_weights <- sweep((para$alpha_z - 1), 1, rowSums(para$alpha_z -1), FUN="/")
#   originial_weights <- data.frame(eid = para$eid, loading = patient_weights)
#   corr_topics <- c()
#   for(i in 1:para$K){
#     combined_weights <- new_weights %>%
#       select(eid, paste0("topic_weights.", i)) %>%
#       left_join(select(originial_weights, eid, paste0("loading.", i)), by = c("eid"))
#     corr_topics[i] <- cor(combined_weights[,2], combined_weights[,3])
#   }
#   expect_gt(mean(corr_topics), 0.7)
# })

