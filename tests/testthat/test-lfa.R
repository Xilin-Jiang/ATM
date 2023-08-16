test_that("test lfa inference", {
  para <- topic_init_lfa_cvb(rec_data = HES_age_example, ds_list=UKB_349_disease, topic_num=10)
  # update the zn; each step should increase the lower bound
  lb1 <- CVB_lb_lfa(para)
  para <- CVB0_E_zn(para)
  lb2 <- CVB_lb_lfa(para)
  expect_gt(lb2, lb1)
  # update the beta
  para <- update_beta_lfa(para)
  lb3 <- CVB_lb_lfa(para)
  expect_gt(lb3, lb2)
})
