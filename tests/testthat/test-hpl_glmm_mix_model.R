cmdstan_version <- try(cmdstanr::cmdstan_version(), silent = TRUE)
found_cmdstan <- inherits(cmdstan_version, "try-error")
test_that("test skips if no cmdstan install on system", {
  skip_if_not(found_cmdstan)
})

# setup args ----

# test that simulated response has expected structure ----
test_that("hpl_glmm_mix simulates data with the expected structure", {
  skip_if_not(found_cmdstan)
  expect_equal(dim(y_sim_single), c(args$G, nrow(args$X_g)))
})

# fit the model to the simulated data set ----

# test that model fits simulated data without error
test_that("model fits without error", {
  skip_if_not(found_cmdstan)
  expect_no_error(do.call(run_hpl_glmm_mix_model, args = sample_args))
})
