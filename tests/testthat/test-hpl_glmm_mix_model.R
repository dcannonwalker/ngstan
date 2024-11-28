cmdstan_version <- try(cmdstanr::cmdstan_version(), silent = TRUE)
found_cmdstan <- inherits(cmdstan_version, "try-error")
test_that("test skips if no cmdstan install on system", {
  skip_if_not(found_cmdstan)
})

# setup args ----
args <- list(
  G = 5,
  comps = rbind(1, c(1, 1, 0)),
  prob = c(0.2, 0.8),
  X_g = cbind(
    rep(c(0, 1), each = 4),
    rep(c(0, 1, 0, 1), each = 2), c(rep(0, 6), rep(1, 2))
  ),
  Z_g = rbind(diag(1, 4), diag(1, 4)),
  run_estimation = 0,
  iter_warmup = 10,
  iter_sampling = 10,
  chains = 1,
  seed = 1
)

# test that simulated response has expected structure ----
test_that("hpl_glmm_mix simulates data with the expected structure", {
  skip_if_not(found_cmdstan)
  expect_equal(dim(y_sim_single), c(args$G, nrow(args$X_g)))
})

# fit the model to the simulated data set ----
sample_args <- args
sample_args$y <- c(t(y_sim_single))
sample_args$run_estimation <- 1
sample_args$iter_warmup <- 10
sample_args$iter_sampling <- 10
sample_args$chains <- 1
sample_args$parallel_chains <- 1
sample_args$seed <- 2

# test that model fits simulated data without error
testthat("model fits without error", {
  skip_if_not(found_cmdstan)
  expect_no_error(do.call(run_hpl_glmm_mix_model, args = sample_args))
})
