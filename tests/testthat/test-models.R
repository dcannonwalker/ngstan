# test that the sampler works
args <- list(
  G = 10,
  mixture_probabilities = c(1, 1, 0.2),
  X_g = cbind(
    rep(c(0, 1), each = 16),
    rep(c(0, 1, 0, 1), each = 8), c(rep(0, 24), rep(1, 8))
  ),
  Z_g = rbind(diag(1, 8), diag(1, 8)),
  run_estimation = 0,
  iter_warmup = 100,
  iter_sampling = 100,
  chains = 1,
  seed = 1,
  normfactors_known = FALSE,
  use_neg_binomial_response = TRUE
)
sim <- do.call(run_hpl_glmm_mix_model, args)
test_that("simulation fails with mismatched X and Z", {
  expect_error(sim$draws())
})

args$Z_g <- rbind(diag(1, 16), diag(1, 16))

sim <- do.call(run_hpl_glmm_mix_model, args)
test_that("simulation runs for negative binomial response model", {
  expect_no_error(sim$draws())
})

test_that("multithread gives a speed advantage", {
  args$use_multithread <- FALSE
  start <- Sys.time()
  sim1 <- do.call(run_hpl_glmm_mix_model, args)
  stop <- Sys.time()
  diff1 <- stop - start
  args$grainsize <- 12
  args$use_multithread <- TRUE
  start <- Sys.time()
  sim2 <- do.call(run_hpl_glmm_mix_model, args)
  stop <- Sys.time()
  diff2 <- stop - start
  expect_true(as.numeric(diff2) / as.numeric(diff1) < 0.8)
})

test_that("default normfactors_known works", {
  args$normfactors_known <- TRUE
  sim <- do.call(run_hpl_glmm_mix_model, args)
  draws <- sim$draws()
  S <- posterior::extract_variable_array(draws, "S")
  expect_equal(sum(S), expected = 0)
  expect_error(draws[, , "S_PARAM[1]"])
})

test_that("NULL y is not allowed if run_estimation == TRUE", {
  args$run_estimation <- TRUE
  expect_error(do.call(run_hpl_glmm_mix_model, args))
})
