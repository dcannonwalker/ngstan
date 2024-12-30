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
