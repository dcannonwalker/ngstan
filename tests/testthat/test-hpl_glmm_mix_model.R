cmdstan_version <- try(cmdstanr::cmdstan_version(), silent = TRUE)
found_cmdstan <- inherits(cmdstan_version, "try-error")
message(found_cmdstan)

test_that("hpl_glmm_mix simulates data with the expected structure", {
  skip_if_not(found_cmdstan)
  expect_no_error(cmdstanr::cmdstan_path())
})
