test_that("cmdstan", {
  cmdstan_version <- try(cmdstanr::cmdstan_version(), silent = TRUE)
  found_cmdstan <- !inherits(cmdstan_version, "try-error")
  message(found_cmdstan)
  message(cmdstan_version)
  skip_if_not(found_cmdstan)
  path <- cmdstanr::cmdstan_default_path() %||%
    cmdstanr::cmdstan_default_path(old = TRUE)
  cmdstanr::set_cmdstan_path()
  expect_no_error(cmdstanr::cmdstan_path())
})
