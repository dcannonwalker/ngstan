#' @title Fit the hpl_glmm_mix_multithread model
#' @export
#' @family models
#' @description Fit the hpl_glmm_mix Stan model and return fit object
#' @return An object of class `CmdStanMCMC`
#' @param G Number of groups
#' @param X_g Fixed effects design matrix for a single group
#' @param Z_g Random effects design matrix for a single group
#' @param y Numeric vector of output values, sorted by group
#' @param normfactors_known Use fixed normalization factors extrinsic to
#' the model?
#' @param S_DATA Numeric vector of normalization factors
#' for each sample; optional, if `NULL` and `normfactors_known == TRUE`,
#' normalization factors will be estimated by the `TMM` method as
#' described in the `{edgeR}` package
#' @param comps Matrix encoding bernoulli mixture priors
#' @param prob Vector giving the prior probability of each combination of
#' bernoulli components, i.e. of each row in `comps`
#' @param run_estimation one of `c(0, 1)`; if 0, samples from the prior only
#' and ignores the data
#' @param a_sig2 vector of shape parameters
#' for Inverse Gamma priors on `sig2`
#' @param b_sig2 vector of scale parameters
#' for Inverse Gamma priors on `sig2`
#' @param a_sig2_mu vector of shape parameters
#' for Inverse Gamma priors on `sig2_mu`
#' @param b_sig2_mu vector of scale parameters
#' for Inverse Gamma priors on `sig2_mu`
#' @param a_mu_offset vector of shape parameters
#' for Inverse Gamma priors on `mu_offset`
#' @param b_mu_offset vector of scale parameters
#' for Inverse Gamma priors on `mu_offset`
#' @param a_sig2_offset vector of shape parameters
#' for Inverse Gamma priors on `sig2_offset`
#' @param b_sig2_offset vector of scale parameters
#' for Inverse Gamma priors on `sig2_offset`
#' @param a_sig2_u vector of shape parameters
#' for Inverse Gamma priors on `sig2_u`
#' @param b_sig2_u vector of scale parameters
#' for Inverse Gamma priors on `sig2_u`
#' @param method One of `c("sample", "vb", "pathfinder")`
#' @param grainsize Grainsize for multithread
#' @param ... Named arguments to the `sample()` method of CmdStan model
#'   objects: <https://mc-stan.org/cmdstanr/reference/model-method-sample.html>
run_hpl_glmm_mix_mt_model <- function(method = c("sample", "vb", "pathfinder"),
                                      run_estimation = 0,
                                      grainsize = 1,
                                      G, X_g, Z_g, comps, prob,
                                      y = NULL,
                                      a_sig2 = NULL, b_sig2 = NULL,
                                      a_sig2_mu = NULL, b_sig2_mu = NULL,
                                      a_mu_offset = NULL, b_mu_offset = NULL,
                                      a_sig2_offset = NULL,
                                      b_sig2_offset = NULL,
                                      a_sig2_u = NULL, b_sig2_u = NULL,
                                      normfactors_known = FALSE,
                                      S_DATA = NULL, ...) {
  method <- match.arg(method)
  N_g <- nrow(X_g)
  K <- ncol(X_g)
  U <- ncol(Z_g)
  N_comps <- nrow(comps)
  which_mix <- sort(unique(unlist(apply(comps, 1, function(r) which(r == 0)))))
  N_mix <- length(which_mix)
  comps_per_mix <- N_comps / 2
  mix_idx <- matrix(nrow = N_mix, ncol = comps_per_mix)
  for (i in seq(1, N_mix)) {
    mix_idx[i, ] <- which(apply(comps, 1, function(r) r[which_mix[i]] == 1))
  }
  if (run_estimation == 0) {
    y <- y %||% matrix( # nolint
      nrow = G, ncol = N_g, byrow = TRUE,
      rpois(G * N_g, 10)
    ) # allow null y if run_estimation == 0
  } else {
    if (is.null(y)) {
      stop("y cannot be NULL if run_estimation != 0")
    }
  }
  if (normfactors_known) {
    # TODO: use calc_norm_factors() instead
    S_DATA <- S_DATA %||% rep(0, N_g) # nolint
  } else {
    S_DATA <- numeric(0)
  }

  standata <- list(
    G = G, X_g = X_g, Z_g = Z_g,
    N_g = N_g, K = K, U = U,
    comps = comps, N_comps = N_comps, N_mix = N_mix,
    comps_per_mix = N_comps / 2,
    mix_idx = mix_idx, prob = prob,
    run_estimation = run_estimation,
    grainsize = grainsize,
    y = y,
    a_sig2 = a_sig2 %||% rep(10, K), # nolint
    b_sig2 = b_sig2 %||% rep(1, K), # nolint
    a_sig2_mu = a_sig2_mu %||% rep(10, K), # nolint
    b_sig2_mu = b_sig2_mu %||% rep(10, K), # nolint
    a_mu_offset = a_mu_offset %||% 5, # nolint
    b_mu_offset = b_mu_offset %||% 1, # nolint
    a_sig2_offset = a_sig2_offset %||% 10, # nolint
    b_sig2_offset = b_sig2_offset %||% 1, # nolint
    a_sig2_u = a_sig2_u %||% rep(10, U), # nolint
    b_sig2_u = b_sig2_u %||% rep(1, U), # nolint
    normfactors_known = normfactors_known,
    S_DATA = S_DATA
  )

  model <- instantiate::stan_package_model(
    name = "hpl_glmm_mix_multithread",
    package = "ngstan",
    compile = TRUE,
    cpp_options = list(stan_threads = TRUE)
  )

  fit <- switch(method,
    sample = model$sample(data = standata, ...),
    vb = model$variational(data = standata, ...),
    pathfinder = model$pathfinder(data = standata, ...)
  )
  return(fit)
}
