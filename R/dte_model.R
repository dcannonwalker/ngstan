#' @title Fit the DTE model
#' @export
#' @family models
#' @description Fit the DTE Stan model and return fit object
#' @return An object of class `CmdStanMCMC`
#' @param G Number of groups
#' @param X_g Fixed effects design matrix for a single group
#' @param Z_g Random effects design matrix for a single group
#' @param y Numeric vector of output values, sorted by group
#' @param norm_factors Numeric vector of normalization factors
#' for each sample
#' @param a_p scale parameter for Beta prior on `p`
#' @param b_p shape parameter for Beta prior on `p`
#' @param a_sig2 vector of shape parameters for Inverse Gamma priors on `sig2`
#' @param b_sig2 vector of scale parameters for Inverse Gamma priors on `sig2`
#' @param sig2_mu vector of variances for the Gaussian prior on `mu`
#' @param sig2_u vector of variances for the Gaussian prior on random effects
#' @param method One of `c("sampling", "vb", "pathfinder")`
#' @param ... Named arguments to the `sample()` method of CmdStan model
#'   objects: <https://mc-stan.org/cmdstanr/reference/model-method-sample.html>
run_dte_model <- function(G, X_g, Z_g, y,
           norm_factors,
           a_p = 1, b_p = 2,
           a_sig2 = NULL, b_sig2 = NULL,
           sig2_mu = NULL, sig2_u = NULL,
           method = c("sampling", "vb", "pathfinder"), ...) {
    # construct standata
    K <- ncol(X_g)
    U <- ncol(Z_g)
    N_g <- length(y) / G
    a_sig2 <- a_sig2 %||% rep(10, K)
    b_sig2 <- b_sig2 %||% rep(0.1, K)
    sig2_mu <- sig2_mu %||% rep(1, K)
    sig2_u <- sig2_u %||% rep(0.1, U)
    standata <- list(G = G, N_g = N_g, K = K, U = U,
                     X_g = X_g, Z_g = Z_g, y = y,
                     a_p = a_p, b_p = b_p,
                     a_sig2 = a_sig2, b_sig2 = b_sig2,
                     sig2_mu = sig2_mu,
                     sig2_u = sig2_u, norm_factors = norm_factors)

    model <- instantiate::stan_package_model(
      name = "dte",
      package = "ngstan"
    )

    fit <- switch(method,
                  sampling = model$sample(data = standata, ...),
                  vb = model$variational(data = standata, ...),
                  pathfinder = model$pathfinder(data = standata, ...)
    )
    return(fit)
}
