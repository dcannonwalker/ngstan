#' Bayesian hierarchical Poisson GLMM with discrete mixture prior for
#' detection of differential translational efficiency in RiboSeq studies
#'
#' @export
#' @param G Number of groups
#' @param X_g Fixed effects design matrix for a single group
#' @param Z_g Random effects design matrix for a single group
#' @param y Numeric vector of output values, sorted by group
#' @param a_p scale parameter for Beta prior on `p`
#' @param b_p shape parameter for Beta prior on `p`
#' @param a_sig2 vector of shape parameters for Inverse Gamma priors on `sig2`
#' @param b_sig2 vector of scale parameters for Inverse Gamma priors on `sig2`
#' @param sig2_mu vector of variances for the Gaussian prior on `mu`
#' @param sig2_u vector of variances for the Gaussian prior on random effects
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
dte_stan <- function(G, X_g, Z_g, y,
                     a_p = 1, b_p = 2,
                     a_sig2 = NULL, b_sig2 = NULL,
                     sig2_mu = NULL, sig2_u = NULL, ...) {
    K <- ncol(X_g)
    U <- ncol(Z_g)
    N_g <- length(y) / G
    a_sig2 <- a_sig2 %||% rep(10, K)
    b_sig2 <- b_sig2 %||% rep(0.1, K)
    sig2_mu <- sig2_mu %||% rep(1, K)
    sig2_u <- sig2_u %||% rep(0.1, K)
    standata <- list(G = G, N_g = N_g, K = K, U = U,
                     X_g = X_g, Z_g = Z_g, y = y,
                     a_p = a_p, b_p = b_p,
                     a_sig2 = a_sig2, b_sig2 = b_sig2,
                     sig2_mu = sig2_mu,
                     sig2_u = sig2_u)
    out <- rstan::sampling(stanmodels$dte, data = standata, ...)
    return(out)
}
