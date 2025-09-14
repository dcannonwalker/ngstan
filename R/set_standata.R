#' Set the value of one of the modifiable elements of `standata`
#' @param standata The data list to be updated
#' @param beta_phi_prior 2-vector giving location and scale parameters for
#' the prior distribution of `beta_phi`
#' @param S_DATA Numeric vector of normalization factors
#' for each sample; optional, if `NULL` and `normfactors_known == TRUE`,
#' normalization factors will be estimated by the `TMM` method as
#' described in the `{edgeR}` package
#' @param A_S location parameter for `S_PARAM` if `!normfactors_known`
#' @param B_S scale parameter for `S_PARAM` if `!normfactors_known`
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
#' @param normfactors_known Use fixed normalization factors extrinsic to
#' the model?
#' @export
set_standata <- function(
    standata,
    a_sig2 = NULL,
    b_sig2 = NULL,
    a_sig2_mu = NULL,
    b_sig2_mu = NULL,
    a_mu_offset = NULL,
    b_mu_offset = NULL,
    a_sig2_offset = NULL,
    b_sig2_offset = NULL,
    a_sig2_u = NULL,
    b_sig2_u = NULL,
    beta_phi_prior = NULL,
    A_S = NULL,
    B_S = NULL,
    S_DATA = NULL,
    normfactors_known = NULL
    ) {
  if (!is.null(a_sig2)) {
    standata[["a_sig2"]] <- a_sig2
  }

  if (!is.null(b_sig2)) {
    standata[["b_sig2"]] <- b_sig2
  }

  if (!is.null(a_sig2_mu)) {
    standata[["a_sig2_mu"]] <- a_sig2_mu
  }

  if (!is.null(b_sig2_mu)) {
    standata[["b_sig2_mu"]] <- b_sig2_mu
  }

  if (!is.null(a_mu_offset)) {
    standata[["a_mu_offset"]] <- a_mu_offset
  }

  if (!is.null(b_mu_offset)) {
    standata[["b_mu_offset"]] <- b_mu_offset
  }

  if (!is.null(a_sig2_offset)) {
    standata[["a_sig2_offset"]] <- a_sig2_offset
  }

  if (!is.null(b_sig2_offset)) {
    standata[["b_sig2_offset"]] <- b_sig2_offset
  }

  if (!is.null(a_sig2_u)) {
    standata[["a_sig2_u"]] <- a_sig2_u
  }

  if (!is.null(b_sig2_u)) {
    standata[["b_sig2_u"]] <- b_sig2_u
  }

  if (!is.null(beta_phi_prior)) {
    standata[["beta_phi_prior"]] <- beta_phi_prior
  }

  if (!is.null(A_S)) {
    standata[["A_S"]] <- A_S
  }

  if (!is.null(B_S)) {
    standata[["B_S"]] <- B_S
  }

  if (!is.null(S_DATA)) {
    standata[["S_DATA"]] <- S_DATA
  }

  if (!is.null(normfactors_known)) {
    standata[["normfactors_known"]] <- normfactors_known
  }

  return(standata)
}
