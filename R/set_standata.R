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
#' @param use_neg_binomial_response Use the negative binomial distribution
#' instead of Poisson for the errors?
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
    normfactors_known = NULL,
    use_neg_binomial_response = NULL) {
  if (!is.null(a_sig2)) {
    if (length(a_sig2) !=
      length(standata[["a_sig2"]])) {
      warning(
        "The new value of 'a_sig2' does not have the same dimension"
      )
    }
    standata[["a_sig2"]] <- a_sig2
  }

  if (!is.null(b_sig2)) {
    if (length(b_sig2) !=
      length(standata[["b_sig2"]])) {
      warning(
        "The new value of 'b_sig2' does not have the same dimension"
      )
    }
    standata[["b_sig2"]] <- b_sig2
  }

  if (!is.null(a_sig2_mu)) {
    if (length(a_sig2_mu) !=
      length(standata[["a_sig2_mu"]])) {
      warning(
        "The new value of 'a_sig2_mu' does not have the same dimension"
      )
    }
    standata[["a_sig2_mu"]] <- a_sig2_mu
  }

  if (!is.null(b_sig2_mu)) {
    if (length(b_sig2_mu) !=
      length(standata[["b_sig2_mu"]])) {
      warning(
        "The new value of 'b_sig2_mu' does not have the same dimension"
      )
    }
    standata[["b_sig2_mu"]] <- b_sig2_mu
  }

  if (!is.null(a_mu_offset)) {
    if (length(a_mu_offset) !=
      length(standata[["a_mu_offset"]])) {
      warning(
        "The new value of 'a_mu_offset' does not have the same dimension"
      )
    }
    standata[["a_mu_offset"]] <- a_mu_offset
  }

  if (!is.null(b_mu_offset)) {
    if (length(b_mu_offset) !=
      length(standata[["b_mu_offset"]])) {
      warning(
        "The new value of 'b_mu_offset' does not have the same dimension"
      )
    }
    standata[["b_mu_offset"]] <- b_mu_offset
  }

  if (!is.null(a_sig2_offset)) {
    if (length(a_sig2_offset) !=
      length(standata[["a_sig2_offset"]])) {
      warning(
        "The new value of 'a_sig2_offset' does not have the same dimension"
      )
    }
    standata[["a_sig2_offset"]] <- a_sig2_offset
  }

  if (!is.null(b_sig2_offset)) {
    if (length(b_sig2_offset) !=
      length(standata[["b_sig2_offset"]])) {
      warning(
        "The new value of 'b_sig2_offset' does not have the same dimension"
      )
    }
    standata[["b_sig2_offset"]] <- b_sig2_offset
  }

  if (!is.null(a_sig2_u)) {
    if (length(a_sig2_u) !=
      length(standata[["a_sig2_u"]])) {
      warning("The new value of 'a_sig2_u' does not have the same dimension")
    }
    standata[["a_sig2_u"]] <- a_sig2_u
  }

  if (!is.null(b_sig2_u)) {
    if (length(b_sig2_u) !=
      length(standata[["b_sig2_u"]])) {
      warning("The new value of 'b_sig2_u' does not have the same dimension")
    }
    standata[["b_sig2_u"]] <- b_sig2_u
  }

  if (!is.null(beta_phi_prior)) {
    if (length(beta_phi_prior) !=
      length(standata[["beta_phi_prior"]])) {
      warning(
        "The new value of 'beta_phi_prior' does not have the same dimension"
      )
    }
    standata[["beta_phi_prior"]] <- beta_phi_prior
  }

  if (!is.null(A_S)) {
    if (length(A_S) !=
      length(standata[["A_S"]])) {
      warning("The new value of 'A_S' does not have the same dimension")
    }
    standata[["A_S"]] <- A_S
  }

  if (!is.null(B_S)) {
    if (length(B_S) !=
      length(standata[["B_S"]])) {
      warning("The new value of 'B_S' does not have the same dimension")
    }
    standata[["B_S"]] <- B_S
  }

  if (!is.null(S_DATA)) {
    if (length(S_DATA) !=
      length(standata[["S_DATA"]])) {
      warning(
        "The new value of 'S_DATA' does not have the same dimension"
      )
    }
    standata[["S_DATA"]] <- S_DATA
  }

  if (!is.null(normfactors_known)) {
    if (length(normfactors_known) !=
      length(standata[["normfactors_known"]])) {
      warning(
        "The new value of 'normfactors_known' does not have the same dimension"
      )
    }
    standata[["normfactors_known"]] <- normfactors_known
  }

  if (!is.null(use_neg_binomial_response)) {
    if (length(use_neg_binomial_response) !=
      length(standata[["use_neg_binomial_response"]])) {
      warning(
        "The new value of 'use_neg_binomial_response' does not have the same dimension" # nolint
      )
    }
  }
  return(standata)
}
