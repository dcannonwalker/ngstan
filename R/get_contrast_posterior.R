#' Get per-draw probability that a given contrast is 0
#'
#' @param contrast numeric vector giving the contrast
#' @param comps taken from `standata` list
#' @param d_pmf posterior draws for `d_mpf`
#' @param beta posterior draws for `beta`
#'
#' @export
get_contrast_posterior <- function(contrast, comps, d_pmf, beta) {
  # TBD: the posterior probability that any contrast is 0 while
  # any beta is drawn from the non-zero component of the mixture
  # is (probably) 0; I'm doing this numerical check to avoid having to
  # demonstrate that analytically at the moment
  contrast_mat <- comps %*% diag(contrast)
  contrast_comp_draws <- aperm(
    apply(beta, 1:3, function(b) contrast_mat %*% b),
    c(2, 3, 4, 1)
  )
  contrast_comp_eval_draws <- contrast_comp_draws == 0

  # TBD: if I'm right, this will always be:
  # d_pmf, for any element of the sample space where all relevent
  # betas are drawn from their zero component
  # 0, otherwise
  weighted_contrast_comp_eval_draws <- contrast_comp_eval_draws * d_pmf

  # `contrast_posterior` is the
  contrast_posterior <- apply(
    weighted_contrast_comp_eval_draws,
    c(1, 2, 3), function(g) sum(g)
  )
}
