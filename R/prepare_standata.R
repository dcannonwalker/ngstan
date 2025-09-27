#' Prepare basic standata for in-development models

#' Prepare standata for `simple_mixture`
#' @param G number of genes
#' @param N_g number of samples
#' @param X_g sample-level design vector
#' @param ... additional values for `standata`
prepare_standata.simple_mixture <- function(G = 10, N_g = 4, X_g = NULL, ...) {
  dots <- list(...)
  X_g <- X_g %||% c(rep(0, ceiling(N_g / 2)), rep(1, floor(N_g / 2)))
  out <- list(
    G = G,
    N_g = N_g,
    X_g = X_g,
  )
}
