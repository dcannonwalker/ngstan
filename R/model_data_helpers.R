#' Build the set of mixture prior indicators
#'
#' @param x An indicator vector identifying which
#' regression parameters should use the mixture prior
build_mixture_indicator <- function(x) {
  x_list <- list()
  for (i in seq(1, length(x))) {
    if (x[i] == 0) {
      x_list[[i]] <- 1
    } else if (x[i] == 1) {
      x_list[[i]] <- c(0, 1)
    }
  }
  as.matrix(expand.grid(x_list))
}
