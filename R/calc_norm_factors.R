#' Calculate fixed sample normalization factors for RNA-Seq
#' or Ribo-Seq data
#'
#' This is a lightweight re-implementation of the normalization factors
#' calculation from the `edgeR` package
#' @param counts Matrix of read counts with one row per gene and one column
#' per sample
#' @param method Method to use
#' @param ref_column Column id of counts to use as reference
#' @param ... Additional arguments passed to `method` function
#' @seealso [https://code.bioconductor.org/browse/edgeR]
#' @export
calc_norm_factors <- function(counts, method = "TMM", ref_column = NULL,
                              ...) {
  lib_size <- colSums(counts)
  f <- switch(method,
    TMM = {
      if (is.null(ref_column)) {
        f75 <- .calc_factor_quantile(
          x = counts,
          lib_size = lib_size,
          p = 0.75
        )
        if (median(f75) < 1e-20) {
          ref_column <- which.max((colSums(sqrt(counts))))
        } else {
          ref_column <- which.min(abs(f75 - mean(f75)))
        }
      }
      f <- rep_len(NA_real_, ncol(counts))
      for (i in seq(1, ncol(counts))) {
        f[i] <- .calc_factor_tmm(
          obs = counts[, i],
          ref = counts[, ref_column],
          libsize_obs = lib_size[i],
          libsize_ref = lib_size[ref_column], ...
        )
      }
      f
    }
  )
}

.calc_factor_quantile <- function(x, lib_size, p = 0.75, cutoff = -1e10) {
  f <- rep_len(1, ncol(x))
  for (j in seq_len(ncol(x))) {
    f[j] <- quantile(x[, j], probs = p)
  }
  if (min(f) == 0) {
    warning("One or more quantiles are zero")
  }
  f / lib_size
}

.calc_factor_tmm <- function(obs, ref, libsize_obs, libsize_ref,
                             cutoff = -1e10,
                             logratio_trim = 0.3,
                             sum_trim = 0.05,
                             use_weights = TRUE) {
  log_r <- log2((obs / libsize_obs) / (ref / libsize_ref))
  abs_e <- (log2(obs / libsize_obs) + log2(ref / libsize_ref)) / 2
  v <- (libsize_obs - obs) / libsize_obs / obs +
    (libsize_ref - ref) / libsize_ref / ref # delta method var. of log_r

  fin <- is.finite(log_r) & is.finite(abs_e) & (abs_e > cutoff)

  log_r <- log_r[fin]
  abs_e <- abs_e[fin]
  v <- v[fin]

  if (max(abs(log_r)) < 1e-6) {
    return(1)
  }

  n <- length(log_r)
  low_logratio <- floor(n * logratio_trim) + 1
  high_logratio <- n + 1 - low_logratio
  low_sum <- floor(n * sum_trim) + 1
  high_sum <- n + 1 - low_sum

  keep <- (rank(log_r) > low_logratio & rank(log_r) < high_logratio) &
    (rank(abs_e) > low_sum & rank(log_r) < high_sum)

  if (use_weights) {
    f <- sum(log_r[keep] / v[keep], na.rm = TRUE) /
      sum(1 / v[keep], na.rm = TRUE)
  } else {
    f <- mean(log_r[keep], na.rm = TRUE)
  }

  if (is.na(f)) f <- 0
  2^f
}
