#' Construct required arguments to [dte_stan()] from
#' a data frame of Ribo-Seq and RNA-Seq read counts
#' @param counts data frame of read counts, with one row per gene,
#'  column geneid, and one column per sample
#'  @param sample_names vector of sample column names in `counts`
#' @param X_g matrix representing the fixed effects design for
#'  counts from a single gene
#' @param Z_g matrix representing the random effects design,
#' i.e. which Ribo-Seq and mRNA samples are paired,
#' for counts from a single gene
#' @export
prep_dte_data <- function(counts, sample_names, X_g, Z_g) {
    G <- nrow(counts)
    y <- c(t(as.matrix(counts[, sample_names])))
    list(
        G = G,
        y = y,
        X_g = X_g,
        Z_g = Z_g
    )
}
