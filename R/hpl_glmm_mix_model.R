#' @title Fit the hpl_glmm_mix model
#' @export
#' @family models
#' @description Fit the hpl_glmm_mix Stan model and return fit object
#' @return An object of class `CmdStanMCMC`
#' @param standata A list of data to be passed to the Stan model,
#' matching the format of the `standata` element of a `seqlist` object
#' @param use_multithread Use the multithread-enabled version of the Stan model?
#' @param grainsize Grainsize for multithread
#' @param run_estimation one of `c(0, 1)`; if 0, samples from the prior only
#' and ignores the data
#' @param method One of `c("sample", "vb", "pathfinder")`
#' @param ... Named arguments to the `sample()` method of CmdStan model
#'   objects: <https://mc-stan.org/cmdstanr/reference/model-method-sample.html>
run_hpl_glmm_mix_model <- function(standata,
                                   method = c("sample", "vb", "pathfinder"),
                                   run_estimation = FALSE,
                                   use_multithread = FALSE,
                                   grainsize = NULL, ...) {
  method <- match.arg(method)
  standata[["run_estimation"]] <- run_estimation

  if (!run_estimation) {
    standata[["y"]] <- standata[["y"]] %||% matrix( # nolint
      nrow = standata$G, ncol = standata$N_g, byrow = TRUE,
      rpois(standata$G * standata$N_g, 10)
    ) # allow null counts if run_estimation == 0
  } else {
    if (is.null(standata[["y"]])) {
      stop("'y' cannot be NULL if run_estimation != 0")
    }
  }


  if (use_multithread) {
    grainsize <- grainsize %||% 1 # nolint
    standata[["grainsize"]] <- grainsize
    model <- instantiate::stan_package_model(
      name = "hpl_glmm_mix_multithread",
      package = "ngstan"
    )
  } else {
    model <- instantiate::stan_package_model(
      name = "hpl_glmm_mix",
      package = "ngstan"
    )
  }

  fit <- switch(method,
                sample = model$sample(data = standata, ...),
                vb = model$variational(data = standata, ...),
                pathfinder = model$pathfinder(data = standata, ...)
  )
  return(fit)
}
