#' Extract a single simulated data set from a simulation
#' run of the sampler
#' @param sim_draws draws object from a simulation run
#' @param draw which draw to extract
#' @param chain which chain to extract from
get_sim <- function(sim_draws, draw = 1, chain = 1) {
  y_sim <- posterior::extract_variable_array(sim_draws, "y_sim")
  which_comp <- posterior::extract_variable_array(sim_draws, "which_comp")
  all_variables <- sim_draws[draw, chain, ]
  list(
    y = y_sim[draw, chain, , ],
    which_comp = which_comp[draw, chain, ],
    all_variables = all_variables
  )
}

simulation_study <- R6::R6Class(
  "simulation_study",
  public = list(
    #'
    #' @field
    simulation = NULL,
    #'
    #' @field
    fit = NULL,
    #' @description
    #' Create a new `simulation_study`
    #'
    initialize = function(simulation = NULL, fit = NULL) {
      self$simulation <- simulation
      self$fit <- fit
    },
    #' @description
    #' Summarize the comparison between fit and true value
    summarize = function(variable) {
      var <- posterior::extract_variable()
    }

  )
)
dmnm <- dimnames(fitdraws)
unique(stringr::str_remove_all(dmnm$variable, pattern = "\\[.*\\]"))
x = c("abc[1]", "asdc[1,2]")
library(stringr)
str_remove_all(x, pattern = "\\[.*\\]")
