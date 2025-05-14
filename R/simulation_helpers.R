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
