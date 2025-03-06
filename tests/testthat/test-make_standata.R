model <- instantiate::stan_package_model(
  name = "hpl_glmm_mix",
  package = "ngstan"
)

vars <- model$variables()

data_names <- names(vars$data)

