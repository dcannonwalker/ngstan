mixture_probabilities <- c(1, 1, 0.9)
use_mixture_prior <- mixture_probabilities != 1
test_that(
  "comps matrix produced by 'build_mixture_indicator()` is correct",
  {
    comps <- build_mixture_indicator(use_mixture_prior)
    expect_equal(comps, rbind(c(1, 1, 0), 1), ignore_attr = TRUE)
  }
)
test_that(
  "prob vector produced by 'build_mixture_probabilities()` is correct",
  {
    comps <- build_mixture_indicator(use_mixture_prior)
    prob <- build_mixture_probabilities(mixture_probabilities, comps)
    expect_equal(prob, c(0.1, 0.9))
  }
)
