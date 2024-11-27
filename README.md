<!-- badges: start -->
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/dcannonwalker/ngstan/main.svg)](https://results.pre-commit.ci/latest/github/dcannonwalker/ngstan/main)
[![Codecov test coverage](https://codecov.io/gh/dcannonwalker/ngstan/graph/badge.svg)](https://app.codecov.io/gh/dcannonwalker/ngstan)
[![R-CMD-check](https://github.com/dcannonwalker/ngstan/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dcannonwalker/ngstan/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# ngstan: stan models for next generation sequencing data

`ngstan` provides Stan implementations of statistical models for the 
analysis of Next Generation Sequencing (NGS) data. `ngstan` uses the `cmdstanr`
and `instantiate` packages to provide pre-compiled Stan models that can be fit
using a variety of sampling algorithms available with `cmdstanr` (e.g. ADVI and Pathfinder), 
in addition to the NUTS sampler. 

# Documentation

# Installing `ngstan`

The `ngstan` package depends on the R package
[`CmdStanR`](https://mc-stan.org/cmdstanr/) and the command line tool
[`CmdStan`](https://mc-stan.org/users/interfaces/cmdstan), so it is
important to follow these steps in order:

1.  Install the R package [`CmdStanR`](https://mc-stan.org/cmdstanr/).
    [`CmdStanR`](https://mc-stan.org/cmdstanr/) is not on CRAN, so the
    recommended way to install it is
    `install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))`.
2.  Optional: set environment variables `CMDSTAN_INSTALL` and/or
    `CMDSTAN` to manage the
    [`CmdStan`](https://mc-stan.org/users/interfaces/cmdstan)
    installation. 
3.  Install `ngstan` using the R command below.

| Type        | Source     | Command                                                                     |
|-------------|------------|-----------------------------------------------------------------------------|
| Development | GitHub     | `remotes::install_github("dcannonwalker/ngstan")`                            |

<!---
# Citation
-->
