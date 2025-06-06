#' R6 class representing data from an RNA-Seq or Ribo-Seq experiment
#'
#' @description
#' A `seqlist` has data from the experiment and information about
#' the experimental design and model specification.
#'
#' @export
seqlist <- R6::R6Class(
  "seqlist",
  public = list(
    #' @field counts Numeric matrix of output value, where rows represent genes
    #' or tags and columns represent samples
    counts = NULL,
    #' @field tags A dataframe of tag metadata, matching the row-order of `counts`
    tags = NULL,
    #' @field fixed_design Fixed effects design matrix for a single group
    fixed_design = NULL,
    #' @field random_design Random effects design matrix for a single group
    random_design = NULL,
    #' @field mixture_probabilities Vector giving the prior probabilities that
    #' each regression parameter is drawn from the 0-component of the mixture
    #' prior
    mixture_probabilities = NULL,
    #' @field standata List of data to be passed to Stan
    standata = NULL,
    #' @field fit Fit object returned by `run_model()`
    fit = NULL,
    #' @description Create a new `seqlist`
    #'
    #' @param counts Value for `counts` field
    #' @param tags Value for `tags` field
    #' @param fixed_design Value for `fixed_design` field
    #' @param random_design Value for `random_design` field
    #' @param mixture_probabilities Value for `mixture_probabilities` field
    #' @return A new `seqlist` object
    initialize = function(counts = NULL, tags = NULL,
                          fixed_design = NULL,
                          random_design = NULL,
                          mixture_probabilities = NULL) {
      self$counts <- as.matrix(counts)
      self$tags <- tags
      self$fixed_design <- fixed_design
      self$random_design <- random_design %||% matrix(data = numeric(), ncol = 0) # nolint
      if (!is.null(fixed_design)) {
        mixture_probabilities <- mixture_probabilities %||% rep(1, ncol(fixed_design)) # nolint
      }
    },
    #' set the `counts` field
    #' @param counts the counts to use
    set_counts = function(counts) {
      self$counts <- counts
      invisible(self)
    },
    #' filter the `counts` and `tags` fields
    #' @param keep logical vector indicating which rows to keep
    filter = function(keep) {
      self$counts <- self$counts[keep, ]
      self$tags <- self$tags[keep, ]
      invisible(self)
    },
    #' order the `counts` and `tags` fields
    #' @param index numeric vector ordering the rows
    arrange = function(index) {
      self$counts <- self$counts[index, ]
      self$tags <- self$tags[index, ]
      invisible(self)
    },
    #' set the `fixed_design` field
    #' @param fixed_design the fixed design matrix to use
    set_fixed_design = function(fixed_design) {
      self$fixed_design <- fixed_design
      invisible(self)
    },
    #' set the `random_design` matrix field
    #' @param random_design the random design matrix to use
    set_random_design = function(random_design) {
      self$random_design <- random_design
      invisible(self)
    },
    #' set the `mixture_probabilities` field
    #' @param mixture_probabilities the vector of mixture probabilities
    set_mixture_probabilities = function(mixture_probabilities) {
      self$mixture_probabilities <- mixture_probabilities
      invisible(self)
    },
    #' set the `standata` field
    #' @param ... arguments passed to private methods
    initialize_standata = function(...) {
      if (is.null(self$fixed_design)) {
        message("Not intitializing 'standata' because 'fixed_design' not set")
      } else {
        self$standata <- private$make_standata(...)
        invisible(self)
      }
    },
    #' Set the value of one of the modifiable elements of `standata`
    #' @param ... arguments passed to `ngstan::set_standata()`
    set_standata = function(...) {
      self$standata <- set_standata(self$standata, ...)
      invisible(self)
    },
    #' run the model
    #' @param method the sampling algorithm to use
    #' @param run_estimation include the likelihood in the objective?
    #' @param use_multithread use the multithread-enabled model?
    #' @param grainsize grainsize for multithread
    #' @param modify_in_place add the fit to the `seqlist` object? otherwise
    #' return it
    #' @param ... arguments passed to the `{cmdstanr}` function
    #' identified by `method`
    run_model = function(method = c("sample", "vb", "pathfinder"),
                         run_estimation = FALSE,
                         use_multithread = FALSE,
                         grainsize = NULL,
                         modify_in_place = TRUE,
                         ...) {
      fit <- run_hpl_glmm_mix_model(
        standata = self$standata, method = method,
        run_estimation = run_estimation,
        use_multithread = use_multithread,
        grainsize = grainsize,
        ...
      )
      if (modify_in_place) {
        self$fit <- fit
        invisible(self)
      } else {
        return(fit)
      }
    }
  ),
  private = list(
    make_standata = function(...) {
      standata <- private$make_mandatory_standata()
      standata <- c(standata, private$make_optional_priors(standata, ...))
      return(standata)
    },
    make_mandatory_standata = function() {
      use_mixture_prior <- self$mixture_probabilities != 1
      comps <- private$build_mixture_indicator(use_mixture_prior) # nolint
      prob <- private$build_mixture_probabilities(self$mixture_probabilities, comps) # nolint
      N_comps <- nrow(comps)
      which_mix <- which(use_mixture_prior)
      N_mix <- length(which_mix)
      comps_per_mix <- floor(N_comps / 2)
      mix_idx <- matrix(nrow = N_mix, ncol = comps_per_mix)
      if (N_mix > 0) {
        for (i in seq(1, N_mix)) {
          mix_idx[i, ] <- which(apply(comps, 1, function(r) r[which_mix[i]] == 1)) # nolint
        }
      }

      G <- nrow(self$counts)
      N_g <- nrow(self$fixed_design)
      K <- ncol(self$fixed_design)
      U <- ncol(self$random_design)
      mandatory_standata <- list(
        y = self$counts,
        G = G,
        X_g = self$fixed_design,
        Z_g = self$random_design,
        N_g = N_g,
        K = K,
        U = U,
        comps = comps,
        N_comps = N_comps,
        N_mix = N_mix,
        comps_per_mix = N_comps / 2,
        mix_idx = mix_idx,
        prob = prob
      )
      return(mandatory_standata)
    },
    make_optional_priors = function(mandatory_standata,
                                    a_sig2 = NULL,
                                    b_sig2 = NULL,
                                    a_sig2_mu = NULL,
                                    b_sig2_mu = NULL,
                                    a_mu_offset = NULL,
                                    b_mu_offset = NULL,
                                    a_sig2_offset = NULL,
                                    b_sig2_offset = NULL,
                                    a_sig2_u = NULL,
                                    b_sig2_u = NULL,
                                    beta_phi_prior = NULL,
                                    A_S = NULL,
                                    B_S = NULL,
                                    S_DATA = NULL,
                                    normfactors_known = NULL,
                                    use_neg_binomial_response = NULL) {
      K <- mandatory_standata$K
      U <- mandatory_standata$U
      N_g <- mandatory_standata$N_g

      normfactors_known <- normfactors_known %||% FALSE # nolint
      use_neg_binomial_response <- use_neg_binomial_response %||% FALSE # nolint
      if (normfactors_known) {
        # TODO: use calc_norm_factors() instead
        S_DATA <- S_DATA %||% rep(0, N_g) # nolint
        A_S <- numeric(0)
        B_S <- numeric(0)
      } else {
        S_DATA <- numeric(0)
        A_S <- A_S %||% rep(0, N_g) # nolint
        B_S <- B_S %||% rep(0.1, N_g) # nolint
      }
      if (use_neg_binomial_response) {
        beta_phi_prior <- beta_phi_prior %||% c(0, 1) # nolint
      } else {
        beta_phi_prior <- numeric(0)
      }

      optional_priors <- list(
        a_sig2 = a_sig2 %||% rep(10, K), # nolint
        b_sig2 = b_sig2 %||% rep(1, K), # nolint
        a_sig2_mu = a_sig2_mu %||% rep(10, K), # nolint
        b_sig2_mu = b_sig2_mu %||% rep(10, K), # nolint
        a_mu_offset = a_mu_offset %||% 5, # nolint
        b_mu_offset = b_mu_offset %||% 1, # nolint
        a_sig2_offset = a_sig2_offset %||% 10, # nolint
        b_sig2_offset = b_sig2_offset %||% 1, # nolint
        a_sig2_u = a_sig2_u %||% rep(10, U), # nolint
        b_sig2_u = b_sig2_u %||% rep(1, U), # nolint
        normfactors_known = normfactors_known,
        S_DATA = S_DATA,
        A_S = A_S,
        B_S = B_S,
        use_neg_binomial_response = use_neg_binomial_response,
        beta_phi_prior = beta_phi_prior
      )
      return(optional_priors)
    },
    # @description
    # Build the set of mixture prior indicators
    #
    # @param x An indicator vector identifying which
    # regression parameters should use the mixture prior
    build_mixture_indicator = function(x) {
      x_list <- list()
      for (i in seq(1, length(x))) {
        if (x[i] == 0) {
          x_list[[i]] <- 1
        } else if (x[i] == 1) {
          x_list[[i]] <- c(0, 1)
        }
      }
      as.matrix(expand.grid(x_list))
    },

    # @description
    # Build the set of mixture combination prior probabilities
    #
    # @param x A vector of probabilities
    # @param y A matrix with each row representing a combination of
    # draws from binary random variables
    build_mixture_probabilities = function(x, y) {
      apply(
        y, 1, function(r) {
          prod(r * x + abs(1 - r) * (1 - x))
        }
      )
    }
  )
)
