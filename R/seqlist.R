#' @param fixed_design Fixed effects design matrix for a single group
#' @param random_design Random effects design matrix for a single group
#' @param counts Numeric matrix of output values
#' @param mixture_probabilities Vector giving the prior probabilities that
#' each regression parameter is drawn from the 0-component of the mixture
#' prior
#' @export
seqlist <- R6::R6Class(
  "seqlist",
  public = list(
    counts = NULL,
    tags = NULL,
    fixed_design = NULL,
    random_design = NULL,
    mixture_probabilities = NULL,
    standata = NULL,
    initialize = function(counts, tags = NULL,
                          fixed_design = NULL,
                          random_design = NULL,
                          mixture_probabilities = NULL) {
      self$counts <- as.matrix(counts)
      self$tags <- tags
      self$fixed_design <- fixed_design
      self$random_design <- random_design %||% matrix(data = numeric(), ncol = 0)
      if (!is.null(fixed_design)) {
        mixture_probabilities <- mixture_probabilities %||% rep(1, ncol(fixed_design))
      }
    },
    set_fixed_design = function(fixed_design) {
      self$fixed_design <- fixed_design
    },
    set_random_design = function(fixed_design) {
      self$random_design <- fixed_design
    },
    set_mixture_probabilities = function(mixture_probabilities) {
      self$mixture_probabilities <- mixture_probabilities
    },
    initialize_standata = function(...) {
      if (is.null(self$fixed_design)) {
        message("Not intitializing 'standata' because 'fixed_design' not set")
      }
      self$standata <- make_standata(...)
    },
    #' Set the value of one of the modifiable elements of `standata`
    #' @param beta_phi_prior 2-vector giving location and scale parameters for
    #' the prior distribution of `beta_phi`
    #' @param S_DATA Numeric vector of normalization factors
    #' for each sample; optional, if `NULL` and `normfactors_known == TRUE`,
    #' normalization factors will be estimated by the `TMM` method as
    #' described in the `{edgeR}` package
    #' @param A_S location parameter for `S_PARAM` if `!normfactors_known`
    #' @param B_S scale parameter for `S_PARAM` if `!normfactors_known`
    #' @param a_sig2 vector of shape parameters
    #' for Inverse Gamma priors on `sig2`
    #' @param b_sig2 vector of scale parameters
    #' for Inverse Gamma priors on `sig2`
    #' @param a_sig2_mu vector of shape parameters
    #' for Inverse Gamma priors on `sig2_mu`
    #' @param b_sig2_mu vector of scale parameters
    #' for Inverse Gamma priors on `sig2_mu`
    #' @param a_mu_offset vector of shape parameters
    #' for Inverse Gamma priors on `mu_offset`
    #' @param b_mu_offset vector of scale parameters
    #' for Inverse Gamma priors on `mu_offset`
    #' @param a_sig2_offset vector of shape parameters
    #' for Inverse Gamma priors on `sig2_offset`
    #' @param b_sig2_offset vector of scale parameters
    #' for Inverse Gamma priors on `sig2_offset`
    #' @param a_sig2_u vector of shape parameters
    #' for Inverse Gamma priors on `sig2_u`
    #' @param b_sig2_u vector of scale parameters
    #' for Inverse Gamma priors on `sig2_u`
    #' @param normfactors_known Use fixed normalization factors extrinsic to
    #' the model?
    #' @param use_neg_binomial_response Use the negative binomial distribution
    #' instead of Poisson for the errors?
    set_standata <- function(a_sig2 = NULL,
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

    }
  ),
  private = list(
    make_standata = function(...) {
      standata <- make_mandatory_standata()
      standata <- c(standata, make_optional_priors(standata, ...))
      return(standata)
    },
    make_mandatory_standata = function() {
      use_mixture_prior <- self$mixture_probabilities != 1
      comps <- ngstan:::build_mixture_indicator(use_mixture_prior) # nolint
      prob <- ngstan:::build_mixture_probabilities(mixture_probabilities, comps) # nolint
      N_comps <- nrow(comps)
      which_mix <- which(use_mixture_prior)
      N_mix <- length(which_mix)
      comps_per_mix <- floor(N_comps / 2)
      mix_idx <- matrix(nrow = N_mix, ncol = comps_per_mix)
      if (N_mix > 0) {
        for (i in seq(1, N_mix)) {
          mix_idx[i, ] <- which(apply(comps, 1, function(r) r[which_mix[i]] == 1))
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
        prob = prob,
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

      normfactors_known <- normfactors_known %||% FALSE # nolint
      use_neg_binomial_response <- use_neg %||% FALSE # nolint
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
    }
  )
)
