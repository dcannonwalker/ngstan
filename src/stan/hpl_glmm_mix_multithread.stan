functions {
  vector poisson_log_lpmf_mixture(
    array[] int y, // length N_g
    array[] matrix design_comps, // length N_comps, N_g x (K + U)
    vector all_glmm_pars, // length K + U
    real log_offset,
    vector S,
    array[] real prob, // length N_comps
    int run_estimation
    ) {
      int n_comps = size(prob);
      vector[n_comps] comps;
      for (i in 1:n_comps) {
        comps[i] = log(prob[i]);
        if (run_estimation == 1) {
          // comps[i] += poisson_log_lpmf(y | log_lambda[i]);
          // comps[i] += poisson_log_glm_lupmf(y |
          //  design_comps[i], log_offset, all_glmm_pars)
          comps[i] += poisson_log_glm_lpmf(y |
          design_comps[i],
          S + log_offset,
          all_glmm_pars);
        }
      }
      return(comps);
    }
    vector neg_binomial_2_log_lpmf_mixture(
      array[] int y, // length N_g
      array[] matrix design_comps, // length N_comps, N_g x (K + U)
      vector all_glmm_pars, // length K + U
      real log_offset,
      vector S,
      array[] real prob, // length N_comps
      int run_estimation,
      int use_neg_binomial_response,
      real phi
      ) {
        int n_comps = size(prob);
        vector[n_comps] comps;
        for (i in 1:n_comps) {
          comps[i] = log(prob[i]);
          if (run_estimation == 1) {
            comps[i] += neg_binomial_2_log_glm_lpmf(y |
            design_comps[i],
            S + log_offset,
            all_glmm_pars,
            phi);
          }
        }
        return(comps);
      }
      real partial_sum(
        array[,] int y_slice, // G x N_g
        int start, int end,
        // int N_comps,
        array[] real log_offset, // G
        array[] matrix design_comps, // N_comps x N_g x (K + U)
        array[] vector all_glmm_pars, // K + U
        array[] real prob,
        vector S, // N_g
        int run_estimation,
        int use_neg_binomial_response,
        array[] real phi
        ) {
          int L = end - start + 1;
          vector[L] partial_vec;
          array[L] vector[size(prob)] lp;
          array[L] real lse;
          for (g in start:end) {
            int idx = g - start + 1;
            if (use_neg_binomial_response) {
              lp[idx] = neg_binomial_2_log_lpmf_mixture(
                y_slice[idx],
                design_comps, all_glmm_pars[g],
                log_offset[g], S, prob, run_estimation,
                use_neg_binomial_response, phi[g]);
                lse[idx] = log_sum_exp(lp[idx]);
            } else {
              lp[idx] = poisson_log_lpmf_mixture(
                y_slice[idx],
                design_comps, all_glmm_pars[g],
                log_offset[g], S, prob, run_estimation);
                lse[idx] = log_sum_exp(lp[idx]);
            }
          }
          return(sum(lse));
        }
}
data {
  int<lower=1> G;
  int<lower=1> N_g;
  int<lower=0> K;
  int<lower=0> U;
  matrix[N_g, K] X_g;
  matrix[N_g, U] Z_g;
  array[G, N_g] int<lower=0> y;
  int<lower=0, upper=1> normfactors_known; // fixed normalization factors?
  vector[normfactors_known ? N_g : 0] S_DATA;
  array[normfactors_known ? 0 : N_g] real A_S;
  array[normfactors_known ? 0 : N_g] real<lower=0> B_S;
  int<lower=0> N_mix;
  int<lower=1> N_comps;
  array[N_comps] row_vector[K] comps;
  array[N_comps] real<lower=0, upper=1> prob;
  int<lower=1> comps_per_mix;
  array[N_mix, comps_per_mix] int mix_idx;
  array[K] real<lower=0> a_sig2; // shape
  array[K] real<lower=0> b_sig2; // scale
  array[K] real<lower=0> a_sig2_mu; // shape
  array[K] real<lower=0> b_sig2_mu; // scale
  array[U] real<lower=0> a_sig2_u; // shape
  array[U] real<lower=0> b_sig2_u; // scale
  real a_mu_offset; // location; should probably be positive
  real<lower=0> b_mu_offset; // scale
  real<lower=0> a_sig2_offset; // shape
  real<lower=0> b_sig2_offset; // scale
  int<lower=0, upper=1> run_estimation;
  int<lower=0, upper=1> use_neg_binomial_response;
  vector[use_neg_binomial_response ? 2 : 0] beta_phi_prior;
  int grainsize;
}
transformed data {
  array[N_comps] matrix[N_g, K] X_g_comps;
  array[N_comps] matrix[N_g, K + U] design_comps;
  for (i in 1:N_comps) {
    matrix[N_g, K] rbind_comps;
    for (r in 1:N_g) {
      rbind_comps[r, ] = comps[i];
    }
    X_g_comps[i] = X_g .* rbind_comps;
    design_comps[i] = append_col(X_g, Z_g);
  }
}
parameters {
  array[G] vector[K] beta;
  array[G] vector[U] u;
  array[G] real log_offset; // log scale offset or intercept
  vector[K] mu;
  vector<lower=0>[K] sig2;
  real mu_offset;
  real<lower=0> sig2_offset;
  vector<lower=0>[U] sig2_u;
  vector<lower=0>[K] sig2_mu;
  vector[normfactors_known ? 0 : N_g] S_PARAM;
  vector<lower=0>[use_neg_binomial_response ? 2 : 0] beta_phi;
}
transformed parameters {
  array[G] vector[K + U] all_glmm_pars;
  array[use_neg_binomial_response ? G : 0] real phi;
  vector[N_g] S;
  if (normfactors_known) {
    S = S_DATA;
  } else {
    S = S_PARAM;
  }
  for (g in 1:G) {
    if (use_neg_binomial_response) {

      phi[g] = beta_phi[1] + beta_phi[2] * log_offset[g];
    }
    all_glmm_pars[g] = append_row(beta[g], u[g]);
  }
}
model {
  mu ~ normal(0, sig2_mu);
  sig2 ~ inv_gamma(a_sig2, b_sig2);
  log_offset ~ normal(mu_offset, sig2_offset);
  mu_offset ~ normal(a_mu_offset, b_mu_offset);
  sig2_offset ~ inv_gamma(a_sig2_offset, b_sig2_offset);
  sig2_mu ~ inv_gamma(a_sig2_mu, b_sig2_mu);
  sig2_u ~ inv_gamma(a_sig2_u, b_sig2_u);
  if (!normfactors_known) {
    S_PARAM ~ normal(A_S, B_S);
  }
  if (use_neg_binomial_response) {
    beta_phi ~ normal(beta_phi_prior[1], beta_phi_prior[2]);
  }
  for (g in 1:G) {
    beta[g] ~ normal(mu, sig2);
    u[g] ~ normal(0, sig2_u);
  }
  target += reduce_sum(partial_sum, y, grainsize, log_offset,
  design_comps, all_glmm_pars, prob, S, run_estimation,
  use_neg_binomial_response, phi);
}
generated quantities {
  array[G] vector[N_comps] lp;
  array[G, N_comps] vector[N_g] log_lambda;
  for (g in 1:G) {
    for (i in 1:N_comps) {
      log_lambda[g, i] = log_offset[g] + X_g_comps[i] * beta[g] + Z_g * u[g] + S;
      lp[g, i] = log(prob[i]) + poisson_log_lpmf(y[g] | log_lambda[g, i]);
    }
  }
  array[G] vector[N_mix] p_dg;
  array[G] vector[N_comps] d_pmf;
  array[G] vector[N_comps] numerator;
  array[G] real denominator;
  for (i in 1:N_mix) {
    for (g in 1:G) {
      numerator[g] = exp(lp[g]);
      denominator[g] = sum(exp(lp[g]));
      p_dg[g][i] = sum(numerator[g][mix_idx[i]]) / denominator[g];
      d_pmf[g] = numerator[g] / denominator[g];
    }
  }
  array[G, N_g] int y_sim;
  array[G] int which_comp;
  for (g in 1:G) {
    which_comp[g] = categorical_rng(to_vector(prob));
    if (use_neg_binomial_response) {
      y_sim[g] = neg_binomial_2_rng(exp(log_lambda[g, which_comp[g]]),
      phi[g]);
    } else {
      y_sim[g] = poisson_log_rng(log_lambda[g, which_comp[g]]);
    }
  }
  array[N_comps] matrix[N_g, K] X_g_comps_out = X_g_comps;
  matrix[N_g, U] Z_g_out = Z_g;
}
