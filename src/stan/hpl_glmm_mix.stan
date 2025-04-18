functions {
  vector poisson_log_lpmf_mixture(
    array[] int y,
    array[] vector log_lambda,
    array[] real prob,
    int run_estimation
    ) {
      int n_comps = size(log_lambda);
      if (size(prob) != n_comps) {
        fatal_error("size(prob) must match size(log_lambda)");
      }
      vector[n_comps] comps;
      for (i in 1:n_comps) {
        if (size(log_lambda[i]) != size(y)) {
          fatal_error("size(log_lambda[i]) must match size(y)");
        }
        comps[i] = prob[i];
        if (run_estimation == 1) {
          comps[i] += poisson_log_lpmf(y | log_lambda[i]);
        }
      }
      return(comps);
    }
    vector neg_binomial_2_log_lpmf_mixture(
      array[] int y,
      array[] vector log_lambda,
      array[] real prob,
      int run_estimation,
      real phi
      ) {
        int n_comps = size(log_lambda);
        if (size(prob) != n_comps) {
          fatal_error("size(prob) must match size(log_lambda)");
        }
        vector[n_comps] comps;
        for (i in 1:n_comps) {
          if (size(log_lambda[i]) != size(y)) {
            fatal_error("size(log_lambda[i]) must match size(y)");
          }
          comps[i] = prob[i];
          if (run_estimation == 1) {
            comps[i] += neg_binomial_2_log_lpmf(y | log_lambda[i], phi);
          }
        }
        return(comps);
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
}
transformed data {
  array[N_comps] matrix[N_g, K] X_g_comps;
  for (i in 1:N_comps) {
    matrix[N_g, K] rbind_comps;
    for (r in 1:N_g) {
      rbind_comps[r, ] = comps[i];
    }
    X_g_comps[i] = X_g .* rbind_comps;
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
  array[G] vector[N_comps] lp;
  array[G] real lse; // log sum of exponents for each group
  array[G, N_comps] vector[N_g] log_lambda;
  array[use_neg_binomial_response ? G : 0] real phi;
  vector[N_g] S;
  if (normfactors_known) {
    S = S_DATA;
  } else {
    S = S_PARAM;
  }
  for (g in 1:G) {
    for (i in 1:N_comps) {
      log_lambda[g, i] = log_offset[g] +
        X_g_comps[i] * beta[g] +
        Z_g * u[g] +
        S;
    }
    if (use_neg_binomial_response) {
      phi[g] = beta_phi[1] + beta_phi[2] * log_offset[g];
      lp[g] = neg_binomial_2_log_lpmf_mixture(
        y[g], log_lambda[g],
        prob,
        run_estimation,
        phi[g]);
        lse[g] = log_sum_exp(lp[g]);
    } else {
      lp[g] = poisson_log_lpmf_mixture(
        y[g], log_lambda[g],
        prob,
        run_estimation);
        lse[g] = log_sum_exp(lp[g]);
    }

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
  if (use_neg_binomial_response) {
    beta_phi ~ normal(beta_phi_prior[1], beta_phi_prior[2]);
  }
  if (!normfactors_known) {
    S_PARAM ~ normal(A_S, B_S);
  }
  for (g in 1:G) {
    beta[g] ~ normal(mu, sig2);
    u[g] ~ normal(0, sig2_u);
  }
  target += sum(lse);
}
generated quantities {
  array[G] vector[N_mix] p_dg;
  array[G] vector[N_comps] d_pmf;
  for (i in 1:N_mix) {
    for (g in 1:G) {
      vector[N_comps] numerator = exp(lp[g]);
      real denominator = sum(exp(lp[g]));
      p_dg[g][i] = sum(numerator[mix_idx[i]]) / denominator;
      d_pmf[g] = numerator / denominator;
    }
  }
  array[G, N_g] int y_sim;
  array[G] int which_comp;
  for (g in 1:G) {
    which_comp[g] = categorical_rng(to_vector(prob));
    if (use_neg_binomial_response) {
      y_sim[g] = neg_binomial_2_log_rng(log_lambda[g, which_comp[g]],
        phi[g]);
    } else {
      y_sim[g] = poisson_log_rng(log_lambda[g, which_comp[g]]);
    }
  }
}
