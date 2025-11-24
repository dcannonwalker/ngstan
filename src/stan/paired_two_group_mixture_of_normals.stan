data {
  int<lower=1> G; // number of groups
  int<lower=1> N_g; // observations per group
  vector[N_g] X_g;
  int<lower=1> U; // number of pairs
  matrix[N_g, U] Z_g;
  array[G, N_g] int<lower=0> y;
  int<lower=0, upper=1> run_estimation;
  int<lower=0, upper=1> sim_data;
  real<lower=0> a_sig; // shape for ALL inv_gamma priors
  real<lower=0> b_sig; // scale for ALL inv_gamma priors
  real<lower=0> a_mu; // shape for ALL inv_gamma priors
  real<lower=0> b_mu; // scale for ALL inv_gamma priors
  real<lower=0> a_u; // shape for ALL inv_gamma priors
  real<lower=0> b_u; // scale for ALL inv_gamma priors
  real<lower=0> a_offset; // shape for ALL inv_gamma priors
  real<lower=0> b_offset; // scale for ALL inv_gamma priors
  real m; // prior mean for mu_offset
  real M1; // prior mean for first mixture component
  real M2; // prior mean for second mixture component
}
transformed data {
  array[G, 2] real w;
  w[, 1] = rep_array(1, G);
  w[, 2] = rep_array(0, G);
}
parameters {
  array[G] real beta;
  array[G] vector[U] u;
  array[G] real log_offset; // log scale offset or intercept
  real mu1; // mean for the first normal mixture component of beta prior
  real mu2; // mean for the second normal mixture component of beta prior
  real<lower=0> sig; // sd for beta prior
  real mu_offset;
  real<lower=0> sig_offset;
  real<lower=0> sig_mu; // sd for hierarchical prior on mu1 & mu2
  real<lower=0> sig_u;
}
transformed parameters {
  // array[G] vector[2] w;
  // w[, 1] = beta;
  // w[, 2] = rep_array(0, G);
  array[G, 2] vector[N_g] log_lambda;
  array[G] vector[2] lp;
  array[G] real lse;
  array[G] real beta_contr;
  array[G] real u_contr;
  for (g in 1:G) {
    for (i in 1:2) {
      // if i is 1, then w[g, i] is 1 and the treatment effect is included
      log_lambda[g, i] = log_offset[g] + X_g * w[g, i] * beta[g] + Z_g * u[g];
      lp[g, i] = log(0.5);
      if (run_estimation == 1) {
        lp[g, i] += poisson_log_lpmf(y[g] | log_lambda[g, i]);
      }
    }
    lse[g] = log_sum_exp(lp[g]);
    beta_contr[g] = log_sum_exp(normal_lpdf(beta[g] | mu1, sig),
                                normal_lpdf(beta[g] | mu2, sig));
    u_contr[g] = normal_lpdf(u[g] | 0, sig_u);
  }
}
model {
  mu1 ~ normal(M1, sig_mu);
  mu2 ~ normal(M2, sig_mu);
  sig ~ inv_gamma(a_sig, b_sig);
  log_offset ~ normal(mu_offset, sig_offset);
  mu_offset ~ normal(m, sig_mu);
  sig_offset ~ inv_gamma(a_offset, b_offset);
  sig_mu ~ inv_gamma(a_mu, b_mu);
  sig_u ~ inv_gamma(a_u, b_u);
  target += sum(lse);
  target += sum(beta_contr);
  target += sum(u_contr);
}
generated quantities {
  array[G] real p_dg; // mean parameter for bernoulli mixture; the probability of null
  if (sim_data == 1) { // sometimes poisson_log_rng() is a problem
    array[G, N_g] int y_sim;
    array[G] int which_comp;
    for (g in 1:G) {
      which_comp[g] = categorical_rng(rep_vector(0.5, 2));
      y_sim[g] = poisson_log_rng(log_lambda[g, which_comp[g]]);
    }
  }
  for (g in 1:G) {
    vector[2] logdiffs;
    logdiffs[2] = lp[g][1] - lp[g][2];
    p_dg[g] = 1 / (1 + exp(logdiffs[2]));
  }
}
