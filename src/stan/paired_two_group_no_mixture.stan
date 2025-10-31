data {
  int<lower=1> G; // number of groups
  int<lower=1> N_g; // observations per group
  vector[N_g] X_g;
  int<lower=1> U; // number of pairs
  matrix[N_g, U] Z_g;
  array[G, N_g] int<lower=0> y;
  int<lower=0, upper=1> run_estimation;
  real<lower=0> a_sig2; // shape for ALL inv_gamma priors
  real<lower=0> b_sig2; // scale for ALL inv_gamma priors
  real<lower=0> a_mu; // shape for ALL inv_gamma priors
  real<lower=0> b_mu; // scale for ALL inv_gamma priors
  real<lower=0> a_u; // shape for ALL inv_gamma priors
  real<lower=0> b_u; // scale for ALL inv_gamma priors
  real<lower=0> a_offset; // shape for ALL inv_gamma priors
  real<lower=0> b_offset; // scale for ALL inv_gamma priors
  real m; // prior mean for mu_offset
}
parameters {
  array[G] real beta;
  array[G] vector[U] u;
  array[G] real log_offset; // log scale offset or intercept
  real mu;
  real<lower=0> sig2;
  real mu_offset;
  real<lower=0> sig2_offset;
  real<lower=0> sig2_mu;
  real<lower=0> sig2_u;
}
transformed parameters {
  array[G] vector[N_g] log_lambda;
  array[G] real lp;
  array[G] real beta_contr;
  array[G] real u_contr;
  for (g in 1:G) {
    log_lambda[g] = log_offset[g] + X_g * beta[g] + Z_g * u[g];
    lp[g] = 0;
    if (run_estimation == 1) {
      lp[g] += poisson_log_lpmf(y[g] | log_lambda[g]);
    }
    beta_contr[g] = normal_lpdf(beta[g] | mu, sig2);
    u_contr[g] = normal_lpdf(u[g] | 0, sig2_u);
  }
}
model {
  mu ~ normal(0, sig2_mu);
  sig2 ~ inv_gamma(a_sig2, b_sig2);
  log_offset ~ normal(mu_offset, sig2_offset);
  mu_offset ~ normal(m, sig2_mu);
  sig2_offset ~ inv_gamma(a_offset, b_offset);
  sig2_mu ~ inv_gamma(a_mu, b_mu);
  sig2_u ~ inv_gamma(a_u, b_u);
  target += sum(lp);
  target += sum(beta_contr);
  target += sum(u_contr);
}
generated quantities {
  array[G, N_g] int y_sim;
  for (g in 1:G) {
    y_sim[g] = poisson_log_rng(log_lambda[g]);
  }
}
