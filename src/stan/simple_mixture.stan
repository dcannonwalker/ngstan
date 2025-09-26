data {
  int<lower=1> G;
  int<lower=1> N_g;
  vector[N_g] X_g;
  array[G, N_g] int<lower=0> y;
  int<lower=0, upper=1> run_estimation;
  real<lower=0> a;
  real<lower=0> b;
  real m;
}
transformed data {
  array[G, 2] real w;
  w[, 1] = rep_array(1, G);
  w[, 2] = rep_array(0, G);
}
parameters {
  array[G] real beta;
  array[G] real log_offset; // log scale offset or intercept
  real mu;
  real<lower=0> sig2;
  real mu_offset;
  real<lower=0> sig2_offset;
  real<lower=0> sig2_mu;
}
transformed parameters {
  // array[G] vector[2] w;
  // w[, 1] = beta;
  // w[, 2] = rep_array(0, G);
  array[G, 2] vector[N_g] log_lambda;
  array[G] vector[2] lp;
  array[G] real lse;
  for (g in 1:G) {
    for (i in 1:2) {
      // real wgi = w[g, i];
      log_lambda[g, i] = log_offset[g] + X_g * w[g, i] * beta[g];
      lp[g, i] = log(0.5);
      if (run_estimation == 1) {
        lp[g, i] += poisson_log_lpmf(y[g] | log_lambda[g, i]);
      }
    }
    lse[g] = log_sum_exp(lp[g]);
  }
}
model {
  mu ~ normal(0, sig2_mu);
  sig2 ~ inv_gamma(a, b);
  log_offset ~ normal(mu_offset, sig2_offset);
  mu_offset ~ normal(m, 1);
  sig2_offset ~ inv_gamma(a, b);
  sig2_mu ~ inv_gamma(a, b);
  for (g in 1:G) {
    beta[g] ~ normal(mu, sig2);
  }
  target += sum(lse);
}
generated quantities {
  array[G] real p_dg;
  array[G] vector[2] d_pmf;
  array[G] real p_dg2;
  array[G] vector[2] d_pmf2;
  array[G, N_g] int y_sim;
  array[G] int which_comp;
  for (g in 1:G) {
    vector[2] numerator = exp(lp[g]);
    vector[2] logdiffs;
    logdiffs[1] = lp[g][2] - lp[g][1];
    logdiffs[2] = lp[g][1] - lp[g][2];
    real denominator = sum(exp(lp[g]));
    p_dg[g] = numerator[2] / denominator;
    d_pmf[g] = numerator / denominator;
    p_dg2[g] = 1 / (1 + exp(logdiffs[2]));
    d_pmf2[g][1] = 1 / (1 + exp(logdiffs[1]));
    d_pmf2[g][2] = p_dg2[g];
    which_comp[g] = categorical_rng(rep_vector(0.5, 2));
    y_sim[g] = poisson_log_rng(log_lambda[g, which_comp[g]]);
  }
}
