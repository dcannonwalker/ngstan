data {
    int<lower=1> G; //number of groups
    int<lower=1> N_g; // observations per group
    int<lower=1> K; //number of columns in the model matrix X_g
    int<lower=0> U; //number of columns in the random effects matrix Z_g
    int<lower=0> y[N_g * G]; //response
    matrix[N_g, K] X_g; // per group fixed effects design
    matrix[N_g, U] Z_g; // per group random effects design
    real<lower=0> a_p;
    real<lower=0> b_p;
    real<lower=0> a_sig2[K];
    real<lower=0> b_sig2[K];
    real<lower=0> sig2_mu[K];
    real<lower=0> sig2_u[U];
    // int<lower=0, upper=1> model_normalization;
    vector[N_g] norm_factors; // fixed normalization factors
}
transformed data {
    array[G, N_g] int<lower=0> y_g;
    for (g in 1:G) {
        y_g[g] = y[N_g * (g - 1) + 1:N_g * g];
    }
    vector[N_g] S = norm_factors;
}
parameters {
    array[G] vector[K] beta;
    vector[K] mu;
    vector[K] sig2;
    real<lower=0, upper=1> p;
    array[G] vector[U] u;
    // vector[N_g * (1 - model_normalization)] S; // normalization factors
}
transformed parameters {
    array[G] vector[2] lp;
    for (g in 1:G) {
        for (D in 0:1) {
            vector[K] beta_star = beta[g];
            beta_star[K] = beta_star[K] * D;
            lp[g][D + 1] = bernoulli_lpmf(D | p);
            for (i in 1:N_g) {
                real log_lambda;
                log_lambda = X_g[i, ] * beta_star + Z_g[i, ] * u[g] + S[i];
                lp[g][D + 1] = lp[g][D + 1] +
                poisson_log_lpmf(y_g[g][i] | log_lambda);
            }
        }
    }

}
model {
    // if (model_normalization == 1) {
    //     S ~ normal(0, 1);
    // }
    mu ~ normal(0, sig2_mu);
    sig2 ~ inv_gamma(a_sig2, b_sig2);
    p ~ beta(a_p, b_p);
    for (g in 1:G) {
        u[g] ~ normal(0, sig2_u);
        beta[g] ~ normal(mu, sig2);
        target += log_sum_exp(lp[g]);
    }
}
generated quantities {
    array[G] real<lower=0, upper=1> p_g;
    for (g in 1:G) {
        p_g[g] = (1 + exp(lp[g][1] - lp[g][2]) * (1 - p) / p)^-1;
    }
}
