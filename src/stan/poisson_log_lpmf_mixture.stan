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
