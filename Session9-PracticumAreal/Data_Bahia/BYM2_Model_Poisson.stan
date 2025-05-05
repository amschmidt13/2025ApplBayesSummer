data {
  int<lower=1> N;
  int<lower=1> N_edges;
  int<lower=1> p;
  matrix[N, p] X;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] E;           // exposure 
  real<lower=0> scaling_factor; // scales the variance of the spatial effects
}

transformed data {
  vector[N] log_E = log(E);
}

parameters {
  real beta0;                     // intercept
  vector[p] beta;                 // coefficients for fixed effects
  real<lower=0> sigma;            // overall standard deviation
  real<lower=0, upper=1> lambda;  // mixing parameter
  vector[N] theta;                // heterogeneous effects
  vector[N] s;                    // spatial effects
}

transformed parameters {
  vector[N] convolved_re; // variance of each component should be approximately equal to 1
  convolved_re = sqrt(1 - lambda) * theta + sqrt(lambda / scaling_factor) * s;
}

model {
  // Likelihood
  y ~ poisson_log(log_E + beta0 + X * beta + convolved_re * sigma);

  // ICAR prior for spatial effects
  target += -0.5 * dot_self(s[node1] - s[node2]);

  // Hard sum-to-zero constraint
  sum(s) ~ normal(0, 0.001 * N);

  // Priors
  for (j in 1:p) {
    beta[j] ~ normal(0.0, 1.0);
  }
  beta0 ~ normal(0.0, 1.0);
  theta ~ normal(0.0, 1.0);
  sigma ~ normal(0.0, 1.0); 
  lambda ~ uniform(0.0, 1.0);
}

generated quantities {
  vector[N] mu = exp(log_E + beta0 + X * beta + convolved_re * sigma);
  vector[N] log_lik;

  for (i in 1:N) {
    log_lik[i] = poisson_lpmf(y[i] | mu[i]);
  }
}