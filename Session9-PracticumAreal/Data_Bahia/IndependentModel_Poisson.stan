data {
  int<lower=1> N;
  int<lower=1> p;
  matrix[N, p] X;
  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] E;           // exposure 
}

transformed data {
  vector[N] log_E = log(E);
}

parameters {
  real beta0;                     // intercept
  vector[p] beta;                 // coefficients for fixed effects
  real<lower=0> sigma_theta;      // Standard deviation of unstructured random effects
  vector[N] theta_raw;            // Standardized random effects
}

transformed parameters {
  vector[N] theta = log_E + beta0 + X * beta + sigma_theta * theta_raw;
}

model {
  y ~ poisson_log(theta); 

  theta_raw ~ normal(0.0, 1.0); 

  for (j in 1:p) {
    beta[j] ~ normal(0.0, 1.0);
  }
  
  beta0 ~ normal(0.0, 1.0);
  sigma_theta ~ cauchy(0.0, 1.0);
}

generated quantities {
  vector[N] mu = exp(theta);
  vector[N] log_lik;

  for (i in 1:N) {
    log_lik[i] = poisson_lpmf(y[i] | mu[i]);
  }
}