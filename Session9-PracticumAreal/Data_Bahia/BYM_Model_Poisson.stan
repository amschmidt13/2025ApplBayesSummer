data {
  int<lower=1> N;
  int<lower=1> N_edges;
  int<lower=1> p;
  matrix[N, p] X;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] E;           // exposure 
}

transformed data {
  vector[N] log_E = log(E);
}

parameters {
  real beta0;                     // intercept
  vector[p] beta;                 // coefficients for fixed effects
  real<lower=0> sigma_s;          // standard deviation of spatial effects
  real<lower=0> sigma_theta;      // standard deviation of unstructured random effects
  vector[N] s;                    // spatial effects
}

transformed parameters {
  vector[N] theta = log_E + beta0 + X * beta + s * sigma_s;
}

model {
  // Likelihood
  y ~ poisson_log(theta);

  // Prior for theta
  theta ~ normal(0, sigma_theta);

  // ICAR prior for spatial effects
  target += -0.5 * dot_self(s[node1] - s[node2]);

  // Hard sum-to-zero constraint
   sum(s) ~ normal(0, 0.001 * N);

  // Priors
  for (j in 1:p) {
    beta[j] ~ normal(0.0, 1.0);
  }
  beta0 ~ normal(0.0, 10.0);
  sigma_s ~ cauchy(0.0, 1.0);  //
  sigma_theta ~ cauchy(0.0, 1.0);  // 
}

generated quantities {
  vector[N] mu = exp(theta);
  vector[N] log_lik;

  for (i in 1:N) {
    log_lik[i] = poisson_lpmf(y[i] | mu[i]);
  }
}