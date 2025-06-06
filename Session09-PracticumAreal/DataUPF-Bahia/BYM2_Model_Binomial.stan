data {
  int<lower=1> N;  // Number of observations
  int<lower=1> N_edges;  // Number of edges
  int<lower=1> p;  // Number of predictors
  matrix[N, p] X;  // Covariate matrix
  array[N_edges] int<lower=1, upper=N> node1;  // node1[i] adjacent to node2[i]
  array[N_edges] int<lower=1, upper=N> node2;  // and node1[i] < node2[i]
  array[N] int<lower=0> y;  // Number of individuals who consumed UPF
  array[N] int<lower=0> n;  // Number of accompanied individuals
  real<lower=0> scaling_factor;  // Scales the variance of the spatial effects
}

parameters {
  real alpha;  // Intercept
  vector[p] beta;  // Coefficients for predictors
  real<lower=0> sigma;  // Overall standard deviation
  real<lower=0, upper=1> lambda;  // Mixing parameter
  vector[N] theta;  // Heterogeneous effects
  vector[N] s;  // Spatial effects
}

transformed parameters {
  vector[N] convolved_re;  // Convolved random effects
  convolved_re = sqrt(1 - lambda) * theta + sqrt(lambda / scaling_factor) * s;
}

model {
  vector[N] mu;

  // Priors
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  sigma ~ cauchy(0, 1);
  lambda ~ beta(1, 1);
  theta ~ normal(0, 1);

  // CAR prior for structured random effects
  target += -0.5 * dot_self(s[node1] - s[node2]);
  sum(s) ~ normal(0, 0.001 * N);

  // Linear predictor
  mu = alpha + X * beta + convolved_re * sigma;

  // Likelihood
  y ~ binomial_logit(n, mu);
}

generated quantities {
  vector[N] mu_pred;
  vector[N] log_lik;

  mu_pred = alpha + X * beta + convolved_re * sigma;
  
  for (i in 1:N) {
    log_lik[i] = binomial_logit_lpmf(y[i] | n[i], mu_pred[i]);
  }
}
