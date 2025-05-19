data {
  int<lower=1> N;  // Number of observations
  int<lower=1> N_edges;  // Number of edges (pairs of neighbours)
  int<lower=1> p;  // Number of predictors
  matrix[N, p] X;  // Covariate matrix
  array[N_edges] int<lower=1, upper=N> node1;  // node1[i] adjacent to node2[i]
  array[N_edges] int<lower=1, upper=N> node2;  // and node1[i] < node2[i]
  array[N] int<lower=0> y;  // Number of individuals who consumed UPF
  array[N] int<lower=0> n;  // Number of accompanied individuals
}

parameters {
  real alpha;  // Intercept
  vector[p] beta;  // Coefficients for predictors
  real<lower=0> sigma_s;  // Overall standard deviation of spatial effects
  vector[N] s;  // Spatial effects
}

transformed parameters {
  array[N] real<lower=0,upper=1> prob;  // Probabilities for binomial distribution

  for (i in 1:N) {
    prob[i] = inv_logit(alpha + dot_product(X[i], beta) + s[i] * sigma_s);
  }
}

model {
  // Priors
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  sigma_s ~ cauchy(0, 1);

  // Likelihood
  for (i in 1:N) {
    y[i] ~ binomial(n[i], prob[i]);
  }

  // CAR prior for spatial effects
 
    target += -0.5 * dot_self(s[node1] - s[node2]);
  sum(s) ~ normal(0, 0.001 * N);  // Soft sum-to-zero constraint
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = binomial_logit_lpmf(y[i] | n[i], alpha + dot_product(X[i], beta) + s[i]*sigma_s);
  }
}
