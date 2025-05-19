data {
  int<lower=1> N;  // Number of observations
  int<lower=1> p;  // Number of predictors
  matrix[N, p] X;  // Covariate matrix
  array[N] int<lower=0> y;  // Number of individuals who consumed UPF
  array[N] int<lower=0> n;  // Number of accompanied individuals
}

parameters {
  real alpha;  // Intercept
  vector[p] beta;  // Coefficients for predictors
}

transformed parameters {
  vector[N] mu;  // Linear predictor
  array[N] real<lower=0,upper=1> prob;  // Probabilities for binomial distribution

  for (i in 1:N) {
    mu[i] = alpha + dot_product(X[i], beta);
    prob[i] = inv_logit(mu[i]);
  }
}

model {
  // Priors
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);

  // Likelihood
  for (i in 1:N) {
    y[i] ~ binomial(n[i], prob[i]);
  }
}

generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    log_lik[i] = binomial_lpmf(y[i] | n[i], prob[i]);
  }
}
