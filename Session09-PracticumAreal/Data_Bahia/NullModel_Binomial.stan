data {
  int<lower=1> N;  // Number of observations
  int<lower=1> p;  // Number of predictors
  matrix[N, p] X;  // Covariate matrix
  int<lower=0> y[N];  // Number of individuals who consumed UPF
  int<lower=0> n[N];  // Number of accompanied individuals
}

parameters {
  real alpha;  // Intercept
  vector[p] beta;  // Coefficients for predictors
}

transformed parameters {
  vector[N] mu;  // Linear predictor
  vector[N] prob;  // Probabilities for binomial distribution

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
