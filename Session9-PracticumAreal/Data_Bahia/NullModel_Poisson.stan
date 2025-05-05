data {
  int<lower=1> N;
  int<lower=1> p;
  matrix[N,p] X;
  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] E;           // exposure 
}

transformed data {
  vector[N] log_E = log(E);
}

parameters {
  real beta0;            // intercept
  vector[p] beta; // coefficients for fixed effects
}


model {
  y ~ poisson_log(log_E + beta0 + X*beta); 
  
  for(j in 1:p){
    beta[j] ~ normal(0.0, 10.0);
  }
  
  beta0 ~ normal(0.0, 10.0);
}

generated quantities {
  vector[N] mu = exp(log_E + beta0 + X * beta);
  vector[N] log_lik;

  for (i in 1:N) {
    log_lik[i] = poisson_lpmf(y[i] | mu[i]);
  }
}
