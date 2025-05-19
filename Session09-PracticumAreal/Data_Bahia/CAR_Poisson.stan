data {
  int<lower=1> N;
  int<lower=1> N_edges;
  int<lower=1> p;
  matrix[N,p] X;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] E;           // exposure 
}

transformed data {
  vector[N] log_E = log(E);
}

parameters {
  real beta0;            // intercept
  vector[p] beta;
  real<lower=0> sigma_s;        // standard deviation of spatial effects
  vector[N] s;         // spatial effects
}

model {
  y ~ poisson_log(log_E + beta0 + X * beta + s * sigma_s); 
  // This is the prior for s! (up to proportionality)
  target += -0.5 * dot_self(s[node1] - s[node2]);
  sum(s) ~ normal(0, 0.001 * N);
  
  for(j in 1:p){
    beta[j] ~ normal(0.0, 1);
  }
  
  beta0 ~ normal(0.0, 1);
  sigma_s ~ cauchy(0.0, 1);
}

generated quantities {
  vector[N] mu=exp(log_E + beta0 + X * beta + s * sigma_s);
  vector[N] log_lik;
  
  for(i in 1:N){
    log_lik[i] = poisson_lpmf(y[i] | mu[i]);
  }
}
