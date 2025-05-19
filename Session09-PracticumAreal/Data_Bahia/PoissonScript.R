library(dplyr)
library(ggplot2)
library(geobr)
library(corrplot)
library(spdep)
library(SpatialEpi)
library(MCMCvis)
library(readxl)
library(rstan)
library(loo)
library(INLA)
#---------------- Spatial effects --------------

df_ba <- readRDS("DataframeBahia_Poisson.rds")

W <- readRDS("W_Matrix_Bahia.rds")
W <- as.matrix(W)
N = length(unique(df_ba$code_micro))
adjlist = function(W,N){ 
  adj=0
  for(i in 1:N){  
    for(j in 1:N){  
      if(W[i,j]==1){adj = append(adj,j)}
    }
  }
  adj = adj[-1]
  return(adj)
}

mungeCARdata4stan = function(adjBUGS,numBUGS) {
  N = length(numBUGS);
  nn = numBUGS;
  N_edges = length(adjBUGS) / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  iAdj = 0;
  iEdge = 0;
  for (i in 1:N) {
    for (j in 1:nn[i]) {
      iAdj = iAdj + 1;
      if (i < adjBUGS[iAdj]) {
        iEdge = iEdge + 1;
        node1[iEdge] = i;
        node2[iEdge] = adjBUGS[iAdj];
      }
    }
  }
  return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2));
}

#----- Adjusting the spatial structure to use in Stan ------

neigh = adjlist(W, N)
numneigh = apply(W,2,sum)
nbs = mungeCARdata4stan(neigh, numneigh)
N = nbs$N ; node1 = nbs$node1 ; node2 = nbs$node2 ; N_edges = nbs$N_edges

#------- Scaling factor based on graph -------
Q_INLA = Diagonal(N, rowSums(W)) - W # precision Q
Q_pert = Q_INLA + Diagonal(N) * max(diag(Q_INLA)) * sqrt(.Machine$double.eps)
Q_inv_INLA = inla.qinv(Q_pert, constr=list(A = matrix(1,1,N),e=0)) # Q^- : general inverse of Q, under linear constraint
scaling_factor = exp(mean(log(diag(Q_inv_INLA))))

#------ Coding on Stan ------

#-------------- Null Model ----------------------
y = df_ba$obesity ; E = df_ba$offset ; X=matrix(c(df_ba$bfp, df_ba$health_teams, df_ba$illit,
                                                  df_ba$percap_income),N, 4) ; p=4

modelNULL = stan("./NullModel_Poisson.stan", 
                 data=list(N=N,y=y,E=E, X=X, p=p),warmup=5000, iter=15000, 
                 chains=2, thin = 5)

# Check convergence:
#pairs(SSB_fit, pars = c("beta", "sigma", "lambda"))
rstan::traceplot(modelNULL, pars = c("beta"))
samplesNULL=extract(modelNULL)

# WAIC
log_lik_matrixNULL = extract_log_lik(modelNULL)
waic(log_lik_matrixNULL)

# Betas
betasNULL = samplesNULL$beta
beta_meansNULL=exp(apply(betasNULL, 2, mean));beta_quantilesNULL=t(exp(apply(betasNULL, 2, quantile, probs = c(0.025, 0.975))))

(cbind(beta_meansNULL, beta_quantilesNULL))

#-------- Independent model ----------------
y = df_ba$obesity ; E = df_ba$offset ; X=matrix(c(df_ba$bfp, df_ba$health_teams, df_ba$illit,
                                                  df_ba$percap_income),N, 4) ; p=4

modelIND = stan("./IndependentModel_Poisson.stan", 
                data=list(N=N,y=y,E=E, X=X, p=p),warmup=5000, iter=15000, 
                chains=2, thin = 5)

# Check convergence:
#pairs(SSB_fit, pars = c("beta", "sigma", "lambda"))
rstan::traceplot(modelIND, pars = c("beta", "theta[1]", "theta[6]", "theta[32]",
                                    "sigma_theta"))
samplesIND=extract(modelIND)

# WAIC
log_lik_matrixIND = extract_log_lik(modelIND)
waic(log_lik_matrixIND)

# Betas
betasIND = samplesIND$beta
beta_meansIND=exp(apply(betasIND, 2, mean));beta_quantilesIND=t(exp(apply(betasIND, 2, quantile, probs = c(0.025, 0.975))))

(cbind(beta_meansIND, beta_quantilesIND))
#-------------- CAR ---------------
N=N ; N_edges=N_edges ; node1=node1 ; node2=node2
y = df_ba$obesity ; E = df_ba$offset ; X=matrix(c(df_ba$bfp, df_ba$health_teams, df_ba$illit,
                                                  df_ba$percap_income),N, 4) ; p=4

modelCAR = stan("./CAR_Poisson.stan", 
                data=list(N=N,N_edges=N_edges,node1=node1,node2=node2,y=y,E=E, 
                          X=X, p=p),warmup=5000, iter=15000, chains=2, thin = 5)

# Check convergence:
#pairs(SSB_fit, pars = c("beta", "sigma", "lambda"))
rstan::traceplot(modelCAR, pars = c("beta", "s[1]", "s[14]", "s[31]",
                                    "sigma_s"))
samplesCAR=extract(modelCAR)

# WAIC
log_lik_matrixCAR = extract_log_lik(modelCAR)
waic(log_lik_matrixCAR)

# Betas
betasCAR = samplesCAR$beta
beta_meansCAR=exp(apply(betasCAR, 2, mean));beta_quantilesCAR=t(exp(apply(betasCAR, 2, quantile, probs = c(0.025, 0.975))))

(cbind(beta_meansCAR, beta_quantilesCAR))

#------------------- BYM ------------------
N=N ; N_edges=N_edges ; node1=node1 ; node2=node2
y = df_ba$obesity ; E = df_ba$offset ; X=matrix(c(df_ba$bfp, df_ba$health_teams, df_ba$illit,
                                                  df_ba$percap_income),N, 4) ; p=4

modelBYM = stan("./BYM_Model_Poisson.stan", 
                data=list(N=N,N_edges=N_edges,node1=node1,node2=node2,y=y,E=E, 
                          X=X, p=p),warmup=2500, iter=15000, chains=2, thin = 5)

# Check convergence:
#pairs(SSB_fit, pars = c("beta", "sigma", "lambda"))
rstan::traceplot(modelBYM, pars = c("beta", "s[1]", "s[12]", "s[26]",
                                    "sigma_s"))
samplesBYM=extract(modelBYM)

# WAIC
log_lik_matrixBYM = extract_log_lik(modelBYM)
waic(log_lik_matrixBYM)

# Betas
betasBYM = samplesBYM$beta
beta_meansBYM=exp(apply(betasBYM, 2, mean));beta_quantilesBYM=t(exp(apply(betasBYM, 2, quantile, probs = c(0.025, 0.975))))

(cbind(beta_meansBYM, beta_quantilesBYM))
#------------------- BYM2 ------------------
N=N ; N_edges=N_edges ; node1=node1 ; node2=node2 ; scaling_factor = scaling_factor
y = df_ba$obesity ; E = df_ba$offset ; X=matrix(c(df_ba$bfp, df_ba$health_teams, df_ba$illit,
                                                  df_ba$percap_income),N, 4) ; p=4

modelBYM2 = stan("./BYM2_Model_Poisson.stan", 
                 data=list(N=N,N_edges=N_edges,node1=node1,node2=node2,y=y,E=E, 
                           X=X, p=p, scaling_factor=scaling_factor), warmup=2500, 
                 iter=15000, chains=2, thin = 5)

# Check convergence:
#pairs(SSB_fit, pars = c("beta", "sigma", "lambda"))
rstan::traceplot(modelBYM2, pars = c("beta", "s[3]", "s[6]", "s[18]"))
samplesBYM2=extract(modelBYM2)

# WAIC
log_lik_matrixBYM2 = extract_log_lik(modelBYM2)
waic(log_lik_matrixBYM2)

# Betas
betasBYM2 = samplesBYM2$beta
beta_meansBYM2=exp(apply(betasBYM2, 2, mean));beta_quantilesBYM2=t(exp(apply(betasBYM2, 2, quantile, probs = c(0.025, 0.975))))

(cbind(beta_meansBYM2, beta_quantilesBYM2))

# Mixing parameter

lambdaBYM2 = samplesBYM2$lambda
summary(lambdaBYM2)