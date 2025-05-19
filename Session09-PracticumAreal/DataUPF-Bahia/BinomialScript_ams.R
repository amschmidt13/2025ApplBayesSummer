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

df_ba <- readRDS("DataframeBahia_binomial.rds")

#df_ba<- load("DataframeBahia_Binomial.Rdata")

code_micro<-df_ba$code_micro
obese<-as.numeric(df_ba$obesity)
followed<-df_ba$accompanieds
bfp<-df_ba$bfp
prop_illiterate<-df_ba$illit
avg_percapita_income<-(df_ba$percap_income)
health_teams<-df_ba$health_teams


par(mfrow=c(2,2))
hist(obese/followed,prob=1)
hist(bfp,prob=1)
hist(prop_illiterate,prob=1)
hist(avg_percapita_income,prob=1)

obes_bahia<-data.frame(code_micro=code_micro,
                       obese=obese,followed=followed,
                       bfp=bfp,prop_illiterate=prop_illiterate,
                       avg_percapita_income=avg_percapita_income,
                       health_teams=health_teams)

#write.csv(obes_bahia,"child_obesity_Bahia.txt",row.names=FALSE)



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
y = df_ba$obesity ; n = df_ba$accompanieds ; 
X=matrix(c(df_ba$bfp, df_ba$health_teams, 
           df_ba$illit,df_ba$percap_income),N, 4) ; p=4

modelNULL = stan("NullModel_Binomial.stan", 
                 data=list(N=N,y=y,n=n, X=X, p=p),warmup=5000, iter=15000, 
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
y = df_ba$obesity ; n = df_ba$accompanieds ; X=matrix(c(df_ba$bfp, df_ba$health_teams, df_ba$illit,
                                                  df_ba$percap_income),N, 4) ; p=4

modelIND = stan("./IndependentModel_Binomial.stan", 
                data=list(N=N,y=y,n=n, X=X, p=p),warmup=2000, iter=10000, 
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
y = df_ba$obesity ; n = df_ba$accompanieds ; X=matrix(c(df_ba$bfp, df_ba$health_teams, df_ba$illit,
                                                  df_ba$percap_income),N, 4) ; p=4

modelCAR = stan("./CAR_Binomial.stan", 
                data=list(N=N,N_edges=N_edges,node1=node1,node2=node2,y=y,n=n, 
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
y = df_ba$obesity ; n = df_ba$accompanieds ; X=matrix(c(df_ba$bfp, df_ba$health_teams, df_ba$illit,
                                                  df_ba$percap_income),N, 4) ; p=4

modelBYM = stan("./BYM_Model_Binomial.stan", 
                data=list(N=N,N_edges=N_edges,node1=node1,node2=node2,y=y,n=n, 
                          X=X, p=p, scaling_factor = scaling_factor),
                warmup=2500, iter=10000, chains=2, thin = 5)

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

# Some maps
# Predicted probabilities

prob_pred <- data.frame(prob_pred = apply(samplesBYM$prob_pred, 2, mean))
prob_pred$code_micro <- df_ba$code_micro

prob_pred <- merge(prob_pred, df_ba, by = "code_micro")

(mapProbPred <- 
    ggplot() +
    geom_sf(data = prob_pred, aes(fill = prob_pred, geometry = geom)) +
    scale_fill_viridis_c(name = "Predicted Probability", na.value = "grey") +
    labs(title = NULL) + 
    theme_bw())
#------------------- BYM2 ------------------
N=N ; N_edges=N_edges ; node1=node1 ; node2=node2 ; scaling_factor = scaling_factor
y = df_ba$obesity ; n = df_ba$accompanieds ; X=matrix(c(df_ba$bfp, df_ba$health_teams, df_ba$illit,
                                                  df_ba$percap_income),N, 4) ; p=4

modelBYM2 = stan("./BYM2_Model_Binomial.stan", 
                 data=list(N=N,N_edges=N_edges,node1=node1,node2=node2,y=y,n=n, 
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