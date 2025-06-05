setwd('Session11-PracticumDLM')

rm(list=ls())
library(tidyverse)
library(Rcpp)
library(RcppEigen)
library(R.matlab)
library(tictoc)
library(dlm)
library(coda)
library(bsts)
library(abind)


source("main_function.R")
source("postprocessing.R")


# Define parameters
t <- 500


# generate covariates
set.seed(1)
x1 <- rnorm(t) * 3
set.seed(10)
x2 <- rnorm(t) * 5 + 2
set.seed(11)
z <- rnorm(t) * 2

# simulate beta
tempo <- (1:t)
beta1 <- sin(tempo/t*2*pi) + 3
plot(beta1, type = 'l', main = 'Beta 1')


beta2 = 0.5*sin(tempo/t*2*pi) + 0.5*cos(tempo/t*4*pi)
plot(beta2, type = 'l', main = 'Beta 2')

beta0 <- 0*(seq(0,1,length.out = t)-0.5) - 100*(seq(0,1,length.out = t)-0.5)^2 + 30# time varying intercept (quadratic trend)
plot(beta0, type = 'l', main = 'Beta 0')

# Concatenate
beta_all = array(0, dim = c(1, t, 3))
beta_all[,,1] <- beta0
beta_all[,,2] <- beta1
beta_all[,,3] <- beta2


# generate Y and plot
gamma <- 5
set.seed(42)
eps <- rnorm(t)*2
Y <- beta0 + beta1*x1 + beta2*x2 + gamma*z + eps
plot(Y, type = 'l', main = 'Y', xlab='Time')



# Concatenate arrays
X <- cbind(1, x1, x2)
Z <- cbind(z)

nrep <- 2000
nburn <- 5000
thin <- 1


# package "dlm" ######
myMod <- dlmModReg(cbind(X, Z), addInt = FALSE, 
                   dW = c(rep(1, 3), 0))
tic()
mod0 = dlmGibbsDIG(Y, myMod, 
                       shape.y = 0.01, rate.y = 0.01, 
                       shape.theta = 0.01, rate.theta = 0.01,
                       n.sample = nrep+nburn, thin = thin)
toc() # took 4-5 min on my laptop!

# dynamic coefficients
dc = mod0$theta[-1,,-(1:nburn)] #(t x n. covariates x iter)
for (k in 1:dim(X)[2]) {
  dc_ = data.frame(
    Time = 1:t,
    Bmean = rowMeans(dc[,k,]),
    B025 = apply(dc[,k,], 1, quantile, probs = 0.025),
    B975 = apply(dc[,k,], 1, quantile, probs = 0.975),
    Btrue = beta_all[1,,k]
  )
  gg = ggplot(dc_, aes(x=Time)) +
    geom_line(aes(y= Bmean)) +
    geom_line(aes(y= B025), linetype='dashed') +
    geom_line(aes(y= B975), linetype='dashed') +
    ylab(paste('beta', k-1)) +
    geom_line(aes(y=Btrue), col='red')
  print(gg)
}
# the MCMC needs to run much more



# fitted values
Yfitted = matrix(NA, nrep, t)
for (i in 1:nrep) {
  coeff = mod0$theta[-1,,nburn+i]
  mu_y = rowSums(cbind(1, x1, x2, z) * coeff)
  sd_y = sqrt(dV[i])
  Yfitted[i,] = rnorm(t, mu_y, sd_y)
}
plot(Y, t='l')
lines(colMeans(Yfitted), col=2)
legend('topleft', legend = c('Observed', 'Fitted'), col=1:2, lty=1)
mean((Y-colMeans(Yfitted))^2, na.rm = T)




# package "bsts" #####
ss <- AddLocalLevel(list(), Y)
ss = AddDynamicRegression(ss,Y ~ x1+x2)
tic()
mod1 <- bsts(Y~z-1, state.specification = ss, niter = nrep+nburn, ping=100,
             family = 'gaussian')
toc() # took ~20 sec
plot(mod1, "dynamic", burn = nburn)
plot(mod1, "components", burn = nburn) # contribution of each state component
plot(mod1, "coefficients", burn = nburn)
plot(mod1, "predictors", burn = nburn)
plot(mod1, "residuals", burn = nburn) # posterior distribution of the residuals given complete data (i.e. looking forward and backward in time)

# dynamic coefficients
dc = abind(
  mod1$state.contributions[-(1:nburn), "trend",],
  mod1$dynamic.regression.coefficients[-(1:nburn),,],
  along = 2) #(iter x n. covariates x t)
for (k in 1:dim(X)[2]) {
  dc_ = data.frame(
    Time = 1:t,
    Bmean = colMeans(dc[,k,]),
    B025 = apply(dc[,k,], 2, quantile, probs = 0.025),
    B975 = apply(dc[,k,], 2, quantile, probs = 0.975),
    Btrue = beta_all[1,,k]
  )
  gg = ggplot(dc_, aes(x=Time)) +
    geom_line(aes(y= Bmean)) +
    geom_line(aes(y= B025), linetype='dashed') +
    geom_line(aes(y= B975), linetype='dashed') +
    ylab(paste('beta', k-1)) +
    geom_line(aes(y=Btrue), col='red')
  print(gg)
}

# posterior for gamma
summary(mcmc(mod1$coefficients[-(1:nburn),]))
plot(mcmc(mod1$coefficients[-(1:nburn),]))

# fitted values
Yfitted = Y - colMeans(mod1$one.step.prediction.errors[-(1:nburn),])
plot(Y, t='l')
lines(Yfitted, col=2)
legend('topleft', legend = c('Observed', 'Fitted'), col=1:2, lty=1)
mean((Y-Yfitted)^2, na.rm = T)



# Chan & Jeliazkov (2009) implementation #######

prior_list = list(
  V_beta_0 = 1e4, # Prior variance of initial state
  V_gamma = 1e6,  # Prior variance of constant coefficients
  a1 = 0.01,      # Prior shape for temporal variance
  b1 = 0.01,      # Prior rate for temporal variance
  s2_a = 0.01,    # Prior shape for measurement error variance
  s2_b = 0.01     # Prior rate for measurement error variance
)
print.interval = 100

tic()
mod2 <- DLM(Y, X, Z, nrep, nburn, thin, print.interval, prior_list)
toc() # took 3 sec!
ave <- mod2$ave
out <- mod2$out



plot.tvc(mod2, beta_all)
plot.fitted(mod2, Y)
mean((Y-ave$Ypred_mean[1,])^2, na.rm = T)

chains = posterior.chains(mod2)
summary(chains)
plot(chains)

log_like = out$store_llike
dim(log_like) = c(t, nrep)
loo::waic(t(log_like))


save(mod0, mod1, mod2, file='simulation_dlm_uv.RData')

rm(mod0, mod1, mod2)
gc()


##############################
# Application: Pollution #####

data <- read.delim("pollution_milan.txt")
data$Time = seq(ymd('2018-01-01'), ymd('2022-12-31'), by='1 day')
head(data)


t = nrow(data)

Y = scale(log(data$pm2p5))[,1]
plot(data$Time, Y, t='l')

x1 = scale(data$temp)[,1]
x2 = scale(data$rh)[,1]
X <- cbind(1, x1, x2)

# very simple OLS fit
modOLS = lm(Y~-1+x1+x2)
summary(modOLS)
plot(modOLS)

# increase if necessary!
nrep <- 10000
nburn <- 15000
thin <- 1


# package "bsts" #####
ss <- AddLocalLevel(list(), Y)
ss = AddDynamicRegression(ss,Y ~ x1+x2)
tic()
mod1 <- bsts(Y, state.specification = ss, niter = nrep+nburn, ping=1000,
             family = 'gaussian')
toc()

# some trace plot inspections
plot(mcmc(mod1$sigma.obs[-(1:nburn)]^2))
plot(mcmc(mod1$sigma.level[-(1:nburn)]^2))
plot(mcmc(mod1$x1.sigma[-(1:nburn)]^2))
plot(mcmc(mod1$x2.sigma[-(1:nburn)]^2))


# dynamic coefficients
dc = abind(
  mod1$state.contributions[-(1:nburn), "trend",],
  mod1$dynamic.regression.coefficients[-(1:nburn),,],
  along = 2) #(iter x n. covariates x t)
for (k in 1:dim(X)[2]) {
  dc_ = data.frame(
    Time = data$Time,
    Bmean = colMeans(dc[,k,]),
    B025 = apply(dc[,k,], 2, quantile, probs = 0.025),
    B975 = apply(dc[,k,], 2, quantile, probs = 0.975)
  )
  gg = ggplot(dc_, aes(x=Time)) +
    geom_line(aes(y= Bmean)) +
    geom_line(aes(y= B025), linetype='dashed') +
    geom_line(aes(y= B975), linetype='dashed') +
    ylab(paste('beta', k-1))
  print(gg)
}


# Chan & Jeliazkov (2009) implementation #######
prior_list = list(
  V_beta_0 = 1e4, # Prior variance of initial state
  V_gamma = 1e6,  # Prior variance of constant coefficients
  a1 = 0.01,      # Prior shape for temporal variance
  b1 = 0.01,      # Prior rate for temporal variance
  s2_a = 0.01,    # Prior shape for measurement error variance
  s2_b = 0.01     # Prior rate for measurement error variance
)

tic()
mod2 <- DLM(Y, X, NULL, nrep, nburn, thin, 1000, prior_list)
toc()
ave <- mod2$ave
out <- mod2$out


plot.tvc(mod2)
plot.fitted(mod2, Y)
mean((Y-ave$Ypred_mean[1,])^2, na.rm = T)

chains = posterior.chains(mod2)
summary(chains)
plot(chains)

log_like = out$store_llike
dim(log_like) = c(t, floor(nrep/thin))
loo::waic(t(log_like))
