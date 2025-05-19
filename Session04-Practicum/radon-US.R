library(tidyverse)
library(nible)
library(coda)

radon = read.csv("radon.txt", sep=' ')
ran = radon %>% group_by(county) %>% summarise(n=n()) %>% select(n) %>% range()
nh = radon %>% group_by(county) %>% summarise(n_i=n())
little = nh %>% filter(n_i<=5) %>% nrow()

# number of counties where $n_i$ measurements were collected, for $n_i=1,2, \dots$
ggplot(nh, aes(x=n_i)) + geom_histogram(color = 1, fill = "white", bins = 50) +
  xlab("Num. of measurements within county") + ylab("Number of counties")


# response variable
ggplot(radon, aes(x=log_radon)) +
  geom_histogram(color = 1, fill = "white", bins = 50)


# we can create a new variable in our
# dataset, called `Category`, which assumes values "low", "medium", or "high", 
# and it is such that the counties are grouped with respect to the uranium levels.
r = range(radon$log_uranium)
br = seq(r[1], r[2], length.out = 4)
baseline_cat = rep(0,nrow(radon))
for(i in 1:(length(br)-1)) baseline_cat[radon$log_uranium >= br[i]] <- i
baseline_cat[baseline_cat==1]="Low"
baseline_cat[baseline_cat==2]="Medium"
baseline_cat[baseline_cat==3]="High"
radon$Category = factor(baseline_cat, levels = c("Low", "Medium", "High"))

appo = radon %>% group_by(Category) %>%
  summarise(num_counties=length(unique(county)))
appo$lower_bound = br[1:(length(br)-1)]
appo$upper_bound = br[2:length(br)]

ggplot(radon, aes(x=floor , y=log_radon)) + geom_point() +
  geom_smooth(method=lm, se=F) +
  facet_wrap(~Category, nrow = 1, labeller = label_both) +
  scale_x_continuous(breaks = 0:1)
# increase in the radon level as the uranium level increases
# however, levels of radon are higher in the basement on average, 
# given a certain level of uranium



#####################################
##  linear mixed-effects model (LMM)
##  model with varying intercept
#####################################
codeM1 = nimbleCode({
  ## Specify likelihood
  for(i in 1:n){
    y[i] ~ dnorm(mu[i],tau_epsilon)
    mu[i] <- b0[group[i]] + beta1*floor[i] + beta2*logu[i]
    fitted[i] ~ dnorm(mu[i],tau_epsilon)
  }
  ## Specify priors
  beta0 ~ dnorm(0, 0.0001)
  beta1 ~ dnorm(0, 0.0001)
  beta2 ~ dnorm(0, 0.0001)
  for(j in 1:m){
    b0[j] ~ dnorm(beta0, tau0)
  }
  tau0~dgamma(0.1,0.1)
  sigma20<-1/tau0
  tau_epsilon~dgamma(0.1,0.1)
  sigma2_epsilon<-1/tau_epsilon
})

      


data<-list(y = radon$log_radon,
           floor = radon$floor, logu = radon$log_uranium)
constants = list(group = as.numeric(as.factor(radon$county)),
                 n = nrow(radon), m = length(unique(radon$county)))
initials<-function(){
  list(beta0 = rnorm(1), beta1 = rnorm(1), beta2 = rnorm(1),
       b0=rnorm(constants$m),
       tau0=rgamma(1,2,.1), tau_epsilon=rgamma(1,2,.1),
       fitted = rnorm(constants$n, 0, sd(na.omit(data$y))))
}

params<-c("beta0","beta1","beta2", "b0","sigma20","sigma2_epsilon",
          "fitted")


modelM1 <- nimbleModel(codeM1, constants = constants, data = data, inits = initials())
modelM1$initializeInfo()
cModelM1 <- compileNimble(modelM1)
confM1 <- configureMCMC(modelM1, monitors = params, enableWAIC = T)
MCMC.M1 <- buildMCMC(confM1)
cMCMC.M1 <- compileNimble(MCMC.M1, project = cModelM1)
samplesM1 <- runMCMC(cMCMC.M1, niter = 10000, nburnin = 5000, thin = 5,
                     nchains = 2, 
                     WAIC = TRUE, samplesAsCodaMCMC = TRUE)

samps = samplesM1$samples


para_extended = colnames(samps[[1]])
summary(samps[, grep('beta.*', para_extended)])
summary(samps[, 'sigma2_epsilon'])
summary(samps[, 'sigma20'])

sumM1<-summary(samps[, grep('b0.*', para_extended)])
index<-seq(1,constants$m,1)

ymin<-min(sumM1$quantiles[,1])
ymax<-max(sumM1$quantiles[,5])
plot(index,sumM1$quantiles[,3],bty="n",ylim=c(ymin,ymax),pch=19,cex=0.7,xlab="County")
for(i in 1:constants$m)
  segments(index[i],sumM1$quantiles[i,1],index[i],sumM1$quantiles[i,5])  



#####################################
##  linear mixed-effects model (LMM)
##  model with varying intercept and slope
#####################################
codeM2 = nimbleCode({
  ## Specify likelihood
  for(i in 1:n){
    y[i] ~ dnorm(mu[i],tau_epsilon)
    mu[i] <- b0[group[i]] + b1[group[i]]*floor[i] + beta2*logu[i]
    fitted[i] ~ dnorm(mu[i],tau_epsilon)
  }
  ## Specify priors
  beta0 ~ dnorm(0, 0.0001)
  beta1 ~ dnorm(0, 0.0001)
  beta2 ~ dnorm(0, 0.0001)
  for(j in 1:m){
    b0[j] ~ dnorm(beta0, tau0)
    b1[j] ~ dnorm(beta1, tau1)
  }
  tau0~dgamma(0.1,0.1)
  sigma20<-1/tau0
  tau1~dgamma(0.1,0.1)
  sigma21<-1/tau1
  tau_epsilon~dgamma(0.1,0.1)
  sigma2_epsilon<-1/tau_epsilon
})




data<-list(y = radon$log_radon,
           floor = radon$floor, logu = radon$log_uranium)
constants = list(group = as.numeric(as.factor(radon$county)),
                 n = nrow(radon), m = length(unique(radon$county)))
initials<-function(){
  list(beta0 = rnorm(1), beta1 = rnorm(1), beta2 = rnorm(1),
       b0=rnorm(constants$m), b1=rnorm(constants$m),
       tau0=rgamma(1,2,.1), tau1=rgamma(1,2,.1), 
       tau_epsilon=rgamma(1,2,.1),
       fitted = rnorm(constants$n, 0, sd(na.omit(data$y))))
}

params<-c("beta0","beta1","beta2", "b0","b1","sigma20","sigma21","sigma2_epsilon",
          "fitted")


modelM2 <- nimbleModel(codeM2, constants = constants, data = data, inits = initials())
modelM2$initializeInfo()
cModelM2 <- compileNimble(modelM2)
confM2 <- configureMCMC(modelM2, monitors = params, enableWAIC = T)
MCMC.M2 <- buildMCMC(confM2)
cMCMC.M2 <- compileNimble(MCMC.M2, project = cModelM2)
samplesM2 <- runMCMC(cMCMC.M2, niter = 10000, nburnin = 5000, thin = 5,
                     nchains = 2, 
                     WAIC = TRUE, samplesAsCodaMCMC = TRUE)

samps = samplesM2$samples


para_extended = colnames(samps[[1]])
summary(samps[, grep('beta.*', para_extended)])
summary(samps[, 'sigma2_epsilon'])
summary(samps[, 'sigma20'])
summary(samps[, 'sigma21'])

sumM2<-summary(samps[, grep('b0.*', para_extended)])
index<-seq(1,constants$m,1)

ymin<-min(sumM2$quantiles[,1])
ymax<-max(sumM2$quantiles[,5])
plot(index,sumM2$quantiles[,3],bty="n",ylim=c(ymin,ymax),pch=19,cex=0.7,xlab="County",
     main = "Random intercepts")
for(i in 1:constants$m)
  segments(index[i],sumM2$quantiles[i,1],index[i],sumM2$quantiles[i,5])  


sumM2<-summary(samps[, grep('b1.*', para_extended)])
index<-seq(1,constants$m,1)

ymin<-min(sumM2$quantiles[,1])
ymax<-max(sumM2$quantiles[,5])
plot(index,sumM2$quantiles[,3],bty="n",ylim=c(ymin,ymax),pch=19,cex=0.7,xlab="County",
     main = "Random slopes")
for(i in 1:constants$m)
  segments(index[i],sumM2$quantiles[i,1],index[i],sumM2$quantiles[i,5])  



# compare models
samplesM1$WAIC
samplesM2$WAIC




