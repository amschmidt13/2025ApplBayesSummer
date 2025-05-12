rm(list=ls())
#dataset
y <- c(1,4,9,13,18,20,0,2,6,10,12,16)
n <- rep(20,12)
prop <- y/n
dose <- rep(c(1,2,4,8,16,32),2)
ldose <- logb(dose,2)
male <- c(rep(1,6),rep(0,6))



#packages to run under the Bayesian framework
library("nimble") 
library("igraph")
library("coda")
#code in NIMBLE
budwormCode<-nimbleCode({
  ## Specify likelihood
  for(i in 1:12){
    y[i] ~ dbin(p[i], n[i])
    logit(p[i]) <- b0 + b1*log(dose[i])+b2*male[i]
    fitted[i] ~ dbinom(p[i],n[i])
  }
  ## Specify priors
  b0 ~ dnorm(0, sd=100)
  b1 ~ dnorm(0, sd=100)
  b2 ~ dnorm(0, sd=100)
})

#defining constants, data and initial values

budworm.data <- list(y=y,n=n,dose=dose,male=male)#"y","n","dose","male")
budworm.params<-c("b0","b1","b2","fitted[]","p[]")
budworm.inits <-list("b0"=rnorm(1),"b1"=rnorm(1),"b2"=rnorm(1),"fitted"=rbinom(12,20,0.5))

#creates an object representing a Bayesian graphical model, specified with a 
#BUGS-language description of the prior distribution, and a set of data.

budworm <- nimbleModel(code = budwormCode, name = 'budworm',
                    data = budworm.data, inits = budworm.inits)


#graphical model representation
par(mfrow=c(1,1))
plot(budworm$graph)

#showing the stochastic nodes
budworm$getNodeNames(stochOnly=TRUE)
#The default MCMC configuration includes monitors on all top-level stochastic nodes of the
#model.
mcmcConf <- configureMCMC(budworm, monitors = c("b0", "b1","b2","fitted"))
budwormMCMC <- buildMCMC(mcmcConf)


#Running within R (takes much longer than the code below)
# budwormMCMC <- buildMCMC(budworm)
# budwormMCMC$run(500)
# MCMCsamples <- as.matrix(budwormMCMC$mvSamples)
# 
# par(mfrow=c(3,2))
# plot(MCMCsamples[ , 'b0'], type = 'l', xlab = 'iteration',  ylab = "b0")
# hist(MCMCsamples[ , 'b0'],prob=1,ylab='b0')
# plot(MCMCsamples[ , 'b1'], type = 'l', xlab = 'iteration',  ylab = "b1")
# hist(MCMCsamples[ , 'b1'],prob=1,ylab='b1')
# plot(MCMCsamples[ , 'b2'], type = 'l', xlab = 'iteration',  ylab = "b2")
# hist(MCMCsamples[ , 'b2'],prob=1,ylab='b2')


#
#compiling the codes into C++
#

L<-10000
Cbudworm <- compileNimble(budworm)
CbudwormMCMC <- compileNimble(budwormMCMC, project = budworm)
CbudwormMCMC$run(L)
MCMCsamples <- as.matrix(CbudwormMCMC$mvSamples)

summary(MCMCsamples)
sta<-2000
thinning<-8
coda_samples <- mcmc(MCMCsamples)
coda_samples <- window(coda_samples,start=sta,end=L,thin=thinning)
plot(coda_samples)

varnames<-c("fitted[1]","fitted[2]","fitted[3]","fitted[4]","fitted[5]","fitted[6]","fitted[7]","fitted[8]",
            "fitted[9]","fitted[10]","fitted[11]","fitted[12]")

summary_coda<-summary(coda_samples, quantiles = c(0.025, 0.5, 0.975))
summaryfitted<-summary_coda$statistics[varnames,"Mean"]
q025fitted<-summary_coda$quantiles[varnames,"2.5%"]
q975fitted<-summary_coda$quantiles[varnames,"97.5%"]



par(mfrow=c(3,2))
plot(MCMCsamples[ , 'b0'], type = 'l', xlab = 'iteration',  ylab = "b0")
hist(MCMCsamples[ , 'b0'],prob=1,ylab='b0')
plot(MCMCsamples[ , 'b1'], type = 'l', xlab = 'iteration',  ylab = "b1")
hist(MCMCsamples[ , 'b1'],prob=1,ylab='b1')
plot(MCMCsamples[ , 'b2'], type = 'l', xlab = 'iteration',  ylab = "b2")
hist(MCMCsamples[ , 'b2'],prob=1,ylab='b2')

par(mfrow=c(3,2))
plot(MCMCsamples[ , 'fitted[1]'], type = 'l', xlab = 'iteration',  ylab = "b0")
hist(MCMCsamples[ , 'fitted[1]'],prob=1,ylab='b0')
plot(MCMCsamples[ , 'fitted[2]'], type = 'l', xlab = 'iteration',  ylab = "b1")
hist(MCMCsamples[ , 'fitted[2]'],prob=1,ylab='b1')
plot(MCMCsamples[ , 'fitted[3]'], type = 'l', xlab = 'iteration',  ylab = "b2")
hist(MCMCsamples[ , 'fitted[3]'],prob=1,ylab='b2')



#MLE's
fitMLE<-glm(formula = cbind(y, n - y) ~ log(dose) + male,family = binomial)
##### Calculate SEs
X <- cbind(1,log(dose),male)
mu <- fitted(fitMLE)
V.mu <- n*mu*(1-mu)
est.fcn <- X*(y - n*mu)

# Fisher info:
WtVX <- t(X) %*% (V.mu*X)
Iobs.inv <- solve(WtVX)
se.model <- sqrt(diag(Iobs.inv))
round(se.model,4)
Iobs.inv
#
predMLE<-predict.glm(fitMLE,se.fit=TRUE,type="response")

fit.inf1<-n*(predMLE$fit-1.96*predMLE$se.fit)
fit.inf1

fit.sup1<-n*(predMLE$fit+1.96*predMLE$se.fit)
fit.sup1

###############################################
#
#comparing the fitted values
#
###############################################
#predictions <- data.frame(dose,summaryfitted,q025fitted,q975fitted)
#predictions

prds <- data.frame(dose, summaryfitted,q025fitted,q975fitted)

par(mfrow=c(2,1),pty="m")
plot(dose[1:6], y[1:6], cex = 1, col = 1, pch = 19, ylab = "# of Individuals", 
     xlab = "Dose",bty="n",ylim=c(0,20))
lines(prds[1:6, 1], prds[1:6, 2], lwd = 2,lty=2,col="gray70")
lines(prds[1:6, 1], prds[1:6, 3], lwd = 2, col = "gray70")
lines(prds[1:6, 1], prds[1:6, 4], lwd = 2,lty=2,col="gray70")
lines(prds[1:6,1],fitMLE$fitted.values[1:6]*n[1:6],col="red",lty=3,lwd=2)
lines(prds[1:6,1],fit.inf1[1:6],col="red",lty=4,lwd=2)
lines(prds[1:6,1],fit.sup1[1:6],col="red",lty=4,lwd=2)
legend(25,18, legend = c("Male  Bayes","B-95% Male","MLE","MLE-95% Male"), 
       col = c("gray70","gray70","red","red"),
       lty=c(1,2,3,4),lwd = c(2, 2,2,2),bty="n")


plot(dose[7:12], y[7:12], cex = 1, col = 1, pch = 19, ylab = "# of Individuals", 
     xlab = "Dose",bty="n",ylim=c(0,20))
lines(prds[7:12, 1], prds[7:12, 2], lwd = 2,lty=2,col="gray70")
lines(prds[7:12, 1], prds[7:12, 3], lwd = 2, col = "gray70")
lines(prds[7:12, 1], prds[7:12, 4], lwd = 2,lty=2,col="gray70")
lines(prds[7:12,1],fitMLE$fitted.values[7:12]*n[7:12],col="blue",lty=3,lwd=2)
lines(prds[7:12,1],fit.inf1[7:12],col="blue",lty=4,lwd=2)
lines(prds[7:12,1],fit.sup1[7:12],col="blue",lty=4,lwd=2)
legend(25,15, legend = c("Fem  Bayes","B-95% Fem","MLE","MLE-95% Fem"), 
       col = c("gray70","gray70","blue","blue"),
       lty=c(1,2,3,4),lwd = c(2, 2,2,2),bty="n")



