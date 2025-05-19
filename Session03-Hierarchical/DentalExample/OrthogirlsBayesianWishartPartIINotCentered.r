library("rjags") 
library("coda")
library("nlme")

data(Orthodont)
Orthgirl <- Orthodont[Orthodont$Sex=="Female",]
attach(Orthgirl)
i <- rep(4,11); x <- 1:11
Subject <- as.numeric(Subject)

zero.u<-c(0,0)
R.u<-matrix(ncol=2,nrow=2)
R.u[1,1] <- 1
R.u[2,2] <- 1
R.u[1,2] <- 0
R.u[2,1] <- 0


jags.data <- list(y=distance,x=I(age-11),group=Subject,n=length(distance),
                  ngirls=11,zero.u=zero.u,R.u=R.u)
jags.params<-c("beta0","beta1","b0","b1","sigma_epsilon",
               "sigma00","sigma11","rhob")

jags.inits <- function(){
  list("beta0"=rnorm(1),"beta1"=rnorm(1),"b0"=rnorm(11,1),
       "b1"=rnorm(11,1),tau0=0.5,tau_epsilon=0.5)
}


#saving the code for the model in a file
cat("
    model{
    
    ## Specify likelihood
    for(i in 1:n){
     y[i] ~ dnorm(mu[i],tau_epsilon)
     mu[i] <- beta0+(beta1+b1[group[i]])*x[i]+b0[group[i]]
    }
    ## Specify priors
    beta0 ~ dnorm(0, 0.0001)
    beta1 ~ dnorm(0, 0.0001)
    for( j in 1:ngirls )
    {
      u[j,1:2] ~ dmnorm(zero.u,invSigma.u)
      b0[j]<-u[j,1]
      b1[j]<-u[j,2]
    }
    invSigma.u ~ dwish(R.u,5)
    Sigma.u <- inverse(invSigma.u)
    sigma00<-Sigma.u[1,1]
    sigma11<-Sigma.u[2,2]
    rhob<-Sigma.u[2,1]/sqrt(Sigma.u[1,1]*Sigma.u[2,2])
    tau0~dunif(0,100)#dgamma(2,0.1)
    sigma20<-sqrt(1/tau0)
    tau_epsilon~dunif(0,100)#dgamma(2,0.1)
    sigma_epsilon<-sqrt(1/tau_epsilon)
    }
    ",file="ModelGirlsBayesianWishartNotCentered.txt")

#creates an object representing a Bayesian graphical model, specified with a 
#BUGS-language description of the prior distribution, and a set of data.

jagsfit<-jags.model("ModelGirlsBayesianWishartNotCentered.txt", data = jags.data, n.chains = 3, n.adapt = 2000)

samps <- coda.samples(jagsfit, jags.params, n.iter = 15000,thin=15)
burn.in <- 5000
summary(window(samps, start = burn.in,thin=10))
plot(window(samps, start = burn.in,thin=10))
autocorr.plot(window(samps, start = burn.in,thin=10))

par(mfrow=c(1,1))
crosscorr.plot(window(samps, start = burn.in,thin=10))
gelman.diag(window(samps, start = burn.in,thin=10))

m <- as.matrix(window(samps,start = burn.in,thin=10))
par(mfrow=c(1,1))
boxplot(m[,1:11])

print(summary(window(samps, start = burn.in,thin=10))



#frequentist approach
lmefit<-summary(lme( distance ~ I(age-11), data = Orthgirl, random = ~1+I(age-11) | Subject ))
#comparison between the fits
m <- as.matrix(window(samps, start = burn.in,thin=10))
par(mfrow=c(1,1))
boxplot(m[,1:11])
points(1:11,lmefit$coefficients$random$Subject[,1],col=2,pch=19,lwd=2)

boxplot(m[,12:21])
points(1:11,lmefit$coefficients$random$Subject[,2],col=2,pch=19,lwd=2)


boxplot(m[,23],ylim=c(19,27))
points(lmefit$coefficients$fixed[1],col=2,lwd=2,pch=19)

boxplot(m[,24])
points(lmefit$coefficients$fixed[2],col=2,lwd=2,pch=19)
