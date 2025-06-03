rm(list=ls())

heart=read.table(file="HeartData.txt", header=T)
head(heart)
regression.out= lm(mphr ~ hr + bp + pkhr + sbp + age + baseef + gender, data=heart)
summary(regression.out)


######################################
#
# Bayesian regression using nimble
#
######################################

library("nimble") 
library("igraph")
library("coda")
#code in NIMBLE
heartCode<-nimbleCode({
  ## Specify likelihood
  for(i in 1:n){
    y[i] ~ dnorm(mu[i], sd=sigma)
    mu[i] <- b0 + b1*hr[i] + b2*bp[i] + b3*pkhr[i] + b4*sbp[i] + 
      b5*age[i] + b6*baseef[i] + b7*gender[i]
    fitted[i] ~ dnorm(mu[i],sd=sigma)
  }
  ## Specify priors
  b0 ~ dnorm(0, sd=100)
  b1 ~ dnorm(0, sd=100)
  b2 ~ dnorm(0, sd=100)
  b3 ~ dnorm(0, sd=100)
  b4 ~ dnorm(0, sd=100)
  b5 ~ dnorm(0, sd=100)
  b6 ~ dnorm(0, sd=100)
  b7 ~ dnorm(0, sd=100)
  sigma~dunif(0,100)
  
})

#attaching the data.frame
#do the analysis with and without standardizing 
#the covariates
#notice that the chains of the coefficients do 
# not converge when the covariates are not standardized

sheart<-heart
sheart[,1:7]<-apply(heart[,1:7],2,scale)
names(sheart)
attach(sheart)



#defining constants, data and initial values
heart.const<- list(n=558,hr=hr,bp=bp,pkhr=pkhr,sbp=sbp,age=age,baseef=baseef,gender=gender)
heart.data <- list(y=mphr)
heart.params<-c("b0","b1","b2","b3","b4","b5","b6","b7","sigma")
heart.inits <-list("b0"=rnorm(1),"b1"=rnorm(1),"b2"=rnorm(1),
                   "b3"=rnorm(1),"b4"=rnorm(1),"b5"=rnorm(1),
                   "b6"=rnorm(1),"b7"=rnorm(1),sigma=1)

#creates an object representing a Bayesian graphical model, specified with a 
#BUGS-language description of the prior distribution, and a set of data.

heartB <- nimbleModel(code = heartCode, name = 'heartBayes',
                      data = heart.data, constants = heart.const,
                      inits = heart.inits)

Cheart<-compileNimble(heartB)
heartConf <- configureMCMC(heartB, print = TRUE)
heartMCMC <- buildMCMC(heartConf)
CheartMCMC <- compileNimble(heartMCMC, project = heartB)


niter <- 10000
set.seed(1)
CheartMCMC$run(niter)
## NULL
samples <- as.matrix(CheartMCMC$mvSamples)

par(mfrow = c(2, 2), mai = c(.6, .4, .1, .2))
plot(samples[ , "b0"], type = "l", xlab = "iteration",
     ylab = "b0")
plot(samples[ , "b1"], type = "l", xlab = "iteration",
     ylab = "b1")
plot(samples[ , "b2"], type = "l", xlab = "iteration",
     ylab = "b2")
plot(samples[ , "b3"], type = "l", xlab = "iteration",
     ylab = "b3")
plot(samples[ , "b4"], type = "l", xlab = "iteration",
     ylab = "b4")
par(mfrow=c(2,2))
plot(samples[ , "b5"], type = "l", xlab = "iteration",
     ylab = "b5")
plot(samples[ , "b6"], type = "l", xlab = "iteration",
     ylab = "b6")
plot(samples[ , "b7"], type = "l", xlab = "iteration",
     ylab = "b7")
plot(samples[ , "sigma"], type = "l", xlab = "iteration",
     ylab = "sigma")

heartConf$addSampler(target = c("b0","b1","b2","b3","b4","b5","b6","b7"), type = "RW_block",
                    control = list(adaptInterval = 100))


summary(samples)

######################################
#
# Bayesian regression using jags
#
######################################

library("rjags") 

jags.data <- list(y=mphr,n=558,hr=hr,bp=bp,pkhr=pkhr,sbp=sbp,age=age,baseef=baseef,gender=gender,
                  hr.new=0.048,bp.new=-1.42,pkhr.new=0.25,sbp.new=-.6,age.new=1.56,baseef.new=-.33,gender.new=0)

jags.params<-c("b0","b1","b2","b3","b4","b5","b6","b7","sigma","y.new")
jags.inits <- function(){
  list("b0"=rnorm(1),"b1"=rnorm(1),"b2"=rnorm(1),
       "b3"=rnorm(1),"b4"=rnorm(1),"b5"=rnorm(1),
       "b6"=rnorm(1),"b7"=rnorm(1),sigma=1)
}

#saving the code for the model in a file
cat("
model{

   for(i in 1:n){
    y[i] ~ dnorm(mu[i], prec)
    mu[i] <- b0 + b1*hr[i] + b2*bp[i] + b3*pkhr[i] + b4*sbp[i] + 
      b5*age[i] + b6*baseef[i] + b7*gender[i]
   }
    y.new~dnorm(mu.new,prec)
    mu.new<-b0 + b1*hr.new + b2*bp.new + b3*pkhr.new + b4*sbp.new + b5*age.new + b6*baseef.new + b7*gender.new
  ## Specify priors
  b0 ~ dnorm(0,0.001)
  b1 ~ dnorm(0, 0.001)
  b2 ~ dnorm(0,0.001)
  b3 ~ dnorm(0,0.001)
  b4 ~ dnorm(0,0.001)
  b5 ~ dnorm(0,0.001)
  b6 ~ dnorm(0,0.001)
  b7 ~ dnorm(0,0.001)
  prec<-1/pow(sigma,2)
  #sigma~dunif(0,100)
  sigma~dgamma(2,0.001)
  

}
",file="modelheartBayesian.txt")

#creates an object representing a Bayesian graphical model, specified with a 
#BUGS-language description of the prior distribution, and a set of data.

jagsfit<-jags.model("modelheartBayesian.txt", data = jags.data, 
                    n.chains = 3, n.adapt = 1000)

#You usually need to throw away the initial samples ("burn-in"):

update(jagsfit, n.iter = 1000)

#function to extract random samples from the posterior distribution of a jags model

jm.sample <- jags.samples(jagsfit, variable.names = jags.params, n.iter = 1000, thin = 1)

plot(as.mcmc.list(jm.sample$b0), main = "Beta_0")
plot(as.mcmc.list(jm.sample$b1), main = "Beta_1")
plot(as.mcmc.list(jm.sample$b2), main = "Beta_2")

plot(as.mcmc.list(jm.sample$b3), main = "Beta_3")
plot(as.mcmc.list(jm.sample$b4), main = "Beta_4")
plot(as.mcmc.list(jm.sample$b5), main = "Beta_5")

plot(as.mcmc.list(jm.sample$b4), main = "Beta_6")
plot(as.mcmc.list(jm.sample$b5), main = "Beta_7")

plot(as.mcmc.list(jm.sample$y.new), main = "y.new")
