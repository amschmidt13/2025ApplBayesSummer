library(nimble)
library(coda)



setwd('Session04-Practicum')

bangladesh <- read.csv("bangladesh.csv", sep=";")

sort(unique(bangladesh$district))
bangladesh$district_id <- as.integer(as.factor(bangladesh$district))

sort(unique(bangladesh$district_id))
dist<-as.factor(bangladesh$district)

table(bangladesh$district,bangladesh$use.contraception) |> plot()
x<-model.matrix(use.contraception~as.factor(district_id),data=bangladesh)

m0<-glm(use.contraception~-1+dist,family="binomial",data=bangladesh)
summ0<-summary(m0)
round(summ0$coefficients,4)

MLEfunction<-function(ni,sumy,grid){
  
  loglik<-c()
  prob<-c()
  for(i in 1:length(grid)){
    prob[i]<-exp(grid[i])/(1+exp(grid[i]))
    loglik[i]<-sumy*log(prob[i])+(ni-sumy)*log(1-prob[i])
  }
  loglik
  
  
}

bangladesh$use.contraception[bangladesh$district_id==3]


#district 3
grid<-seq(-3,20,length=1000)
MLEbeta<-MLEfunction(2,2,grid)
plot(grid,MLEbeta,type="l",bty="n")

#district 11
length(bangladesh$use.contraception[bangladesh$district_id==11])
grid<-seq(-20,3,length=1000)
MLEbeta<-MLEfunction(21,0,grid)
plot(grid,MLEbeta,type="l",bty="n")



#district 49
bangladesh$use.contraception[bangladesh$district_id==49]
grid<-seq(-20,3,length=1000)
MLEbeta<-MLEfunction(4,0,grid)
plot(grid,MLEbeta,type="l",bty="n")


####################################
#
# Nimble fixed effects (M1)
#
####################################
codeM1 = nimbleCode({
  for (i in 1:N){
    contracept[i]~dbern(p[i])
    logit(p[i])<-inprod(beta[], x[i,])
  }
  for(j in 1:60){
    beta[j]~dnorm(0,0.1)
    prob[j]<-expit(beta[j])
  }
})


#
# Initial estimates     
#          
initials<-list(beta=rnorm(60,0,0.1))

#
# Data
#
datcontracept<-list(contracept=bangladesh$use.contraception,
                    x=x)
constants = list(N=1934)

params<-c("beta","prob")


modelM1 <- nimbleModel(codeM1, constants = constants, data = datcontracept, inits = initials)
modelM1$initializeInfo()
cModelM1 <- compileNimble(modelM1)
confM1 <- configureMCMC(modelM1, monitors = params, enableWAIC = T)
MCMC.M1 <- buildMCMC(confM1)
cMCMC.M1 <- compileNimble(MCMC.M1, project = cModelM1)
samplesM1 <- runMCMC(cMCMC.M1, niter = 10000, nburnin = 5000, thin = 5,
                     nchains = 2, 
                   WAIC = TRUE, samplesAsCodaMCMC = TRUE)

samps = samplesM1$samples
sumM1<-summary(samps)


index<-seq(1,60,1)

ymin<-min(sumM1$quantiles[61:120,1])
ymax<-max(sumM1$quantiles[61:120,5])
plot(index,sumM1$quantiles[61:120,3],bty="n",ylim=c(ymin,ymax),pch=19,cex=0.7,xlab="District")
for(i in 1:60)
  segments(index[i],sumM1$quantiles[60+i,1],index[i],sumM1$quantiles[60+i,5])  

#posteiror full conditional of a fixed effect
fullConditionalbeta<-function(ni,m0,c0,sumy,grid){
  
  fullcond<-c()
  prob<-c()
  for(i in 1:length(grid)){
    prob[i]<-exp(grid[i])/(1+exp(grid[i]))
    fullcond[i]<-sumy*log(prob[i])+(ni-sumy)*log(1-prob[i])+dnorm(grid[i],m0,c0,log=TRUE)
  }
  fullcond
  
}

bangladesh$use.contraception[bangladesh$district_id==3]
grid<-seq(-5,15,length=1000)
fullbeta<-fullConditionalbeta(2,0,10,2,grid)
plot(grid,fullbeta,type="l",bty="n")


####################################
#
# Nimble random effects (M2)
#
####################################
codeM2 = nimbleCode({
  for (i in 1:N){
    logit(p[i])<-beta0+b[district_id[i]]
    contracept[i]~dbern(p[i])
  }
  beta0~dnorm(0,0.1)
  for(j in 1:60){
    b[j]~dnorm(0,tau)
    pdist[j]<-expit(beta0+b[j])
  }
  tau~dgamma(0.1,0.1)
  sigma<-1/pow(tau,2)
})


initials<-list(beta0=0.5,b=rnorm(60,0,1),tau=0.1)

#
# Data
#
datcontracept<-list(contracept=bangladesh$use.contraception)
constants = list(N=1934, district_id=bangladesh$district_id)

params<-c("beta0","b","pdist","tau")


modelM2 <- nimbleModel(codeM2, constants = constants, data = datcontracept, inits = initials)
modelM2$initializeInfo()
cModelM2 <- compileNimble(modelM2)
confM2 <- configureMCMC(modelM2, monitors = params, enableWAIC = T)
MCMC.M2 <- buildMCMC(confM2)
cMCMC.M2 <- compileNimble(MCMC.M2, project = cModelM2)
samplesM2 <- runMCMC(cMCMC.M2, niter = 20000, nburnin = 5000, thin = 10,
                     nchains = 2, 
                     WAIC = TRUE, samplesAsCodaMCMC = TRUE)

samps = samplesM2$samples
sumM2<-summary(samps)

index<-seq(1,60,1)

ymin<-min(sumM2$quantiles[62:121,1])
ymax<-max(sumM2$quantiles[62:121,5])
plot(index,sumM2$quantiles[62:121,3],bty="n",ylim=c(ymin,ymax),pch=19,cex=0.7,xlab="District")
for(i in 1:60)
 segments(index[i],sumM2$quantiles[61+i,1],index[i],sumM2$quantiles[61+i,5])  



# just to recap...
# Frequentist estimates
ci0 = confint.default(m0)
ymin<-min(ci0[,1])
ymax<-max(ci0[,2])
plot(index,summ0$coefficients[,1],bty="n",ylim=c(ymin,ymax),pch=19,cex=0.7,xlab="District")
for(i in 1:60)
  segments(index[i],ci0[i,1],index[i],ci0[i,2])  

plot(index,summ0$coefficients[,1],bty="n",pch=19,cex=0.7,xlab="District")


# model M1
ymin<-min(sumM1$quantiles[61:120,1])
ymax<-max(sumM1$quantiles[61:120,5])
plot(index,sumM1$quantiles[61:120,3],bty="n",ylim=c(ymin,ymax),pch=19,cex=0.7,xlab="District")
for(i in 1:60)
  segments(index[i],sumM1$quantiles[60+i,1],index[i],sumM1$quantiles[60+i,5])  


# model M2
ymin<-min(sumM2$quantiles[62:121,1])
ymax<-max(sumM2$quantiles[62:121,5])
plot(index,sumM2$quantiles[62:121,3],bty="n",ylim=c(ymin,ymax),pch=19,cex=0.7,xlab="District")
for(i in 1:60)
  segments(index[i],sumM2$quantiles[61+i,1],index[i],sumM2$quantiles[61+i,5])  

