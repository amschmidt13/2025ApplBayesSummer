#######################################################
#									                                    #
#		MAXIMUM LIKELIHOOD ESTIMATOR			                #
#			(NORMAL CASE)		                              	#
#######################################################

L=100
n1=10
n2=25
n3=50
n4=1000
x1=matrix(ncol=L,nrow=n1)
x2=matrix(ncol=L,nrow=n2)
x3=matrix(ncol=L,nrow=n3)
x4=matrix(ncol=L,nrow=n4)
xbar=matrix(ncol=1,nrow=L)
xvar=matrix(ncol=1,nrow=L)
xbar1=xbar
xbar2=xbar
xbar3=xbar
xbar4=xbar
xvar1=xvar
xvar2=xvar
xvar3=xvar
xvar4=xvar
for(i in 1:L){
	x1[,i]=rnorm(n1,0,1)
	xbar1[i]=mean(x1[,i])
	xvar1[i]=var(x1[,i])
	x2[,i]=rnorm(n2,0,1)
	xbar2[i]=mean(x2[,i])
	xvar2[i]=var(x2[,i])
	x3[,i]=rnorm(n3,0,1)
	xbar3[i]=mean(x3[,i])
	xvar3[i]=var(x3[,i])
	x4[,i]=rnorm(n4,0,1)
	xbar4[i]=mean(x4[,i])
	xvar4[i]=var(x4[,i])
}

#distribution of xbar under each sample size
par(mfrow=c(2,2))
hist(xbar1,prob=1,xlim=c(-2,2))
hist(xbar2,prob=1,xlim=c(-2,2))
hist(xbar3,prob=1,xlim=c(-2,2))
hist(xbar4,prob=1,xlim=c(-2,2))

# distribution of sample variance under each sample size
par(mfrow=c(2,2))
hist(xvar1,prob=1,xlim=c(0,2))
hist(xvar2,prob=1,xlim=c(0,2))
hist(xvar3,prob=1,xlim=c(0,2))
hist(xvar4,prob=1,xlim=c(0,2))


###############################################################################
#												                                                    	#
#			Example with a Normal sample 	                                					#
#		Prior mean=0 prior variance=1	                                						#
#												                                                    	#
###############################################################################
#grid with possible values of theta
theta<-seq(-2.5,2.5,0.01)
mu<-0
v2<-1
sigma2<-sqrt(25)
xbarra<- -0.95
n<-20
#prior density
prior1<-dnorm(theta,mu,v2)
#likelihood function
likelihood1<- dnorm(theta,xbarra,sigma2/n)
#parameters of the posterior according to conjugacy
mu11<- ((sigma2*mu)+(n*v2*xbarra))/(sigma2+n*v2)
v121<- sigma2*v2/(sigma2+n*v2)
#posterior distribution
posterior1<-dnorm(theta,mu11,v121)
par(mfrow=c(1,1))
plot(theta,likelihood1,type="l",lty=1,ylim=c(min(likelihood1,prior1,posterior1),max(likelihood1,prior1,posterior1)),ylab="Density",xlab="theta",bty="n")
lines(theta,prior1,type="l",lty=2)
lines(theta,posterior1,type="l",lty=3)
legend("topright",legend=c("likelihood","prior","posterior"),lty=c(1,2,3),bty="n")

#### increasing the sample size ########
#### n=100 ###
theta<-seq(-2.5,2.5,0.01)
#sampe prior as before
mu<-0
v2<-1
sigma2<-sqrt(25)
xbar<- -0.95
n<-100
prior<-dnorm(theta,mu,v2)
likelihood<- dnorm(theta,xbar,sigma2/n)
mu1<- ((sigma2*mu)+(n*v2*xbar))/(sigma2+n*v2)
v12<- sigma2*v2/(sigma2+n*v2)
posterior<-dnorm(theta,mu1,v12)

par(mfrow=c(1,1))
plot(theta,likelihood,type="l",lty=1,ylim=c(min(likelihood,prior,posterior),max(likelihood,prior,posterior)),ylab="Density",xlab="theta",bty="n")
lines(theta,prior,type="l",lty=2)
lines(theta,posterior,type="l",lty=3)
legend("topright",legend=c("likelihood","prior","posterior"),lty=c(1,2,3),bty="n")

#comparing the posteriors when n=20 and n=100

par(mfrow=c(2,1),pty="m")
plot(theta,likelihood1,type="l",lty=1,ylim=c(min(likelihood1,prior1,posterior1),max(likelihood1,prior1,posterior1)),ylab="Density",xlab="theta",lwd=1.5,bty="n")
lines(theta,prior1,type="l",lty=2,lwd=1.2)
lines(theta,posterior1,type="l",lty=3,lwd=1.2)
abline(v=xbarra,col=3,lwd=1.2,lty=2)
abline(v=mu11,lwd=1.2,col=2,lty=3)
title("n=20")
legend("topright",legend=c("likelihood","prior","posterior"),lty=c(1,2,3),bty="n")

plot(theta,likelihood,type="l",lty=1,ylim=c(min(likelihood,prior,posterior),max(likelihood,prior,posterior)),ylab="Density",xlab="theta",lwd=1.5,bty="n")
lines(theta,prior,type="l",lty=2,lwd=1.2)
lines(theta,posterior,type="l",lty=3,lwd=1.2)
abline(v=xbarra,col=3,lwd=1.5,lty=2)
abline(v=mu1,lwd=1.2,col=2,lty=3)
title("n=100")
legend("topright",legend=c("likelihood","prior","posterior"),lty=c(1,2,3),bty="n")


################################################################################

####### changing the value of v2, the prior variance #####
theta<-seq(-2.5,2.5,0.01)
mu<-0
#"big" variance
v2<-sqrt(225)
sigma2<-sqrt(25)
xbarra<- -0.95
n<-20
prior<-dnorm(theta,mu,v2)
likelihood<- dnorm(theta,xbarra,sigma2/n)
mu1<- ((sigma2*mu)+(n*v2*xbarra))/(sigma2+n*v2)
v12<- sigma2*v2/(sigma2+n*v2)
posterior<-dnorm(theta,mu1,v12)
### plotting all the functions together, prior, likelihood and posterior
par(mfrow=c(1,1))
plot(theta,likelihood,type="l",lty=1,ylim=c(min(likelihood,prior,posterior),max(likelihood,prior,posterior)),ylab="Density",xlab="theta",bty="n")
lines(theta,prior,type="l",lty=2)
lines(theta,posterior,type="l",lty=3)
abline(v=xbarra,col=3,lwd=1.5,lty=2)
abline(v=mu1,lwd=1.2,col=2,lty=3)
legend("topright",legend=c("likelihood","prior","posterior"),lty=c(1,2,3),bty="n")


###########################################################################################
#															                                                            #
#	Exercise with Bernoulli likelihood and a beta prior			                      	        #
# theta proportion of defective items							                                    		#
# number of defective items 3 and n=100								                                  	#
# prior: beta(2,200)												                                              #
#														                                                             	#
###########################################################################################

p=3
n=100
alpha=2
beta=200
theta=seq(0,1,by=0.001)
priori=dbeta(theta,alpha,beta)
veros=theta^p*(1-theta)^(n-p)
post=dbeta(theta,p+alpha,n-p+beta)
# focusing on the region with highest probability density
theta=seq(0,0.1,by=0.001)
priori1=dbeta(theta,alpha,beta)
#veros=theta^p*(1-theta)^(n-p)
veros=dbeta(theta,p+1,n-p+1)
post1=dbeta(theta,p+alpha,n-p+beta)

#changing the prior
alpha=2
beta=49
theta=seq(0,0.1,by=0.001)
priori2=dbeta(theta,alpha,beta)
#veros=theta^p*(1-theta)^(n-p)
veros=dbeta(theta,p+1,n-p+1)
post2=dbeta(theta,p+alpha,n-p+beta)

#comparing priors and posteriors
par(mfrow=c(1,2))
plot(theta,priori1,type="l",lty=1,ylim=c(min(priori1,veros,post1),max(priori1,veros,post1)),bty="n")
points(theta,veros,type="l",lty=2)
points(theta,post1,type="l",lty=3)
legend("topright",legend=c("likelihood","prior","posterior"),lty=c(1,2,3),bty="n")

plot(theta,priori2,type="l",lty=1,ylim=c(min(priori2,veros,post2),max(priori2,veros,post2)),bty="n")
points(theta,veros,type="l",lty=2)
points(theta,post2,type="l",lty=3)
legend("topright",legend=c("likelihood","prior","posterior"),lty=c(1,2,3),bty="n")

#############################################################
#								                                        		#
#   			Exponential - Gamma (Exercise 1)	              	#
# theta: time in minutes to assist a client	            		#
# Observed xbar=3.8 and n=20 				                    		#
# Prior gamma(0.04,0.2)					                        		#
#									                                        	#
#############################################################

theta=seq(0.001,1,by=0.001)
alpha=0.04
beta=0.2
n=20
xbarra=3.8
somax=xbarra*n
priori1=dgamma(theta,alpha,beta)
veros=dgamma(theta,n+1,somax)
post1=dgamma(theta,n+alpha,somax+beta)


par(mfrow=c(1,1),pty="s")
plot(theta,priori1,type="l",lty=1,lwd=1.2,bty="n",ylab="Density")
points(theta,veros,type="l",lty=2,lwd=1.2)
points(theta,post1,type="l",lty=3,lwd=1.2)
abline(v=1/xbarra,col=2,lty=2,lwd=1.4)
abline(v=(n+alpha)/(beta+n*xbarra),col=3,lty=3,lwd=1.4)
legend("topright",legend=c("prior","likelihood","posterior"),lty=c(1,2,3),bty="n")



alpha=2
beta=4
priori3=dgamma(theta,alpha,beta)
veros=dgamma(theta,n+1,somax)
post3=dgamma(theta,n+alpha,somax+beta)
plot(theta,priori3,type="l",lty=1,lwd=1.2,ylim=c(min(priori3,veros,post3),max(priori3,veros,post3)),bty="n")
points(theta,veros,type="l",lty=2,lwd=1.2)
points(theta,post3,type="l",lty=3,lwd=1.2)
legend("topright",legend=c("prior","likelihood","posterior"),lty=c(1,2,3),bty="n")


alpha=20
beta=20
priori4=dgamma(theta,alpha,beta)
veros=dgamma(theta,n+1,somax)
post4=dgamma(theta,n+alpha,somax+beta)
plot(theta,priori4,type="l",lty=1,lwd=1.2,ylim=c(min(priori4,veros,post4),max(priori4,veros,post4)))
points(theta,veros,type="l",lty=2,lwd=1.2)
points(theta,post4,type="l",lty=3,lwd=1.2)
legend("topright",legend=c("prior","likelihood","posterior"),lty=c(1,2,3))




par(mfrow=c(2,2),pty="s")
plot(theta,priori1,type="l",lty=1,lwd=1.2,bty="n")
points(theta,veros,type="l",lty=2,lwd=1.2)
points(theta,post1,type="l",lty=3,lwd=1.2)
legend("topright",legend=c("prior","likelihood","posterior"),lty=c(1,2,3),bty="n")


plot(theta,priori3,type="l",lty=1,lwd=1.2,ylim=c(min(priori3,veros,post3),max(priori3,veros,post3)),bty="n")
points(theta,veros,type="l",lty=2,lwd=1.2)
points(theta,post3,type="l",lty=3,lwd=1.2)
legend("topright",legend=c("prior","likelihood","posterior"),lty=c(1,2,3),bty="n")

plot(theta,priori4,type="l",lty=1,lwd=1.2,ylim=c(min(priori4,veros,post4),max(priori4,veros,post4)),bty="n")
points(theta,veros,type="l",lty=2,lwd=1.2)
points(theta,post4,type="l",lty=3,lwd=1.2)
legend("topright",legend=c("prior","likelihood","posterior"),lty=c(1,2,3),bty="n")

#############################################################
#										                                        #
#   			Exponential - Gamma 		                          #
# Predictive Distribution						                        #
# 										                                      #
#############################################################

# Assume that, conditional on beta, the lifetime of fluorescent 
# lamps are independent exponential random variables with parameter 
# beta. Suppose that beta has a prior distribution that is a gamma
# distribution with parameters 5 and 20000. After we observe
# five lamps with lifetimes 3150, 3450, 3230, 3509 and 3110 (in hours),
# we want to predict the lifetime X_6 of the next lamp.
# We need to find the distribution of f(x_6 | x_1,x_2,x_3,x_4,x_5)
#

x=c(3150,3450,3230,3509,3110)
n=5
alphapri=5
gammapri=20000
sumx=sum(x)
gammapost=gammapri+sumx
alphapos=alphapri+n
x6=seq(0.01,10000,by=1)
predf=((gammapost^alphapos)*alphapos)/((gammapost+x6)^(alphapos+1))
par(mfrow=c(1,1))
plot(x6,predf,type="l")
abline(v=3100,col=2,lwd=1.3)
