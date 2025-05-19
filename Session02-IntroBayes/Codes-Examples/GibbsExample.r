#################################################
#								                                #
# Codes Written by Hedibert F. Lopes        		#
#	INSPER - Brazil                               #
#								                                #
#################################################



# Simple Example : sampling from N(0,1)
# try different values for n

n=1000
x=rnorm(n,0,1)
xx=seq(-3,3,by=0.01)
hist(x,prob=1)
lines(xx,dnorm(xx,0,1))
mean(x)
var(x)





#############################################
#
# Sampling from the Multivariate Normal Distribution
#
#############################################

rmvnorm = function(am,mu,S) {
  s = chol(S)
  am.star = dim(S)[[1]]*am
  aux = matrix(rnorm(am.star),dim(S)[[1]],am)
  aux2 = mu + t(s)%*%aux
  t(aux2)
}




###################################################
# Exemplo 1:
# x=(x1,x2) ~ N(m,S)  where m=(m1,m2)' and
# 
#           |  s1 s12 |
#        S =|         |
#           | s21  s2 |
#
#  Sabe-se tambem que:
#
#    x1|x2 ~ N(m1+(s12/s2)*(x2-m2);s1-s12*s21/s2)
#    x2|x1 ~ N(m2+(s21/s1)*(x1-m1);s2-s21*s12/s1)
#
###################################################
m<-c(2,1)
S<-matrix(c(1,0.7,0.7,1),2,2)

# True distribution
x1<-seq(m[1]-3*sqrt(S[1,1]),m[1]+3*sqrt(S[1,1]),by=0.25)
x2<-seq(m[2]-3*sqrt(S[2,2]),m[2]+3*sqrt(S[2,2]),by=0.25)
fnorm<-matrix(0,length(x1),length(x2))
iS<-solve(S)
const<-1/sqrt(2*pi*abs(S[1,1]*S[2,2]-S[1,2]^2))
for (i in 1:length(x1))
	for (j in 1:length(x2)){
		fnorm[i,j]<-const*exp(-0.5*(t(c(x1[i],x2[j])-m)%*%iS%*%(c(x1[i],x2[j])-m)))
   }

par(mfrow=c(1,1))
contour(x1,x2,fnorm,drawlabels=F)
persp(x1,x2,fnorm, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
           ltheta = 120, shade = 0.75, ticktype = "detailed")

# Gibbs Sampler
M<-1000
x<-c(-1,3)
xs<-x
for (i in 1:M){
 x[1]<-rnorm(1,m[1]+S[1,2]/S[2,2]*(x[2]-m[2]),sqrt(S[1,1]-S[1,2]^2/S[2,2]))
 x[2]<-rnorm(1,m[2]+S[2,1]/S[1,1]*(x[1]-m[1]),sqrt(S[2,2]-S[2,1]^2/S[1,2]))
 xs<-rbind(xs,x)
}

# Gibbs Sampler paths
xss<-matrix(0,2*M,2)
xss[1,]<-xs[1,1:2]
xss[2,2]<-xs[1,2]
for (i in 2:M){
  xss[2*(i-1),1]<-xs[i,1]
  xss[2*i-1  ,1]<-xs[i,1]
  xss[2*i-1,2]<-xs[i,2]
  xss[2*i  ,2]<-xs[i,2]
}

par(mfrow=c(1,1))
contour(x1,x2,fnorm,xlim=c(min(xs[,1]),max(xs[,1])),ylim=c(min(xs[,2]),max(xs[,2])),drawlabels=F)
lines(xss[1:10,],col=2,lwd=3)
title("10 iterations")

contour(x1,x2,fnorm,xlim=c(min(xs[,1]),max(xs[,1])),ylim=c(min(xs[,2]),max(xs[,2])),drawlabels=F)
lines(xss[1:20,],col=2,lwd=3)
title("20 iterations")

contour(x1,x2,fnorm,xlim=c(min(xs[,1]),max(xs[,1])),ylim=c(min(xs[,2]),max(xs[,2])),drawlabels=F)
lines(xss[1:100,],col=2,lwd=3)
title("100 iterations")

contour(x1,x2,fnorm,xlim=c(min(xs[,1]),max(xs[,1])),ylim=c(min(xs[,2]),max(xs[,2])),drawlabels=F)
lines(xss,col=2,lwd=3)
title("1000 iterations")


plot(xss,xlim=c(min(xs[,1]),max(xs[,1])),ylim=c(min(xs[,2]),max(xs[,2])))
contour(x1,x2,fnorm,drawlabels=F,add=T)
title("1000 iterations")


# Marginal posterior densities
par(mfrow=c(2,1))
hist(xs[101:M,1],prob=T,col=0,xlab="x1",main="")
lines(x1,dnorm(x1,m[1],sqrt(S[1,1])))
abline(v=m[1],lwd=5)
hist(xs[101:M,2],prob=T,col=0,xlab="x2",main="")
lines(x2,dnorm(x2,m[2],sqrt(S[2,2])))
abline(v=m[2],lwd=5)

# Ergodic means
par(mfrow=c(2,1))
plot(cumsum(xs[,1])/1:(M+1),type="l",xlab="itera??o",ylab="x1")
abline(h=m[1])
plot(cumsum(xs[,2])/1:(M+1),type="l",xlab="itera??o",ylab="x2")
abline(h=m[2])

#############################################################################
#
#          EXAMPLE: Poisson Process with a Change Point:  Originally 
#          analysed by Carlin, Gelfand and Smith (1992) and used in 
#          Tanner's Book - Tools for Statistical Inference (page 147)
#	    (Yearly (1851-1962)) Coalmining Disaster Data	
#
#
#############################################################################
rgama<-function(a,b){
  return(rgamma(1,a)/b)
}
rinvgama<-function(a,b){
  return(b/rgamma(1,a))
}
dgama<-function(x,a,b){
	b*dgamma(b*x,a)
}

y<-c(4,5,4,1,0,4,3,4,0,6,3,3,4,0,2,6,3,3,5,4,5,3,1,4,4,1,5,5,3,
    4,2,5,2,2,3,4,2,1,3,2,2,1,1,1,1,3,0,0,1,0,1,1,0,0,3,1,0,3,
    2,2,0,1,1,1,0,1,0,1,0,0,0,2,1,0,0,0,1,1,0,2,3,3,1,1,2,1,1,
    1,1,2,4,2,0,0,0,1,4,0,0,0,1,0,0,0,0,0,1,0,0,1,0,1)
n<-length(y)

# Likelihood function
lik<-function(k,lambda,phi,y,n){
	t1<-ifelse(k>1,sum(y[1:k]),0)
	t2<-ifelse(k<n,sum(y[(k+1):n]),0)
   return(lambda^t1*exp(-k*lambda)*phi^t2*exp(-(n-k)*phi))
}

# Time-series plot
par(mfrow=c(2,2))
plot(y,ylab="Counts",xlab="Time")
lines(y)
abline(v=41)
count<-rep(0,7);for (i in 0:6){count[i+1]<-sum(y==i)}
plot(0:6,count,xlab="counts",ylab="density",type="h")
title("112 observations")
count<-rep(0,7);for (i in 0:6){count[i+1]<-sum(y[1:41]==i)}
plot(0:6,count,xlab="counts",ylab="density",type="h")
title("First 41 observations")
count<-rep(0,7);for (i in 0:6){count[i+1]<-sum(y[42:n]==i)}
plot(0:6,count,xlab="counts",ylab="density",type="h")
title("Last 71 observations")


# Prior's hyperparameters
# -------------------------
alpha<-1;beta<-1
gamma<-1;delta<-1
par(mfrow=c(2,1))
x<-seq(0.001,10,0.1)
plot(x,dgama(x,alpha,beta),type="l",ylab="",xlab="gamma")
plot(x,dgama(x,gamma,delta),type="l",ylab="",xlab="phi")

# Marginal Posterior of k
ks<-seq(1,n)
postk<-NULL
for (k in ks){
  t1<-ifelse(k>1,sum(y[1:k]),0)
  t2<-ifelse(k<n,sum(y[(k+1):n]),0)
  postk<-c(postk,lgamma(t2+gamma)+lgamma(t1+alpha)-
        (t2+gamma)*log(n-k+delta)-(t1+alpha)*log(k+beta))
}
par(mfrow=c(1,1))
plot(ks,postk,type="h",xlab="k",ylab="log-posterior")


# Joint Posterior for lambda and phi
N<-20
lambdas<-seq(2,4,length=N)
phis<-seq(0.5,1.5,length=N)
post<-matrix(0,N,N)
for (i in 1:N){
  print(i)
  for (j in 1:N)
    for (l in 1:n) # l here represents the possible values of k
      post[i,j]<-post[i,j]+lik(l,lambdas[i],phis[j],y,n)*
                (lambdas[i]^(alpha-1)*exp(-beta*lambdas[i]))*
                (phis[j]^(gamma-1)*exp(-delta*phis[j]))
}
post1<-post/sum(post)
par(mfrow=c(2,1))
contour(lambdas,phis,100*post1,xlab="lambda",ylab="phi")
points(mean(y[1:41]),mean(y[42:n]),type="p",lwd=10,cex=2)
persp(lambdas,phis,100*post1,xlab="lambda",ylab="phi",zlab="densidade")


par(mfrow=c(2,1))
plot(lambdas,apply(100*post1,1,sum)/max(apply(100*post1,1,sum)),
     xlab="lambda",ylab="posteriori",type="l",axes=F,xlim=c(0,4))
axis(1);axis(2)
abline(v=mean(y[1:41]))
plot(phis,apply(100*post1,2,sum)/max(apply(100*post1,2,sum))*3.3,
     xlab="phi",ylab="posteriori",type="l",axes=F,xlim=c(0,4))
axis(1);axis(2)
abline(v=mean(y[42:n]))

###########################################
#
# Gibbs Sampler
#
###########################################
# initial value for k
# -------------------
k<-41

# Gibbs Sampler
# -------------
M<-1000
draws<-NULL
# Gibbs Sampler
for (i in 1:M){
	print(i)
	t1<-ifelse(k>1,sum(y[1:k]),0)
	t2<-ifelse(k<n,sum(y[(k+1):n]),0)
	lambda<-rgama(t1+alpha,k+beta)
	phi<-rgama(t2+gamma,n-k+delta)
        liks<-matrix(ncol=1,nrow=n)
        for (j in 1:n)
	        liks[j]<-lik(j,lambda,phi,y,n)
        k<-sample(1:n,size=1,prob=liks)
        draws<-rbind(draws,c(lambda,phi,k))
}
nome<-c("lambda","phi","k")
par(mfrow=c(3,3))
for (i in 1:3){
	plot(draws[,i],xlab="iteration",ylab="",main=nome[i],type="l")
	acf(draws[,i])
	hist(draws[,i],xlab=nome[i],prob=T,col=0)
}

par(mfrow=c(2,1))
hist(draws[,1],xlim=c(0,4),prob=T,col=0,xlab="lambda",nclass=20)
lines(lambdas,apply(100*post1,1,sum)/max(apply(100*post1,1,sum))*1.4)
axis(1);axis(2)
abline(v=mean(y[1:41]))
hist(draws[,2],xlim=c(0,4),prob=T,col=0,xlab="phi")
lines(phis,apply(100*post1,2,sum)/max(apply(100*post1,2,sum))*3.3)
axis(1);axis(2)
abline(v=mean(y[42:n]))

# Marginal Posterior of k
postk1<-NULL
for (k in 1:n) postk1<-c(postk1,sum(draws[,3]==k))
par(mfrow=c(1,1))
plot(ks,postk1/M,type="h",xlab="k",ylab="posteriori")


##################################################################################
#
# Mixture of Normals using Metropolis-Hastings Algorithm
#
##################################################################################

dnormm<-function(x,m,P){
	exp(-0.5*(t(x-m)%*%P%*%(x-m)))
}
dmixtnormm<-function(x,m1,m2,P1,P2,dP1,dP2){
	0.7*dP1*exp(-0.5*(t(x-m1)%*%P1%*%(x-m1)))+
	0.3*dP2*exp(-0.5*(t(x-m2)%*%P2%*%(x-m2)))
}
m1<-c(4,5)
m2<-c(1,4)
S1<-matrix(c(1,0.7,0.7,1),2,2)
iS1<-solve(S1)
iS1<-(iS1+t(iS1))/2
S2<-matrix(c(1,-0.7,-0.7,1),2,2)
iS2<-solve(S2)
iS2<-(iS2+t(iS2))/2
#determinant of covariance matrix
dP1<-prod(diag(chol(iS1)))
dP2<-prod(diag(chol(iS2)))

#"True" Distribution
x1<-seq(-1,7,by=0.25)
x2<-seq(2,8,by=0.25)
fnorm<-matrix(0,length(x1),length(x2))
for(i in 1:length(x1))
for(j in 1:length(x2))
fnorm[i,j]<-dmixtnormm(c(x1[i],x2[j]),m1,m2,iS1,iS2,dP1,dP2)

par(mfrow=c(2,1))
contour(x1,x2,fnorm)
persp(x1,x2,fnorm)

#Metropolis-Hastings using independent proposals V=diag(1,1)
V<-diag(1,2)
#Metropolis-Hastings using proposal with variance equals  S2
V<-S2

cV<-t(chol(V))
iV<-solve(V)
M<-1000
x<-c(1,4)
xs<-x
accepted<-0
for(i in 1:M){
	if (i %%1==0) print(i)
	xnew<-x+cV%*%rnorm(2)
	acceptance<-min(1,(dmixtnormm(xnew,m1,m2,iS1,iS2,dP1,dP2)/dmixtnormm(x,m1,m2,iS1,iS2,dP1,dP2))*(dnormm(x,xnew,iV)/dnormm(xnew,x,iV)))
	u<-runif(1)
	if (u<acceptance){
		x<-xnew
		accepted<-accepted+1
	}
	xs<-c(xs,x)
}
xs<-matrix(xs,(M+1),2,byrow=T)

xs1<-xs    #V=S2

xs2<-xs    #V=I_2

#Metropolis Paths
xss<-matrix(0,2*M,2)
xss[1,]<-xs[1,1:2]
xss[2,2]<-xs[1,2]
for(i in 2:M){
	xss[2*(i-1),1]<-xs[i,1]
	xss[2*i-1,1]<-xs[i,1]
	xss[2*i-1,2]<-xs[i,2]
	xss[2*i,2]<-xs[i,2]
}

minx<-min(c(xs[,1],x1,xss[,1]))
maxx<-max(c(xs[,1],x1,xss[,1]))
miny<-min(c(xs[,2],x2,xss[,2]))
maxy<-max(c(xs[,2],x2,xss[,2]))

# Markov chain 
par(mfrow=c(2,2))
for(i in c(20,100,300,2*M)){
 contour(x1,x2,fnorm,xlim=c(minx,maxx),ylim=c(miny,maxy),main=paste(i/2, " iterations"))
 lines(xss[1:i,])
}

#comparing the trace plots of both chains
par(mfrow=c(2,2))
plot(xs1[,1],type="l")
plot(xs1[,2],type="l")
plot(xs2[,1],type="l")
plot(xs2[,2],type="l")

# autocorrelations amongst the sampled values
par(mfrow=c(2,2))
acf(xs1[,1])
acf(xs1[,2])
acf(xs2[,1])
acf(xs2[,2])

# Metropolis-Hastings (single moves)
M<-1000
x<-c(1,4)
xs<-x
for(i in 1:M){
	print(i)
	m1new<-x[1]+rnorm(1)
	acceptance<-min(1,(dmixtnormm(c(m1new,x[2]),m1,m2,iS1,iS2,dP1,dP2)/dmixtnormm(x,m1,m2,iS1,iS2,dP1,dP2))*(dnormm(x,c(m1new,x[2]),iV)/dnormm(c(m1new,x[2]),x,iV)))
	u<-runif(1)
	if (u<acceptance){
		x[1]<-m1new
	}
	m2new<-x[2]+rnorm(1)
	acceptance<-min(1,(dmixtnormm(c(x[1],m2new),m1,m2,iS1,iS2,dP1,dP2)/dmixtnormm(x,m1,m2,iS1,iS2,dP1,dP2))*(dnormm(x,c(x[1],m2new),iV)/dnormm(c(x[1],m2new),x,iV)))
	u<-runif(1)
	if (u<acceptance){
		x[2]<-m2new
	}
	xs<-c(xs,x)
}
xs<-matrix(xs,(M+1),2,byrow=T)
xs3<-xs

xs<-xs3
#Metropolis Paths
xss<-matrix(0,2*M,2)
xss[1,]<-xs[1,1:2]
xss[2,2]<-xs[1,2]
for(i in 2:M){
	xss[2*(i-1),1]<-xs[i,1]
	xss[2*i-1,1]<-xs[i,1]
	xss[2*i-1,2]<-xs[i,2]
	xss[2*i,2]<-xs[i,2]
}

minx<-min(c(xs[,1],x1,xss[,1]))
maxx<-max(c(xs[,1],x1,xss[,1]))
miny<-min(c(xs[,2],x2,xss[,2]))
maxy<-max(c(xs[,2],x2,xss[,2]))

# Markov chain
par(mfrow=c(2,2))
for(i in c(20,100,300,2*M)){
 contour(x1,x2,fnorm,xlim=c(minx,maxx),ylim=c(miny,maxy),main=paste(i/2, " iteracoes"))
 lines(xss[1:i,])
}

#comparing the trace plots of all chains
par(mfrow=c(3,2))
plot(xs1[,1],type="l")
plot(xs1[,2],type="l")
plot(xs2[,1],type="l")
plot(xs2[,2],type="l")
plot(xs3[,1],type="l")
plot(xs3[,2],type="l")

# autocorrelations amongst the sampled values
par(mfrow=c(3,2))
acf(xs1[,1])
acf(xs1[,2])
acf(xs2[,1])
acf(xs2[,2])
acf(xs3[,1])
acf(xs3[,2])


