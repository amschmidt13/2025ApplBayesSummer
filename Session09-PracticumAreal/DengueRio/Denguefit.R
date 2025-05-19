library(spdep)
library(maptools)
library(rstan)   
library(sp)
library(RColorBrewer)
library(classInt)

options(mc.cores = parallel::detectCores())  

#shapefile of Rio de Janeiro
shape_rj=readShapePoly('DataShapeRJ') 

#Districts of Rio de Janeiro
#pay attention to the islands


plot(shape_rj)

#plotting the distribution of HDI across the city
colors <- brewer.pal(9, "YlOrRd") #set breaks for the 9 colors 
brks<-classIntervals(shape_rj$IDH, n=9, style="quantile")
brks<- brks$brks #plot the map
plot(shape_rj, col=colors[findInterval(shape_rj$IDH, brks,all.inside=TRUE)], axes=FALSE,xlim=c(-44,-42.8),main="Distribution of HDI")
legend(x=-43.1, y=-22.9, legend=leglabs(round(brks,2)), fill=colors, bty="n",x.intersp = .5, y.intersp = .5)

#computing the SMR
SMR <- shape_rj$Ncases/shape_rj$Ecases

#plotting the distribution of SMR across the city
colors <- brewer.pal(9, "YlOrRd") #set breaks for the 9 colors 
brks<-classIntervals(SMR, n=9, style="quantile")
brks<- brks$brks #plot the map
plot(shape_rj, col=colors[findInterval(SMR, brks,all.inside=TRUE)], axes=FALSE,xlim=c(-44,-42.8),main="Distribution of SMR")
legend(x=-43.1, y=-22.9, legend=leglabs(round(brks,2)), fill=colors, bty="n",x.intersp = .5, y.intersp = .5)

########################################
#
# Computing the Moran test statistic
#
########################################

# poly2nb builds a neighbours list based on regions with contiguous boundaries, 
#that is sharing one or more boundary point.
ccNb = poly2nb(shape_rj,row.names=shape_rj$NOME)

# converting the neighbor structure into a neighbor list 
W<-nb2listw(ccNb,zero.policy=TRUE,style="W")

#testing for positive spatial correlation, areas with no 
#neighbors are set to zero
moran.test(shape_rj$Ncases, W,zero.policy=TRUE)

#A plot of spatial data against its spatially lagged values, augmented by 
#reporting the summary of influence measures for the linear relationship 
#between the data and the lag. If zero policy is TRUE, 
#such observations are also marked if they occur.

moran.plot(shape_rj$Ncases,W,zero.policy=TRUE)

#getting the neighborhood structure to fit the BYM model
#poly2nb is within package spdep
ccNb = poly2nb(shape_rj,row.names=shape_rj$NOME)
#spdep function to create a weight matrix for a neighbours list with spatial weights
ccMat=nb2mat(ccNb,style='B', zero.policy=TRUE)  


# Connecting the disconnected districts into the main land
ccMat2 <- ccMat
# Cid Univ
which(rownames(ccMat)=="Cidade Universitária")
which(rownames(ccMat)=="Maré")
which(rownames(ccMat)=="Caju")
which(rownames(ccMat)=="Galeão")

ccMat2[42,which(rownames(ccMat)=="Maré")] <- 1 
ccMat2[42,which(rownames(ccMat)=="Caju")] <- 1 
ccMat2[42,which(rownames(ccMat)=="Galeão")] <- 1 

ccMat2[which(rownames(ccMat)=="Maré"),42] <- 1 
ccMat2[which(rownames(ccMat)=="Caju"),42] <- 1 
ccMat2[which(rownames(ccMat)=="Galeão"),42] <- 1 

# Galeao
which(rownames(ccMat)=="Penha Circular")
which(rownames(ccMat)=="Penha")
which(rownames(ccMat)=="Cordovil")

ccMat2[3,which(rownames(ccMat)=="Penha Circular")] <- 1 
ccMat2[3,which(rownames(ccMat)=="Penha")] <- 1 
ccMat2[3,which(rownames(ccMat)=="Cordovil")] <- 1 
ccMat2[3,which(rownames(ccMat)=="Maré")] <- 1 

ccMat2[which(rownames(ccMat)=="Penha Circular"),3] <- 1 
ccMat2[which(rownames(ccMat)=="Penha"),3] <- 1 
ccMat2[which(rownames(ccMat)=="Cordovil"),3] <- 1 
ccMat2[which(rownames(ccMat)=="Maré"),3] <- 1 



#Building the information about Neighborhood structure 
#as required by WinBugs/OpenBugs

Neigh<-c()
sumNeigh<-matrix(ncol=1,nrow=160)
for(i in 1:ncol(ccMat2)){
  Neigh<-c(Neigh,which(ccMat2[i,]==1))
  sumNeigh[i]<-sum(ccMat2[i,])
}
sumNeigh

#transforming the neighborhood structure to be read into the stan code
source("mungeCARdata4stan.R")  
nbs = mungeCARdata4stan(Neigh, sumNeigh);
N = nbs$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;


#assigning the data to the variables that go into the stan code

#number of cases of dengue fever
y <- shape_rj$Ncases
#standardized HDI
x <- (shape_rj$IDH - mean(shape_rj$IDH))/sd(shape_rj$IDH)
#Expected number of cases based on the population of each district
E <- shape_rj$Ecases

##########################################
#
#Poisson regression - a classical approach
#
#########################################

summary(glm(y~x,offset=log(E),family="poisson"))

#Accounting for  overdispersion
summary(glm(y~x,offset=log(E),family="quasipoisson"))

#fitting a BYM model using Stan
bym_dengue_stanfit = stan("bym_predictor_plus_offset.stan", 
                          data=list(N=N,N_edges=N_edges,node1=node1,node2=node2,y=y,x=x,E=E), 
                        control=list(adapt_delta = 0.97, stepsize = 0.1), chains=3, warmup=9000, iter=10000, 
                        save_warmup=FALSE);

#printing the posterior summary of some parameters
print(bym_dengue_stanfit, pars=c("beta0", "beta1", "sigma_phi", "tau_phi", "sigma_theta", 
                               "tau_theta", "rr[5]", "phi[5]", "theta[5]"), probs=c(0.025, 0.5, 0.975));

#storing the chains of all the parameters in the model
mcmc_chain<-as.matrix(bym_dengue_stanfit)

plot(mcmc_chain[,"beta0"],type="l")
hist(mcmc_chain[,"beta0"],prob=1)


plot(mcmc_chain[,"beta1"],type="l")
hist(mcmc_chain[,"beta1"],prob=1)

#summary of phi (spatial random effect)
sumtheta<-summary(bym_dengue_stanfit, pars=c("phi"))$summary
#summary of the relative risk 
sumrr<-summary(bym_dengue_stanfit, pars=c("rr"))$summary



#trace plots of the chains of some parameters
traceplot(bym_dengue_stanfit, pars = c("beta0", "beta1", "sigma_phi", "tau_phi", "sigma_theta", 
                            "tau_theta", "rr[5]", "phi[5]", "theta[5]"))

#plotting the estimated relative risk
shape_rj$spatial<-sumtheta[,1]
shape_rj$rr<-sumrr[,1]


colors <- brewer.pal(9, "YlOrRd") #set breaks for the 9 colors 
brks<-classIntervals(shape_rj$rr, n=9, style="quantile")
brks<- brks$brks #plot the map
plot(shape_rj, col=colors[findInterval(shape_rj$rr, brks,all.inside=TRUE)], axes=FALSE,xlim=c(-44,-42.8),main="Estimated relative risk")
legend(x=-43.1, y=-22.9, legend=leglabs(round(brks,2)), fill=colors, bty="n",x.intersp = .5, y.intersp = .5)


colors <- brewer.pal(9, "YlOrRd") #set breaks for the 9 colors 
brks<-classIntervals(shape_rj$spatial, n=9, style="quantile")
brks<- brks$brks #plot the map
plot(shape_rj, col=colors[findInterval(shape_rj$spatial, brks,all.inside=TRUE)], axes=FALSE,xlim=c(-44,-42.8),main="Posterior mean spatial effect")
legend(x=-43.1, y=-22.9, legend=leglabs(round(brks,2)), fill=colors, bty="n",x.intersp = .5, y.intersp = .5)


