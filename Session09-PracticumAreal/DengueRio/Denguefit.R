library(spdep)
library(sf)
library(rstan)   
library(sp)
library(RColorBrewer)
library(classInt)
library(ggplot2)

setwd("Session09-PracticumAreal/DengueRio")
options(mc.cores = parallel::detectCores()-1)


#shapefile of Rio de Janeiro
shape_rj=st_read('DataShapeRJ.shp', 'DataShapeRJ') 

#Districts of Rio de Janeiro
#pay attention to the islands


plot(shape_rj)

#plotting the distribution of HDI across the city
colors <- brewer.pal(9, "YlOrRd") #set breaks for the 9 colors 
brks<-classIntervals(shape_rj$IDH, n=9, style="quantile")
brks<- brks$brks #plot the map
ggplot(shape_rj) + geom_sf(aes(fill=IDH)) + 
  ggtitle("Distribution of HDI") +
  scale_fill_binned(breaks = brks, labels = \(x) round(x,2))

#computing the SMR
SMR <- shape_rj$Ncases/shape_rj$Ecases
shape_rj$SMR=SMR

#plotting the distribution of SMR across the city
colors <- brewer.pal(9, "YlOrRd") #set breaks for the 9 colors 
brks<-classIntervals(SMR, n=9, style="quantile")
brks<- brks$brks #plot the map
ggplot(shape_rj) + geom_sf(aes(fill=SMR)) + 
  ggtitle("Distribution of SMR") +
  scale_fill_binned(breaks = brks, labels = \(x) round(x,2))


########################################
#
# Computing the Moran test statistic
#
########################################

# poly2nb builds a neighbours list based on regions with contiguous boundaries, 
#that is sharing one or more boundary point.
ccNb = poly2nb(st_make_valid(shape_rj),row.names=shape_rj$NOME)

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


#spdep function to create a weight matrix for a neighbours list with spatial weights
ccMat=nb2mat(ccNb,style='B', zero.policy=TRUE)  


# Connecting the disconnected districts into the main land
ccMat2 <- ccMat
# Cid Univ
which(rownames(ccMat)=="Cidade Universitária")
which(rownames(ccMat)=="Maré")
which(rownames(ccMat)=="Caju")
# Galeao
which(rownames(ccMat)=="Galeão")

ccMat2[42,which(rownames(ccMat)=="Maré")] <- 1 
ccMat2[42,which(rownames(ccMat)=="Caju")] <- 1 
ccMat2[42,which(rownames(ccMat)=="Galeão")] <- 1 

ccMat2[which(rownames(ccMat)=="Maré"),42] <- 1 
ccMat2[which(rownames(ccMat)=="Caju"),42] <- 1 
ccMat2[which(rownames(ccMat)=="Galeão"),42] <- 1 

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
                        control=list(adapt_delta = 0.97, stepsize = 0.1), chains=3, warmup=900, iter=1500)

#printing the posterior summary of some parameters
print(bym_dengue_stanfit, pars=c("beta0", "beta1", "sigma_phi", "tau_phi", "sigma_theta", 
                               "tau_theta", "rr[5]", "phi[5]", "theta[5]"), probs=c(0.025, 0.5, 0.975));

#storing the chains of all the parameters in the model
mcmc_chain<-as.matrix(bym_dengue_stanfit)


#summary of phi (spatial random effect)
sumtheta<-summary(bym_dengue_stanfit, pars=c("phi"))$summary
#summary of the relative risk 
sumrr<-summary(bym_dengue_stanfit, pars=c("rr"))$summary



#trace plots of the chains of some parameters
rstan::traceplot(bym_dengue_stanfit, pars = c("beta0", "beta1", "sigma_phi", "tau_phi", "sigma_theta", 
                            "tau_theta", "rr[5]", "phi[5]", "theta[5]"))

#plotting the estimated relative risk
shape_rj$spatial<-sumtheta[,1]
shape_rj$rr<-sumrr[,1]


colors <- brewer.pal(9, "YlOrRd") #set breaks for the 9 colors 
brks<-classIntervals(shape_rj$rr, n=9, style="quantile")
brks<- brks$brks #plot the map
ggplot(shape_rj) + geom_sf(aes(fill=rr)) + 
  ggtitle("Estimated relative risk") +
  scale_fill_binned(breaks = brks, labels = \(x) round(x,2))


colors <- brewer.pal(9, "YlOrRd") #set breaks for the 9 colors 
brks<-classIntervals(shape_rj$spatial, n=9, style="quantile")
brks<- brks$brks #plot the map
ggplot(shape_rj) + geom_sf(aes(fill=spatial)) + 
  ggtitle("Posterior mean spatial effect") +
  scale_fill_binned(breaks = brks, labels = \(x) round(x,2))

