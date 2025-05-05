######################################
#                                    #
#           Lab 9 codes              #
#                                    #
######################################
library(maps)
library(spdep)
library(fields)
library(R2WinBUGS)

# read the csv file.
dat<-read.csv("Obama2012.csv",header=TRUE)

## before we do all the analysis, check that the county labels are the same 
#as in the plot.nc function
county_map <- map("county","north carolina",plot=FALSE)
county_label <- rep(0,length(county_map$names[-(28:29)]))                       # 28 and 29 are replicates with 27
for (i in 1:length(county_map$names[-(28:29)])){
county_label[i] <- toupper(strsplit(county_map$names[-(28:29)],",")[[i]][2])
}
county_label[27] <- "CURRITUCK"

## now county_label is the county order we use in plot.nc function
## county in our file
county <- toupper(as.character(dat[,1]))
sum(county==county_label)

## so, it's not 100, means there are counties' oreders not the same with plot.nc
## adjustment
new_dat <- dat
for (ii in 1:length(county_label)){
new_dat[ii,] <- dat[county==county_label[ii],]
}                    # end of ii loop

dat <- new_dat

# response
Y <- dat[,2]

# 1. choose covariates
# a. look at the correlation plot
# b. fit linear model to all combination and pick the one by AIC.
# c. Reasons to pick some of the covariates
# d. Others

# These are I randomly picked....so don't pick exactly the same as mine in the HW.....
pctBlack <- dat[,8]
pctAme <- dat[,9]
pctAsi <- dat[,10]
pctHis <- dat[,11]
pctunemployed <- dat[,18]

# 2. plot response and residuals

plot.nc <- function(data,main="",breaks=5,sig_digits=0){
   d<-data[c(1:26,27,27,27:100)] #county 27 has three polygons!
   if(length(breaks)==1){                                                       # means if we give the number of breaks only
     breaks<-quantile(d,seq(0,1,length=breaks+1),na.rm=T)
   }
   nbreaks<-length(breaks)-1                                                    # number of intervals
   low<-breaks[-(nbreaks+1)]
   high<-breaks[-1]
   col <- rep(1,length(d))
   for(j in 2:nbreaks){
     col<-ifelse(d>low[j],j,col)
   }
   shade <- 1-col/nbreaks                                                       # shade needs to between 0 to 1
   shades<- 1-1:nbreaks/nbreaks                                                 # colors in the labels

   low<-round(low,sig_digits)
   high<-round(high,sig_digits)
   map("county","north carolina",fill=T,col=gray(shade))
   title(main)
   legend("bottomleft",paste(low,"-",high),lwd=5,ncol=2,col=gray(shades))
}

lm1 <- lm(Y~pctBlack + pctAme + pctAsi + pctHis + pctunemployed)
resid1 <- lm1$resid

X11()
plot.nc(Y,main="Response",sig_digits=3)
X11()
plot.nc(resid1,main="Residual",sig_digits=3)


# 3. Moran test
## create 2 weight objects

## first: adjacency matrix
NCcentroids <- read.csv("NCcentroids.csv", header=TRUE)
s <- NCcentroids[,4:3]                                                          
# longitude then latitude
NCADJ<-read.csv("NCADJ.csv", header=TRUE, row.names=1)
NCADJ<-as.matrix(NCADJ)
######   Create and plot the nb class object based on "adjacency"   #######:
nb1 <- mat2listw(NCADJ, row.names = NULL, style="M")                            
# we do not have row name in the matrix, M means matrix style
summary(nb1)
X11()
map("county","north carolina")
plot.listw(nb1,coords=s,add=TRUE,col=4,lwd=2)                                   
# add=TRUE, add to the existing plot!!

## second:

################################################################################
################################################################################

# Now, look at the moran's I test for response and first type of weight object
test1 <- moran.test(Y,listw=nb1,randomisation=TRUE,alternative="two.sided")                             
# randomisation here represents two different ways in calculate the sd of moran I.
                                                                                
# then the p-value is the Z-score result
                                                                                
# moran.mc
# test2
# test3
# test4


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# CAR model

############################################################
###########           Set up for BUGS           ############
############################################################

n<-ncol(NCADJ)
num<-colSums(NCADJ)
adj<-NULL
for(j in 1:n){
  adj<-c(adj,which(NCADJ[j,]==1))
}
adj<-as.vector(adj)
num<-as.vector(num)
weights<-1+0*adj

# num[i] = number of regions that adjacent to region i
# The first num[1] elements of adj are region 1's neighbors,
# the next num[2] elements of adj are region 2's neighbors, etc.
# weights being all one specify that all neighbors are weighted equally

############################################################
#########     Specify the intrinsic CAR model    ###########
############################################################

sink("D:/ST 733/NCSU Site/docs/Lab 9/icar_model.txt")                           # way to output the following model into a text file.
cat("

model{
   for(i in 1:n){  #Likelihood
       Y[i] ~  dnorm(mu[i],tau[1])
        mu[i] <- b[1]*pctBlack[i] + b[2]*pctAme[i] + b[3]*pctAsi[i] + b[4]*pctHis[i] + b[5]*pctunemployed[i] + S[i]
         S[i] <- int + R[i]  #Add in an intercept to get the unrestricted CAR.
    }

   #CAR model (in BUGS, the R are forced to sum to zero)
   R[1:n] ~ car.normal(adj[], weights[], num[], tau[2])

   #Priors
   int~dflat()    # Must be flat to get the ICAR
   for(j in 1:5){
   b[j] ~ dnorm(0,0.1)
   }
   for(j in 1:2){
      tau[j] ~ dgamma(0.1,0.1)
      sigma[j]<-pow(tau[j],-0.5)
   }
}

}", fill = TRUE)
sink()

############################################################
#########              Call BUGS                 ###########
############################################################


data<-list(pctBlack=pctBlack, pctAme=pctAme, Y=Y, pctAsi=pctAsi, pctHis=pctHis, pctunemployed=pctunemployed, n=n, num=num, adj=adj, weights=weights)
inits<-function(){
    list(int=0,b=rep(0,5),tau=rep(10,2))
}
keepers <- c("b","sigma","S","mu","int")

out <- bugs(
	data=data,
	inits=inits,
	parameters.to.save=keepers,
	model.file="D:/ST 733/NCSU Site/docs/Lab 9/icar_model.txt",
	n.iter=10000,
	n.chains=3,
	n.burnin=1000,
	n.thin=1,
	debug=TRUE,
#Change this to the location of the openBUGS executable on your computer
        bugs.directory = "c:/Program Files (x86)/WinBUGS14/",
        #OpenBUGS.pgm="C:\\Program Files (x86)\\OpenBUGS\\OpenBUGS321\\OpenBUGS.exe"  # using OpenBUGS need this line and also use library(R2OpenBUGS)
                                                                                 
	)

############################################################
#########       Summarize the results            ###########
############################################################

print(round(out$summary[1:7,c(1,2,3,7)],3))

X11()
plot.nc(out$mean$mu,main="Posterior mean mu",sig_digits=3)

X11()
plot.nc(out$mean$S,main="Posterior mean S",sig_digits=3)

int_result <- round(out$summary[208,c(1,2,3,7)],3)

R <- out$mean$S-int_result[1]

X11()
plot.nc(R,main="Posterior mean R",sig_digits=3)

out$DIC

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

############## CAR model ###################
model1 <- spautolm(Y~pctBlack + pctAme + pctAsi + pctHis + pctunemployed,listw=nb1,family="CAR")

#?spautolm
summary(model1)
# look at the LR test value
LL<-model1$LL
LL0<-model1$LL0
LR <- -2*(LL0-LL)
pvalue <- 1-pchisq(LR, 1)

X11()
plot.nc(model1$fit$signal_trend,main="Non-spatial part",sig_digits=3)
X11()
plot.nc(model1$fit$signal_stochastic,main="Spatial part",sig_digits=3)

model1$lambda

# compare with the OLS results
summary(lm1)


