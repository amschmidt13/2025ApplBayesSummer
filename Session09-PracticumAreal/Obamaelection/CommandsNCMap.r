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


plot.nc(Y,main="Response",sig_digits=3)

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
map("county","north carolina")
plot.listw(nb1,coords=s,add=TRUE,col=4,lwd=2)  

## second: based on distances
coords <- matrix(c(s$LONGITUDE,s$LATITUDE),byrow=FALSE,ncol=2)
k1 <- knn2nb(knearneigh(coords,longlat=TRUE))
all.linked <- max(unlist(nbdists(k1, coords)))
col.nb.0.all <- dnearneigh(coords, 0, all.linked)
summary(col.nb.0.all, coords)
map("county","north carolina")
plot(col.nb.0.all, coords, add=TRUE)
title(main=paste("Distance based neighbours 0-",  format(all.linked),
                 " distance units", sep=""))
