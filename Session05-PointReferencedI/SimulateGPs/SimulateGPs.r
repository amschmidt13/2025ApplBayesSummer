#############################################################
#
#
# Analysis of spatially structured data
# Summer School in Health Statistics
# Alexandra M. Schmidt - 01/06/2018
#
#############################################################

library(geoR)


#
# Variogram models with the same "practical" range:
#
v.f <- function(x, ...){1-cov.spatial(x, ...)}
#
curve(v.f(x, cov.pars=c(1, .2)), from = 0, to = 1,
      xlab = "distance", ylab = expression(gamma(h)),
      main = "variograms with equivalent \"practical range\"")
curve(v.f(x, cov.pars = c(1, .6), cov.model = "sph"), 0, 1,
      add = TRUE, lty = 2)
curve(v.f(x, cov.pars = c(1, .6/sqrt(3)), cov.model = "gau"),
      0, 1, add = TRUE, lwd = 2)
legend(0.5,.3, c("exponential", "spherical", "gaussian"),
       lty=c(1,2,1), lwd=c(1,1,2))


#
# Matern models with equivalent "practical range"
# and varying smoothness parameter
#
par(mfrow=c(1,1))
curve(v.f(x, cov.pars = c(1, 0.25), kappa = 0.5),from = 0, to = 1,
      xlab = "distance", ylab = expression(gamma(h)), lty = 2,
      main = "models with equivalent \"practical\" range")
curve(v.f(x, cov.pars = c(1, 0.188), kappa = 1),from = 0, to = 1,
      add = TRUE)
curve(v.f(x, cov.pars = c(1, 0.14), kappa = 2),from = 0, to = 1,
      add = TRUE, lwd=2, lty=2)
curve(v.f(x, cov.pars = c(1, 0.117), kappa = 2),from = 0, to = 1,
      add = TRUE, lwd=2)
legend(0.4,.4, c(expression(paste(kappa == 0.5, " and ",
                                  phi == 0.250)),
                 expression(paste(kappa == 1, " and ", phi == 0.188)),
                 expression(paste(kappa == 2, " and ", phi == 0.140)),
                 expression(paste(kappa == 3, " and ", phi == 0.117))),
       lty=c(2,1,2,1), lwd=c(1,1,2,2))

####################################################################################################
#
# simulating different GPs with different values of phi, sigma2 in the correlation function 
#
###############################################################################################3



## Different values of phi
##
for(i in 0:9){
  #  jpeg(paste("phi",i, ".jpg",  sep=""), wid=600, hei=600)
  par(mfrow=c(2,2), mar=c(1.5,.5,1.5,0), mgp=c(1, .5, 0))
  set.seed(234+i)
  ap1 <- grf(961, grid="reg", cov.pars=c(1, 0))
  set.seed(234+i)
  ap2 <- grf(961, grid="reg", cov.pars=c(1, .1))
  set.seed(234+i)
  ap3 <- grf(961, grid="reg", cov.pars=c(1, .25))
  set.seed(234+i)
  ap4 <- grf(961, grid="reg", cov.pars=c(1, .75))
  iis <- range(c(ap1$data, ap2$data, ap3$data, ap4$data))
  image(ap1, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(phi==0), cex=1.5)
  image(ap2, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(phi==0.10), cex=1.5)
  image(ap3, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(phi==0.25), cex=1.5)
  image(ap4, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(phi==0.75), cex=1.5)
  #dev.off()
}

##
## Different values of sigma^2
##
for(i in 0:9){
#  jpeg(paste("sigma",i, ".jpg",  sep=""), wid=600, hei=600)
  par(mfrow=c(2,2), mar=c(1.5,.5,1.5,0), mgp=c(1, .5, 0))
  set.seed(234+i)
  ap1 <- grf(961, grid="reg", cov.pars=c(1, 0.3))
  set.seed(234+i)
  ap2 <- grf(961, grid="reg", cov.pars=c(2, .3))
  set.seed(234+i)
  ap3 <- grf(961, grid="reg", cov.pars=c(3, .3))
  set.seed(234+i)
  ap4 <- grf(961, grid="reg", cov.pars=c(5, .3))
  iis <- range(c(ap1$data, ap2$data, ap3$data, ap4$data))
  image(ap1, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(sigma^2==1), cex=1.5)
  image(ap2, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(sigma^2==2), cex=1.5)
  image(ap3, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(sigma^2==3), cex=1.5)
  image(ap4, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(sigma^2==5), cex=1.5)
 # dev.off()
}

##
## Different nugget effects
##
for(i in 0:9){
  # jpeg(paste("tau",i, ".jpg",  sep=""), wid=600, hei=600)
  par(mfrow=c(2,2), mar=c(1.5,.5,1.5,0), mgp=c(1, .5, 0))
  set.seed(234+i)
  ap1 <- grf(961, grid="reg", cov.pars=c(1, 0.3), nug=0)
  set.seed(234+i)
  ap2 <- grf(961, grid="reg", cov.pars=c(.75, .3), nug=0.25)
  set.seed(234+i)
  ap3 <- grf(961, grid="reg", cov.pars=c(.5, .3), nug=0.5)
  set.seed(234+i)
  ap4 <- grf(961, grid="reg", cov.pars=c(.1, .3), nug=.9)
  iis <- range(c(ap1$data, ap2$data, ap3$data, ap4$data))
  image(ap1, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(paste(sigma^2==1, " e ", tau^2 == 0), cex=1.5))
  image(ap2, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(paste(sigma^2==0.75, " e ", tau^2 == 0.25), cex=1.5))
  image(ap3, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(paste(sigma^2==0.5, " e ", tau^2 == 0.5), cex=1.5))
  image(ap4, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(paste(sigma^2==0.1, " e ", tau^2 == 0.9), cex=1.5))
  # dev.off()
}

##
## Different correlation functions
##
for(i in 0:9){
  # jpeg(paste("model",i, ".jpg",  sep=""), wid=600, hei=600)
  par(mfrow=c(2,2), mar=c(1.5,.5,1.5,0), mgp=c(1, .5, 0))
  set.seed(234+i)
  ap1 <- grf(961, grid="reg", cov.pars=c(1, .25))
  set.seed(234+i)
  ap2 <- grf(961, grid="reg", cov.pars=c(1, .75), cov.model="sph")
  set.seed(234+i)
  ap3 <- grf(961, grid="reg", cov.pars=c(1, .14), kappa=2)
  set.seed(234+i)
  ap4 <- grf(961, grid="reg", cov.pars=c(1, .095), kappa=5)
  iis <- range(c(ap1$data, ap2$data, ap3$data, ap4$data))
  image(ap1, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext("exponencial", cex=1.5)
  image(ap2, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext("esferico", cex=1.5)
  image(ap3, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(paste("Matern com ", kappa==2)), cex=1.5)
  image(ap4, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(paste("Matern com ", kappa==5)), cex=1.5)
  # dev.off()
}

##
## Different anisotropy
##
for(i in 0:9){
  #jpeg(paste("aniso",i, ".jpg",  sep=""), wid=600, hei=600)
  par(mfrow=c(2,2), mar=c(1.5,.5,1.5,0), mgp=c(1, .5, 0))
  set.seed(234+i)
  ap1 <- grf(961, grid="reg", cov.pars=c(1, .25), aniso.pars=c(pi/4, 2))
  set.seed(234+i)
  ap2 <- grf(961, grid="reg", cov.pars=c(1, .25), aniso.pars=c(pi/4, 4))
  set.seed(234+i)
  ap3 <- grf(961, grid="reg", cov.pars=c(1, .25), aniso.pars=c(2*pi/3, 2))
  set.seed(234+i)
  ap4 <- grf(961, grid="reg", cov.pars=c(1, .25), aniso.pars=c(2*pi/3, 4))
  iis <- range(c(ap1$data, ap2$data, ap3$data, ap4$data))
  image(ap1, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(paste(phi[a]==pi/4, " \ ,\ ", phi[r] == 2), cex=1.5))
  image(ap2, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(paste(phi[a]==pi/4, " \ ,\ ", phi[r] == 4), cex=1.5))
  image(ap3, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(paste(phi[a]== 2*pi/3, " \ ,\ ", phi[r] == 2), cex=1.5))
  image(ap4, xlab="", ylab="", col=gray(seq(1,0,l=21)), zlim=iis)
  mtext(expression(paste(phi[a]==2*pi/3, " \ ,\ ", phi[r] == 4), cex=1.5))
  #dev.off()
}
