
require('INLA')
load("Navarre.RData")
data.navarra <- data.frame(ZBS=brainnav$ZBS,
                           Y=brainnav$OBSERVED,E=brainnav$EXPECTED)

navarra.graph <- inla.read.graph("Navarra.graph")


formula.zip <- Y~1 + f(ZBS, model="bym", graph=navarra.graph,
                         hyper=list(prec.unstruct=list(prior="gaussian",param=c(0,1)),
                                    prec.spatial=list(prior="gaussian",param=c(0,1))))


mod.zip1 <- inla(formula.zip,family="zeroinflatedpoisson1",
                   data=data.navarra, offset = log(E),
                   control.predictor=list(compute=TRUE))

round(mod.zip1$summary.hyperpar,3)

mod.zip0 <- inla(formula.zip,family="zeroinflatedpoisson0",
                 data=data.navarra, E = E,
                 control.predictor=list(compute=TRUE))

round(mod.zip0$summary.hyperpar,3)

Nareas <- nrow(data.navarra)
# Extract the random effects
zeta.navarra0 <- data.frame(zeta=unlist(lapply(mod.zip0$marginals.random$ZBS[1:Nareas],function(x)inla.emarginal(exp,x))))
zeta.navarra1 <- data.frame(zeta=unlist(lapply(mod.zip1$marginals.random$ZBS[1:Nareas],function(x)inla.emarginal(exp,x))))
# Create factor variables
RR.cutoff<- c(0.5, 0.8, 0.95,  1.05,  1.2, 2.5)
RR.navarra0 <- cut(zeta.navarra0$zeta,breaks=RR.cutoff,include.lowest=TRUE)
RR.navarra1 <- cut(zeta.navarra1$zeta,breaks=RR.cutoff,include.lowest=TRUE)

results <- data.frame(ZBS=data.navarra$ZBS,RR.navarra0, RR.navarra1)
data.navarra.shp <- attr(brainnav, "data")
attr(brainnav, "data") <- merge(data.navarra.shp, results, by="ZBS")

require('lattice')
trellis.par.set(axis.line=list(col=NA))
spplot(obj=brainnav, zcol="RR.navarra1", col.regions=gray(4.5:0.5/5),main="")
# ***