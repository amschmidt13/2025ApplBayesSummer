library(sf)
f <- file.path("https://www.paulamoraga.com/book-spatial/",
               "data/PM25USA2022.csv")
d <- read.csv(f)
d <- st_as_sf(d, coords = c("longitude", "latitude"))
st_crs(d) <- "EPSG:4326"

library(rnaturalearth)
map <- ne_countries(type = "countries",
                    country = "United States of America",
                    scale = "medium", returnclass = "sf")
map <- st_crop(map, xmin = -130, xmax = -60, ymin = 18, ymax = 72)

d <- st_filter(d, map)
nrow(d)


library(ggplot2)
library(viridis)
ggplot() + geom_sf(data = map) +
  geom_sf(data = d, aes(col = value)) +
  scale_color_viridis()

library(sf)
library(terra)

# raster grid covering map
grid <- terra::rast(map, nrows = 100, ncols = 100)
# coordinates of all cells
xy <- terra::xyFromCell(grid, 1:ncell(grid))

dp <- st_as_sf(as.data.frame(xy), coords = c("x", "y"),
               crs = st_crs(map))

# indices points within the map
indicespointswithin <- which(st_intersects(dp, map,
                                           sparse = FALSE))

# points within the map
dp <- st_filter(dp, map)

# plot
ggplot() + geom_sf(data = map) +
  geom_sf(data = dp)

#covariates
library(geodata)
covtemp <- worldclim_global(var = "tavg", res = 10,
                            path = tempdir())
covprec <- worldclim_global(var = "prec", res = 10,
                            path = tempdir())

#compute the averages over months and extract the values at the observation and prediction 
#locations with the extract() function of terra.
d$covtemp <- extract(mean(covtemp), st_coordinates(d))[, 1]
d$covprec <- extract(mean(covprec), st_coordinates(d))[, 1]
# Extract at prediction locations
dp$covtemp <- extract(mean(covtemp), st_coordinates(dp))[, 1]
dp$covprec <- extract(mean(covprec), st_coordinates(dp))[, 1]

#temperature and precipitation
library("patchwork")
p1 <- ggplot() + geom_sf(data = map) +
  geom_sf(data = d, aes(col = covtemp)) +
  scale_color_viridis()
p2 <- ggplot() + geom_sf(data = map) +
  geom_sf(data = d, aes(col = covprec)) +
  scale_color_viridis()
p1 / p2

#The data we are dealing with have a geographic CRS that references locations using longitude and latitude values. 
#In order to work with kilometers instead of degrees, we use st_transform() to transform the CRS of 
#the sf objects with the data corresponding to the observed (d) and the prediction (dp) locations from geographic 
#to a projected CRS.

st_crs("EPSG:3857")$proj4string
projMercator<-"+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0
+x_0=0 +y_0=0 +k=1 +units=km +nadgrids=@null +wktext +no_defs"
d <- st_transform(d, crs = projMercator)
dp <- st_transform(dp, crs = projMercator)


# Observed coordinates
coo <- st_coordinates(d)

# Predicted coordinates
coop <- st_coordinates(dp)

#fitting a geostatistical model in INLA using the SPDE approach
library(INLA)

summary(dist(coo)) # summary of distances between locations

mesh <- inla.mesh.2d(loc = coo, max.edge = c(200, 500),
                     cutoff = 1)
mesh$n

plot(mesh)
points(coo, col = "red")
axis(1)
axis(2)

#Use the inla.spde2.matern() function to build the SPDE model. This function has parameters mesh with the 
#triangulated mesh constructed and constr = TRUE to impose an integrate-to-zero constraint. Moreover,  set 
#the smoothness parameter ν equal to 1. In the spatial case d=2 and α=ν+d/2=2

spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)


#This function generates a list with the vector s ranging from 1 to spde$n.spde. Additionally, it creates 
#two vectors, s.group and s.repl, containing all elements set to 1 and lengths equal 
#to the number of mesh vertices.

indexs <- inla.spde.make.index("s", spde$n.spde)
lengths(indexs)

#We use the inla.spde.make.A() function of R-INLA passing the mesh (mesh) and the coordinates (coo) 
#to easily construct a projection matrix A that projects the spatially continuous Gaussian random field from the observations to the mesh nodes.

A <- inla.spde.make.A(mesh = mesh, loc = coo)
dim(A)
nrow(coo)
#We also create a projection matrix for the prediction locations.
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

# stack for estimation stk.e
stk.e <- inla.stack(tag = "est",
                    data = list(y = d$value), A = list(1, A),
                    effects = list(data.frame(b0 = rep(1, nrow(A)),
                                              covtemp = d$covtemp, covprec = d$covprec),
                                   s = indexs))

# stack for prediction stk.p
stk.p <- inla.stack(tag = "pred",
                    data = list(y = NA), A = list(1, Ap),
                    effects = list(data.frame(b0 = rep(1, nrow(Ap)),
                                              covtemp = dp$covtemp, covprec = dp$covprec),
                                   s = indexs))

# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

formula <- y ~ 0 + b0 + covtemp + covprec + f(s, model = spde)

res <- inla(formula, family = "gaussian",
            data = inla.stack.data(stk.full),
            control.predictor = list(compute = TRUE,
                                     A = inla.stack.A(stk.full)),
            control.compute = list(return.marginals.predictor = TRUE),verbose = TRUE)
            
res$summary.fixed

#Mapping predicted PM2.5 values
index <- inla.stack.index(stack = stk.full, tag = "pred")$data
pred_mean <- res$summary.fitted.values[index, "mean"]
pred_ll <- res$summary.fitted.values[index, "0.025quant"]
pred_ul <- res$summary.fitted.values[index, "0.975quant"]

grid$mean <- NA
grid$ll <- NA
grid$ul <- NA

grid$mean[indicespointswithin] <- pred_mean
grid$ll[indicespointswithin] <- pred_ll
grid$ul[indicespointswithin] <- pred_ul

summary(grid) # negative values for the lower limit

library(rasterVis)
levelplot(grid, layout = c(1, 3),
          names.attr = c("Mean", "2.5 percentile", "97.5 percentile"))
