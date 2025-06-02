library(sf)
library(ggplot2)
library(tidyr)
library(tidyverse)

pts1 <- tidyrpts1 <- st_point(c(10,40)); pts2 <- st_point(c(40,30))
pts3 <- st_point(c(20,20)); pts4 <- st_point(c(30,10))
st_geometry(pts1)
#points

pts <- st_multipoint(cbind(c(10,40,20,30),c(40,30,20,10)))
st_geometry(pts)

ggplot(pts) + geom_sf(col = "blue", size = 2) + theme_bw()

#linestrings

coords1 <- cbind(c(10,20,10), c(10,20,40))
ls1 <- st_linestring(coords1)

coords2 <- cbind(c(40,30,40,30), c(40,30,20,10))
ls2 <- st_linestring(coords2)
ls <- st_multilinestring(list(coords1, coords2))
st_geometry(ls)

st_geometry(ls)
ggplot(ls) + geom_sf(col = "blue") + theme_bw()


#polygons
cor1 <- list(cbind(c(30,45,10,30),c(20,40,40,20)))
pol1 <- st_polygon(cor1)
cor2 <- list(cbind(c(15,40,10,5,15),c(5,10,20,10,5)))
pol2 <- st_polygon(cor2)
st_geometry(pol2)

#sf objects
pol <- st_multipolygon(list(cor1, cor2))
st_geometry(pol)

ggplot(pol) +
  geom_sf(aes(geometry = geometry), col = "blue",
          linewidth = 1, fill = "gray90") + theme_bw()


cor1 <- list(cbind(c(40,20,45,40),c(40,45,30,40)))
cor21 <- list(cbind(c(20,10,10,30,45,20),c(35,30,10,5,20,35)))
cor22 <- list(cbind(c(30,20,20,30),c(20,15,25,20)))
pol <- st_multipolygon(list(cor1, c(cor21,cor22)))
st_geometry(pol)

ggplot(pol) +
  geom_sf(aes(geometry = geometry), col = "blue",
          linewidth = 1, fill = "gray90") + theme_bw()

#Create simple feature object
Toronto <- st_point(c(-79.3819, 43.6525))
Montreal <- st_point(c(-73.5833, 45.50))
Ottawa <- st_point(c(-75.7003, 45.4201))

sfc <- st_sfc(Toronto, Montreal, Ottawa)
sfc

dt <- data.frame(Name = c("Toronto", "Montreal","Ottawa"),
                 geom = sfc)
dt


cities <- st_sf(dt)
cities

cities <- cities %>% mutate(ID = 1:3)
cities

ggplot() + geom_sf(data=cities, shape=1, size=5, col=2) +
  theme_bw()

#retrieve coordinate system
st_crs(cities)
st_crs(4326)

st_crs(2958)

cities_sf <- st_set_crs(cities, value = st_crs(4326))
ggplot(cities_sf) + geom_sf() + theme_bw()

#representing these points on a map
library(leaflet)
leaflet(cities_sf) %>% addTiles() %>% addMarkers()

#transforming the coordinate system to UTM
cities_sf <- st_set_crs(cities, value = st_crs(4326))
cities_utm <- st_transform(cities_sf, crs = st_crs(2958))
cities_utm

#compute relative distances
sqrt(sum((st_coordinates(cities_utm)[1,] -
            st_coordinates(cities_utm)[2,])^2))

st_distance(cities_utm, which = "Euclidean")

#Return bounding of a simple feature or simple feature set
st_bbox(cities_sf)

st_bbox(cities_utm)

st_coordinates(cities_sf)

st_coordinates(cities_utm)

#create sf point object from data.frame

stuff <- data.frame(
  lon = c(-79.3819, -73.5833, -75.7003),
  lat = c(43.6525, 45.5, 45.4201))
stuff_sf <- st_as_sf(x = stuff, coords = c("lon", "lat"),
                     crs = st_crs(4326))
stuff_sf

#some operations with sf object
stuff_sf <- stuff_sf %>%
  mutate(val1 = runif(3), val2 = runif(3)) %>%
  mutate(ss = val1 + val2)
stuff_sf

pol1 <- st_polygon(list(cbind(c(30,20,10,30), c(10,40,20,10))))
pol2 <- st_polygon(list(cbind(c(30,40,20,30), c(10,40,40,10))))
pol_sfc <- st_sfc(pol1,pol2)
pol_sf <- st_sf(data.frame(Name = c("A","B"), geom = pol_sfc))
ggplot(pol_sf) + geom_sf(aes(fill = Name)) + theme_bw()

pol_sf <- pol_sf %>% mutate(Value = runif(2))
ggplot(pol_sf) + geom_sf(aes(fill = Value)) + theme_bw()

pol_centroid <- st_centroid(pol_sf)
ggplot(pol_sf) + geom_sf(fill = "gray90") +
  geom_sf(data = pol_centroid, col = "blue") + theme_bw()

pol_A_sf <- pol_sf %>% filter(Name == "A")
ggplot(pol_A_sf) + geom_sf(fill = "gray90") + theme_bw()

pol_union_sf <- st_union(pol_sf)
ggplot(pol_union_sf) + geom_sf(fill = "gray90") + theme_bw()

pol_buffer_sf <- st_buffer(pol_union_sf, dist = 1)
ggplot() + geom_sf(data = pol_union_sf, fill = "gray90") +
  geom_sf(data = pol_buffer_sf, fill = NA, col = "blue",
          linewidth = 1) + theme_bw()

pol_grid_sf <- st_make_grid(pol_union_sf, n = c(10,10),
                            what = "polygons")
corner_grid_sf <- st_make_grid(pol_union_sf, n = c(10,10),
                               what = "corners")
center_grid_sf <- st_make_grid(pol_union_sf, n = c(10,10),
                               what = "centers")

ggplot() + geom_sf(data = pol_union_sf, fill = "gray90") +
  geom_sf(data = pol_grid_sf, fill = NA, col = "green") +
  geom_sf(data = corner_grid_sf, col = "red") +
  geom_sf(data = center_grid_sf, col = "blue") + theme_bw()


center_in_pol_sf <- st_intersection(pol_union_sf, center_grid_sf)
corner_in_pol_sf <- st_intersection(pol_union_sf, corner_grid_sf)


ggplot() + geom_sf(data = pol_union_sf, fill = "gray90") +
  geom_sf(data = center_in_pol_sf, col = "red") +
  geom_sf(data = corner_in_pol_sf, col = "blue") + theme_bw()


#Maps in R
library(rnaturalearth)
world <- ne_countries(scale = "medium", type = "countries",
                      returnclass = "sf")
class(world)

#retrieve coordinates of object world
st_crs(world)

ggplot(world) + geom_sf(fill = NA) +
  coord_sf() + theme_bw()

canada <- world %>% filter(admin == "Canada")
class(canada)

ggplot(canada) + geom_sf(fill = NA) + coord_sf() + theme_bw()

#Italy
italy <- world %>% filter(admin == "Italy")
class(italy)

ggplot(italy) + geom_sf(fill = NA) + coord_sf() + theme_bw()

#maps

library(leaflet)
leaflet(canada) %>% addTiles() %>%
  addPolygons(fillColor = "red", weight = 1, color = "black")

#obtaining the boundaries for Ontario
on <- read_sf("./Province/Province.shp")
st_crs(on)

# obtaining coordinates with different projections
on_sf <- st_transform(on, crs = st_crs(4326))
on_utm <- st_transform(on, crs = st_crs(6836))

library(leaflet)
leaflet(on_sf) %>% addTiles() %>%
  addPolygons(fillColor = "red", weight = 1, color = "blue")

#Raster data
library(raster)
r <- raster(ncol = 10, nrow = 10, xmn = 0, xmx = 1,
            ymn = 0, ymx = 1)
values(r) <- runif(ncell(r))
plot(r)

#Maps of raster data
library(terra)
fpath <- system.file("ex/elev.tif", package = "terra")
spr <- rast(fpath)
spr

spr_df <- as.data.frame(spr, xy = TRUE)
str(spr_df)

#Maps of raster data
terra_elev_sf<-st_as_sf(spr_df,coords =c("x","y"), crs=st_crs(4326))
ggplot(terra_elev_sf) + geom_sf(aes(col = elevation)) +
                            scale_color_viridis_c() + xlab("")+ 
                       ylab("") + theme_bw()


library(leaflet)
rb <- brick(spr)
pal <- colorNumeric("YlOrRd", values(spr),
                    na.color = "transparent")
leaflet() %>% addTiles() %>%
  addRasterImage(spr, colors = pal, opacity = 0.8) %>%
  addLegend(pal = pal,
            values = values(spr),
            title = "elevation")

#elevation data download
library(elevatr)
on_elev <- get_elev_raster(on_sf, z = 1)
names(on_elev) <- "elev_meter"
plot(on_elev)


on_elev_sf<-st_as_sf(as.data.frame(on_elev,xy =TRUE),
                     coords=c("x","y"),crs=st_crs(4326))
library(colorspace)                       
ggplot(on_elev_sf)+geom_sf(aes(col =elev_meter))+
scale_color_viridis_c() +
geom_sf(data =on_sf,fill =NA,col ="blue") +
theme_bw()+xlab("")+ylab("")+theme(legend.title =element_blank())


# elevation data on grid
on_grid <- st_make_grid(on_sf, n = c(30,30), crs = st_crs(4326),
                        what = "centers")
on_grid_sf <- st_as_sf(on_grid)
all.equal(st_crs(on_sf), st_crs(on_grid_sf))


within_on_grid <- st_intersection(on_sf, on_grid_sf)
on_elev_grid <- get_elev_point(locations = within_on_grid, prj =4326)

ggplot(on_elev_grid) + geom_sf(aes(col = elevation)) +
  geom_sf(data = on_sf, fill = NA, col = "gray80", linewidth = 0.7)+
          scale_color_viridis_c() + theme_bw()
