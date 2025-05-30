library(tidyverse)
library(sf)
library(rnaturalearth)
library(abind)
library(Rcpp)
library(RcppEigen)
library(tictoc)
library(coda)
library(ggpubr)

setwd('Session13-Practicum')

##### AGRIMONIA DATASET #######
# https://zenodo.org/records/7956006

load('Agrimonia_Dataset_v_3_0_0.Rdata')
rawdata = AgrImOnIA_Dataset_v_3_0_0
range(rawdata$Time)


stations = rawdata[!duplicated(rawdata$IDStations), c(1:3, 5)] %>% 
  st_as_sf(coords = c('Longitude', 'Latitude'), crs=4326)
metadata = read_csv('Metadata_monitoring_network_registry_v_2_0_1.csv')
metadata
metadata = metadata %>% filter(Pollutant=='NO2')
stations_NO2 = stations %>% filter(IDStations %in% metadata$IDStation)

# import study region and plot stations
lombardy = ne_states(country = 'Italy') %>% 
  filter(region == 'Lombardia')
ggplot() + geom_sf(data = lombardy, fill='white') + 
  geom_sf(data = stations) +
  geom_sf(data = stations_NO2, col = 'red')


# exploratory analysis
data_filtered = rawdata %>% 
  filter(IDStations %in% metadata$IDStation,
         Time >= '2021-01-01') # restricting the time just for illustration
prop.NaN = data_filtered %>% 
  select(IDStations, AQ_no2) %>% 
  group_by(IDStations) %>% 
  summarise(prop.NaN = mean(is.nan(AQ_no2)), .groups = 'drop')
prop.NaN %>% 
  ggplot(aes(IDStations, prop.NaN)) + geom_point()

# drop stations with >15% missing data
prop.NaN = prop.NaN %>% filter(prop.NaN <= 0.15)
stations_NO2 = stations_NO2 %>% filter(IDStations %in% prop.NaN$IDStations)
ggplot() + geom_sf(data = lombardy, fill='white') + 
  geom_sf(data = stations) +
  geom_sf(data = stations_NO2, col = 'red')

# time series interpolation
data_filtered = rawdata %>% 
  filter(IDStations %in% stations_NO2$IDStations,
         Time >= '2021-01-01')
summary(data_filtered)
min(data_filtered$AQ_no2, na.rm = T)



# Preparation
p = nrow(stations_NO2)
t = length(unique(data_filtered$Time))

Y = matrix(log(data_filtered$AQ_no2+.5), p, t, byrow = T)
Y[is.nan(Y)] = NA
matplot(t(Y), t='l')
ggplot() + geom_sf(data = lombardy, fill='white') + 
  geom_sf(aes(col=Y[,100]), data = stations_NO2, size=4)

blh = matrix(scale(data_filtered$WE_blh_layer_max), p, t, byrow = T)
temp2m = matrix(scale(data_filtered$WE_temp_2m), p, t, byrow = T)
ws10m = matrix(scale(data_filtered$WE_wind_speed_10m_mean), p, t, byrow = T)
prec = matrix(scale(data_filtered$WE_tot_precipitation), p, t, byrow = T)
spr = matrix(scale(data_filtered$WE_surface_pressure), p, t, byrow = T)
X = abind(matrix(1, p, t), blh, temp2m, prec, spr, ws10m, along=3)
which(is.na(X)) # NA's are not allowed here!



dist_mat = fields::rdist(st_coordinates(stations_NO2))


# Model Fitting
source('main_function.R')

prior_list = list(
  V_beta_0 = 1e4, # Prior variance of initial state
  V_gamma = 1e6,  # Prior variance of constant coefficients
  a1 = 0.01,      # Prior shape for temporal variance
  b1 = 0.01,      # Prior rate for temporal variance
  s2_a = 0.01,    # Prior shape for measurement error variance
  s2_b = 0.01     # Prior rate for measurement error variance
)

nrep <- 2
nburn <- 2
thin <- 1
point.referenced = TRUE
print.interval = 1

mod = list()
tic()
mod <- DLM.st(Y, X, NULL, dist_mat, point.referenced, nrep, nburn, thin,
              print.interval, prior_list, mod$out)
toc()
saveRDS(mod, "dlm_lombardy.RDS")
str(mod$ave)
str(mod$out)


# some trace plot inspections
plot(mcmc(t(mod$out$S2_err_mis_))) # measurement error variance
plot(mcmc(t(mod$out$Q1inv_time^-1))) # temporal evolution variances
plot(mcmc(t(mod$out$RHO1_space_))) # partial sill of spatial effects
plot(mcmc(t(mod$out$RHO2_space_))) # range of spatial effects
plot(mcmc(t(mod$out$RHO1_space_time_))) # partial sill of spatio-temporal effects
plot(mcmc(t(mod$out$RHO2_space_time_))) # range of spatio-temporal effects

source('postprocessing.R')
b_overall = as.vector(mod$ave$B_postmean)
sdb_overall = sqrt(as.vector(mod$ave$B2_postmean) - b_overall^2)
print(data.frame(Mean=b_overall, ci_low=b_overall -qnorm(.975)*sdb_overall, 
                 ci_high=b_overall+qnorm(.975)*sdb_overall))
                 
plot.tvc(mod)
plot.svc(mod, stations_NO2, lombardy)
plot.stvc(mod, c(1, 10, 20, 30, 40, 50, 60, 70, 80))

plot.fitted(mod, Y, 10) # pick a number between 1 and p

log_like = mod$out$store_llike
dim(log_like) = c(p*t, nrep)
loo::waic(t(log_like))




# changing the hyperparameters #############
prior_list = list(
  V_beta_0 = 1e4, # Prior variance of initial state
  V_gamma = 1e6,  # Prior variance of constant coefficients
  a1 = 100,      # Prior shape for temporal variance
  b1 = 0.01,      # Prior rate for temporal variance
  s2_a = 0.01,    # Prior shape for measurement error variance
  s2_b = 0.01     # Prior rate for measurement error variance
)


mod2 = list()
tic()
mod2 <- DLM.st(Y, X, NULL, dist_mat, point.referenced, nrep, nburn, thin,
              print.interval, prior_list, mod2$out)
toc()
saveRDS(mod2, "dlm_lombardy_2.RDS")
str(mod2$ave)
str(mod2$out)


# some trace plot inspections
plot(mcmc(t(mod2$out$S2_err_mis_))) # measurement error variance
plot(mcmc(t(mod2$out$Q1inv_time^-1))) # temporal evolution variances
plot(mcmc(t(mod2$out$RHO1_space_))) # partial sill of spatial effects
plot(mcmc(t(mod2$out$RHO2_space_))) # range of spatial effects
plot(mcmc(t(mod2$out$RHO1_space_time_))) # partial sill of spatio-temporal effects
plot(mcmc(t(mod2$out$RHO2_space_time_))) # range of spatio-temporal effects


b_overall = as.vector(mod2$ave$B_postmean)
sdb_overall = sqrt(as.vector(mod2$ave$B2_postmean) - b_overall^2)
print(data.frame(Mean=b_overall, ci_low=b_overall -qnorm(.975)*sdb_overall, 
                 ci_high=b_overall+qnorm(.975)*sdb_overall))

plot.tvc(mod2)
plot.svc(mod2, stations_NO2, lombardy)
plot.stvc(mod2, c(1, 10, 20, 30, 40, 50, 60, 70, 80))

plot.fitted(mod2, Y, 10) # pick a number between 1 and p