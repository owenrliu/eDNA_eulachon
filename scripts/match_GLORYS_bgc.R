## Match modeled GLORYS daily temperature/salinity to eDNA data
library(tidyverse)
library(tidync)
library(here)
library(sf)
library(nngeo)

# this version of the GLORYS data was manually downloaded (in two chunks) because the R package isn't working
ncdf4::nc_open(here('data','glorys','cmems_mod_glo_bgc_my_0.25_P1M-m_1709245354119.nc')) # for some better metadata

nc <- tidync(here('data','glorys','cmems_mod_glo_bgc_my_0.25_P1M-m_1709245354119.nc'))
nc

ncvars <- hyper_vars(nc)

# for now, the variables we want are chl, no3, nppv, o2, and phyc
glorys <- nc %>% hyper_tibble()

# time units is 'seconds since 1970-01-01
glorys <- glorys %>% 
  mutate(date=as_datetime(time),
         month=month(date),
         year=year(date))

# match to the eDNA by nearest neighbors (nngeo package)
glorys_locs <- glorys %>% distinct(longitude,latitude) %>% 
  st_as_sf(coords=c("longitude","latitude"),crs=4326,remove=F) %>% 
  rename(lon_glorys=longitude,lat_glorys=latitude)

eDNA_glorys_locs <- read_rds(here('data','qPCR','eulachon qPCR 2019 and 2021 samples clean.rds')) %>% 
  st_as_sf(coords=c("lon","lat"),crs=4326,remove=F) %>% 
  # join nearest neighbors
  st_join(glorys_locs,st_nn) %>% st_set_geometry(NULL)

glorys_locs_needed <- eDNA_glorys_locs %>% 
  distinct(lon_glorys,lat_glorys,.keep_all = T) %>% 
  dplyr::select(lon_glorys,lat_glorys) %>% 
  mutate(need=1L)

# last, we need to make sure depths match up
# we make a matching key to associate glorys depth cells with depth categories from the eDNA collection
# however, we have to make sure that the bathymetry lines up: if the eDNA sample was taken
# from "deeper" than the bottom of a GLORYS cell, then we need to match the bottom cell

# first, filter by x-y locations we know we will need
glorys_filt1 <- glorys %>% 
  left_join(glorys_locs_needed,by=c('longitude'='lon_glorys','latitude'='lat_glorys')) %>% 
  filter(!is.na(need))

glorys_filt2 <- map_df(unique(eDNA_glorys_locs$depth_cat),~(mutate(glorys_filt1,depth_cat=.))) %>% 
  mutate(diff=abs(depth_cat-depth)) %>% 
  filter(!is.na(chl)) %>% 
  group_by(date,year,month,latitude,longitude,depth_cat) %>% 
  slice_min(diff) %>% 
  ungroup() %>% 
  dplyr::select(year,month,latitude,longitude,chl,no3,nppv,o2,phyc,depth_cat) %>% 
  mutate(across(c(latitude,longitude),~round(.,4)),.keep = "unused")

# match to the eDNA data
d_obs_glorys <- eDNA_glorys_locs %>% 
  mutate(across(c(lat_glorys,lon_glorys),~round(.,4)),.keep = "unused") %>% 
  left_join(glorys_filt2,by=c("lat_glorys"="latitude","lon_glorys"="longitude","month","year","depth_cat"))

d_obs_glorys %>% 
  filter(year==2019,depth_cat==0) %>% 
  ggplot(aes(lon,lat,size=no3))+
  geom_point()+
  coord_equal()

# this looked like it matched correctly (1 to 1, and the locations and depths look right)
write_rds(d_obs_glorys,here('data','eDNA_glorys_bgc_matched.rds'))
