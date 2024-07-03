library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(here)
library(nngeo)

# Blake's raster with just grid cell ID numbers
dat_raster=rast(here('Data','raster_grid_blake','fivekm_grid.tif'))
dat_samps <- read_rds(here('data','qPCR','eulachon qPCR 2019 and 2021 samples clean.rds'))

# Depth associated with grid cell IDs
raster_depth <- read_csv(here('data','raster_grid_blake','weighted_mean_NGDC_depths_for_5km_gridcells.csv')) %>% 
  rename(depth_m=WM_depth_m)

# join depths to cell numbers
cellnums <- tibble(rastID=as.numeric(values(dat_raster))) %>% 
  mutate(rastcell=row_number()) %>% 
  left_join(raster_depth,by=c("rastID"="Gridcell_ID"))

# make bathy raster
rast_depth <- setValues(dat_raster,cellnums$depth_m)
plot(rast_depth)

# reduce to a bounding box that matches the mesh
locs <- dat_samps %>%
  # mutate(utm.lat.km=utm.lat.m/1000,utm.lon.km=utm.lon.m/1000) %>% 
  distinct(year,station,utm.lon.m,utm.lat.m) %>% 
  st_as_sf(coords=c("utm.lon.m", "utm.lat.m"))
domain <-fmesher::fm_nonconvex_hull(locs,
                                    concave = -0.025,
                                    convex = -0.025) %>% 
  st_as_sf(crs=crs(dat_raster))

# mask raster
rast_depth_mask <- rast_depth %>% mask(domain,inverse=FALSE)
plot(rast_depth_mask)

# as a tibble for prediction
grid.pred <- as.data.frame(rast_depth_mask,xy=T) %>% 
  mutate(utm.lat.km=y/1000,utm.lon.km=x/1000,bathy.bottom.depth=-fivekm_grid)

ggplot(grid.pred)+
  geom_point(aes(utm.lon.km,utm.lat.km),size=0.5)+
  coord_equal()

write_rds(grid.pred,here('data','prediction_grid_5km_sdmTMB_no_covars.rds'))

# add depth categories
grid.pred0 <- grid.pred %>% mutate(depth_cat=0)
grid.pred50 <- grid.pred %>% filter(bathy.bottom.depth >= 50) %>% mutate(depth_cat=50)
grid.pred150 <- grid.pred %>% filter(bathy.bottom.depth >= 100) %>% mutate(depth_cat=150)

grid.predz <- bind_rows(grid.pred0,grid.pred50,grid.pred150) %>% 
  mutate(depth_fct=as.factor(depth_cat))

write_rds(grid.predz,here('data','prediction_grid_5km_3depths_no_covars.rds'))

#background map
coast <- ne_states(country='United States of America',returnclass = 'sf') %>% 
  filter(name %in% c('California','Oregon','Washington','Nevada')) %>%
  st_transform(crs = crs(dat_raster))

# plot by depth
bb <- st_bbox(rast_depth_mask)
ggplot()+
  geom_point(data=grid.predz,aes(utm.lon.km*1000,utm.lat.km*1000),size=0.3)+
  geom_sf(data=coast)+
  xlim(bb[1],bb[3])+ylim(bb[2],bb[4])+
  facet_wrap(~depth_cat)

# write_rds(grid.predz,here('data','prediction_grid_5km_sdmTMB.rds'))

## JOIN COVARIATES
library(tidync)
# for now, the variables we want are temperature (thetao), salinity (so) in psu, and bottom temperature (bottomT)
glorys <- purrr::map_df(list.files(here('data','glorys'),full.names = T),function(k){
  vars4D <- tidync(k) %>% hyper_tibble()
  vars3D <- tidync(k) %>% activate("D2,D1,D3") %>% hyper_tibble()
  left_join(vars4D,vars3D)
})

# time units is 'seconds since 1970-01-01
glorys <- glorys %>% 
  mutate(date=as_datetime(time),
         month=month(date),
         year=year(date))

# match to the grid by nearest neighbors (nngeo package)
glorys_locs <- glorys %>% 
  distinct(longitude,latitude) %>% 
  st_as_sf(coords=c("longitude","latitude"),crs=4326,remove=F) %>%
  st_transform(crs(dat_raster)) %>% 
  rename(lon_glorys=longitude,lat_glorys=latitude)

grid_glorys_locs <- grid.predz %>% 
  st_as_sf(coords=c("x","y"),crs=crs(dat_raster),remove=F) %>% 
  # join nearest neighbors
  st_join(glorys_locs,st_nn) %>% st_set_geometry(NULL)

glorys_locs_needed <- grid_glorys_locs %>% 
  distinct(lon_glorys,lat_glorys,.keep_all = T) %>% 
  dplyr::select(lon_glorys,lat_glorys) %>% 
  mutate(need=1L)

# last, we need to make sure depths match up
# we make a matching key to associate glorys depth cells with depth categories from the prediction grid

# first, filter by x-y locations we know we will need
glorys_filt1 <- glorys %>% 
  left_join(glorys_locs_needed,by=c('longitude'='lon_glorys','latitude'='lat_glorys')) %>% 
  filter(!is.na(need))

glorys_filt2 <- map_df(unique(grid.predz$depth_cat),~(mutate(glorys_filt1,depth_cat=.))) %>% 
  # make a copy of the filtered glorys data for each depth category,
  # find the difference between the depth category and the actual glorys depth
  mutate(diff=abs(depth_cat-depth)) %>% 
  filter(!is.na(thetao)) %>% 
  # then keep the matched row that has the smallest depth difference, but is not NA in GLORYS
  group_by(date,year,month,latitude,longitude,depth_cat) %>% 
  slice_min(diff) %>% 
  ungroup() %>% 
  dplyr::select(year,month,latitude,longitude,so,thetao,uo,vo,bottomT,depth_cat) %>% 
  mutate(across(c(latitude,longitude),~round(.,4)),.keep = "unused")

# match to the prediction grid.
# we will keep only values for AUGUST
grid_glorys_join <- map_df(c(2019,2021),~mutate(grid_glorys_locs,year=.)) %>% 
  mutate(month=8)
grid_glorys_join <- grid_glorys_join %>% 
  # copy for years 2019 and 2021,
  mutate(across(c(lat_glorys,lon_glorys),~round(.,4)),.keep = "unused") %>% 
  left_join(glorys_filt2,by=c("lat_glorys"="latitude","lon_glorys"="longitude","month","year","depth_cat"))

# the last value we have to join is for krill
# see krill matching script for details
# krill_preds <- read_rds(here('model output','krill_sdm_matched_to_pred_grid.rds'))
k4 <- read_rds(here('model output','krill_k2_k4_k5_matched_to_pred_grid.rds')) %>% 
  dplyr::select(-depth_cat) %>% 
  rename(year=Year)
grid_glorys_krill_join <- grid_glorys_join %>% 
  # left_join(krill_preds,by=join_by(x, y, fivekm_grid, bathy.bottom.depth, year)) %>% 
  left_join(k4)
grid_glorys_krill_join <- grid_glorys_krill_join %>% 
  filter(!is.na(thetao))

# this looked like it matched correctly (1 to 1, and the locations and depths look right)
write_rds(grid_glorys_krill_join,here('data','prediction_grid_5km_sdmTMB_with_covars.rds'))
