# river discharge
#RC4USCoast: A river chemistry dataset for regional ocean model application in the U.S. East, Gulf of Mexico, and West Coasts 
#from 1950-01-01 to 2022-12-31 (NCEI Accession 0260455)
#https://www.ncei.noaa.gov/data/oceans/ncei/ocads/metadata/0260455.html

library(tidyverse)
library(tidync)
library(here)
library(sf)
library(viridis)
library(nngeo)

# coordinate reference system
pred.crs <- terra::rast(here('data','raster_grid_blake','fivekm_grid.tif')) %>% st_crs()

x <- tidync(here('data','river discharge','mclim_19902022_disc.nc'))
ncdf4::nc_open(here('data','river discharge','mclim_19902022_disc.nc'))
x
locs <- x %>% activate("D1") %>% hyper_tibble() %>% filter(region=="West Coast")
disc <- x %>% hyper_tibble()

climat <- locs %>% left_join(disc,by = join_by(RC4USCoast_ID)) %>% 
  mutate(time=as_date("1950-01-01")+days(time)) %>% 
  ggplot(aes(month(time),disc,color=river_name))+
  geom_line()

locs_disc <- locs %>% left_join(disc,by = join_by(RC4USCoast_ID)) %>% 
  mutate(time=as_date("1950-01-01")+days(time)) %>%
  mutate(month=month(time)) %>% 
  filter(month %in% c(1,2,3,4,5)) %>% 
  group_by(mouth_lat,mouth_lon,river_name) %>% 
  summarise(discharge_mean=mean(disc))

latmap <- locs_disc %>% 
  ggplot(aes(mouth_lat,1))+
  coord_flip()+
  geom_rug(sides='b')+
  theme(axis.text.x=element_blank(),axis.title = element_blank())
latmap

# locs_disc_filt <- locs_disc %>% ungroup() %>% slice_max(n=10,order_by=discharge_mean)
# 
# abun_index_plot <- ggplot()+
#   geom_line(data=abun_yr_lat,aes(lat_glorys,log10(abun_idx),color=factor(depth_cat)))+
#   geom_vline(data=locs_disc_filt,aes(xintercept=lat),color='red')+
#   coord_flip()+
#   scale_color_manual(values=PNWColors::pnw_palette("Starfish",3))+
#   facet_wrap(~year)+
#   labs(y="Log10(Total eDNA)",x="Latitude",color="Depth")+
#   theme(text=element_text(size=14))
# abun_index_plot

# Calculate a metric for eulachon SDMs
# We'll use an inverse distance-weighted surface of "river influence". We take our prediction grid, then use the locations
# of river mouths to create a "surface" of each river's influence weighted by distance. We can then stack these across all rivers
# and multiply them by each river's March-May mean discharge to make our derived metric

gr <- read_rds(here('data','prediction_grid_5km_sdmTMB_no_covars.rds')) %>% st_as_sf(coords=c('x','y'),crs=pred.crs)
glimpse(gr)
river_locs <- locs %>% dplyr::select(river_name,mouth_lat,mouth_lon) %>% 
  st_as_sf(coords=c('mouth_lon','mouth_lat'),crs=4326) %>% 
  st_transform(pred.crs) %>% 
  mutate(utm.lon.km=st_coordinates(.)[,1]/1000,utm.lat.km=st_coordinates(.)[,2]/1000)

# weighting factors will be 1 divided by the square of the distance to each river
make_weights_surface <- function(name_of_river,decay_term=10){ # decay term here controls the strength of decay over distance
  river_info <- locs_disc %>% filter(river_name==name_of_river)
  river_xy <- river_info %>% dplyr::select(mouth_lon,mouth_lat) %>% 
    st_as_sf(coords=c('mouth_lon','mouth_lat'),crs=4326) %>% 
    st_transform(pred.crs)
  nns <- st_nn(river_xy,gr,k=nrow(gr),returnDist = T)
  idx <- nns$nn[[1]]
  dists=nns$dist[[1]]/1000 %>% as.numeric
  weighted_discharge<-river_info$discharge_mean*exp(-dists/decay_term) # this is the key exponential decay relationship!
  discharge_df <- tibble(gid=idx,distance=dists,weighted_disc=weighted_discharge)
  out <- gr %>% mutate(gid=row_number()) %>% mutate(river_name=name_of_river,
                                                    disc_mean=river_info$discharge_mean) %>% 
    left_join(discharge_df,by="gid")
  out
}

riverine_influence_long <- purrr::map_df(unique(river_locs$river_name),make_weights_surface,decay_term=10) %>% 
  st_set_geometry(NULL)

# test
riverine_influence_long %>%  ggplot(aes(distance,weighted_disc,color=river_name))+
  geom_point(alpha=0.5)+
  xlim(0,400)

riverine_influence <- riverine_influence_long %>% 
  group_by(gid,utm.lat.km,utm.lon.km,bathy.bottom.depth) %>% 
  summarise(river_input=sum(weighted_disc)) %>% 
  ungroup()

# test map for one river
ggplot()+ 
  geom_raster(data=riverine_influence_long %>% filter(river_name=='Klamath'), aes(utm.lon.km, utm.lat.km, fill = weighted_disc)) +
  theme_minimal()+
  labs(x="",y="")+
  coord_equal()+
  scale_fill_viridis(option='magma')+
  theme(panel.border = element_rect(color='black',fill=NA),
        axis.text.x=element_blank())


ggplot()+ 
  geom_raster(data=riverine_influence, aes(utm.lon.km, utm.lat.km, fill = river_input)) +
  theme_minimal()+
  labs(x="",y="",title="River Influence Metric")+
  coord_equal()+
  scale_fill_viridis(option='magma',trans='log10')+
  theme(panel.border = element_rect(color='black',fill=NA),
        axis.text.x=element_blank())

# save
write_rds(riverine_influence,here('data','river_influence_grid_matched.rds'))
