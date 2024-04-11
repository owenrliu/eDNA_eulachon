# Map of sampling locations, 2019 and 2021
library(tidyverse)
library(sf)
library(rnaturalearth)

pred.crs <- terra::rast(here('data','raster_grid_blake','fivekm_grid.tif')) %>% st_crs()
coast <- ne_states(country='United States of America',returnclass = 'sf') %>%
  filter(name %in% c('California','Oregon','Washington','Nevada')) %>%
  st_transform(crs = pred.crs)
samps <- read_rds(here('data','eDNA_glorys_matched.rds')) %>% 
  mutate(Ct=replace_na(Ct,0)) %>% # delta-models in sdmTMB need this
  mutate(utm.lon.km=utm.lon.m/1000,
         utm.lat.km=utm.lat.m/1000) %>% 
  st_as_sf(coords=c("utm.lon.m","utm.lat.m"),crs=pred.crs)
bbox <- st_bbox(samps)
sampssfc <- samps %>% filter(depth_cat==0)

samples_plot <- ggplot()+
  geom_sf(data=coast,fill='gray50')+
  geom_sf(data=sampssfc,color='gray30',size=1)+
  facet_wrap(~year)+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
  coord_sf(datum=NA)

samples_plot
ggsave(here::here('plots','samples_map_basic.png'),samples_plot,h=6,w=5,bg='white')

# number of positive samples
numsamps <- samps %>% 
  st_set_geometry(NULL) %>% 
  filter(depth_cat %in% c(0,50,150)) %>% 
  mutate(ispos=Ct>0) %>% 
  count(year, depth_cat,ispos)
glimpse(numsamps)

numsamps_plot <- numsamps %>% 
  mutate(posd=ifelse(ispos,"Amplified","Did not\namplify")) %>% 
  ggplot(aes(factor(depth_cat),n,fill=posd))+
  geom_col()+
  facet_wrap(~year)+
  scale_fill_manual(values=c('gray20','#BD3786FF'))+
  labs(x="Depth of Sample",y="Number of Samples",fill="")
numsamps_plot
