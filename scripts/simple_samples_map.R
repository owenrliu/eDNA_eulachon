# Map of sampling locations, 2019 and 2021
library(tidyverse)
library(sf)
library(rnaturalearth)
library(here)

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
sampssfc <- samps %>% 
  filter(depth_cat %in% c(0,50,150)) %>%
  mutate(ispos=Ct>0) %>% 
  mutate(ispos=ifelse(ispos==T,"Amplified","Did not\namplify"))

samples_plot <- ggplot()+
  geom_sf(data=coast,fill='gray50')+
  geom_sf(data=sampssfc,color='gray30',size=1)+
  facet_wrap(~year)+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
  coord_sf(datum=NA)

samples_plot
ggsave(here::here('plots','samples_map_basic.png'),samples_plot,h=6,w=5,bg='white')

#presence/absence
samples_pa_plot <- ggplot()+
  geom_sf(data=coast,fill='gray50')+
  geom_sf(data=sampssfc,aes(color=factor(ispos),size=factor(ispos)),alpha=0.7)+
  scale_color_manual(values=c('#5DC863FF','gray20'))+
  # scale_color_manual(values=c('#BD3786FF','gray20'))+
  scale_size_manual(values=c(2,1))+
  facet_wrap(~year)+
  labs(color="",size="")+
  theme(legend.position=c(0.85,0.6),legend.background = element_rect(fill="white",color='black'))+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])+
  coord_sf(datum=NA)

samples_pa_plot
ggsave(here::here('plots','samples_map_eul_pa.png'),samples_pa_plot,h=6,w=5,bg='white')

# number of positive samples
numsamps <- samps %>% 
  st_set_geometry(NULL) %>% 
  filter(depth_cat %in% c(0,50,150)) %>% 
  mutate(ispos=Ct>0) %>% 
  count(year, depth_cat,ispos)
glimpse(numsamps)

numsamps_plot <- numsamps %>% 
  mutate(depth_cat=factor(depth_cat,levels=c("150","50","0"))) %>%  
  mutate(posd=ifelse(ispos,"Amplified","Did not\namplify")) %>% 
  ggplot(aes(depth_cat,n,fill=posd))+
  geom_col()+
  facet_wrap(~year,nrow=2)+
  coord_flip()+
  scale_fill_manual(values=c('#5DC863FF','gray20'))+
  # scale_fill_manual(values=c('gray20','#BD3786FF'))+
  labs(x="Depth of Sample",y="Number of Samples",fill="")
numsamps_plot

# what about cumulative amplications with latitude

cumulative_samps <- sampssfc %>% 
  mutate(ispos=if_else(ispos=="Amplified",1,0)) %>% 
  dplyr::select(year,lat,depth_cat,ispos) %>% 
  group_by(year,depth_cat) %>% 
  arrange(lat) %>% 
  mutate(cumpos=cumsum(ispos),cumtot=row_number()) %>% 
  ungroup()
p_cumulative_samps<- cumulative_samps %>% 
  ggplot(aes(lat,cumpos,color=factor(depth_cat)))+
  facet_wrap(~year)+
  geom_line(size=1.5)+
  geom_rug(sides='b')+
  geom_vline(xintercept=41.0,linetype=2)+
  scale_color_manual(values=PNWColors::pnw_palette("Starfish",3))+
  theme_minimal()+
  labs(x="Latitude",y="Cumulative samples successfully amplified",color="Sample\nDepth")+
  theme(panel.border = element_rect(color='black',fill=NA))

p_cumulative_samps  
ggsave(here::here('plots','cumulative_amplifications_latitude.png'),p_cumulative_samps,h=5,w=7,bg='white')
    