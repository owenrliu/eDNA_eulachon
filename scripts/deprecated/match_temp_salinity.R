# Explore MOCHA Data for Temperature and Salinity
# https://www.nodc.noaa.gov/archive/arc0212/0277984/2.2/data/0-data/
# https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0277984

library(tidyverse)
library(here)
library(magrittr)
library(sf)

dat <- read_csv(here('data','aggregated_daily_dataset.csv'),
                col_types ='Ddddddddddddddddddddddddddddddddcccccccccc')
dat %<>%
  filter(year(time_utc)>2018)

dat2019 <- dat %>% 
  filter(year(time_utc)==2019) %>% 
  dplyr::select(time_utc:sal_pss)

locs_temperature <- dat2019 %>% 
  filter(!is.na(t_C)) %>% 
  distinct(latitude,longitude) %>% 
  ggplot(aes(longitude,latitude))+
  geom_point()+
  coord_equal()

locs_temperature

locs_sal <- dat2019 %>% 
  filter(!is.na(sal_pss)) %>% 
  distinct(latitude,longitude) %>% 
  ggplot(aes(longitude,latitude))+
  geom_point()+
  coord_equal()

locs_sal

## GLORYS DATA!!
# Script to download GLORYS data directly from Copernicus Servers
# Requires a Copernicus marine service account
# https://data.marine.copernicus.eu/products

# install.packages("CopernicusMarine")
library(CopernicusMarine)
library(tidyverse)
library(tidync)
library(here)


