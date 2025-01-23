# JOIN MULTIPLE YEARS OF PROCESSED QPCR DATA, INCLUDING UNKNOWN (FIELD) SAMPLES, STANDARDS
# ALSO DO SOME CALCULATIONS OF IMPORTANT OFFSETS AND KEEP TRACK OF WHAT GETS FILTERED OUT

library(tidyverse)
library(here)
library(PNWColors)

# import samples data
dat19 <- read_rds(here('data','qPCR','eulachon qPCR 2019 joined cleaned 11_15_2023.rds'))
dat21 <- read_rds(here('data','qPCR','eulachon qPCR 2021 joined cleaned.rds'))

# import standards data
stand19 <- read_rds(here('data','qPCR','eulachon qPCR 2019 standards cleaned 11_20_2023.rds'))
stand21 <- read_rds(here('data','qPCR','eulachon qPCR 2021 standards cleaned 11_20_2023.rds'))

## Thin, clean, filter, bind samples data
## For now, we include only the necessary identifiers and minimal covariates

# 2019 data
glimpse(dat19)
# samples by depth category
dat19 %>% count(depth_cat)
stations_depths_count_2019 <- dat19 %>% count(station,depth_cat) # there should be a lot of 6s (2 bio reps x 3 tech. reps)
# Some depth/station combinations have more, because they went through dilution tests etc. For now, we will keep all uninhibited samples (se filtering below)

dat19_thin <- dat19 %>%
  filter(is.na(Zymo)) %>% 
  dplyr::select(date,year,month,day,time,station,lat,lon,utm.lon.m,utm.lat.m,qPCR,inhibition_rate,drop.sample,dilution,volume,Ct,copies_ul,depth_cat,bathy.bottom.depth,transect_dist_m)
# how many samples are diluted?
dat19_thin %>% count(dilution)

# 2021 data
glimpse(dat21)
# samples by depth category
dat21 %>% count(depth_cat)
stations_depths_count_2021 <- dat21 %>% count(station,depth_cat)
dat21_thin <- dat21 %>% 
  dplyr::select(date,year,month,day,time,station,lat,lon,utm.lon.m,utm.lat.m,qPCR,inhibition_rate,drop.sample,dilution,volume,Ct,copies_ul,depth_cat,bathy.bottom.depth,transect_dist_m) %>% 
  mutate(inhibition_rate=as.numeric(inhibition_rate),
         volume=as.numeric(volume))

datjoin <- dat19_thin %>% 
  bind_rows(dat21_thin)
# 16886 samples

# count of samples by year, station, and depth category. Note that this still includes inhibited samples and various controls and tests
stations_depths_count_all_raw <- datjoin %>% count(year,station,depth_cat)

# Filter (but keep track of what we filter out for later discussion)
# first filter= only keep non-inhibited samples
datjoin <- datjoin %>%
  filter(inhibition_rate<0.5)
# 14498 samples

# count NAs across vars
NA_tracker <- datjoin %>% summarise(across(everything(), ~sum(is.na(.))))
NA_tracker
# most vars have 2084 rows with NA- look at these
NADate <- datjoin %>% filter(is.na(date))# 2084 obs
# a lot of these are non-template controls or other types of methodological controls, but not 100% sure
datjoin <- datjoin %>% 
  filter(!is.na(date))
# 12414 samples

NA_tracker2 <- datjoin %>% summarise(across(everything(), ~sum(is.na(.))))
NA_tracker2

# depth category is NA for some other controls
NAdepthcat <- datjoin %>% filter(is.na(depth_cat))

datjoin <- datjoin %>% 
  filter(!is.na(inhibition_rate),!is.na(depth_cat))
# 12293 samples

NA_tracker3 <- datjoin %>% summarise(across(everything(), ~sum(is.na(.))))
NA_tracker3

# dilution is NA for many samples, but NAs should be 1. But keep a record of them in case we need for later QA/QC
NAdilut <- datjoin %>% filter(is.na(dilution))

datjoin <- datjoin %>% 
  mutate(dilution=ifelse(is.na(dilution),1,dilution))

NA_tracker4 <- datjoin %>% summarise(across(everything(), ~sum(is.na(.))))
NA_tracker4

# volume is NA in 9 other samples
NAvolume <- datjoin %>% filter(is.na(volume))
datjoin <- datjoin %>% filter(!is.na(volume))
#12284 samples

stations_depths_count_all_filt<- datjoin %>% count(year,station,depth_cat)

# Thin, clean, bind standards data

standjoin <- mutate(stand19,year=2019)%>% 
  bind_rows(mutate(stand21,year=2021)) %>% 
  rename(plate=qPCR) %>% 
  mutate(Ct=na_if(Ct,-99)) %>% 
  # make a dummy row identifier so we can more easily make substitutions
  mutate(sid=row_number())

### IMPORTANT DATA QA/QC STEP ###
# some of the standards seem to have been labeled as unknowns
# find them and fix them
prbs <- standjoin %>% filter(hake_task=="STANDARD",task!="STANDARD") #mislabeled as unknowns
# find rows where the first character of known copies is NOT 1 or 5
prbs2 <- standjoin %>%filter(!(as.numeric(substr(copies_ul, 1, 1))%in%c(1,5))) 
prbs <- bind_rows(prbs,prbs2) %>% distinct()
# fix copies_ul for these samples based on the following rule
# if sample== E00, that should be 1 copy
# likewise, E01 should be 10 copies, E02 is 100 copies, etc.
# "5" should actually be the weird E00/01 or E0.5, which is 5 copies
standjoin$copies_ul[standjoin$sid==117] <- 5
standjoin$copies_ul[standjoin$sid==287] <- 1
standjoin$copies_ul[standjoin$sid==531] <- 1
standjoin$copies_ul[standjoin$sid==687] <- 100
standjoin$copies_ul[standjoin$sid==697] <- 1
standjoin$copies_ul[standjoin$sid==723] <- 1
standjoin$copies_ul[standjoin$sid==762] <- 10
standjoin$copies_ul[standjoin$sid==771] <- 1
standjoin$copies_ul[standjoin$sid==791] <- 1
standjoin$copies_ul[standjoin$sid==73] <- 1
standjoin$copies_ul[standjoin$sid==784] <- 10

#SHOULD CHECK SID #1342- ONLY STANDARD AT 5 COPIES THAT DID NOT AMPLIFY?

###

# select only the essential columns for output
standjoin_thin <- standjoin %>% 
  dplyr::select(plate,Ct,known_conc_ul=copies_ul,year)

glimpse(standjoin_thin)
unique(standjoin_thin$known_conc_ul) # looks way better!

#### OFFSETS AND CORRECTIONS ####

# ADD CORRECTION FOR 2uL ADDITION INSTEAD OF 1
standjoin_thin <- standjoin_thin %>% mutate(known_copies = known_conc_ul*2)

# ADD OFFSETS FOR DILUTION AND VOLUME FOR UNKNOWN SAMPLES
datjoin <- datjoin %>% 
  # volume offset
  mutate(ln_vol_offset=log(volume/2.5),
         # dilution/inhibition offset       
         ln_dil_offset=log(dilution))
# expansion factor to copies/L
datjoin <- datjoin %>% 
  mutate(ln_expand_offset=log(1/20))

# we can add the offsets together to create just one offset vectors since they are additive in log space
datjoin <- datjoin %>% 
  mutate(offsets_all=ln_vol_offset+ln_dil_offset+ln_expand_offset)
  
# for the 2019 wash error
# find samples that were washed, then find their pairs
datjoin <- datjoin %>% 
  mutate(washed=ifelse(drop.sample %in% c("30EtOH","30EtOHpaired"),1,0))
#

# count up the distribution of these offsets
datjoin %>% 
  count(volume,ln_vol_offset)
datjoin %>% 
  count(dilution,ln_dil_offset)
datjoin %>% 
  count(washed)

# Set some final variable names
datjoin_final <- datjoin %>% 
  rename(plate=qPCR) %>% 
  mutate(Ct=na_if(Ct,-99)) %>% 
  dplyr::select(-inhibition_rate)

# Save!
write_rds(standjoin_thin,here('data','qPCR','eulachon qPCR 2019 and 2021 standards clean.rds'))
write_rds(datjoin_final,here('data','qPCR','eulachon qPCR 2019 and 2021 samples clean.rds'))
