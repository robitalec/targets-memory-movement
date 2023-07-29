# Script to demonstrate how to set up SSF from Merkle et al. 2019 (Ecol Lett)
# WRitten by Jerod Merkle, January 2023

require(tidyverse)
require(sf)
require(raster)
require(lubridate)
library(amt)

# location where the spatial data are stored
setwd("G:/Shared drives/wyo-coop-merkle/Memory_Movement_review/data")

#----------------------------#
# Read in positional data ####
#----------------------------#

# This is 1 individual mule deer collared for nearly 2 years
# in NW Wyoming. These data were collected by Hall Sawyer (WEST Inc)
dat_all <- readRDS("GPSdata_muldeer445.rds")

#-----------------------#
# Prep movement data ####
#-----------------------#

# import or develop table of migration start and stop dates for each year and season
migtable <- data.frame(id=445, 
                       year=c(2014, 2014, 2015, 2015),
                       season=c("spring","fall","spring","fall"),
                       date_start=c("2014-05-04","2014-11-02","2015-04-22","2015-11-13"),
                       date_end=c("2014-05-17","2014-11-21","2015-05-13","2015-11-24")) %>% 
  mutate(date_start=ymd(date_start, tz="America/Denver"),
         date_end=ymd(date_end, tz="America/Denver"))

# add column to denote year of monitoring
migtable <- migtable %>% 
  group_by(id) %>% 
  arrange(year, desc(season)) %>% 
  mutate(year_cumm=year-min(year)+1)

# Isolate data to be analyzed 
# (i.e., analyze second year of data, while using the first year to index of memory)
# Here, we grab the spring migration for the second year of monitoring
dat2keep <- migtable %>% 
  filter(year_cumm==2, season=="spring")

dat <- do.call(rbind, lapply(1:nrow(dat2keep), function(i){
  dat_all %>% 
    filter(id==dat2keep$id[i], 
           date >= dat2keep$date_start[i],
           date <= dat2keep$date_end[i]) %>% 
    mutate(season="spring") %>% 
    return()
}))

#---------------------------#
# Prep indices of memory ####
#---------------------------#

# 1) Isolate previous year's migration routes
# Loop through each id/season identified in dat2keep
prev.routes <- do.call(rbind, lapply(1:nrow(dat2keep), function(i){
  
  # get migtable organized for the previous year of data but for the id in question
  mig.i <- migtable %>% 
    filter(id==dat2keep$id[i], year_cumm==1, season=="spring")
  
  # isolate the points for the previous year's migration in the same season
  pts <- dat_all %>% 
    filter(id==mig.i$id, 
           date >= mig.i$date_start,
           date <= mig.i$date_end)
  
  # turn into a sf line dataframe
  ln <- pts %>% 
    st_combine() %>% 
    st_cast("LINESTRING") %>% 
    data.frame(id=mig.i$id, 
               season="spring",
               year=mig.i$year,
               geometry=.) %>% 
    st_as_sf(sf_column_name = "geometry")
  
  return(ln)
}))


# 2) Isolate centroid of previous year's summer range
# Again, loop through each id/season identified in dat2keep
prev.ranges <- do.call(rbind, lapply(1:nrow(dat2keep), function(i){
  
  # get migtable organized for the previous year of data but for the id in question
  mig.i <- migtable %>% 
    filter(id==dat2keep$id[i], year_cumm==1)
  
  # isolate the points for the period bewteen the end of the spring migration
  # and the start of the fall migration (i.e., summer)
  pts <- dat_all %>% 
    filter(id==mig.i$id, 
           date >= mig.i$date_end[mig.i$season=="spring"],  
           date <= mig.i$date_start[mig.i$season=="fall"])
  
  # calculate centroid and turn into a sf point dataframe
  centr <- pts %>% 
    st_union() %>% 
    st_centroid() %>% 
    data.frame(id=mig.i$id[1], 
               season="summer",
               year=mig.i$year[1],
               geometry=.) %>% 
    st_as_sf(sf_column_name = "geometry")
  
  return(centr)
}))

rm(dat2keep, migtable)

#--------------#
# Build SSF ####
#--------------#

# translate location data to tracks
# Note, from here on out, there is no code to handle multiple individuals.
trk <- dat %>% 
  mutate(x = st_coordinates(.)[,1], # add x and y as columns
         y = st_coordinates(.)[,2]) %>% 
  st_drop_geometry() %>%   #need to remove the geometry column
  make_track(x, y, date, id = id, season = season, 
             crs = st_crs(dat))
trk

# transform into a crs with meters as unit
trk <- trk %>% 
  transform_coords(5072)

# have a look at step lengths
trk <- trk %>% 
  mutate(sl_ = step_lengths(.))

hist(trk$sl_)
median(trk$sl_, na.rm=TRUE)

# sampling rate is very consistent so don't need to worry about resampling  
summarize_sampling_rate(trk) 

# create steps
ssf1 <- trk %>% 
  steps()
ssf1

# check to see if some directional persistence is present
hist(ssf1$ta_) # indeed there is, so the SSF is properly set up!

# draw random steps
ssf1 <- ssf1 %>% 
  random_steps()
ssf1

#-------------------------------#
# extract habitat covariates ####
#-------------------------------#

# read in variables first
# slope derived from USGS NED at approx 30m resolution
slope <- raster("./Slope30m.tif")

# Terrain Position Index derived from USGS NED at approx 30m resolution
tpi <- raster("./TPI_9x9window_30m.tif")

# TRASP derived from USGS NED at approx 30m resolution calculated based on
# trasp = ( 1 - cos((pi/180)(a-30) ) / 2 where; a = aspect in degrees
# values near 1 are SSW slopes and values near 0 are NNE. 
# Note: to save space, values were multiplied by 100 and saved as integers
trasp <- raster("./TRASP30m.tif")

# Tree cover from the National Landcover Dataset (2019): https://www.mrlc.gov/
treecov <- raster("./TreeCover_nlcd2019_30m.tif")

# Biomass of Forbs and grasses from the RangeLand Analysis Platform
# year 2015, units are lbs/acre (https://rangelands.app/)
herbs <- raster("./RAP_2015_Biomass_ForbsGrasses.tif")

# Need to build separate sf object for extraction
# because amt cannot deal with different projections 
# of the enviro data
step_ends <- ssf1 %>% 
  as.data.frame() %>% 
  select(x2_, y2_) %>% 
  st_as_sf(coords=c("x2_","y2_"), dim="XY", crs=5072)

ssf1$tpi <- raster::extract(tpi, st_transform(step_ends, projection(tpi)))
ssf1$trasp <- raster::extract(trasp, st_transform(step_ends, projection(trasp)))
ssf1$slope <- raster::extract(slope, st_transform(step_ends, projection(slope)))
ssf1$treecov <- raster::extract(treecov, st_transform(step_ends, projection(treecov)))
ssf1$herbs <- raster::extract(herbs, st_transform(step_ends, projection(herbs)))

rm(tpi, trasp, slope, treecov, herbs)

#------------------------------#
# extract memory components ####
#------------------------------#

# 1) distance to previous route
prev.routes <- prev.routes %>% 
  st_transform(5072)  # transform to projection of the tracks

ssf1$dist2prevroute <- as.numeric(st_distance(step_ends, prev.routes))

# 2) direction to previous summer range
prev.ranges <- prev.ranges %>% 
  st_transform(5072) %>% # transform to projection of the tracks
  mutate(x = st_coordinates(.)[,1], # add x and y as columns
         y = st_coordinates(.)[,2])

# calculate direction from source location to previous range (in radians)
direct2prev_rng_source <- atan2((prev.ranges$y-ssf1$y1_), (prev.ranges$x-ssf1$x1_))
# calculate direction from source location to target location (in radians)
direct2target <- atan2((ssf1$y2_-ssf1$y1_), (ssf1$x2_-ssf1$x1_)) # change this to direct2target
# subtract the two and take the cosine
ssf1$direct2prev_rng <- cos(direct2prev_rng_source-direct2target)
# Values of 1 represent a step directly towards previous range
# Values of 0 represent a step directly away from previous range

#-------------#
# Fit SSFs ####
#-------------#

# add log of step length
ssf1 <- ssf1 %>% 
  mutate(sl_log = log(sl_))

# fit environmental model
m_envir <- ssf1 %>% 
  fit_clogit(case_ ~ sl_log+tpi+slope+trasp+herbs+treecov+
               strata(step_id_))

summary(m_envir)
# Interpretation: This deer select steps that go to
# more ridgy areas, warmer aspects with higher biomass of herbs


# fit environmental model + memory
m_envir_memor <- ssf1 %>% 
  fit_clogit(case_ ~ sl_log+tpi+slope+trasp+herbs+treecov+
               dist2prevroute+direct2prev_rng+
               strata(step_id_))

summary(m_envir_memor)

# compare relative empirical support for each model
AIC(m_envir)
AIC(m_envir_memor) 

# Interpretation: This deer select steps that go to
# more ridgy areas, warmer aspects with higher biomass of herbs
# and they select steps closer to their previous route and 
# in the direction of their previous summer range
