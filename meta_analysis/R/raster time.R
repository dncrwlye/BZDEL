library(tidyverse)
library(stringi)
library(readxl)
library(binom)
library(lubridate)
library(plotly)
library(ggmap)
library(raster)
library(rgdal)
library(broom)
library(maptools)
library(rgeos)
library(sp)
library(googleway)

#raster fun
#.................load data from seroprevalence_clean_script.R........................

key <- "AIzaSyD40FLStnvp085UB-FKXbyONuxV3ke4umY"
register_google(key = key)

setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
load("Data/seroprevalence.Rdata")

seroprevalence_search <- as.data.frame(unique(seroprevalence$sampling_location)) %>%
  rename(sampling_location = `unique(seroprevalence$sampling_location)`) %>%
  mutate(sampling_location = as.character(sampling_location)) %>%
  mutate(north=NA) %>%
  mutate(south=NA) %>%
  mutate(west=NA) %>%
  mutate(east=NA) %>%
  mutate(address = NA) 

for(i in 1:nrow(seroprevalence_search))
{
  query <- seroprevalence_search$sampling_location[i] 
  rd <- geocode(query, output = 'more', source = 'google')
  Sys.sleep(4)
  if(is.atomic(rd)==FALSE)
  {
    if (is.na(rd$lon) == FALSE) 
    {
      seroprevalence_search$north[i] <- rd$north
      seroprevalence_search$south[i] <- rd$south
      seroprevalence_search$east[i]  <- rd$east
      seroprevalence_search$west[i]  <- rd$west
      seroprevalence_search$address[i]  <- rd$address
    }
  }
}

filter(seroprevalence_search, is.na(address))
seroprevalence_search_2 <- as.data.frame(unique(seroprevalence$sampling_location_two)) %>%
  rename(sampling_location_two = `unique(seroprevalence$sampling_location_two)`) %>%
  mutate(sampling_location_two = as.character(sampling_location_two)) %>%
  filter(!is.na(sampling_location_two)) %>%
  mutate(north_two=NA) %>%
  mutate(south_two=NA) %>%
  mutate(west_two=NA) %>%
  mutate(east_two=NA) %>%
  mutate(address_two = NA) 


for(i in 1:nrow(seroprevalence_search_2))
{
  query <- seroprevalence_search_2$sampling_location[i] 
  rd <- geocode(query, output = 'more', source = 'google')
  Sys.sleep(10)
  if(is.atomic(rd)==FALSE)
  {
    if (is.na(rd$lon) == FALSE) 
    {
      seroprevalence_search_2$north_two[i] <- rd$north
      seroprevalence_search_2$south_two[i] <- rd$south
      seroprevalence_search_2$east_two[i]  <- rd$east
      seroprevalence_search_2$west_two[i]  <- rd$west
      seroprevalence_search_2$address_two[i]  <- rd$address
    }
  }
}

filter(seroprevalence_search_2, is.na(address_two))

seroprevalence <- full_join(seroprevalence, seroprevalence_search, by=c("sampling_location"))
seroprevalence <- full_join(seroprevalence, seroprevalence_search_2, by=c('sampling_location_two'))

seroprevalence_geolocation <- seroprevalence %>%
  mutate(north_final = pmax(north, north_two, na.rm=TRUE)) %>%
  mutate(south_final = pmin(south, south_two, na.rm=TRUE)) %>%
  mutate(west_final = pmin(west, west_two, na.rm=TRUE)) %>%
  mutate(east_final = pmax(east, east_two, na.rm=TRUE)) 


#geocode('ghana, tanoboase', output = 'more', source = 'google')

save(seroprevalence_geolocation, file='data/seroprevalence_geolocation.Rdata')


#..............moving onto the ecoregions analyses.................................

load(file = 'data/seroprevalence_geolocation.Rdata')
ecos <- shapefile('data/official_teow/wwf_terr_ecos.shp')

#we should probably go back and use polygons over polygons, but that is proving way too hard rn

seroprevalence_x <- seroprevalence_geolocation%>%
  dplyr::select(-c(north, south, west, east,  north_two, south_two, west_two, east_two)) %>%
  mutate(coordinate_box = NA)

seroprevalence_x_unique <- seroprevalence_x %>%
  dplyr::ungroup() %>%
  dplyr::select(north_final, south_final, west_final, east_final) %>%
  unique()

eco_regions_placeholder <- as.data.frame(matrix(ncol=5))
colnames(eco_regions_placeholder) <- (c("ECO_NAME","north_final","south_final","west_final","east_final"))

coordinate_box <- list()

for(i in 1:nrow(seroprevalence_x_unique))
{
  if(is.na(seroprevalence_x_unique[i,'north_final'])) next 
  x <- as(raster::extent(as.numeric(seroprevalence_x_unique[i,3]), as.numeric(seroprevalence_x_unique[i,4]), 
                         as.numeric(seroprevalence_x_unique[i,2]), as.numeric(seroprevalence_x_unique[i,1])), "SpatialPolygons")
  
  #proj4string(x) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  proj4string(x) <- proj4string(ecos)
  coordinate_box[[i]] <- x
  
  m <- as.data.frame(over(coordinate_box[[i]], ecos, returnList = FALSE)) %>%
    unique() %>%
    mutate(north_final = as.numeric(seroprevalence_x_unique[i, "north_final"])) %>%
    mutate(south_final = as.numeric(seroprevalence_x_unique[i, "south_final"])) %>%  
    mutate(west_final = as.numeric(seroprevalence_x_unique[i, "west_final"])) %>%
    mutate(east_final = as.numeric(seroprevalence_x_unique[i, "east_final"])) 
  
  eco_regions_placeholder <- bind_rows(eco_regions_placeholder, m)
  print(i)
}

seroprevalence_geolocation_final <- full_join(seroprevalence_x, eco_regions_placeholder)

y <- seroprevalence_x_final %>%
  filter(is.na(ECO_NAME))

save(seroprevalence_geolocation_final, file ="data/seroprevalence_ecoregions_final.Rdata")

