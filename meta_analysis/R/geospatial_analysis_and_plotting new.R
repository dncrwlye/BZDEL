### mapping samples to raster polygons

# Libraries --------------------------------------------------------------
library(raster)
library(rgdal)
library(tidyverse)
library(ggmap)
library(broom)
library(maptools)
library(rgeos)
library(plotly)
library(sp)
setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis")


rm(list=ls())
load('data/seroprevalence.Rdata')

by_country_sampling_effort_for_map <- seroprevalence %>%
  mutate(sample_size = ifelse(is.na(sampling_location_two), sample_size,
                       ifelse(is.na(sampling_location_three) & !is.na(sampling_location_two), sample_size/2,
                       ifelse(is.na(sampling_location_four) & !is.na(sampling_location_three) & !is.na(sampling_location_two), sample_size/3,
                       sample_size/4)))) %>%
  gather(key = 'location_repeat', 'location',
         'sampling_location', 
         'sampling_location_two', 
         'sampling_location_three', 
         'sampling_location_four' 
         ) %>%
  filter(!is.na(location)) %>%
  mutate(country = gsub("(.*),.*", "\\1", location)) %>%
  mutate(country = gsub("(.*),.*", "\\1", country)) %>%
  mutate(prev.cat = ifelse(seroprevalence_percentage == 0, 0, 1)) %>%
  mutate(methodology.pcr = ifelse(methodology == 'PCR based method', 'PCR based method', 'Ab')) %>%
  group_by(country, virus, methodology.pcr) %>%
  summarise(n = n())

countries <- shapefile('data/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp')
countries.points = fortify(countries)


countries.data <- countries@data
countries.data$id <- rownames(countries@data)

countries.join <- full_join(countries.data, countries.points, by = 'id')

countries.join <- countries.join %>%
  mutate(country = tolower(ADMIN)) %>%
  mutate(fill.background = 1)

countries.join <- full_join(countries.join, by_country_sampling_effort_for_map, by = 'country')
countries.join <- countries.join %>%
  filter(!is.na(virus))

lat <- mean(countries.join$lat, na.rm =TRUE)
lon <- mean(countries.join$long, na.rm = TRUE)
map <- get_map(location = c(lon = mean(lon), lat = mean(lat)), zoom = 1,
               maptype = "satellite", source = "google")

sampling.map <- ggmap(map)+  
  scale_x_continuous(limits = c(-200,200)) +
  scale_y_continuous(limits = c(-100,100)) +
  geom_polygon(data=countries.points,aes(x=long,y=lat, group=group)) +
  geom_polygon(data=countries.join,aes(x=long,y=lat,fill = n, group=group)) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_continuous(low="green4", high="aquamarine1") +
  facet_grid(virus~methodology.pcr) 



sampling.map
ggsave("figures/sampling.map.png", bg = "transparent", height = 18, width = 18)

#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................
#...............................Positives Positive Positives...................................

by_country_positive_for_map <- seroprevalence %>%
  mutate(sample_size = ifelse(is.na(sampling_location_two), sample_size,
                              ifelse(is.na(sampling_location_three) & !is.na(sampling_location_two), sample_size/2,
                                     ifelse(is.na(sampling_location_four) & !is.na(sampling_location_three) & !is.na(sampling_location_two), sample_size/3,
                                            sample_size/4)))) %>%
  gather(key = 'location_repeat', 'location',
         'sampling_location', 
         'sampling_location_two', 
         'sampling_location_three', 
         'sampling_location_four' 
  ) %>%
  filter(!is.na(location)) %>%
  mutate(country = gsub("(.*),.*", "\\1", location)) %>%
  mutate(country = gsub("(.*),.*", "\\1", country)) %>%
  mutate(prev.cat = ifelse(seroprevalence_percentage == 0, 0, 1)) %>%
  filter(prev.cat == 1) %>%
  #filter(virus  == 'Filovirus') %>%
  #filter(methodology == 'PCR based method')
  mutate(methodology.pcr = ifelse(methodology == 'PCR based method', 'PCR based method', 'Ab')) %>%
  group_by(country, virus, methodology.pcr, prev.cat) %>%
  summarise(n = n()) %>%
  filter(prev.cat == 1)

#countries <- shapefile('data/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp')
#countries.points = fortify(countries)

#countries.data <- countries@data
#countries.data$id <- rownames(countries@data)

countries.join <- full_join(countries.data, countries.points, by = 'id')

countries.join <- countries.join %>%
  mutate(country = tolower(ADMIN)) %>%
  mutate(fill.background = 1)

countries.join.pos <- full_join(countries.join, by_country_positive_for_map, by = 'country')
countries.join.pos <- countries.join.pos %>%
  filter(!is.na(virus))

positive.map <- ggmap(map)+
  scale_x_continuous(limits = c(-200,200)) +
  scale_y_continuous(limits = c(-100,100)) +
 # geom_polygon(data=countries.points,aes(x=long,y=lat, group=group)) +
  geom_polygon(data=countries.join.pos,aes(x=long,y=lat,fill = 1, group=group)) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_continuous(low="green4", high="aquamarine1") +
  facet_grid(virus~methodology.pcr)

positive.map
ggsave("figures/positive.map.png", bg = "transparent", height = 18, width = 18)








git rm --C:\Users\r83c996\Documents\BZDEL\meta_analysis\data\map_data