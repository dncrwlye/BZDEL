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

# Which polygon is our point in? ------------------------------------------
setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis")
#setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
load("data/seroprevalence_ecoregions_final.Rdata")
ecos <- shapefile('data/official_teow/wwf_terr_ecos.shp')

#we should probably go back and use polygons over polygons, but that is proving way too hard rn

seroprevalence_x <- seroprevalence %>%
  select(-c(north, south, west, east,  north_two, south_two, west_two, east_two)) %>%
  mutate(coordinate_box = NA)

seroprevalence_x_unique <- seroprevalence_x %>%
  select(north_final, south_final, west_final, east_final) %>%
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
  
  m <- as.data.frame(over(coordinate_box[[i]], ecos[,'ECO_NAME'], returnList = TRUE)) %>%
    unique() %>%
    mutate(north_final = as.numeric(seroprevalence_x_unique[i, "north_final"])) %>%
    mutate(south_final = as.numeric(seroprevalence_x_unique[i, "south_final"])) %>%  
    mutate(west_final = as.numeric(seroprevalence_x_unique[i, "west_final"])) %>%
    mutate(east_final = as.numeric(seroprevalence_x_unique[i, "east_final"])) 
  
  eco_regions_placeholder <- bind_rows(eco_regions_placeholder, m)
  print(i)
}


seroprevalence_x_final <- full_join(seroprevalence_x, eco_regions_placeholder)

y <- seroprevalence_x_final %>%
  filter(is.na(coordinate_box))

save(seroprevalence_x_final, file ="~/BZDEL/Data/MetaAnalysis/seroprevalence_ecoregions_alternative.Rdata")


load(file ="data/seroprevalence_ecoregions_alternative.Rdata")
#######################################################################################

#............plotting data ....................................


lat <- c(min(seroprevalence_geolocation_final$south_final, na.rm=TRUE) ,max(seroprevalence_geolocation_final$north_final, na.rm=TRUE))
lon <- c(min(seroprevalence_geolocation_final$west_final, na.rm=TRUE),max(seroprevalence_geolocation_final$east_final, na.rm=TRUE))

map <- get_map(location = c(lon = mean(lon), lat = mean(lat)), zoom = 1,
               maptype = "satellite", source = "google")

m <- do.call(bind, coordinate_box)
coordinate_box_fortified <- fortify(m)

eco.points = fortify(ecos)

#.............start here if you just gonna load.................
rm(list=ls())
load(file='data/map_data')
load(file ="data/seroprevalence_ecoregions_final.Rdata")

ecos.data <- ecos@data
ecos.data$id <- rownames(ecos@data)

ecos.data.filter <- ecos.data%>%
  filter(ECO_NAME %in% seroprevalence_geolocation_final$ECO_NAME) %>%
  select(ECO_NAME, ECO_NUM, id)

eco.points.join <- eco.points %>%
  mutate(fill = ifelse(id %in% ecos.data.filter$id, 1, 0))
  #filter(id.2 %in% ecos.data.filter$ECO_NUM)

eco.points.join <- eco.points %>%
  #mutate(fill = ifelse(id.2 %in% ecos.data.filter$ECO_NUM, 1, 0))
   filter(id.2 %in% ecos.data.filter$ECO_NUM)


### When you draw a figure, you limit lon and lat.      
ggmap(map)+
  scale_x_continuous(limits = c(min(seroprevalence_geolocation_final$west_final),max(seroprevalence_geolocation_final$east_final))) +
  scale_y_continuous(limits = c(min(seroprevalence_geolocation_final$west_final),max(seroprevalence_geolocation_final$east_final))) +
  #geom_polygon(aes(x=long, y=lat, group=group), fill='grey', size=.2,color='green', data=coordinate_box_fortified, alpha=.3) +
  geom_polygon(data=eco.points.join,aes(x=long,y=lat,fill = fill, group=group))

save(list = ls(), file = 'data/map_data')


head(ecos@data)



ecos <- shapefile('data/official_teow/wwf_terr_ecos.shp')
countries <- shapefile('data/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp')

load(file='data/map_data')

ggmap(map)+
  #scale_x_continuous(limits = c(min(seroprevalence_geolocation_final$west_final),max(seroprevalence_geolocation_final$east_final))) +
  #scale_y_continuous(limits = c(min(seroprevalence_geolocation_final$west_final),max(seroprevalence_geolocation_final$east_final))) +
  #geom_polygon(aes(x=long, y=lat, group=group), fill='grey', size=.2,color='green', data=coordinate_box_fortified, alpha=.3) +
  geom_polygon(data=countries.points,aes(x=long,y=lat))


countries.points = fortify(countries)

mapdata <- map_data("world")

ggmap(map)+
  geom_polygon(data=mapdata, aes(x=long,y=lat, group=group, fill=region), alpha=0.5) + 
  guides(fill=FALSE)

