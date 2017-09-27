### mapping samples to raster polygons
library(raster)
library(rgdal)
library(magrittr)
library(ggmap)
library(broom)
library(maptools)
library(rgeos)
# install.packages("gpclib", type="source")
# gpclibPermit()

ecos <- shapefile('official_teow/wwf_terr_ecos.shp')
cols <- rainbow(length(unique(ecos$BIOME))) %>% sample(.)

plot(ecos)
i=0
for (biom in unique(ecos$BIOME)){
  i=i+1
  plot(ecos[ecos$BIOME==biom,],col=cols[i],add=TRUE)
}


### alternatively, using maptools & ggmap
# ecos <- readShapePoly('official_teow/wwf_terr_ecos.shp')
world <- map_data('world')

ecos2 <- readOGR('official_teow')
# eco.points = tidy(ecos2,region = 'BIOME')
epts <- tidy(ecos2)
# biomes <- tidy(ecos2,region='BIOME') #this takes too long.

ggplot()+
  geom_polygon(data=world,aes(x=long,y=lat,group=group))+
  geom_polygon(data=eco.points,aes(x=long,y=lat,group=group,fill=group))+
  theme(legend.position = 'none')+
  coord_fixed(1.3)

mp <- ggplot(eco.points)+aes(long,lat,group=group)

# 
# worldmap <- borders('world',colour = 'grey50',fill='grey50')
# 
# mp <- ggplot() + worldmap
# mp+geom_polygon(data=ecos,mapping(group=BIOME))



# Which polygon is our point in? ------------------------------------------
library(tidyverse)
load("~/Desktop/BDEL/BZDEL/Data/MetaAnalysis/seroprevalence.Rdata")

ecos <- shapefile('~/Desktop/BDEL/BZDEL/Data/MetaAnalysis/official_teow/wwf_terr_ecos.shp')

# xy <- seroprevalence[,c('north','east')]
# colnames(xy) <- c('lat','long')
# xy <- xy[(!(is.na(xy$lat)|is.na(xy$long))),] %>% as.data.frame

xy = seroprevalence %>%
  mutate(lat = ave(north_final, south_final, rm.na=TRUE)) %>%
  mutate(lon = ave(west_final, east_final, rm.na=TRUE)) %>%
  select(c(lon, lat)) %>%
  #rename(lat=north,long=east) %>% 
  #filter(!(is.na(lat))) %>%
  #filter(!(is.na(lon))) %>%
  mutate(lon = ifelse(is.na(lon), 10, lon)) %>%
  mutate(lat = ifelse(is.na(lat), 10, lat)) %>%
  as.data.frame

#we should probably go back and use polygons over polygons, but that is proving way too hard rn

spdf <- SpatialPointsDataFrame(coords = xy,data=xy,
                               proj4string = CRS(proj4string(ecos)))

pps <- over(spdf,ecos[,'ECO_NAME'])

seroprevalence_x <- as.data.frame(c(seroprevalence, pps)) %>%
  select(-c(north, south, west, east, administrative_area_level_1, north_two, south_two, west_two, east_two))

