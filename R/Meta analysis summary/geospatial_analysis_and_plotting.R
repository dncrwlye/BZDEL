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

#ecos <- shapefile('official_teow/wwf_terr_ecos.shp')
#cols <- rainbow(length(unique(ecos$BIOME))) %>% sample(.)

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
load("~/BZDEL/Data/MetaAnalysis/seroprevalence.Rdata")

ecos <- shapefile('~/BZDEL/Data/MetaAnalysis/official_teow/wwf_terr_ecos.shp')

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

spdf <- SpatialPointsDataFrame(coords = xy,data=xy,
                               proj4string = CRS(proj4string(ecos)))

pps <- over(spdf,ecos[,'ECO_NAME'])

seroprevalence_x <- as.data.frame(c(seroprevalence, pps)) %>%
  select(-c(north, south, west, east, administrative_area_level_1, north_two, south_two, west_two, east_two))

#we should probably go back and use polygons over polygons, but that is proving way too hard rn

seroprevalence_x <- seroprevalence%>%
  select(-c(north, south, west, east,  north_two, south_two, west_two, east_two)) %>%
  mutate(coordinate_box = NA)

seroprevalence_x_unique <- seroprevalence_x %>%
  select(north_final, south_final, west_final, east_final) %>%
  unique()

coordinate_box <- list()


for(i in 1:nrow(seroprevalence_x_unique))
{
  if(is.na(seroprevalence_x_unique[i,'north_final'])) next 
  
  x <- as(raster::extent(as.numeric(seroprevalence_x_unique[i,3]), as.numeric(seroprevalence_x_unique[i,4]), 
                         as.numeric(seroprevalence_x_unique[i,2]), as.numeric(seroprevalence_x_unique[i,1])), "SpatialPolygons")
  
  #proj4string(x) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  proj4string(x) <- proj4string(ecos)
  coordinate_box[[i]] <- x
  
  seroprevalence_x_unique[i,'coordinate_box'] <- coordinate_box[[i]] %over% ecos[,'ECO_NAME']
  
  print(i)
}

seroprevalence_x <- seroprevalence_x %>%
  select(-c(coordinate_box))

seroprevalence_x_final <- full_join(seroprevalence_x, seroprevalence_x_unique)

y <- seroprevalence_x_final %>%
  filter(is.na(coordinate_box))

save(seroprevalence_x_final, file ="~/BZDEL/Data/MetaAnalysis/seroprevalence_ecoregions.Rdata")

#...
