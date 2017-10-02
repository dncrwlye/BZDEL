### mapping samples to raster polygons
library(raster)
library(rgdal)
library(tidyverse)
library(ggmap)
library(broom)
library(maptools)
library(rgeos)
library(plotly)

# install.packages("gpclib", type="source")
# gpclibPermit()

#ecos <- shapefile('official_teow/wwf_terr_ecos.shp')
#cols <- rainbow(length(unique(ecos$BIOME))) %>% sample(.)


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
mp
# 
# worldmap <- borders('world',colour = 'grey50',fill='grey50')
# 
# mp <- ggplot() + worldmap
# mp+geom_polygon(data=ecos,mapping(group=BIOME))



# Which polygon is our point in? ------------------------------------------
library(tidyverse)
load("~/BZDEL/Data/MetaAnalysis/seroprevalence.Rdata")

ecos <- shapefile('~/BZDEL/Data/MetaAnalysis/official_teow/wwf_terr_ecos.shp')

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

m <- do.call(bind, coordinate_box)
plot(m)


lat <- c(min(seroprevalence_x_final$south_final) ,max(seroprevalence_x_final$north_final))
lon <- c(min(seroprevalence_x_final$west_final),max(seroprevalence_x_final$east_final))

map <- get_map(location = c(lon = mean(lon), lat = mean(lat)), zoom = 2,
               maptype = "satellite", source = "google")

coordinate_box_fortified <- fortify(m)
eco.points = fortify(ecos)
ecos2 = tidy(ecos)
eco.points2 = fortify(ecos2)
### When you draw a figure, you limit lon and lat.      
ggmap(map)+
  scale_x_continuous(limits = c(min(seroprevalence_x_final$west_final),max(seroprevalence_x_final$east_final))) +
  scale_y_continuous(limits = c(min(seroprevalence_x_final$west_final),max(seroprevalence_x_final$east_final))) +
  #geom_polygon(aes(x=long, y=lat, group=group), fill='grey', size=.2,color='green', data=coordinate_box_fortified, alpha=.3) #+
  geom_polygon(data=eco.points,aes(x=long,y=lat,group=group,fill=group))

foo <- ggplot()+
  geom_polygon(data=world,aes(x=long,y=lat,group=group))+
  geom_polygon(data=eco.points,aes(x=long,y=lat,group=group,fill=group))+
  theme(legend.position = 'none')+
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(min(seroprevalence_x_final$west_final),max(seroprevalence_x_final$east_final))) +
  scale_y_continuous(limits = c(min(seroprevalence_x_final$west_final),max(seroprevalence_x_final$east_final))) +
  geom_polygon(aes(x=long, y=lat, group=group), fill='grey', size=.5,color='black', data=coordinate_box_fortified, alpha=.5) #+

ggplotly(foo)


y <- seroprevalence_x_final %>%
  dplyr::select(c(title, sampling_location, coordinate_box)) %>%
  unique()
  
