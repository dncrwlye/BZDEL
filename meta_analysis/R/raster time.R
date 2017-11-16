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

#raster fun

MetaAnalysis_Data_New_Version <- read_excel("~/Dropbox_gmail/Dropbox/bat virus meta-analysis/MetaAnalysis Data New Version.xlsx", 
                                            col_types = c("text", "numeric", "text", 
                                                          "text", "text", "text", "text", "text", 
                                                          "text", "text", "text", "text", "text", 
                                                          "text", "text", "text", "text", "numeric", "numeric", 
                                                          "numeric", "numeric", "numeric", 
                                                          "numeric", "numeric", "numeric", 
                                                          "numeric", "text", "numeric", "numeric", 
                                                          "numeric", "numeric", "date", "date", 
                                                          "date", "text", "text", "text", "text", 
                                                          "numeric", "text"))

seroprevalence <- MetaAnalysis_Data_New_Version %>%
  filter(outcome == 'Seroprevalence') %>%
  dplyr::select(title, last_name_of_first_author, virus, study_type, study_design, methodology, species, sex, age_class, sampling_location, sampling_location_two, sample_size, seroprevalence_percentage, single_sampling_point, sampling_date_single_time_point, start_of_sampling, end_of_sampling) %>%
  mutate(virus = ifelse(grepl('Ebola|Marburg|Sudan', virus), 'Filovirus',
                 ifelse(grepl('Henipa|Hendra|Nipah', virus), 'Henipavirus',
                 ifelse(grepl('Tioman', virus), 'Tioman', virus)))) %>%
  mutate(sampling_location = tolower(sampling_location)) %>%
  mutate(country = stri_extract_first_regex(sampling_location, '[a-z]+')) %>%
  mutate(species = tolower(species)) %>%
  mutate(species = trimws(species)) %>%
  mutate(species = stri_extract_first_regex(species, '[a-z]+ [a-z]+')) %>%
  mutate(methodology = ifelse(grepl('PCR', methodology), 'PCR based method',
                          ifelse(grepl('ELISA|Luminex', methodology), 'non nAb based method',
                          ifelse(grepl('VNT|SNT|Neutralizing', methodology), "nAb based method", methodology)))) %>%
  #mutate(methodology_old = ifelse(methodology == 'PCR'| methodology == 'RT-PCR'| methodology == "RT-PCR  (urine)"| methodology == "RT-PCR (Oro-pharangyeal swab)", 'PCR based method', 
  #                            ifelse(methodology == "ELISA" | methodology == "ELISA + WB" | methodology == "Luminex", 'non nAb based method',
  #                                   ifelse(methodology == "VNT"| methodology == "SNT"| methodology == "Unclear, Presumably Neutralizing Antibodies", "nAb based method", methodology)))) %>%
  mutate(successes = round((1/100)*seroprevalence_percentage * sample_size,0)) 
 
table(seroprevalence$virus, seroprevalence$virus.dc)


seroprevalence_search <- as.data.frame(unique(seroprevalence$sampling_location)) %>%
  rename(sampling_location = `unique(seroprevalence$sampling_location)`) %>%
  mutate(sampling_location = as.character(sampling_location)) %>%
  mutate(north=NA) %>%
  mutate(south=NA) %>%
  mutate(west=NA) %>%
  mutate(east=NA) %>%
  mutate(address = NA) 

seroprevalence <- seroprevalence %>% unique()


library(googleway)

key <- "AIzaSyD40FLStnvp085UB-FKXbyONuxV3ke4umY"

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

seroprevalence <- full_join(seroprevalence, seroprevalence_search, by=c("sampling_location"))
seroprevalence <- full_join(seroprevalence, seroprevalence_search_2, by=c('sampling_location_two'))

seroprevalence <- seroprevalence %>%
  mutate(north_final = pmax(north, north_two, na.rm=TRUE)) %>%
  mutate(south_final = pmin(south, south_two, na.rm=TRUE)) %>%
  mutate(west_final = pmin(west, west_two, na.rm=TRUE)) %>%
  mutate(east_final = pmax(east, east_two, na.rm=TRUE)) 
  
x<-seroprevalence %>%
  filter(is.na(west_final))
#geocode('ghana, tanoboase', output = 'more', source = 'google')

setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
save(seroprevalence, file='data/seroprevalence.Rdata')

#..............moving onto the ecoregions analyses.................................

ecos <- shapefile('data/official_teow/wwf_terr_ecos.shp')
#ecos <- shapefile('~/Desktop/BZDEL/meta_analysis/data/official_teow/wwf_terr_ecos.shp')

#we should probably go back and use polygons over polygons, but that is proving way too hard rn

seroprevalence_x <- seroprevalence%>%
  dplyr::select(-c(north, south, west, east,  north_two, south_two, west_two, east_two)) %>%
  mutate(coordinate_box = NA)

seroprevalence_x_unique <- seroprevalence_x %>%
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
  filter(is.na(ECO_NAME))

save(seroprevalence_x_final, file ="data/seroprevalence_ecoregions_alternative.Rdata")


seroprevalence_x_final <- seroprevalence_x_final %>%
  dplyr::select(-c(ECO_NAME)) %>%
  unique()
#######################################################################################

#............plotting data ....................................
# 
# lat <- c(min(seroprevalence_x_final$south_final, na.rm=TRUE) ,max(seroprevalence_x_final$north_final, na.rm=TRUE))
# lon <- c(min(seroprevalence_x_final$west_final, na.rm=TRUE),max(seroprevalence_x_final$east_final, na.rm=TRUE))
# 
# map <- get_map(location = c(lon = mean(lon), lat = mean(lat)), zoom = 2,
#                maptype = "satellite", source = "google")
# 
# m <- do.call(bind, coordinate_box)
# coordinate_box_fortified <- fortify(m)
# 
# eco.points = fortify(ecos)
# 
# 
# ### When you draw a figure, you limit lon and lat.      
# ggmap(map)+
#   scale_x_continuous(limits = c(min(seroprevalence_x_final$west_final),max(seroprevalence_x_final$east_final))) +
#   scale_y_continuous(limits = c(min(seroprevalence_x_final$west_final),max(seroprevalence_x_final$east_final))) +
#   geom_polygon(aes(x=long, y=lat, group=group), fill='grey', size=.2,color='green', data=coordinate_box_fortified, alpha=.3) +
#   geom_polygon(data=eco.points,aes(x=long,y=lat,group=group,fill=group))






#sp<-readOGR('MetaAnalysis/official_teow')
#sp <-tidy(sp)
#plot(sp)
#globe <- get_map('planet earth', zoom= 3)
#ggmap(globe)

#globe <- globe + geom_polygon(aes(x=long, y=lat, group=group), fill='grey', size=.2,color='green', data=sp, alpha=0)
#ggmap(globe)


