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

Bat_Birth_Pulse_Data <- read_excel("~/Dropbox_gmail/Dropbox/bat virus meta-analysis/bat birth pulse/Bat Birth Pulse Data.xlsx")

Bat_Birth_Pulse_Data <- Bat_Birth_Pulse_Data %>%
  mutate(species = tolower(species)) %>%
  mutate(species = trimws(species)) %>%
  mutate(species = gsub('epomorphus' ,"epomophorus", species)) %>%
  #mutate(species = gsub('hipperosiderus' ,"hipposideros", species)) %>%
  #mutate(species = gsub('megarops' ,"megaerops", species)) %>%
  mutate(species = gsub('schreibersi' ,"schreibersii", species)) %>%
  mutate(species = gsub('schreibersiii' ,"schreibersii", species)) %>%
  mutate(species = gsub('minopterus' ,"miniopterus", species)) %>%
  mutate(species = gsub('myonicterus' ,"myonycteris", species)) %>%
  mutate(species = gsub('jagoli' ,"jagori", species)) %>%
  mutate(species = gsub('leschenaulti' ,"leschenaultii", species)) %>%
  mutate(species = gsub('leschenaultiii' ,"leschenaultii", species)) %>%
  mutate(species = gsub('khuli' ,"kuhlii", species)) %>%
  mutate(species = gsub('roussetus' ,"roussettus", species)) %>%
  mutate(species = gsub('lavartus' ,"larvatus", species)) 

write.csv(Bat_Birth_Pulse_Data, "~/Dropbox_gmail/Dropbox/bat virus meta-analysis/bat birth pulse/Bat Birth Pulse Data_final.csv")

Bat_Birth_Pulse_Data <- Bat_Birth_Pulse_Data %>%
  dplyr::select(c(species,location, birth_pulse_1_quant, birth_pulse_2_quant)) %>%
  #rename(country= location) %>%
  mutate(location = tolower(location)) %>%
  mutate(species = tolower(species)) %>%
  mutate(species = trimws(species)) %>%
  mutate(species = stri_extract_first_regex(species, '[a-z]+ [a-z]+')) %>%
  mutate(birth_pulse_1_quant = as.numeric(birth_pulse_1_quant)) %>%
  mutate(birth_pulse_1_quant_flipped = ifelse(birth_pulse_1_quant > 6, birth_pulse_1_quant -12, birth_pulse_1_quant))
  
Bat_Birth_Pulse_Data_unique <- as.data.frame(unique(Bat_Birth_Pulse_Data$location)) %>%
  rename(location = `unique(Bat_Birth_Pulse_Data$location)`) %>%
  filter(!is.na(location)) %>%
  mutate(location = as.character(location)) %>%
  mutate(north=NA) %>%
  mutate(south=NA) %>%
  mutate(west=NA) %>%
  mutate(east=NA) %>%
  mutate(address = NA) %>%
  filter(location != 'na')


for(i in 1:nrow(Bat_Birth_Pulse_Data_unique))
{
  query <- Bat_Birth_Pulse_Data_unique$location[i] 
  rd <- geocode(query, output = 'more', source = 'google', force= TRUE, override_limit = TRUE)
  if(is.atomic(rd)==FALSE)
  {
    if (is.na(rd$lon) == FALSE) 
    {
      Bat_Birth_Pulse_Data_unique$north[i] <- rd$north
      Bat_Birth_Pulse_Data_unique$south[i] <- rd$south
      Bat_Birth_Pulse_Data_unique$east[i]  <- rd$east
      Bat_Birth_Pulse_Data_unique$west[i]  <- rd$west
      Bat_Birth_Pulse_Data_unique$address[i]  <- rd$address
    }
  }
}

Bat_Birth_Pulse_Data <- full_join(Bat_Birth_Pulse_Data, Bat_Birth_Pulse_Data_unique)

save(Bat_Birth_Pulse_Data, file = "~/Desktop/BDEL/BZDEL/meta_analysis/data/bat_birth_pulse_geocode.Rdata")


#..............moving onto the ecoregions analyses.................................

setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
ecos <- shapefile('data/official_teow/wwf_terr_ecos.shp')

#we should probably go back and use polygons over polygons, but that is proving way too hard rn

Bat_Birth_Pulse_Data_unique <- Bat_Birth_Pulse_Data %>%
  dplyr::select(c(north, south,east, west))%>%
  unique()

eco_regions_placeholder <- as.data.frame(matrix(ncol=5))
colnames(eco_regions_placeholder) <- (c("ECO_NAME","north","south","west","east"))
coordinate_box <- list()

for(i in 1:nrow(Bat_Birth_Pulse_Data_unique))
{
  if(is.na(Bat_Birth_Pulse_Data_unique[i,'north'])) next 
  
  x <- as(raster::extent(as.numeric(Bat_Birth_Pulse_Data_unique[i,3]), as.numeric(Bat_Birth_Pulse_Data_unique[i,4]), 
                         as.numeric(Bat_Birth_Pulse_Data_unique[i,2]), as.numeric(Bat_Birth_Pulse_Data_unique[i,1])), "SpatialPolygons")
  
  #proj4string(x) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  proj4string(x) <- proj4string(ecos)
  coordinate_box[[i]] <- x
  
  m <- as.data.frame(over(coordinate_box[[i]], ecos[,'ECO_NAME'], returnList = TRUE)) %>%
    unique() %>%
    mutate(north = as.numeric(Bat_Birth_Pulse_Data_unique[i, "north"])) %>%
    mutate(south = as.numeric(Bat_Birth_Pulse_Data_unique[i, "south"])) %>%  
    mutate(west = as.numeric(Bat_Birth_Pulse_Data_unique[i, "west"])) %>%
    mutate(east = as.numeric(Bat_Birth_Pulse_Data_unique[i, "east"])) 
  
  eco_regions_placeholder <- bind_rows(eco_regions_placeholder, m)
  print(i)
}

Bat_Birth_Pulse_Data_final <- full_join(Bat_Birth_Pulse_Data, eco_regions_placeholder)

save(Bat_Birth_Pulse_Data_final, file ="data/Bat_Birth_Pulse_Data_final_alternative.Rdata")

# 
# coordinate_box <- list()
# 
# for(i in 1:nrow(Bat_Birth_Pulse_Data_x_unique))
# {
#   if(is.na(Bat_Birth_Pulse_Data_x_unique[i,'north'])) next 
#   
#   x <- as(raster::extent(as.numeric(Bat_Birth_Pulse_Data_x_unique[i,3]), as.numeric(Bat_Birth_Pulse_Data_x_unique[i,4]), 
#                          as.numeric(Bat_Birth_Pulse_Data_x_unique[i,2]), as.numeric(Bat_Birth_Pulse_Data_x_unique[i,1])), "SpatialPolygons")
#   
#   #proj4string(x) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#   proj4string(x) <- proj4string(ecos)
#   coordinate_box[[i]] <- x
#   
#   Bat_Birth_Pulse_Data_x_unique[i,'coordinate_box'] <- coordinate_box[[i]] %over% ecos[,'ECO_NAME']
#   
#   print(i)
# }
# 
# Bat_Birth_Pulse_Data <- full_join(Bat_Birth_Pulse_Data, Bat_Birth_Pulse_Data_x_unique)
# 
# save(Bat_Birth_Pulse_Data, file ="data/Bat_Birth_Pulse_Data_ecoregions.Rdata")
# 















