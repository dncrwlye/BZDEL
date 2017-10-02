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
  select(c(species,location, birth_pulse_1_quant, birth_pulse_2_quant)) %>%
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
  mutate(address = NA) 

for(i in 1:nrow(Bat_Birth_Pulse_Data_unique))
{
  query <- Bat_Birth_Pulse_Data_unique$location[i] 
  rd <- geocode(query, output = 'more', source = 'google')
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

