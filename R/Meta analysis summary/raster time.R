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
                                                          "numeric"))

seroprevalence <- MetaAnalysis_Data_New_Version %>%
  filter(outcome == 'Seroprevalence') %>%
  dplyr::select(title, last_name_of_first_author, virus, study_type, study_design, methodology, species, sex, age_class, sampling_location, sampling_location_two, sample_size, seroprevalence_percentage, single_sampling_point, sampling_date_single_time_point, start_of_sampling, end_of_sampling) %>%
  mutate(virus = ifelse((virus == "Ebola" | 
                           virus == "Marburg" | 
                           virus == "Zaire Ebola"|
                           virus == "Sudan virus" | 
                           virus == "Zaire Ebolavirus" | 
                           virus == "Reston Ebola"), "Filovirus", 
                        ifelse(virus  == "Nipah" | virus == "Hendra" | virus  == "Henipavirus", "Henipavirus", virus))) %>%
  mutate(global_location = ifelse(grepl('Thailand|Malaysia|Cambodia|Phillippines', sampling_location), 'South East Asia', 
                                  ifelse(grepl('Brazil', sampling_location), 'South America', 
                                         ifelse(grepl('Bangladesh|India',sampling_location), 'Central Asia', 
                                                ifelse(grepl('Ghana|DRC|Gabon|Annobon|Uganda|Congo|Zambia|Madagascar',sampling_location), "sub_saharan_africa", NA))))) %>%
  mutate(sampling_location = tolower(sampling_location)) %>%
  mutate(country = stri_extract_first_regex(sampling_location, '[a-z]+')) %>%
  mutate(species = tolower(species)) %>%
  mutate(species = trimws(species)) %>%
  mutate(species = stri_extract_first_regex(species, '[a-z]+ [a-z]+')) %>%
  mutate(methodology = ifelse(methodology == 'PCR'| methodology == 'RT-PCR'| methodology == "RT-PCR  (urine)"| methodology == "RT-PCR (Oro-pharangyeal swab)", 'PCR based method', 
                              ifelse(methodology == "ELISA" | methodology == "ELISA + WB" | methodology == "Luminex", 'non nAb based method',
                                     ifelse(methodology == "VNT"| methodology == "SNT"| methodology == "Unclear, Presumably Neutralizing Antibodies", "nAb based method", methodology)))) %>%
  mutate(successes = round((1/100)*seroprevalence_percentage * sample_size,0)) 
 
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

seroprevalence <- full_join(seroprevalence, seroprevalence_search)
seroprevalence <- full_join(seroprevalence, seroprevalence_search_2)

seroprevalence <- seroprevalence %>%
  mutate(north_final = pmax(north, north_two, na.rm=TRUE)) %>%
  mutate(south_final = pmin(south, south_two, na.rm=TRUE)) %>%
  mutate(west_final = pmin(west, west_two, na.rm=TRUE)) %>%
  mutate(east_final = pmax(east, east_two, na.rm=TRUE)) 
  
#x<-seroprevalence %>%
#  filter(!is.na(west_two))

setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/Data/")
save(seroprevalence, file='MetaAnalysis/seroprevalence.Rdata')

#sp<-readOGR('MetaAnalysis/official_teow')
#sp <-tidy(sp)
#plot(sp)
#globe <- get_map('planet earth', zoom= 3)
#ggmap(globe)

#globe <- globe + geom_polygon(aes(x=long, y=lat, group=group), fill='grey', size=.2,color='green', data=sp, alpha=0)
#ggmap(globe)


