library(tidyverse)
library(stringi)
library(readxl)
library(binom)
library(lubridate)
library(plotly)
library(ggmap)
#raster fun

MetaAnalysis_Data_New_Version <- read_excel("~/Dropbox_gmail/Dropbox/bat virus meta-analysis/MetaAnalysis Data New Version.xlsx", 
                                            col_types = c("text", "numeric", "text", 
                                                          "text", "text", "text", "text", "text", 
                                                          "text", "text", "text", "text", "text", 
                                                          "text", "text", "numeric", "numeric", 
                                                          "numeric", "numeric", "numeric", 
                                                          "numeric", "numeric", "numeric", 
                                                          "numeric", "text", "numeric", "numeric", 
                                                          "numeric", "numeric", "date", "date", 
                                                          "date", "text", "text", "text", "text", 
                                                          "numeric"))

seroprevalence <- MetaAnalysis_Data_New_Version %>%
  filter(outcome == 'Seroprevalence') %>%
  select(title, last_name_of_first_author, virus, study_type, study_design, methodology, species, sex, age_class, sampling_location, sample_size, seroprevalence_percentage, single_sampling_point, sampling_date_single_time_point, start_of_sampling, end_of_sampling) %>%
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


