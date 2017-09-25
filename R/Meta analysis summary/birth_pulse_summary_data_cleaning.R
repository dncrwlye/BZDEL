library(tidyverse)
library(stringi)
library(readxl)
library(binom)
library(lubridate)


Bat_Birth_Pulse_Data <- read_excel("~/Dropbox_gmail/Dropbox/bat virus meta-analysis/bat birth pulse/Bat Birth Pulse Data.xlsx")

Bat_Birth_Pulse_Data <- Bat_Birth_Pulse_Data %>%
  select(c(species,location, location_details, birth_pulse_1_quant, birth_pulse_2_quant)) %>%
  rename(country= location) %>%
  mutate(species = tolower(species)) %>%
  mutate(species = trimws(species)) %>%
  mutate(species = stri_extract_first_regex(species, '[a-z]+ [a-z]+')) %>%
  mutate(birth_pulse_1_quant = as.numeric(birth_pulse_1_quant)) %>%
  mutate(birth_pulse_1_quant_flipped = ifelse(birth_pulse_1_quant > 6, birth_pulse_1_quant -12, birth_pulse_1_quant))
  
