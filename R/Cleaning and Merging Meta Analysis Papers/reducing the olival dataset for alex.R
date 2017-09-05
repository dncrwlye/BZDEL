# clean the olival data for phylofactor that will remove semi duplicate rows 
library(tidyverse)
library(fuzzyjoin)

load("~/Desktop/BDEL/BZDEL/Data/Meta Analyses Papers/Olival_et_al_2017.Rdata")

olival_et_al_2017_reduced <- olival_et_al_2017 %>%
  mutate(virus_name_dc = gsub('_[1-9]', ' ', vVirusNameCorrected)) %>%
  ################################################################ the following mutates needs to take their source material from the new variable
  mutate(virus_name_dc = gsub('-[1-9]', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub(' [0-9]', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('_', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('-', ' ', virus_name_dc)) %>% 
  mutate(virus_name_dc = gsub('\\s+$', '', virus_name_dc)) %>% 
  mutate(virus_name_dc = tolower(virus_name_dc)) %>%
  select(-c(vVirusNameCorrected)) %>%
  unique()


olival_et_al_2017_reduced <- olival_et_al_2017_reduced %>%
  group_by(virus_name_dc) %>%
  top_n(n=1, wt = vPubMedCites)

unique(olival_et_al_2017_reduced$virus_name_dc)

save(olival_et_al_2017_reduced, file='/Users/buckcrowley/Desktop/BDEL/BZDEL/Data/Meta Analyses Papers/Olival_Dataset_Reduced.Rdata')

