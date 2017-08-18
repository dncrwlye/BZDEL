# merge olival data with extra information from other papers
library(tidyverse)
library(fuzzyjoin)
#install.packages("fuzzyjoin")

associations <- read_csv("~/Desktop/BDEL/BZDEL/Data/Olival/data/associations.csv")
load("~/Desktop/BDEL/BZDEL/Data/Meta Analyses Papers/Olival_et_al_2017.Rdata")
load("~/Desktop/BDEL/BZDEL/Data/Meta Analyses Papers/Geohegan_Holmes_2015_viral_traits.Rdata")

olival_et_al_2017 <- full_join(olival_et_al_2017, associations)

olival_et_al_2017 <- olival_et_al_2017 %>%
  mutate(virus_name_dc = gsub('_[1-9]', ' ', vVirusNameCorrected)) %>%
  ################################################################ the following mutates needs to take their source material from the new variable
  mutate(virus_name_dc = gsub('-[1-9]', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('_', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('-', ' ', virus_name_dc)) %>% 
  mutate(virus_name_dc = gsub('\\s+$', '', virus_name_dc)) %>% 
  mutate(virus_name_dc = tolower(virus_name_dc)) 
  
Geohegan_Holmes <- Geohegan_Holmes %>%
  mutate(virus_name_dc = gsub('_[1-9]', ' ', Species)) %>%
  ################################################################ the following mutates needs to take their source material from the new variable
  mutate(virus_name_dc = gsub('-[1-9]', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('_', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('-', ' ', virus_name_dc)) %>% 
  mutate(virus_name_dc = tolower(virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('\\*', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('\\s+$', '', virus_name_dc)) %>% 
  #mutate(virus_name_dc = gsub('\\**', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('[0-9]{1,2}', ' ', virus_name_dc))
  

unique(Geohegan_Holmes$virus_name_dc)[1:10]
unique(olival_et_al_2017$virus_name_dc)[1:10]

#merged <- fuzzy_full_join(olival_et_al_2017[virus], Geohegan_Holmes, by=NULL,match_fun = c("virus_name_dc", "virus_name_dc"))
merged <- dplyr::full_join(as.data.frame(olival_et_al_2017[c('virus_name_dc', 'IsZoonotic')]), 
                              as.data.frame(Geohegan_Holmes[,c('virus_name_dc','Humany_to_Human')]), by=c('virus_name_dc','virus_name_dc'))

merged2 <- dplyr::full_join(as.data.frame(olival_et_al_2017), 
                           as.data.frame(Geohegan_Holmes), 
                           by=c('virus_name_dc','virus_name_dc')) %>%
          select(c(vVirusNameCorrected, virus_name_dc, Transmission_Mode)) %>%
          unique()



