# merge olival data with extra information from other papers
library(tidyverse)
library(fuzzyjoin)

Olival_Dataset_with_phyloBin_Variable <- read_excel("~/Desktop/BDEL/BZDEL/Data/Olival/data/Olival_Dataset_with_phyloBin_Variable.xlsx")

associations <- read_csv("~/Desktop/BDEL/BZDEL/Data/Olival/data/associations.csv")
#load("~/Desktop/BDEL/BZDEL/Data/Meta Analyses Papers/Olival_et_al_2017.Rdata")
load("~/Desktop/BDEL/BZDEL/Data/Meta Analyses Papers/Geohegan_Holmes_2015_viral_traits.Rdata")
load("~/Desktop/BDEL/BZDEL/Data/Meta Analyses Papers/Plourde_plos1_2017.rdata")
load("~/Desktop/BDEL/BZDEL/Data/Meta Analyses Papers/jones_et_al.Rdata")

olival_et_al_2017 <- full_join(Olival_Dataset_with_phyloBin_Variable[,2:ncol(Olival_Dataset_with_phyloBin_Variable)], associations)

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
  
Plourde_plos1_2017 <- Plourde_plos1_2017 %>%
  filter(type=='virus') %>%
  mutate(virus_name_dc = gsub('_[1-9]', ' ', pathogen)) %>%
  ################################################################ the following mutates needs to take their source material from the new variable
  mutate(virus_name_dc = gsub('-[1-9]', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('_', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('-', ' ', virus_name_dc)) %>% 
  mutate(virus_name_dc = tolower(virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('\\*', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('\\s+$', '', virus_name_dc)) %>% 
  #mutate(virus_name_dc = gsub('\\**', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('[0-9]{1,2}', ' ', virus_name_dc))

jones_et_al <- jones_et_al %>%
  filter(PathType=='virus') %>%
  mutate(virus_name_dc = gsub('_[1-9]', ' ', `Pathogen responsible for each EID event`)) %>%
  ################################################################ the following mutates needs to take their source material from the new variable
  mutate(virus_name_dc = gsub('-[1-9]', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('_', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('-', ' ', virus_name_dc)) %>% 
  mutate(virus_name_dc = tolower(virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('\\*', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('\\s+$', '', virus_name_dc)) %>% 
  #mutate(virus_name_dc = gsub('\\**', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('[0-9]{1,2}', ' ', virus_name_dc))


#merged <- dplyr::full_join(as.data.frame(olival_et_al_2017[c('virus_name_dc', 'IsZoonotic')]), 
#                              as.data.frame(Geohegan_Holmes[,c('virus_name_dc','Humany_to_Human')]), by=c('virus_name_dc','virus_name_dc'))

olival_et_al_2017_appended <- dplyr::full_join(as.data.frame(olival_et_al_2017), 
                           as.data.frame(Geohegan_Holmes), 
                           by=c('virus_name_dc','virus_name_dc')) %>%
          select(c(phyloBin, vVirusNameCorrected,
                   vOrder, vFamily, vSubfamily, vGenus, ReverseZoonoses, 
                   vGenomeMinLength, vGenomeMaxLength, vGenomeAveLength, 
                   vCytoReplicTF, vSegmentedTF, vVectorYNna, vSSoDS, vDNAoRNA, 
                   vEnvelope, IsZoonotic, IsZoonotic.stringent, 
                   hHostNameFinal, WildDomInReference, DetectionMethod, DetectionQuality,
                   virus_name_dc, 
                   Transmission_Mode, 
                   Mortality_Rate_PCT, 
                   Envelope_Status, 
                   Recombination_Frequency, 
                   Humany_to_Human 
                   )) %>%
          unique()

#100*(1 - sum(is.na(merged2$Transmission_Mode))/nrow(merged2))

olival_et_al_2017_appended <- dplyr::full_join(as.data.frame(olival_et_al_2017_appended), 
                            as.data.frame(Plourde_plos1_2017), 
                            by=c('virus_name_dc','virus_name_dc')) %>%

  select(c(phyloBin, vVirusNameCorrected,
             vOrder, vFamily, vSubfamily, vGenus, ReverseZoonoses, 
             vGenomeMinLength, vGenomeMaxLength, vGenomeAveLength, 
             vCytoReplicTF, vSegmentedTF, vVectorYNna, vSSoDS, vDNAoRNA, 
             vEnvelope, IsZoonotic, IsZoonotic.stringent, 
             hHostNameFinal, WildDomInReference, DetectionMethod, DetectionQuality,
    virus_name_dc, 
    Transmission_Mode, 
    Mortality_Rate_PCT, 
    Envelope_Status, 
    Recombination_Frequency, 
    Humany_to_Human, 
    vGenomeMaxLength, 
    IsZoonotic,
    target,
    transmission,
    vector.id,
    region,
    reservoir)) %>%
  unique()

#100*(1 - sum(is.na(merged3$region))/nrow(merged3))

olival_et_al_2017_appended <- dplyr::full_join(as.data.frame(olival_et_al_2017_appended), 
                            as.data.frame(jones_et_al), 
                            by=c('virus_name_dc','virus_name_dc')) %>%
  select(c(phyloBin, vVirusNameCorrected,
           vOrder, vFamily, vSubfamily, vGenus, ReverseZoonoses, 
           vGenomeMinLength, vGenomeMaxLength, vGenomeAveLength, 
           vCytoReplicTF, vSegmentedTF, vVectorYNna, vSSoDS, vDNAoRNA, 
           vEnvelope, IsZoonotic, IsZoonotic.stringent, 
           hHostNameFinal, WildDomInReference, DetectionMethod, DetectionQuality, 
    virus_name_dc, 
    Transmission_Mode, 
    Mortality_Rate_PCT, 
    Envelope_Status, 
    Recombination_Frequency, 
    Humany_to_Human, 
    vGenomeMaxLength, 
    IsZoonotic,
    target,
    transmission,
    vector.id,
    region,
    reservoir,
    `TranMode Driver`)) %>%
  unique()

#100*(1 - sum(is.na(merged4$`TranMode Driver`))/nrow(merged4))


#bringing in becker's human to human work 

olival_et_al_2017_appended<-olival_et_al_2017_appended %>%
  select(-c(Humany_to_Human))

olival_et_al_2017_appended <- full_join(olival_et_al_2017_appended, human_to_human_zoonosis_data_DB_edit[,7:10], by = c('virus_name_dc'))





save(olival_et_al_2017_appended, file='/Users/buckcrowley/Desktop/BDEL/BZDEL/Data/Meta Analyses Papers/Olival_Appended_Dataset.Rdata')




