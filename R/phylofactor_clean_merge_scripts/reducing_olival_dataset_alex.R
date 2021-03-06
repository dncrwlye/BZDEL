# clean the olival data for phylofactor that will remove semi duplicate rows 
# part 2: bring in Becker's pathogen release data 

library(tidyverse)
library(fuzzyjoin)

load("~/Desktop/BDEL/BZDEL/Data/phylofactor datasets/Olival_et_al_2017.Rdata")

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
  top_n(n=1, wt = vPubMedCites) %>%
  unique()

#save(olival_et_al_2017_reduced, file='/Users/buckcrowley/Desktop/BDEL/BZDEL/Data/Meta Analyses Papers/Olival_Dataset_Reduced.Rdata')

################################
#pathogen_release_data_DB_edit_vector_done <- read_csv("~/Desktop/BDEL/BZDEL/Data/Olival/pathogen release data_DB edit_vector done.csv")

pathogen_release_data_done <- read_csv("~/Desktop/BDEL/BZDEL/Data/Olival/pathogen release data_done.csv")

colnames(pathogen_release_data_done) <- c(colnames(pathogen_release_data_done[,1:3]),paste("becker", colnames(pathogen_release_data_done[,4:11]), sep = "_"))

pathogen_release_data_done <- pathogen_release_data_done %>%
  select(c(virus_name_dc, becker_pr_vector, becker_pr_excrete, becker_pr_slaughter, becker_pr_source, becker_pr_notes))

olival_et_al_2017_reduced_becker_joined <- full_join(olival_et_al_2017_reduced, pathogen_release_data_done)

check <- olival_et_al_2017_reduced_becker_joined %>%
  filter(IsZoonotic == 1) %>%
  filter(is.na(becker_pr_source)) 

save(olival_et_al_2017_reduced_becker_joined, file='/Users/buckcrowley/Desktop/BDEL/BZDEL/Data/phylofactor datasets/olival_et_al_2017_reduced_becker_joined.Rdata')

