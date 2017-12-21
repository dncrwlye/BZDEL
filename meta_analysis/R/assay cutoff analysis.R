setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")

library(tidyverse)
library(readxl)
library(stringi)
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
  filter(outcome != "Experimental ") %>%
  dplyr::select(title, last_name_of_first_author, virus, study_type, study_design, methodology,assay_cutoff, protein_for_elisa, species, sex, age_class, sampling_location, sampling_location_two, sample_size, seroprevalence_percentage, number_positive, single_sampling_point, sampling_date_single_time_point, start_of_sampling, end_of_sampling) %>%
  mutate(virus = ifelse(grepl('Ebola|Marburg|Sudan', virus), 'Filovirus',
                        ifelse(grepl('Henipa|Hendra|Nipah', virus), 'Henipavirus',
                               ifelse(grepl('Tioman', virus), 'Tioman', virus)))) %>%
  mutate(sampling_location = tolower(sampling_location)) %>%
  mutate(country = stri_extract_first_regex(sampling_location, '[a-z]+')) %>%
  filter(species != 'Feral Cats') %>% 
  mutate(species = tolower(species)) %>%
  mutate(species = trimws(species)) %>%
  mutate(species = stri_extract_first_regex(species, '[a-z]+ [a-z]+')) %>%
  mutate(methodology_cleaned = ifelse(grepl('PCR', methodology), 'PCR based method',
                              ifelse(grepl('ELISA|Luminex', methodology), 'non nAb based method',
                                     ifelse(grepl('VNT|SNT|Neutralizing', methodology), "nAb based method", methodology)))) %>%
  mutate(species = gsub('hipperosiderus' ,"hipposideros", species)) %>%
  mutate(species = gsub('megarops' ,"megaerops", species)) %>%
  mutate(species = gsub('schreibersi' ,"schreibersii", species)) %>%
  mutate(species = gsub('schreibersiii' ,"schreibersii", species)) %>%
  mutate(species = gsub('minopterus' ,"miniopterus", species)) %>%
  mutate(species = gsub('miniopertus' ,"miniopterus", species)) %>%
  mutate(species = gsub('myonicterus' ,"myonycteris", species)) %>%
  mutate(species = gsub('jagoli' ,"jagori", species)) %>%
  mutate(species = gsub('leschenaulti' ,"leschenaultii", species)) %>%
  mutate(species = gsub('leschenaultiii' ,"leschenaultii", species)) %>%
  mutate(species = gsub('khuli' ,"kuhlii", species)) %>%
  mutate(species = gsub('roussetus' ,"roussettus", species)) %>%
  mutate(species = gsub('lavartus' ,"larvatus", species)) %>%
  mutate(successes = ifelse(is.na(number_positive), round((1/100)*seroprevalence_percentage * sample_size,0), number_positive)) %>% 
  dplyr::select(-c(number_positive)) %>%
  mutate(seroprevalence_percentage = ifelse(is.na(seroprevalence_percentage), successes/sample_size, seroprevalence_percentage)) %>%
  filter(!(is.na(seroprevalence_percentage)))

#some papers did repeat PCR or ELISA testing, but did like PCR urine, PCR blood...those aren't independent
#im going to group by everything but those values and then take a weighted average of (sero)prevalence 

#look at the x dataframe to see the difference the following lines do 

#x<-seroprevalence%>%
#filter(title=='Large serological survey showing cocirculation of Ebola and Marburg viruses in Gabonese bat populations , and a high seroprevalence of both viruses in Rousettus aegyptiacus')

seroprevalence <- seroprevalence %>%
  group_by_at(vars(-seroprevalence_percentage, -successes, -sample_size)) %>%
  summarise(seroprevalence_percentage.dc = weighted.mean(seroprevalence_percentage, w = sample_size), sample_size.dc = mean(sample_size)) %>%
  dplyr::rename(seroprevalence_percentage = seroprevalence_percentage.dc) %>%
  dplyr::rename(sample_size = sample_size.dc) %>%
  ungroup()

rm(MetaAnalysis_Data_New_Version)

cutoff <- seroprevalence %>%
  filter(methodology_cleaned != 'PCR based method')%>%
  select(methodology_cleaned, assay_cutoff) %>%
  mutate(od = ifelse(grepl('OD', assay_cutoff), as.numeric(stri_extract_first_regex(cutoff$assay_cutoff, '.[0-9]+')), NA)) %>%
  mutate(od = ifelse(od == 3.00, NA, od)) %>% #dont count 'od of 3' from the studies that used OD 3x above the mean
  mutate(od2 = ifelse(grepl('OD', assay_cutoff), TRUE, FALSE)) %>%
  mutate(MFI = ifelse(grepl('MFI', assay_cutoff), as.numeric(stri_extract_first_regex(cutoff$assay_cutoff, '[0-9]+')), NA)) %>%
  mutate(MFI2 = ifelse(grepl('MFI', assay_cutoff), TRUE, FALSE)) %>%
  mutate(titer = ifelse(grepl('[D|d]ilution', assay_cutoff), as.numeric(stri_extract_first_regex(cutoff$assay_cutoff, '(?<=[0-9]:)[0-9]+')),
                 ifelse(grepl('[T|t]iter', assay_cutoff), as.numeric((stri_extract_first_regex(cutoff$assay_cutoff, '[0-9]+'))), NA ))) %>%
  mutate(titer2 = ifelse(grepl('[T|t]iter|[D|d]ilution', assay_cutoff), TRUE, FALSE)) %>%
  mutate(check = ifelse(od2 == FALSE & MFI2 == FALSE & titer2 == FALSE, FALSE, TRUE))

hist(cutoff$od)
hist(cutoff$MFI)
hist(log(cutoff$titer))




