#............................script that will clean prevalence and seroprevalence data...................................................
library(tidyverse)
library(readxl)
library(stringi)

MetaAnalysis_Data_New_Version <- read_excel("C:/Users/r83c996/Dropbox/bat virus meta-analysis/MetaAnalysis Data New Version.xlsx", 
#MetaAnalysis_Data_New_Version <- read_excel("~/Dropbox_gmail/Dropbox/bat virus meta-analysis/MetaAnalysis Data New Version.xlsx", 
                                                 col_types = c("text", "numeric", "text", 
                                                                         "text", "text", "text", "text", "text", 
                                                                         "text", "text", "text", "text", "text", 
                                                                         "text", "text", "text", "text", "text", 
                                                                         "text", "text", "numeric", "text", 
                                                                         "text", "text", "text", "text", "text", 
                                                                         "numeric", "text", "text", "text", 
                                                                         "text", "text", "text", "date", "date", 
                                                                         "date", "text", "text", "text", "text", 
                                                                         "text", "text"))                                      
                                            
## save old virus
MetaAnalysis_Data_New_Version$virusold=MetaAnalysis_Data_New_Version$virus

# Cleaning ----
seroprevalence <- MetaAnalysis_Data_New_Version %>%
  filter(outcome == 'Prevalence_Seroprevalence') %>%
  filter(study_type == "Observational") %>%
  dplyr::select(title, last_name_of_first_author, virus, virusold, study_type, study_design, methodology, species, sex, age_class, sampling_location, sampling_location_two, sampling_location_three, sampling_location_four, sample_size, seroprevalence_percentage, number_positive, single_sampling_point, sampling_date_single_time_point, start_of_sampling, end_of_sampling) %>%
  mutate(virus = ifelse(grepl('Ebola|Marburg|Sudan|Lloviu', virus), 'Filovirus',
                        ifelse(grepl('Henipa|Hendra|Nipah', virus), 'Henipavirus',
                               ifelse(grepl('Tioman', virus), 'Tioman', virus)))) %>%
  mutate(sampling_location       = tolower(sampling_location)) %>%
  mutate(sampling_location_two   = tolower(sampling_location_two)) %>%
  mutate(sampling_location_three = tolower(sampling_location_three)) %>%
  mutate(sampling_location_four  = tolower(sampling_location_four)) %>%
  filter(virus != 'Tioman') %>%
  #mutate(country = stri_extract_first_regex(sampling_location, '[a-z]+')) %>%
  filter(species != 'Feral Cats') %>% 
  mutate(species = tolower(species)) %>%
  mutate(species = trimws(species)) %>%
  #mutate(species = stri_extract_first_regex(species, '[a-z]+ [a-z]+')) %>%
  mutate(methodology = ifelse(grepl('PCR', methodology), 'PCR based method',
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
  mutate(species = gsub('horsfieldi' ,"horsfieldii", species)) %>%
  mutate(species = gsub('horsfieldiii' ,"horsfieldii", species)) %>%
  mutate(species = gsub('veldkampi' ,"veldkampii", species)) %>%
  mutate(species = gsub('veldkampiiii' ,"veldkampii", species)) %>%
  mutate(species = gsub('ferrum-equinum' ,"ferrumequinum", species)) %>%
  mutate(species = gsub('roussettus' ,"rousettus", species)) %>%
  mutate(number_positive = as.numeric(number_positive)) %>%
  mutate(seroprevalence_percentage = as.numeric(seroprevalence_percentage)) %>%
  mutate(sample_size = as.numeric(sample_size)) %>%
  mutate(successes = ifelse(is.na(number_positive), round((1/100)*seroprevalence_percentage * sample_size,0), number_positive)) %>% 
  dplyr::select(-c(number_positive)) %>%
  mutate(seroprevalence_percentage = ifelse(is.na(seroprevalence_percentage), successes/sample_size, seroprevalence_percentage)) %>%
  filter(!(is.na(seroprevalence_percentage)))

#some papers did repeat PCR or ELISA testing, but did like PCR urine, PCR blood...those aren't independent
#im going to group by everything but those values and then take a weighted average of (sero)prevalence 

seroprevalence.x <- seroprevalence %>%
  mutate(row.unique = paste(sample_size, title, study_design, species, sex, methodology, age_class, sampling_location, single_sampling_point, sep = ', '))

seroprevalence.y <- seroprevalence %>%
  group_by_at(vars(-seroprevalence_percentage, -successes, -sample_size)) %>%
  summarise(seroprevalence_percentage.dc = weighted.mean(seroprevalence_percentage, w = sample_size), sample_size.dc = mean(sample_size)) %>%
  dplyr::rename(seroprevalence_percentage = seroprevalence_percentage.dc) %>%
  dplyr::rename(sample_size = sample_size.dc) %>%
  ungroup() %>%
  mutate(row.unique.w = paste(sample_size, title, study_design, species, sex, methodology, age_class, sampling_location, single_sampling_point, sep = ', '))
  
seroprevalence.w <- seroprevalence.x %>%
  filter(!(row.unique %in% seroprevalence.y$row.unique.w))

seroprevalence <- seroprevalence %>%
  group_by_at(vars(-seroprevalence_percentage, -successes, -sample_size)) %>%
  summarise(seroprevalence_percentage.dc = weighted.mean(seroprevalence_percentage, w = sample_size), sample_size.dc = mean(sample_size)) %>%
  dplyr::rename(seroprevalence_percentage = seroprevalence_percentage.dc) %>%
  dplyr::rename(sample_size = sample_size.dc) %>%
  ungroup()

#...........adding on a column for becker, unfortunately many papers did multiple methods.....

seroprevalence <- seroprevalence %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  mutate(date_diff_cat = ifelse(single_sampling_point == 1 | date_diff < 365/12, '<30.4 days', 
                         ifelse(single_sampling_point == 0 | date_diff >= 365/12, '>30.4 days', NA))) %>%
  mutate(substudy_non_annual = paste(title, study_design, species, sex, methodology, age_class, sampling_location, single_sampling_point, sep = ', '))
  
explicit_longitudinal<-seroprevalence %>%
  filter(single_sampling_point == TRUE | date_diff < 365/12) %>%
  dplyr::select(title, substudy_non_annual, sampling_date_single_time_point, start_of_sampling) %>%
  unique() %>%
  group_by(substudy_non_annual, title) %>%
  summarise(counts= n()) %>%
  filter(counts >= 2) %>%
  ungroup() %>%
  #dplyr::select(title) %>%
  unique()

explicit_longitudinal.a<-seroprevalence %>%
  filter(single_sampling_point == TRUE | date_diff < 365/12) %>%
  dplyr::select(title, substudy_non_annual, sampling_date_single_time_point, start_of_sampling) %>%
  unique() %>%
  group_by(substudy_non_annual, title) %>%
  summarise(counts= n()) %>%
  filter(counts == 2| counts == 3) %>%
  ungroup() %>%
  #dplyr::select(title) %>%
  unique()

explicit_longitudinal.b<-seroprevalence %>%
  filter(single_sampling_point == TRUE | date_diff < 365/12) %>%
  dplyr::select(title, substudy_non_annual, sampling_date_single_time_point, start_of_sampling) %>%
  unique() %>%
  group_by(substudy_non_annual, title) %>%
  summarise(counts= n()) %>%
  filter(counts >= 4) %>%
  ungroup() %>%
  unique()

single_time_points_but_decent_range<-seroprevalence %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == TRUE | date_diff < 365/12) %>%
  dplyr::select(title, substudy_non_annual, sampling_date_single_time_point, start_of_sampling) %>%
  unique() %>%
  group_by(title, substudy_non_annual) %>%
  summarise(counts= n()) %>%
  filter(counts == 1) %>%
  ungroup() %>%
  #dplyr::select(title) %>%
  unique()

pooled_estimates_just_horrible<-seroprevalence %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == FALSE & date_diff >= 365/12) %>%
  dplyr::select(title, substudy_non_annual, sampling_date_single_time_point, start_of_sampling, date_diff) %>%
  unique() %>%
  group_by(title, substudy_non_annual, date_diff) %>%
  summarise(counts= n()) %>%
  ungroup() %>%
  unique()

#.....................................improving how we tell which papers fall into multiple categories.......

seroprevalence.compare <- seroprevalence %>%
  mutate(great = ifelse(substudy_non_annual %in% explicit_longitudinal.a$substudy_non_annual, TRUE, FALSE)) %>%
  mutate(good = ifelse(substudy_non_annual %in% explicit_longitudinal.b$substudy_non_annual, TRUE, FALSE)) %>%
  mutate(okay = ifelse(substudy_non_annual %in% single_time_points_but_decent_range$substudy_non_annual, TRUE, FALSE)) %>%
  mutate(bad = ifelse(substudy_non_annual %in% pooled_estimates_just_horrible$substudy_non_annual, TRUE, FALSE)) 
  
explicit_longitudinal$substudy_non_annual %in% single_time_points_but_decent_range$substudy_non_annual
explicit_longitudinal$substudy_non_annual %in% pooled_estimates_just_horrible$substudy_non_annual
single_time_points_but_decent_range$substudy_non_annual %in% pooled_estimates_just_horrible$substudy_non_annual
    
unique(seroprevalence$substudy_non_annual)

seroprevalence <- seroprevalence %>%
  mutate(sampling.strategy = ifelse(substudy_non_annual %in% explicit_longitudinal.a$substudy_non_annual, "> 4 repeated sampling events",
                                ifelse(substudy_non_annual %in% explicit_longitudinal.b$substudy_non_annual, "2-3 repeated sampling events",
                                  ifelse(substudy_non_annual %in% single_time_points_but_decent_range$substudy_non_annual, "1 sampling event (<30 days)",
                                         ifelse(substudy_non_annual %in% pooled_estimates_just_horrible$substudy_non_annual, "pooled multiple sampling events (>30 days)", "unclear sampling strategy")))))
         
save(seroprevalence, file='data/seroprevalence.Rdata')

