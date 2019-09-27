#............................script that will clean prevalence and seroprevalence data...................................................

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(tidyverse)
library(readxl)
library(stringi)

#setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis")
#MetaAnalysis_Data_New_Version <- read_excel("C:/Users/r83c996/Dropbox/bat virus meta-analysis/MetaAnalysis Data New Version.xlsx", 
# MetaAnalysis_Data_New_Version <- read_excel("~/Dropbox/bat virus meta-analysis/MetaAnalysis Data New Version.xlsx",
#                                                  col_types = c("text", "numeric", "text",
#                                                                          "text", "text", "text", "text", "text",
#                                                                          "text", "text", "text", "text", "text",
#                                                                          "text", "text", "text", "text", "text",
#                                                                          "text", "text", "numeric", "text",
#                                                                          "text", "text", "text", "text", "text",
#                                                                          "numeric", "text", "text", "text",
#                                                                          "text", "text", "text", "date", "date",
#                                                                          "date", "text", "text", "text", "text",
#                                                                          "text", "text", "text", "text","text"))
setwd("~/Dropbox (MSU projects)/Spillover postdoc/bat virus meta-analysis")
MetaAnalysis_Data_New_Version <- read_excel("MetaAnalysis Data New Version.xlsx",
                                            col_types = c("text", "numeric", "text",
                                                          "text", "text", "text", "text", "text",
                                                          "text", "text", "text", "text", "text",
                                                          "text", "text", "text", "text", "text",
                                                          "text", "text", "numeric", "text",
                                                          "text", "text", "text", "text", "text",
                                                          "numeric", "text", "text", "text",
                                                          "text", "text", "text", "date", "date",
                                                          "date", "text", "text", "text", "text",
                                                          "text", "text", "text", "text","text"),
                                            trim_ws=T)

## load in new data
new_data <- read_excel("additional_data.xlsx",
                                                     col_types = c("text", "numeric", "text",
                                                                   "text", "text", "text", "text", "text",
                                                                   "text", "text", "text", "text", "text",
                                                                   "text", "text", "text", "text", "text",
                                                                   "text", "text", "numeric", "text",
                                                                   "text", "text", "text", "text", "text",
                                                                   "numeric", "text", "text", "text",
                                                                   "text", "text", "text", "date", "date",
                                                                   "date", "text", "text", "text", "text",
                                                                   "text", "text", "text", "text","text"),
                       trim_ws=T)

## add note
MetaAnalysis_Data_New_Version$search="previous"
new_data$search="new"

## are the names the same
n1=names(MetaAnalysis_Data_New_Version)
n2=names(new_data)
n1%in%n2

## sort correctly
new_data=new_data[names(MetaAnalysis_Data_New_Version)]

## combine
MetaAnalysis_Data_New_Version=rbind.data.frame(MetaAnalysis_Data_New_Version,new_data)

## check species
table(MetaAnalysis_Data_New_Version$species)

## textclean to figure out what's happening
library(textclean)
check_text(unique(MetaAnalysis_Data_New_Version$species))

## replace non-ascii
MetaAnalysis_Data_New_Version$species=replace_non_ascii(MetaAnalysis_Data_New_Version$species)

## check species
table(MetaAnalysis_Data_New_Version$species)

## save old virus
MetaAnalysis_Data_New_Version$virus_specific=MetaAnalysis_Data_New_Version$virus

# Cleaning ----
seroprevalence <- MetaAnalysis_Data_New_Version %>%
  filter(outcome == 'Prevalence_Seroprevalence') %>%
  filter(study_type == "Observational") %>%
  dplyr::select(title, last_name_of_first_author, year_published, virus, virus_specific, study_type, study_design, 
                methodology, assay_cutoff, protein_for_elisa,
                species, sex, age_class, sampling_location, sampling_location_two, sampling_location_three, sampling_location_four, sample_size, seroprevalence_percentage, number_positive, single_sampling_point, sampling_date_single_time_point, start_of_sampling, end_of_sampling, purpose_of_study, secondary_purpose_of_study, search) %>%
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
  mutate(methodology_general = ifelse(grepl('PCR', methodology), 'PCR based method',
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
  #mutate(species = gsub('veldkampi' ,"veldkampii", species)) %>%
  #mutate(species = gsub('veldkampiiii' ,"veldkampii", species)) %>%
  mutate(species = gsub('ferrum-equinum' ,"ferrumequinum", species)) %>%
  mutate(species = gsub('roussettus' ,"rousettus", species))

## fix percentage
seroprevalence$number_positive=as.numeric(as.character(seroprevalence$number_positive))
seroprevalence$seroprevalence_percentage=as.numeric(as.character(seroprevalence$seroprevalence_percentage))
seroprevalence$sample_size=as.numeric(as.character(seroprevalence$sample_size))

## fix wonky percents
set=seroprevalence
set=set[-which(set$seroprevalence_percentage==0),]
set=set[which(set$seroprevalence_percentage<1),]
rm(set)

## remove sample size of zero
seroprevalence=seroprevalence[-which(seroprevalence$sample_size==0),]

## fix for missing
seroprevalence$seroprevalence_percentage=ifelse(is.na(seroprevalence$seroprevalence_percentage),
                                                (seroprevalence$number_positive/seroprevalence$sample_size)*100,
       seroprevalence$seroprevalence_percentage)

## remove fake sample
seroprevalence=seroprevalence[-which(seroprevalence$sample_size==9999999999),]

## check
set=seroprevalence[c("seroprevalence_percentage","sample_size","number_positive")]
set=set[order(set$number_positive),]

## Crowley code
seroprevalence <- seroprevalence %>%
  mutate(number_positive = as.numeric(number_positive)) %>%
  mutate(seroprevalence_percentage = as.numeric(seroprevalence_percentage)) %>%
  mutate(sample_size = as.numeric(sample_size)) %>%
  mutate(successes = (seroprevalence_percentage/100)*sample_size) %>%
  filter(!is.na(sample_size) & !is.na(successes))
 
rm(MetaAnalysis_Data_New_Version)

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
  mutate(sampling.strategy = ifelse(substudy_non_annual %in% explicit_longitudinal.b$substudy_non_annual, "> 4 repeated sampling events",
                                ifelse(substudy_non_annual %in% explicit_longitudinal.a$substudy_non_annual, "2-3 repeated sampling events",
                                  ifelse(substudy_non_annual %in% single_time_points_but_decent_range$substudy_non_annual, "1 sampling event (<30 days)",
                                         ifelse(substudy_non_annual %in% pooled_estimates_just_horrible$substudy_non_annual, "pooled multiple sampling events (>30 days)", "unclear sampling strategy")))))
         

rm(explicit_longitudinal,  explicit_longitudinal.a, explicit_longitudinal.b, pooled_estimates_just_horrible, seroprevalence.compare, single_time_points_but_decent_range)
#rm(list=ls())


#some papers did repeat PCR or ELISA testing, but did like PCR urine, PCR blood...those aren't independent
#additionally, there are some paper that did analyses are different substrains of Filoviruses (marburg, ebola, etc). We now have those all under 'filoviruses, need to account for those)
#im going to group by everything but those values and then take a weighted average of (sero)prevalence

#NOTE OKAY WE NEED TO COME UP WITH A WAY TO DEAL WITH THIS BETTER. RIGHT NOW I CANNOT FIGURE IT OUT


# seroprevalence <- seroprevalence %>%
#   group_by(title,
#            last_name_of_first_author,
#            virus,
#            #not virus_specific, often times these are not independent samples, but could be interesting for future analyes,
#            study_type,
#            study_design,
#            #methodology PCR in urine and PCR in blood not independent probs
#            species, 
#            sex,
#            age_class,
#            sampling_location,
#            sampling_location_two, 
#            sampling_location_three,
#            sampling_location_four,
#            #sample_size
#            #seroprevalence_percentage
#            single_sampling_point,
#            sampling_date_single_time_point,
#            start_of_sampling,
#            end_of_sampling,
#            methodology_general,
#            #successes,
#            date_diff,
#            date_diff_cat,
#            substudy_non_annual,
#            sampling.strategy,
#            purpose_of_study,
#            secondary_purpose_of_study
#            )
# 
# seroprevalence <- seroprevalence %>%
#   #group_by_at(vars(-seroprevalence_percentage, -successes, -sample_size)) %>%
#   dplyr::summarise(seroprevalence_percentage = weighted.mean(seroprevalence_percentage, w = sample_size), 
#                    successes = weighted.mean(successes, w = sample_size),
#                    sample_size = mean(sample_size)) %>%
#   ungroup() %>%
#   mutate(successes.1 = (seroprevalence_percentage/100)*sample_size)

setwd("~/Dropbox (MSU projects)/Spillover postdoc/bat virus meta-analysis")
#save(seroprevalence, file='meta_analysis/data/seroprevalence_revised.Rdata')
save(seroprevalence, file='seroprevalence_revised.Rdata')


#okay clearly messing something up
seroprevalence.weighted.means.1 %>%
  filter(successes < 9999) %>%
  ggplot() +
  geom_point(aes(x = successes, y= successes.1))


