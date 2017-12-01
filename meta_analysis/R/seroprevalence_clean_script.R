#............................script that will clean prevalence and seroprevalence data...................................................

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
  dplyr::select(title, last_name_of_first_author, virus, study_type, study_design, methodology, species, sex, age_class, sampling_location, sampling_location_two, sample_size, seroprevalence_percentage, number_positive, single_sampling_point, sampling_date_single_time_point, start_of_sampling, end_of_sampling) %>%
  mutate(virus = ifelse(grepl('Ebola|Marburg|Sudan', virus), 'Filovirus',
                        ifelse(grepl('Henipa|Hendra|Nipah', virus), 'Henipavirus',
                               ifelse(grepl('Tioman', virus), 'Tioman', virus)))) %>%
  mutate(sampling_location = tolower(sampling_location)) %>%
  mutate(country = stri_extract_first_regex(sampling_location, '[a-z]+')) %>%
  filter(species != 'Feral Cats') %>% 
  mutate(species = tolower(species)) %>%
  mutate(species = trimws(species)) %>%
  mutate(species = stri_extract_first_regex(species, '[a-z]+ [a-z]+')) %>%
  mutate(methodology = ifelse(grepl('PCR', methodology), 'PCR based method',
                              ifelse(grepl('ELISA|Luminex', methodology), 'non nAb based method',
                                     ifelse(grepl('VNT|SNT|Neutralizing', methodology), "nAb based method", methodology)))) %>%
  mutate(species = gsub('hipperosiderus' ,"hipposideros", species)) %>%
  mutate(species = gsub('megarops' ,"megaerops", species)) %>%
  mutate(species = gsub('schreibersi' ,"schreibersii", species)) %>%
  mutate(species = gsub('schreibersiii' ,"schreibersii", species)) %>%
  mutate(species = gsub('minopterus' ,"miniopterus", species)) %>%
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

#...........adding on a column for becker, unfortunately many papers did multiple methods.....

seroprevalence <- seroprevalence %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  mutate(date_diff_cat = ifelse(single_sampling_point == 1 | date_diff < 365/12, 'not horrible', 
                         ifelse(single_sampling_point == 0 | date_diff > 365/12, 'horrible', NA))) %>%
  mutate(substudy_non_annual = paste(title, study_design, species, sex, methodology, age_class, sampling_location, sep = ', '))
  
explicit_longitudinal<-seroprevalence %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 1 | date_diff < 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling) %>%
  unique() %>%
  group_by(substudy_non_annual) %>%
  summarise(counts= n()) %>%
  filter(counts >= 2) %>%
  ungroup() %>%
  #dplyr::select(title) %>%
  unique()

single_time_points_but_decent_range<-seroprevalence %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 1 | date_diff < 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling) %>%
  unique() %>%
  group_by(substudy_non_annual) %>%
  summarise(counts= n()) %>%
  filter(counts == 1) %>%
  ungroup() %>%
  #dplyr::select(title) %>%
  unique()

pooled_estimates_just_horrible<-seroprevalence %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 0 & date_diff >= 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling, date_diff) %>%
  unique() %>%
  group_by(substudy_non_annual, date_diff) %>%
  summarise(counts= n()) %>%
  ungroup() %>%
  unique()

# sum(pooled_estimates_just_horrible$counts) + sum(explicit_longitudinal$counts) + sum(single_time_points_but_decent_range$counts)
# nrow(seroprevalence)
# xy <- seroprevalence %>%
#   filter(!(seroprevalence$substudy_non_annual %in% explicit_longitudinal$substudy_non_annual) &
#          !(seroprevalence$substudy_non_annual %in% single_time_points_but_decent_range$substudy_non_annual) &
#          !(seroprevalence$substudy_non_annual %in% pooled_estimates_just_horrible$substudy_non_annual))
# l <- seroprevalence %>%
#   filter((seroprevalence$substudy_non_annual %in% explicit_longitudinal$substudy_non_annual))
# m <- seroprevalence %>%
#   filter((seroprevalence$substudy_non_annual %in% single_time_points_but_decent_range$substudy_non_annual))
# n <- seroprevalence %>%
#   filter((seroprevalence$substudy_non_annual %in% pooled_estimates_just_horrible$substudy_non_annual))
# l$substudy_non_annual %in% m$substudy_non_annual 
# l$substudy_non_annual %in% n$substudy_non_annual
# m$substudy_non_annual %in% n$substudy_non_annual
# l.x <- l %>%
#   filter(l$substudy_non_annual %in% n$substudy_non_annual)

seroprevalence <- seroprevalence %>%
  mutate(study_type = ifelse(substudy_non_annual %in% explicit_longitudinal$substudy_non_annual, 'explicit_longitudinal',
                      ifelse(substudy_non_annual %in% single_time_points_but_decent_range$substudy_non_annual, 'single_time_points_but_decent_range',
                      ifelse(substudy_non_annual %in% pooled_estimates_just_horrible$substudy_non_annual, 'pooled_estimates_just_horrible', NA )))) %>%
  dplyr::select(-c(substudy_non_annual, date_diff, date_diff_cat)) %>%
  unique()

save(seroprevalence, file='data/seroprevalence.Rdata')

load(file='data/seroprevalence.Rdata')
x <- seroprevalence %>%
  filter(is.na(species))
