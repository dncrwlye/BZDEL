library(tidyverse)
library(stringi)
library(readxl)
library(binom)
library(lubridate)
library(plotly)
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

species_location_table <- seroprevalence %>%
  #filter(virus == 'Henipavirus') %>%
  group_by(species, country) %>%
  summarise(sum = sum(sample_size))# %>%
  #mutate(search_term = paste("(", species, " AND ", country, ")", sep = ""))

#search_term <- paste(species_location_table$search_term[1:178], collapse = " OR ")

#write.csv(species_location_table, '/Users/buckcrowley/Desktop/species_location_table.csv')

filovirus.seroprevalence<-seroprevalence #%>%
  filter(virus=='Henipavirus')

seroprevalence_henipavirus_single_time_point <- filovirus.seroprevalence %>%
    filter(single_sampling_point == 1) %>%
    mutate(month = round_date(sampling_date_single_time_point, unit= 'months')) %>%
    mutate(month = as.numeric(format(month, "%m"))) %>%
    mutate(year = as.numeric(format(sampling_date_single_time_point, "%Y"))) %>%
    #mutate(month = format(as.Date(sampling_date_single_time_point), "%m-%d"))) %>%
    unite(group, c(title, age_class, last_name_of_first_author, year), sep = ", ") %>%
    group_by(month, virus, species, country,  methodology, group) %>%
    summarise(sample_size = sum(sample_size), successes = sum (successes))  %>%
    mutate(seroprevalence_per = 100 * (successes/sample_size)) %>%
    filter(!is.na(seroprevalence_per)) %>%
    mutate(seroprevalence_per_upper_bound = 100*binom.exact(successes, sample_size, conf.level = .95)$upper) %>%
    mutate(seroprevalence_per_lower_bound = 100*binom.exact(successes, sample_size,conf.level = .95)$lower) %>%
    ungroup() 
  
seroprevalence_henipavirus_unclear_time_point <- filovirus.seroprevalence %>%
    filter(single_sampling_point == 0) %>%
    filter(!is.na(start_of_sampling)) %>%
    mutate(difftime = as.numeric(difftime(end_of_sampling, start_of_sampling))) %>%
    mutate(difftime = ifelse(difftime > 1000, 1000, difftime)) %>%
    #mutate(month = round_date(sampling_date_single_time_point, unit= 'months')) %>%
    mutate(month_1 = ifelse(difftime < 364, as.numeric(format(start_of_sampling, "%m")), 1)) %>%
    mutate(month_2 = ifelse(difftime < 364, as.numeric(format(end_of_sampling, "%m")), 12)) %>%
    mutate(year = as.numeric(format(start_of_sampling, "%Y"))) %>%
    unite(group, c(title, age_class, last_name_of_first_author, year), sep = ", ") %>%
    group_by(month_1, month_2, virus, species, country, methodology, group) %>%
    #group_by(month_1, month_2, age_class, species, virus, country, title, last_name_of_first_author, methodology, global_location, sampling_location) %>%
    summarise(sample_size = sum(sample_size), successes = sum (successes),difftime = mean(difftime)) %>%
    mutate(difftime_sqrt = sqrt(difftime)) %>% #going to divide any long sampling period by the square root of the sampling period
    mutate(seroprevalence_per = 100 * (successes/sample_size)) %>%
    filter(!is.na(seroprevalence_per)) %>%
    mutate(seroprevalence_per_upper_bound = 100*binom.exact(successes, sample_size, conf.level = .95)$upper) %>%
    mutate(seroprevalence_per_lower_bound = 100*binom.exact(successes, sample_size,conf.level = .95)$lower) %>%
    mutate(seroprevalence_per = ifelse(difftime > 32, seroprevalence_per/difftime_sqrt, seroprevalence_per)) %>%
    mutate(seroprevalence_per_upper_bound = ifelse(difftime > 32, seroprevalence_per_upper_bound/difftime_sqrt, seroprevalence_per_upper_bound)) %>%
    mutate(seroprevalence_per_lower_bound = ifelse(difftime > 32, seroprevalence_per_lower_bound/difftime_sqrt, seroprevalence_per_lower_bound)) %>%
    ungroup() %>%
    select(-c(difftime, difftime_sqrt)) %>%
    mutate(group = paste(group, "xc"))
  
  x <- seroprevalence_henipavirus_unclear_time_point[,c(2:ncol(seroprevalence_henipavirus_unclear_time_point))]
  y <-seroprevalence_henipavirus_unclear_time_point[,c(1,3:ncol(seroprevalence_henipavirus_unclear_time_point))]
  colnames(y) <- colnames(x)
  seroprevalence_henipavirus_unclear_time_point <- rbind(y,x) %>%
    rename(month = month_2)
  
  
  #colnames(seroprevalence_henipavirus_unclear_time_point) <- colnames(seroprevalence_henipavirus_single_time_point)
seroprevalence_graph <- rbind(seroprevalence_henipavirus_unclear_time_point, seroprevalence_henipavirus_single_time_point)
  
species_location_table_graph <- seroprevalence_graph %>%
  group_by(species, country) %>%
  summarise(sum = sum(sample_size))

#making the adjustments to our time series data

seroprevalence_graph_inner_join <- inner_join(seroprevalence_graph, Bat_Birth_Pulse_Data, by = c('species','country'))

seroprevalence_graph_full_join <- full_join(seroprevalence_graph, Bat_Birth_Pulse_Data, by = c('species','country'))

x<-seroprevalence_graph_full_join %>%
  select(c(species, country, sample_size, birth_pulse_1_quant)) %>%
  group_by(species, country,birth_pulse_1_quant) %>%
  summarise(sample_size = sum(sample_size))
  
#write.csv(x, '/Users/buckcrowley/Desktop/species_location_table.csv')


plot<- seroprevalence_graph_inner_join %>%
  mutate(month.dc.birthpulse = month - birth_pulse_1_quant_flipped) %>%
  unite(group, c(species, group), sep = ": ") %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(x= month.dc.birthpulse, y= seroprevalence_per, colour = group, group = group)) +
  geom_line(lwd=1, aes(x= month.dc.birthpulse, y=seroprevalence_per, colour = group, group = group)) +
  geom_ribbon(alpha = .2, aes(x= month.dc.birthpulse, 
                              ymin=seroprevalence_per_lower_bound, 
                              ymax=seroprevalence_per_upper_bound,
                              fill= group)) +
  theme(legend.position="none") +
  ylab('seroprevalence') +
  ggtitle(paste(c('Filovirus', "Seroprevalence; Data from Meta Analysis")))+
  facet_grid(methodology ~ virus )

ggplotly(plot, tooltip = c("x","y", "group"))

  


