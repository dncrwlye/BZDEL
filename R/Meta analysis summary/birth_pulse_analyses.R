library(tidyverse)
library(stringi)
library(readxl)
library(binom)
library(lubridate)

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

MetaAnalysis_Data_New_Version <- MetaAnalysis_Data_New_Version %>%
  mutate(virus = trimws(virus)) %>%
  mutate(species = trimws(species)) 

seroprevalence <- MetaAnalysis_Data_New_Version %>%
  filter(outcome == 'Seroprevalence') %>%
  select(title, last_name_of_first_author, virus, study_type, study_design, methodology, species, sex, age_class, sampling_location, sample_size, seroprevalence_percentage, single_sampling_point, sampling_date_single_time_point, start_of_sampling, end_of_sampling) %>%
  mutate(virus = ifelse((virus == "Ebola" | 
                           virus == "Marburg" | 
                           virus == "Sudan virus" | 
                           virus == "Zaire Ebolavirus" | 
                           virus == "Reston Ebola"), "Filovirus", virus)) %>%
  mutate(global_location = ifelse(grepl('Thailand|Malaysia|China|Cambodia', sampling_location), 'East Asia', 
                                  ifelse(grepl('Brazil', sampling_location), 'South America', 
                                         ifelse(grepl('Bangladesh|India',sampling_location), 'Central Asia', 
                                                ifelse(grepl('Australia|Indonesia|Phillippines', sampling_location), "Austranesia",
                                                       ifelse(grepl('Ghana|DRC|Gabon|Annobon|Uganda|Congo|Zambia|Madagascar',sampling_location), "sub_saharan_africa", NA)))))) %>%
  mutate(sampling_location = tolower(sampling_location)) %>%
  mutate(country = stri_extract_first_regex(sampling_location, '[a-z]+')) %>%
  mutate(successes = round((1/100)*seroprevalence_percentage * sample_size,0))

species_location_table <- seroprevalence %>%
  group_by(species, country) %>%
  summarise(sum = sum(sample_size))

filovirus.seroprevalence<-seroprevalence #%>%
  #filter(virus=='Filovirus')

seroprevalence_henipavirus_single_time_point <- filovirus.seroprevalence %>%
    filter(single_sampling_point == 1) %>%
    mutate(month = round_date(sampling_date_single_time_point, unit= 'months')) %>%
    mutate(month = as.numeric(format(month, "%m"))) %>%
    #mutate(month = format(as.Date(sampling_date_single_time_point), "%m-%d"))) %>%
    group_by(month, age_class, title, species, virus, country, last_name_of_first_author, methodology, global_location, sampling_location) %>%
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
    group_by(month_1, month_2, age_class, species, virus, country, title, last_name_of_first_author, methodology, global_location, sampling_location) %>%
    summarise(sample_size = sum(sample_size), successes = sum (successes),difftime = mean(difftime))%>%
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
    mutate(title = paste(title, "xc"))
  
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
seroprevalence_graph <- seroprevalence_graph %>%
  mutate(month.dc.birthpulse = ifelse(species == 'Rousettus aegyptiacus ', month - 8.5, #according to amman and towner they give birth in august/september
                               ifelse(species == 'Eidolon helvum', month - 11.5, #this is so effing sketch, the ranges of when this bat are so drastic. I don't have estimates for Zambia (apparently there are 'too many crocodiles', ~*eye roll*~), so I'm going to try Malawi...oh god
                               ifelse(species == 'Epomops franqueti', month - 1, #no estimate for our location (gabon), so using the DRC estimate
                               ifelse(species == "Hypsignathus monstrosus", month -2, #using the gabon estimate for DRC and gabon)
                               ifelse(species == "Pteropus lylei", month - 2, #using Wacha... estimate for Thailand for both Thailand and Cambodia
                               ifelse(species == "Pteropus hypomelanus", month - 4.5, # no estimate for Malaysia or Thailand, so Im using the Phillipines estimate, probably way off      
                               ifelse((species == "Pteropus vampyrus" & country == 'malaysia')|(species == "Pteropus vampyrus" & country == 'thailand'), month - 3.4,
                               ifelse(species == "Pteropus vampyrus" & country == 'indonesia', month - 4.5,
                  
                               ifelse(species == "Myonycteris torquata", month - 2, NA))))))))))

seroprevalence_graph %>%
  filter(!is.na(month.dc.birthpulse)) %>%
  unite(group, c(age_class, species, country, title, methodology), sep = " ") %>%
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
  facet_wrap(~virus)



seroprevalence_graph %>%
    ungroup() %>%
    ggplot() +
    geom_point(aes(x= month, y= seroprevalence_per, colour = title)) +
    geom_line(lwd=1, aes(x= month, y=seroprevalence_per, colour = title)) +
    geom_ribbon(alpha = .2, aes(x= month, 
                                ymin=seroprevalence_per_lower_bound, 
                                ymax=seroprevalence_per_upper_bound,
                                fill= title)) +
    theme(legend.position="none") +
    ylab('seroprevalence') +
    ggtitle(paste(c('Filovirus', "Seroprevalence; Data from Meta Analysis"))) +
    facet_wrap(~global_location)
  
  graph

