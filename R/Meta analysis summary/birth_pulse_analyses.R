library(tidyverse)
library(stringi)
library(readxl)
library(binom)
library(lubridate)
library(plotly)
library(lme4)
load("~/BZDEL/Data/MetaAnalysis/Bat_Birth_Pulse_Data_ecoregions.Rdata")
#load("~/BZDEL/Data/MetaAnalysis/seroprevalence_ecoregions.Rdata")
load("~/BZDEL/Data/MetaAnalysis/seroprevalence_ecoregions_alternative.Rdata")

seroprevalence_single_time_point <- seroprevalence_x_final %>%
    filter(single_sampling_point == 1) %>%
    mutate(month = round_date(sampling_date_single_time_point, unit= 'months')) %>%
    mutate(month = as.numeric(format(month, "%m"))) %>%
    mutate(year = as.numeric(format(sampling_date_single_time_point, "%Y"))) %>%
    #mutate(month = format(as.Date(sampling_date_single_time_point), "%m-%d"))) %>%
    unite(group, c(title, age_class, last_name_of_first_author, year), sep = ", ") %>%
    group_by(month, virus, species, ECO_NAME,  methodology, group) %>%
    summarise(sample_size = sum(sample_size), successes = sum (successes))  %>%
    mutate(seroprevalence_per = 100 * (successes/sample_size)) %>%
    filter(!is.na(seroprevalence_per)) %>%
    mutate(seroprevalence_per_upper_bound = 100*binom.exact(successes, sample_size, conf.level = .95)$upper) %>%
    mutate(seroprevalence_per_lower_bound = 100*binom.exact(successes, sample_size,conf.level = .95)$lower) %>%
    ungroup() 
  
seroprevalence_unclear_time_point <- seroprevalence_x_final %>%
    filter(single_sampling_point == 0) %>%
    filter(!is.na(start_of_sampling)) %>%
    mutate(difftime = as.numeric(difftime(end_of_sampling, start_of_sampling))) %>%
    mutate(difftime = ifelse(difftime > 1000, 1000, difftime)) %>%
    #mutate(month = round_date(sampling_date_single_time_point, unit= 'months')) %>%
    mutate(month_1 = ifelse(difftime < 364, as.numeric(format(start_of_sampling, "%m")), 1)) %>%
    mutate(month_2 = ifelse(difftime < 364, as.numeric(format(end_of_sampling, "%m")), 12)) %>%
    mutate(year = as.numeric(format(start_of_sampling, "%Y"))) %>%
    unite(group, c(title, age_class, last_name_of_first_author, year), sep = ", ") %>%
    group_by(month_1, month_2, virus, species, ECO_NAME, methodology, group) %>%
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
    dplyr::select(-c(difftime, difftime_sqrt)) %>%
    mutate(group = paste(group, "xc"))
  
  x <- seroprevalence_unclear_time_point[,c(2:ncol(seroprevalence_unclear_time_point))]
  y <-seroprevalence_unclear_time_point[,c(1,3:ncol(seroprevalence_unclear_time_point))]
  colnames(y) <- colnames(x)
  seroprevalence_unclear_time_point <- rbind(y,x) %>%
    rename(month = month_2)
  
  
  #colnames(seroprevalence_henipavirus_unclear_time_point) <- colnames(seroprevalence_henipavirus_single_time_point)
seroprevalence_graph <- rbind(seroprevalence_unclear_time_point, seroprevalence_single_time_point)
  
#species_location_table_graph <- seroprevalence_graph %>%
#  group_by(species, country) %>%
#  summarise(sum = sum(sample_size))

#making the adjustments to our time series data


seroprevalence_graph_inner_join <- inner_join(seroprevalence_graph, Bat_Birth_Pulse_Data_final)

seroprevalence_graph_inner_join_filter <- seroprevalence_graph_inner_join %>%
  dplyr::select(-c(ECO_NAME)) %>%
  unique() %>%
  mutate(month.dc.birthpulse = month - birth_pulse_1_quant_flipped) %>%
  mutate(during_birth_pulse = ifelse(month.dc.birthpulse >= 0 & month.dc.birthpulse < 3, 0, 1))



seroprevalence_graph_inner_join %>%
  filter(methodology=='PCR based method') %>%
  mutate(month.dc.birthpulse = month - birth_pulse_1_quant_flipped) %>%
  unite(group, c(species, group), sep = ": ") %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(x= month.dc.birthpulse, y= seroprevalence_per, colour = group, group = group)) +
  geom_line(lwd=1, aes(x= month.dc.birthpulse, y=seroprevalence_per, colour = group, group = group)) +
  #geom_ribbon(alpha = .2, aes(x= month.dc.birthpulse, 
  #                            ymin=seroprevalence_per_lower_bound, 
  #                            ymax=seroprevalence_per_upper_bound,
  #                            fill= group)) +
  theme(legend.position="none") +
  ylab('seroprevalence') +
  ggtitle(paste(c('Filovirus', "Seroprevalence; Data from Meta Analysis")))+
  facet_grid(virus~methodology)

ggplotly(seroprevalence_graph, tooltip = c("x","y", "group"))




  
unique(seroprevalence_graph$species)

