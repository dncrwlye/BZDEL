library(tidyverse)
library(stringi)
library(readxl)
library(binom)
library(lubridate)
library(plotly)
library(lme4)
load("~/BZDEL/Data/MetaAnalysis/Bat_Birth_Pulse_Data_final_alternative.Rdata")
#load("~/BZDEL/Data/MetaAnalysis/seroprevalence_ecoregions.Rdata")
load("~/BZDEL/Data/MetaAnalysis/seroprevalence_ecoregions_alternative.Rdata")

#...........................step 1...... clean the data a little for ones where we have unique time points
seroprevalence_single_time_point <- seroprevalence_x_final %>%
    filter(single_sampling_point == 1) %>%
    mutate(month = round_date(sampling_date_single_time_point, unit= 'months')) %>%
    mutate(month = as.numeric(format(month, "%m"))) %>%
    mutate(year = as.numeric(format(sampling_date_single_time_point, "%Y"))) %>%
    #mutate(month = format(as.Date(sampling_date_single_time_point), "%m-%d"))) %>%
    mutate(substudy = paste(title, year, species, sep = ', ')) %>%
    group_by(title, month, virus, species, ECO_NAME,  methodology, age_class, last_name_of_first_author, year, sampling_location) %>%
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
    #unite(group, c(age_class, last_name_of_first_author, year), sep = ", ") %>%
    group_by(title, month_1, month_2, virus, species, ECO_NAME, methodology, age_class, last_name_of_first_author, year, sampling_location) %>%
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
    mutate(title = paste(title, "xc")) 
   
  
  x <- seroprevalence_unclear_time_point[,c(1:2,4:ncol(seroprevalence_unclear_time_point))]
  y <-seroprevalence_unclear_time_point[,c(1,3:ncol(seroprevalence_unclear_time_point))]
  colnames(y) <- colnames(x)
  seroprevalence_unclear_time_point <- rbind(y,x) %>%
    rename(month = month_1)
  
#colnames(seroprevalence_henipavirus_unclear_time_point) <- colnames(seroprevalence_henipavirus_single_time_point)
seroprevalence_graph <- rbind(seroprevalence_unclear_time_point, seroprevalence_single_time_point)

rm(x,y, seroprevalence_single_time_point, seroprevalence_unclear_time_point)

#.............................................making the adjustments to our time series data

# also join

Bat_Birth_Pulse_Data <- Bat_Birth_Pulse_Data %>%
  rename(ECO_NAME = coordinate_box)

seroprevalence_graph_inner_join <- inner_join(seroprevalence_graph, Bat_Birth_Pulse_Data_final)

#seroprevalence_graph_inner_join_filter.analyses <- seroprevalence_graph_inner_join %>%
#  dplyr::select(species, group,ECO_NAME, month, birth_pulse_1_quant) %>%
#  unique() %>%
#  group_by(species, group) %>%
#  summarise(sd = sd(birth_pulse_1_quant, na.rm=TRUE)) %>%
#  mutate(sd_total = sd(seroprevalence_graph_inner_join$birth_pulse_1_quant, na.rm=TRUE)) %>%
#  mutate(ratio = sd/sd_total)


seroprevalence_graph_inner_join_filter <- seroprevalence_graph_inner_join %>%
  dplyr::select(-c(birth_pulse_1_quant_flipped)) %>%
  group_by_at(vars(-ECO_NAME, birth_pulse_1_quant, birth_pulse_2_quant)) %>%
  summarise(birth_pulse_1_quant_new = mean(birth_pulse_1_quant)) %>%
  unique() %>%
  mutate(month.dc.birthpulse = month - birth_pulse_1_quant_new) 

seroprevalence_graph_inner_join_filter_alex_style <- seroprevalence_graph_inner_join_filter %>%
  #unite(substudy, c(year, species, title), sep = ", ") %>%
  mutate(substudy = paste(title, year, species,sampling_location, sep = ', ')) %>%
  group_by(substudy) %>%
  mutate(value = (seroprevalence_per - mean(seroprevalence_per))/sd(seroprevalence_per, na.rm=TRUE))
  

# ..................................... graphing ...............................
x<-seroprevalence_graph_inner_join_filter_alex_style %>%
  #filter(methodology=='PCR based method') %>%
  #mutate(month.dc.birthpulse = month - birth_pulse_1_quant_flipped) %>%
  #unite(group, c(species, group), sep = ": ") %>%
  ungroup() %>%
  ggplot(aes(x= month.dc.birthpulse, y= value, colour = substudy, group = substudy, text = paste('group', substudy))) +
  geom_point() +
  geom_line() +
  #geom_ribbon(alpha = .2, aes(x= month.dc.birthpulse, 
  #                            ymin=seroprevalence_per_lower_bound, 
  #                            ymax=seroprevalence_per_upper_bound,
  #                            fill= group)) +
  theme(legend.position="none") +
  ylab('seroprevalence') +
  ggtitle(paste(c('Filovirus', "Seroprevalence; Data from Meta Analysis")))+
  facet_grid(virus~methodology)


ggplotly(x, tooltip = 'text')




x<-seroprevalence_x_final %>%
  filter(title == 'Antibodies to Nipah or Nipah-like Viruses in Bats , China' & species == "hipposideros pomona")






seroprevalence_graph_anti_join <- anti_join(seroprevalence_graph, Bat_Birth_Pulse_Data_final)

seroprevalence_graph_anti_join <- seroprevalence_graph_anti_join %>%
  dplyr::select(-c(ECO_NAME)) %>%
  unique()# %>%
  #mutate(month.dc.birthpulse = month - birth_pulse_1_quant_flipped) %>%
  #mutate(during_birth_pulse = ifelse(month.dc.birthpulse >= 0 & month.dc.birthpulse < 3, 0, 1))



ploty_graph <- seroprevalence_graph_anti_join %>%
  #filter(methodology=='PCR based method') %>%
  mutate(month.dc.birthpulse = month) %>%
  unite(group, c(species, group), sep = ": ") %>%
  ungroup() %>%
  ggplot(aes(x= month.dc.birthpulse, y= seroprevalence_per, colour = group, group = group, text = paste('group', group))) +
  geom_point() +
  geom_line() +
  #geom_ribbon(alpha = .2, aes(x= month.dc.birthpulse, 
  #                            ymin=seroprevalence_per_lower_bound, 
  #                            ymax=seroprevalence_per_upper_bound,
  #                            fill= group)) +
  theme(legend.position="none") +
  ylab('seroprevalence') +
  ggtitle(paste(c('Filovirus', "Seroprevalence; Data from Meta Analysis")))+
  facet_grid(virus~methodology)

ggplotly(ploty_graph, tooltip = 'text')


hover_text <- py$get_figure()

str(py$x$layout)

py$x$data[[2]]$text <- c('Laurier', 'Fairmount')
py








