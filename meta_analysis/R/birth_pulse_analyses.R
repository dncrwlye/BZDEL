library(tidyverse)
library(stringi)
library(readxl)
library(binom)
library(lubridate)
library(plotly)
library(lme4)
library(dplyr)
library(mgcv)
library(visreg)

load("~/BZDEL/Data/MetaAnalysis/Bat_Birth_Pulse_Data_final_alternative.Rdata")
load("~/Desktop/BDEL/BZDEL/Data/MetaAnalysis/Bat_Birth_Pulse_Data_final_alternative.Rdata")
#load("~/BZDEL/Data/MetaAnalysis/seroprevalence_ecoregions.Rdata")
load("~/BZDEL/Data/MetaAnalysis/seroprevalence_ecoregions_alternative.Rdata")
load("~/Desktop/BDEL/BZDEL/Data/MetaAnalysis/seroprevalence_ecoregions_alternative.Rdata")

#...........................step 1...... clean the data a little for ones where we have unique time points
seroprevalence_single_time_point <- seroprevalence_x_final %>%
    filter(single_sampling_point == 1) %>%
    mutate(month = round_date(sampling_date_single_time_point, unit= 'months')) %>%
    mutate(month = as.numeric(format(month, "%m"))) %>%
    mutate(day = as.numeric(format(sampling_date_single_time_point, "%d"))) %>%
    mutate(year = as.numeric(format(sampling_date_single_time_point, "%Y"))) %>%
    #mutate(month = format(as.Date(sampling_date_single_time_point), "%m-%d"))) %>%
    #mutate(substudy = paste(title, year, species, sep = ', ')) %>%
    #group_by(title, month, year, day, virus, species, ECO_NAME,  methodology, age_class, last_name_of_first_author, year, sampling_location, address) %>%
    group_by(title, month, year, day, virus, methodology, species, sex, age_class, ECO_NAME, sampling_location, sampling_location_two) %>%
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
    filter(difftime < 365/12) %>% #total arbitrary filter, but we need to come up with some time point
    #mutate(month = round_date(sampling_date_single_time_point, unit= 'months')) %>%
    mutate(month_1 = ifelse(difftime < 364, as.numeric(format(start_of_sampling, "%m")), 1)) %>%
    mutate(month_2 = ifelse(difftime < 364, as.numeric(format(end_of_sampling, "%m")), 12)) %>%
    mutate(year = as.numeric(format(start_of_sampling, "%Y"))) %>%
    mutate(day = as.numeric(format(start_of_sampling, "%d"))) %>%
    #unite(group, c(age_class, last_name_of_first_author, year), sep = ", ") %>%
    #group_by(title, month_1, month_2, virus, species, ECO_NAME, methodology, age_class, last_name_of_first_author, year, sampling_location, north_final, south_final, west_final, east_final) %>%
    group_by(title, month_1, month_2, year, day, virus, methodology, species, sex, age_class, ECO_NAME, sampling_location, sampling_location_two) %>%
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

seroprevalence_graph_inner_join <- left_join(seroprevalence_graph, Bat_Birth_Pulse_Data_final) %>%
  mutate(birth_pulse_1_quant = ifelse(species == 'cynopterus sphinx' & is.na(birth_pulse_1_quant) & sampling_location == 'vietnam',7,birth_pulse_1_quant)) %>%
  #for cynopterus sphinx in thailand we are using the estimates from India, which is pretty funky, and should probably be addressed!
  mutate(birth_pulse_1_quant = ifelse(species == 'hipposideros commersoni' & is.na(birth_pulse_1_quant) & sampling_location == 'democratic republic of the congo, durba',10.5,birth_pulse_1_quant)) %>%
  #for hipposideros commersoni in DRC i am using the zimbabwe estimate
  mutate(birth_pulse_1_quant = ifelse(species == 'hipposideros spp' & is.na(birth_pulse_1_quant) & sampling_location == 'uganda, queen elizabeth national park',2.5,birth_pulse_1_quant)) %>%
  # for the hipposiderod spp in uganda going to use the h. caffer estimate from uganda
  mutate(birth_pulse_1_quant = ifelse(species == 'nycteris hispida spp' & is.na(birth_pulse_1_quant) & sampling_location == 'democratic republic of the congo, durba',0,birth_pulse_1_quant)) %>%
  #going to use the myconecteris estimates from the congo i guess
  mutate(birth_pulse_1_quant = ifelse(species == 'pteropus hypomelanus' & is.na(birth_pulse_1_quant) & sampling_location =='thailand',4,birth_pulse_1_quant)) %>%
  mutate(birth_pulse_1_quant = ifelse(species == 'pteropus hypomelanus' & is.na(birth_pulse_1_quant) & sampling_location =='malaysia, johor',4,birth_pulse_1_quant)) %>%
  mutate(birth_pulse_1_quant = ifelse(species == 'pteropus hypomelanus' & is.na(birth_pulse_1_quant) & sampling_location =='malaysia, perak',4,birth_pulse_1_quant)) %>%
  mutate(birth_pulse_1_quant = ifelse(species == 'pteropus hypomelanus' & is.na(birth_pulse_1_quant) & sampling_location =='malaysia, tioman island',4,birth_pulse_1_quant)) %>%
  #for the pteropus hypomelanus estimates from malaysia and thailand, we are using phillippines estimates
  mutate(birth_pulse_1_quant = ifelse(species == 'pteropus scapulatus' & is.na(birth_pulse_1_quant) & sampling_location =='australia, katherine gorge national park',3.5,birth_pulse_1_quant)) %>%
  mutate(birth_pulse_1_quant = ifelse(species == 'pteropus scapulatus' & is.na(birth_pulse_1_quant) & sampling_location =='australia, flora nature park',3.5,birth_pulse_1_quant)) %>%
  mutate(birth_pulse_1_quant = ifelse(species == 'pteropus scapulatus' & is.na(birth_pulse_1_quant) & sampling_location =='australia, elsey national park',3.5,birth_pulse_1_quant)) %>%
  #our location for 'pteropus scapulatus' is from northern australia, which google thinks is a bar in eastern australia, so have to do this one by hand i guess
  mutate(birth_pulse_1_quant = ifelse(species == 'pteropus vampyrus' & is.na(birth_pulse_1_quant & sampling_location =='indonesia, sumatera'),3.0,birth_pulse_1_quant)) %>%
  #going to use the thailand, malaysia estimates for the missing indonesia estimates
  mutate(birth_pulse_1_quant = ifelse(species == 'rousettus leschenaulti' & is.na(birth_pulse_1_quant & sampling_location =='vietnam'),3.0,birth_pulse_1_quant)) 
  #we have estimate for vietnam, not sure why this isn't working
  
seroprevalence_graph_inner_join_filter <- seroprevalence_graph_inner_join %>%
  dplyr::select(-c(birth_pulse_1_quant_flipped)) %>%
  group_by_at(vars(-ECO_NAME, birth_pulse_1_quant, birth_pulse_2_quant)) %>%
  summarise(birth_pulse_1_quant_new = mean(birth_pulse_1_quant)) %>%
  unique() %>%
  mutate(month.dc.birthpulse = month - birth_pulse_1_quant_new) %>%
  mutate(month.dc.birthpulse = ifelse(month.dc.birthpulse > 6, month.dc.birthpulse - 12,
                               ifelse(month.dc.birthpulse < -6, month.dc.birthpulse + 12, month.dc.birthpulse))) %>%
  ungroup() %>%
  mutate(seroprevalence_per = ifelse(seroprevalence_per == 0, .135, seroprevalence_per)) %>% #need to decide what we want to set 0 seroprevalence to
  mutate(seroprevalence_per_lt = log(seroprevalence_per)) %>%
  mutate(substudy = paste(species, sex, age_class, sampling_location,year, title, sep = ', ')) %>%
  mutate(birth_pulse_2_quant = ifelse(birth_pulse_2_quant=='NA', NA, birth_pulse_2_quant)) %>%
  mutate(birth_pulse_two_cat = ifelse(is.na(birth_pulse_2_quant), 'FALSE', 'TRUE')) %>%
  mutate(birth_pulse_two_cat = as.factor(as.character(birth_pulse_two_cat))) %>%
  mutate(substudy = as.factor(as.character(substudy))) %>%
  mutate(title = as.factor(as.character(title))) %>%
  mutate(virus = as.factor(as.character(virus))) %>%
  mutate(methodology = as.factor(as.character(methodology)))

rm(seroprevalence_graph_inner_join)

#plot(seroprevalence_graph_inner_join_filter_dc_style$value)
#plot(seroprevalence_graph_inner_join_filter_dc_style$seroprevalence_per)
#plot(seroprevalence_graph_inner_join_filter_dc_style$seroprevalence_per_lt, seroprevalence_graph_inner_join_filter_dc_style$value)

#...................................gam..........................................

gam3<-gam(seroprevalence_per_lt ~ s(month.dc.birthpulse, by=birth_pulse_two_cat) + s(substudy, bs='re') + s(title, bs='re') + s(methodology, bs='re'), data = seroprevalence_graph_inner_join_filter, method="REML" ,weights=1/sqrt(sample_size))
summary(gam3)
visreg(gam3, "month.dc.birthpulse", by = "birth_pulse_two_cat", gg=TRUE)
#gam.check(gam3)
#plot(gam3$residuals)

# ..................................... graphing ...............................
seroprevalence_graph_inner_join_filter %>%
  filter(virus=='Henipavirus') %>%
  filter(methodology != "non nAb based method") %>%
  #mutate(month.dc.birthpulse = abs(month.dc.birthpulse)) %>%
  #unite(group, c(species, group), sep = ": ") %>%
  #ungroup() %>%
  #ggplot(aes(x= month.dc.birthpulse, y= seroprevalence_per, colour = substudy, group = substudy, text = paste('group', substudy))) +
  ggplot(aes(x= month.dc.birthpulse, y= value, colour = substudy, group = substudy, text = paste('group', substudy))) +
  geom_point() +
  geom_line(size=1.5) +
  #geom_ribbon(alpha = .2, aes(x= month.dc.birthpulse, 
  #                            ymin=seroprevalence_per_lower_bound, 
  #                            ymax=seroprevalence_per_upper_bound,
  #                            fill= group)) +
  theme(legend.position="none") +
  ylab('change in seroprevalence from mean, adjusted by study') +
  xlab('month since (or until) birth pulse') +
  ggtitle(paste(c("Filovirus & Henipavirus Seroprevalence; Data from Meta Analysis")))+
  facet_grid(~birth_pulse_two_cat)




