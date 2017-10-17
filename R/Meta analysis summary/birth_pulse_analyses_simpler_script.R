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

Bat_Birth_Pulse_Data_final <- Bat_Birth_Pulse_Data_final %>%
  select(species, birth_pulse_1_quant, birth_pulse_2_quant, ECO_NAME)

seroprevalence_join <- left_join(seroprevalence_x_final, Bat_Birth_Pulse_Data_final) %>%
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
  mutate(birth_pulse_1_quant = ifelse(species == 'rousettus leschenaulti' & is.na(birth_pulse_1_quant & sampling_location =='vietnam'),3.0,birth_pulse_1_quant))  %>%
  #we have estimate for vietnam, not sure why this isn't working
  mutate(birth_pulse_1_quant = ifelse(species == 'eidolon helvum' & is.na(birth_pulse_1_quant & sampling_location =='annobon island'),8.5,birth_pulse_1_quant)) 
  #the eidolon helvum results are funky because they are part of the same ecoregion as all the other islands that peel works on, so probably better to just do this by hand


seroprevalence_join2 <- seroprevalence_join %>%
  select(-c(north_final, south_final, west_final, east_final, coordinate_box)) %>%
  unique()

seroprevalence_join3 <- seroprevalence_join2 %>%
  group_by_at(vars(-ECO_NAME, -birth_pulse_1_quant, birth_pulse_2_quant)) %>%
  summarise(birth_pulse_1_quant_new = mean(birth_pulse_1_quant, na.rm=TRUE)) %>%
  unique() %>%
  ungroup() %>%
  mutate(seroprevalence_percentage = ifelse(is.na(seroprevalence_percentage), successes/sample_size, seroprevalence_percentage)) %>%
  filter(!is.na(seroprevalence_percentage))

#x <- seroprevalence_join2 %>%   filter(sampling_location=='annobon island') %>% select(-c(ECO_NAME, birth_pulse_1_quant, birth_pulse_2_quant))%>% unique()
#x <- seroprevalence_join2 %>% filter(sampling_location=='annobon island') %>% group_by_at(vars(-ECO_NAME, -birth_pulse_1_quant, birth_pulse_2_quant)) %>% summarise(x = mean(birth_pulse_1_quant,na.rm=TRUE))

seroprevalence_single_time_point <- seroprevalence_join3 %>%
  filter(single_sampling_point == 1) %>%
  mutate(month = round_date(sampling_date_single_time_point, unit= 'months')) %>%
  mutate(month = as.numeric(format(month, "%m"))) %>%
  mutate(day = as.numeric(format(sampling_date_single_time_point, "%d"))) %>%
  mutate(year = as.numeric(format(sampling_date_single_time_point, "%Y"))) 

seroprevalence_unclear_time_point <- seroprevalence_join3 %>%
  filter(single_sampling_point == 0) %>%
  filter(!is.na(start_of_sampling)) %>%
  mutate(difftime = as.numeric(difftime(end_of_sampling, start_of_sampling))) %>%
  mutate(difftime = ifelse(difftime > 1000, 1000, difftime)) %>%
  filter(difftime < 365/12) %>% #total arbitrary filter, but we need to come up with some time point
  mutate(month_1 = ifelse(difftime < 364, as.numeric(format(start_of_sampling, "%m")), 1)) %>%
  mutate(month_2 = ifelse(difftime < 364, as.numeric(format(end_of_sampling, "%m")), 12)) %>%
  mutate(year = as.numeric(format(start_of_sampling, "%Y"))) %>%
  mutate(day = as.numeric(format(start_of_sampling, "%d"))) %>%
  select(-c(difftime)) %>%
  mutate(title = paste(title, "xc")) 

x <- seroprevalence_unclear_time_point[,c(1:24,26:ncol(seroprevalence_unclear_time_point))]
y <-seroprevalence_unclear_time_point[,c(1:23,25:ncol(seroprevalence_unclear_time_point))]
colnames(y) <- colnames(x)
seroprevalence_unclear_time_point <- rbind(y,x) %>%
  rename(month = month_1)

#colnames(seroprevalence_henipavirus_unclear_time_point) <- colnames(seroprevalence_henipavirus_single_time_point)
seroprevalence_join4 <- rbind(seroprevalence_unclear_time_point, seroprevalence_single_time_point)

rm(x,y, seroprevalence_single_time_point, seroprevalence_unclear_time_point)

#.............................................making the adjustments to our time series data

# also join

seroprevalence_join5 <- seroprevalence_join4 %>%
  mutate(month.dc.birthpulse = month - birth_pulse_1_quant_new) %>%
  mutate(month.dc.birthpulse = ifelse(month.dc.birthpulse > 6, month.dc.birthpulse - 12,
                                      ifelse(month.dc.birthpulse < -6, month.dc.birthpulse + 12, month.dc.birthpulse))) %>%
  mutate(seroprevalence_percentage_lt = ifelse(seroprevalence_percentage == 0, exp(0), seroprevalence_percentage)) %>% #need to decide what we want to set 0 seroprevalence to
  mutate(seroprevalence_per_lt = log(seroprevalence_percentage_lt)) %>%
  mutate(substudy = paste(species, sex, age_class, sampling_location,year, title, sep = ', ')) %>%
  mutate(birth_pulse_2_quant = ifelse(birth_pulse_2_quant=='NA', NA, birth_pulse_2_quant)) %>%
  mutate(birth_pulse_two_cat = ifelse(is.na(birth_pulse_2_quant), 'FALSE', 'TRUE')) %>%
  mutate(birth_pulse_two_cat = as.factor(as.character(birth_pulse_two_cat))) %>%
  mutate(substudy = as.factor(as.character(substudy))) %>%
  mutate(title = as.factor(as.character(title))) %>%
  mutate(virus = as.factor(as.character(virus))) %>%
  mutate(species = as.factor(as.character(species))) %>%
  mutate(methodology = as.factor(as.character(methodology))) 

runif_month<-runif(n=nrow(seroprevalence_join5), -6,6)
seroprevalence_join5 <- cbind(seroprevalence_join5,runif_month) %>%
  mutate(runif_month1= round(as.numeric(runif_month)))

#...................................gam..........................................

seroprevalence_join6 <- seroprevalence_join5 %>%
  mutate(seroprevalence_percentage = seroprevalence_percentage/100) %>%
  mutate(binomial_variance = seroprevalence_percentage*sample_size*(1-seroprevalence_percentage)) %>%
  mutate(inverse_binomial_variance = 1/ binomial_variance)

gam3<-gam(seroprevalence_per_lt ~ s(month.dc.birthpulse, by=birth_pulse_two_cat) + s(species, bs='re') + s(substudy, bs='re') + s(title, bs='re') + s(virus, bs='re')+ s(methodology, bs='re'), data = seroprevalence_join5, method="REML", weights=log(sample_size))
summary(gam3)
#plot(predict.gam(gam3))
visreg(gam3, "month.dc.birthpulse", by = "birth_pulse_two_cat") 
gam.check(gam3)
#plot(gam3$residuals)

gam3<-gam(seroprevalence_per_lt ~ s(runif_month, by=birth_pulse_two_cat) + s(species, bs='re') + s(substudy, bs='re') + s(title, bs='re') + s(virus, bs='re')+ s(methodology, bs='re'), data = seroprevalence_join5, method="REML", weights=log(sample_size))
summary(gam3)
#plot(predict.gam(gam3))
visreg(gam3, "runif_month", by = "birth_pulse_two_cat") 





gam3<-gam(seroprevalence_per_lt ~ s(month.dc.birthpulse, by=birth_pulse_two_cat) + s(substudy, bs='re') + s(title, bs='re') + s(methodology, bs='re'), data = seroprevalence_join6, method="REML" ,weights=inverse_binomial_variance)
summary(gam3)
visreg(gam3, "month.dc.birthpulse", by = "birth_pulse_two_cat", gg=TRUE) 
gam.check(gam3)

seroprevalence_join5 %>%
  group_by(title) %>%
  summarise(n=n())
# ..................................... random scatter ...............................

runif(1,-6,6)




# ..................................... graphing ...............................
x<-seroprevalence_join5 %>%
  filter(last_name_of_first_author == "Amman") %>%
  #filter(virus=='Henipavirus') %>%
  #filter(methodology != "non nAb based method") %>%
  #mutate(month.dc.birthpulse = abs(month.dc.birthpulse)) %>%
  #unite(group, c(species, group), sep = ": ") %>%
  #ungroup() %>%
  #ggplot(aes(x= month.dc.birthpulse, y= seroprevalence_per, colour = substudy, group = substudy, text = paste('group', substudy))) +
  ggplot(aes(x= month.dc.birthpulse, y= seroprevalence_per_lt, colour = substudy, group = substudy, text = paste('group', substudy))) +
  geom_point(aes(size=sample_size)) +
  geom_line(size=1.5) +
  #geom_ribbon(alpha = .2, aes(x= month.dc.birthpulse, 
  #                            ymin=seroprevalence_per_lower_bound, 
  #                            ymax=seroprevalence_per_upper_bound,
  #                            fill= group)) +
  theme(legend.position="none") +
  ylab('change in seroprevalence from mean, adjusted by study') +
  xlab('month since (or until) birth pulse') +
  ggtitle(paste(c("Filovirus & Henipavirus Seroprevalence; Data from Meta Analysis")))+
  facet_grid(methodology~virus~birth_pulse_two_cat) 
ggplotly(x)
