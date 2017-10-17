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
  mutate(month.dc.birthpulse_original = month - birth_pulse_1_quant_new) %>%
  mutate(month.dc.birthpulse = month - birth_pulse_1_quant_new) %>%
  mutate(month.dc.birthpulse = ifelse(month.dc.birthpulse > 6, month.dc.birthpulse - 12,
                                      ifelse(month.dc.birthpulse < -6, month.dc.birthpulse + 12, month.dc.birthpulse))) %>%
  mutate(seroprevalence_percentage = ifelse(seroprevalence_percentage == 0, exp(0), seroprevalence_percentage)) %>% #need to decide what we want to set 0 seroprevalence to
  mutate(seroprevalence_percentage = seroprevalence_percentage/100) %>%
  #mutate(odds = seroprevalence_percentage / (1-seroprevalence_percentage)) %>%
  #mutate(logit = log(odds)) %>%
  mutate(seroprevalence_percentage_lt = log(seroprevalence_percentage)) %>%
  mutate(substudy = paste(species, sex, age_class, sampling_location,year, title, sep = ', ')) %>%
  mutate(birth_pulse_2_quant = ifelse(birth_pulse_2_quant=='NA', NA, birth_pulse_2_quant)) %>%
  mutate(birth_pulse_two_cat = ifelse(is.na(birth_pulse_2_quant), 'FALSE', 'TRUE')) %>%
  mutate(birth_pulse_two_cat = as.factor(as.character(birth_pulse_two_cat))) %>%
  mutate(substudy = as.factor(as.character(substudy))) %>%
  mutate(title = as.factor(as.character(title))) %>%
  mutate(virus = as.factor(as.character(virus))) %>%
  mutate(species = as.factor(as.character(species))) %>%
  mutate(methodology = as.factor(as.character(methodology))) 

hist((seroprevalence_join5$logit))
hist(seroprevalence_join5$seroprevalence_percentage)

# ..................................... random scatter of dates...............................

runif_month<-runif(n=nrow(seroprevalence_join5), -6,6)
seroprevalence_join5 <- cbind(seroprevalence_join5,runif_month) %>%
  mutate(runif_month1= round(as.numeric(runif_month)))

#...................................log binomial variance..........................................

seroprevalence_join6 <- seroprevalence_join5 %>%
  mutate(log_binomial_variance = ((1+(exp(1)^2-1)*seroprevalence_percentage)^sample_size) - ((1+(exp(1)-1)*seroprevalence_percentage)^(2*sample_size))) %>%
  mutate(log_binomial_variance_inverse = 1/log_binomial_variance) %>%
  mutate(inverse_variance = 1/(seroprevalence_percentage*(1-seroprevalence_percentage))/sample_size)



#...................................metafor and vi..........................................

library(metafor)
## convert to effect size (logit transformed)
## xi is whatever your vector of successes are)
## ni  is the vector of sample size
## PLO tells this to calculate logit-transformed proportions (yi) and the corresponding sampling variance (vi)
#data=data.frame(data,escalc(xi=round(data$serop*data$sample,0),ni=data$sample,measure="PLO"))

data=as.data.frame(escalc(xi=round(seroprevalence_join6$seroprevalence_percentage*seroprevalence_join6$sample_size,0),ni=seroprevalence_join6$sample_size,measure="PLO"))

seroprevalence_join6 <- cbind(seroprevalence_join6, data)

gam3<-gam(seroprevalence_percentage_lt ~ s(month.dc.birthpulse, by=birth_pulse_two_cat) + s(species, bs='re') + s(substudy, bs='re') + s(title, bs='re') + s(virus, bs='re')+ s(methodology, bs='re'), data = seroprevalence_join5, method="REML", weights=sqrt(sample_size))
summary(gam3)
#plot(predict.gam(gam3))
visreg(gam3, "month.dc.birthpulse", by = "birth_pulse_two_cat") 
gam.check(gam3)
qqPlot(gam3$residuals)

gam3<-gam(seroprevalence_percentage_lt ~ s(runif_month, by=birth_pulse_two_cat) + s(species, bs='re') + s(substudy, bs='re') + s(title, bs='re') + s(virus, bs='re')+ s(methodology, bs='re'), data = seroprevalence_join5, method="REML", weights=(sample_size))
summary(gam3)
#plot(predict.gam(gam3))
visreg(gam3, "runif_month", by = "birth_pulse_two_cat") 

gam3<-gam(seroprevalence_percentage_lt ~ s(month.dc.birthpulse, by=birth_pulse_two_cat)  + s(substudy, bs='re')+ s(methodology, bs='re'), data = seroprevalence_join6, method="REML", weights=inverse_variance)
summary(gam3)
#plot(predict.gam(gam3))
visreg(gam3, "month.dc.birthpulse", by = "birth_pulse_two_cat") 
gam.check(gam3)


plot(1/(seroprevalence_join6$vi), seroprevalence_join6$sample_size)

#........................becker vi yi thing...................................................

gam3<-gam(yi ~ s(month.dc.birthpulse, by=birth_pulse_two_cat)  + s(substudy, bs='re')+ s(methodology, bs='re'), data = seroprevalence_join6, method="REML", weights=vi)
summary(gam3)
#plot(predict.gam(gam3))
visreg(gam3, "month.dc.birthpulse", by = "birth_pulse_two_cat") 
gam.check(gam3)
qqPlot(gam3$residuals)

gam3<-gam(yi ~ s(month.dc.birthpulse, by=birth_pulse_two_cat)  + s(substudy, bs='re') + s(species, bs='re') + s(methodology, bs='re'), data = seroprevalence_join6, method="REML", weights=(1/vi))
summary(gam3)
visreg(gam3, "month.dc.birthpulse", by = "birth_pulse_two_cat") 
gam.check(gam3)

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

#.....................binomial stuffss for pcr...............................................

seroprevalence_join7 <- seroprevalence_join6 %>%
  mutate(month.dc.birthpulse_original = ifelse(month.dc.birthpulse_original < 0, 12 - abs(month.dc.birthpulse_original), month.dc.birthpulse_original)) %>%
  filter(methodology == 'PCR based method') %>%
  mutate(failures = sample_size-successes) %>%
  data.frame()

df.successes <- seroprevalence_join7[rep(row.names(seroprevalence_join7), seroprevalence_join7$successes), 1:35] %>%
  mutate(successes = 1)
df.failures <- seroprevalence_join7[rep(row.names(seroprevalence_join7), seroprevalence_join7$failures), 1:35] %>%
  mutate(successes = 0)
seroprevalence_join8 <- rbind(df.successes,df.failures)

plot(seroprevalence_join8$month.dc.birthpulse, seroprevalence_join8$month.dc.birthpulse_original)
unique(seroprevalence_join8$month.dc.birthpulse_original)

gam8<-gam(successes ~ s(month.dc.birthpulse_original, by=birth_pulse_two_cat)  + s(substudy, bs='re'), data = seroprevalence_join8, method="REML", family=binomial(link='logit'))
summary(gam8)
visreg(gam8, "month.dc.birthpulse_original", by = "birth_pulse_two_cat") 
gam.check(gam8)

unique(seroprevalence_join7$substudy)

#.....................binomial stuffss for serology...............................................

seroprevalence_join7a <- seroprevalence_join6 %>%
  mutate(month.dc.birthpulse_original = ifelse(month.dc.birthpulse_original < 0, 12 - abs(month.dc.birthpulse_original), month.dc.birthpulse_original)) %>%
  filter(methodology == 'nAb based method' | methodology == "non nAb based method") %>%
  mutate(failures = sample_size-successes) %>%
  data.frame()

df.successes <- seroprevalence_join7a[rep(row.names(seroprevalence_join7a), seroprevalence_join7a$successes), 1:35] %>%
  mutate(successes = 1)
df.failures <- seroprevalence_join7a[rep(row.names(seroprevalence_join7a), seroprevalence_join7a$failures), 1:35] %>%
  mutate(successes = 0)
seroprevalence_join8a <- rbind(df.successes,df.failures)

gam8a<-gam(successes ~ s(month.dc.birthpulse_original, by=birth_pulse_two_cat)  + s(title, bs='re') + s(substudy, bs='re') + s(virus, bs='re') + s(methodology, bs='re'), data = seroprevalence_join8a, method="REML", family=binomial(link='logit'))
summary(gam8a)
visreg(gam8a, "month.dc.birthpulse_original", by = "birth_pulse_two_cat") 
gam.check(gam8a)

#.....................binomial stuffss for all data...............................................

seroprevalence_join7b <- seroprevalence_join6 %>%
  mutate(month.dc.birthpulse_original = ifelse(month.dc.birthpulse_original < 0, 12 - abs(month.dc.birthpulse_original), month.dc.birthpulse_original)) %>%
  #filter(methodology == 'PCR based method') %>%
  mutate(failures = sample_size-successes) %>%
  data.frame()

df.successes <- seroprevalence_join7b[rep(row.names(seroprevalence_join7b), seroprevalence_join7b$successes), 1:35] %>%
  mutate(successes = 1)
df.failures <- seroprevalence_join7b[rep(row.names(seroprevalence_join7b), seroprevalence_join7b$failures), 1:35] %>%
  mutate(successes = 0)
seroprevalence_join8b <- rbind(df.successes,df.failures)

gam8b<-gam(successes ~ s(month.dc.birthpulse_original, by=(birth_pulse_two_cat, methodology))  + s(substudy, bs='re'), data = seroprevalence_join8b, method="REML", family=binomial(link='logit'))
summary(gam8b)
visreg(gam8b, "month.dc.birthpulse_original", by = "methodology") 
gam.check(gam8b)

unique(seroprevalence_join7$substudy)
