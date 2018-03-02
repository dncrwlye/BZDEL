#.......................................devin, if you could figure out how to obtain relative deviance explain (around line 250) that would be awesome....

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
library(stargazer)

load("~/Desktop/BDEL/BZDEL/meta_analysis/data/Bat_Birth_Pulse_Data_final_alternative.Rdata")
load("~/Desktop/BDEL/BZDEL/meta_analysis/data/seroprevalence_ecoregions_alternative.Rdata")

#.................cleaning up some spelling errors, make sure to adjust these at excel sheet level. 

Bat_Birth_Pulse_Data_final <- Bat_Birth_Pulse_Data_final %>%
  mutate(species = gsub('epomorphus' ,"epomophorus", species)) %>%
  #mutate(species = gsub('hipperosiderus' ,"hipposideros", species)) %>%
  #mutate(species = gsub('megarops' ,"megaerops", species)) %>%
  mutate(species = gsub('schreibersi' ,"schreibersii", species)) %>%
  mutate(species = gsub('schreibersiii' ,"schreibersii", species)) %>%
  mutate(species = gsub('minopterus' ,"miniopterus", species)) %>%
  mutate(species = gsub('myonicterus' ,"myonycteris", species)) %>%
  mutate(species = gsub('jagoli' ,"jagori", species)) %>%
  mutate(species = gsub('leschenaulti' ,"leschenaultii", species)) %>%
  mutate(species = gsub('leschenaultiii' ,"leschenaultii", species)) %>%
  mutate(species = gsub('khuli' ,"kuhlii", species)) %>%
  mutate(species = gsub('roussetus' ,"roussettus", species)) %>%
  mutate(species = gsub('lavartus' ,"larvatus", species)) 

seroprevalence_x_final <- seroprevalence_x_final %>%
  ungroup() %>%
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
  mutate(species = gsub('lavartus' ,"larvatus", species)) 


#...........................step 1...... clean the data a little for ones where we have unique time points

Bat_Birth_Pulse_Data_final <- Bat_Birth_Pulse_Data_final %>%
  dplyr::select(species, birth_pulse_1_quant, birth_pulse_2_quant, ECO_NAME)

seroprevalence_join <- left_join(seroprevalence_x_final, Bat_Birth_Pulse_Data_final) %>%
  mutate(birth_pulse_1_quant = ifelse(species == 'cynopterus sphinx' & is.na(birth_pulse_1_quant) & sampling_location == 'vietnam',7,birth_pulse_1_quant)) %>%
  #for cynopterus sphinx in thailand we are using the estimates from India, which is pretty funky, and should probably be addressed!
  mutate(birth_pulse_1_quant = ifelse(species == 'hipposideros commersoni' & is.na(birth_pulse_1_quant) & sampling_location == 'democratic republic of the congo, durba',10.5,birth_pulse_1_quant)) %>%
  #for hipposideros commersoni in DRC i am using the zimbabwe estimate
  mutate(birth_pulse_1_quant = ifelse(species == 'hipposideros spp' & is.na(birth_pulse_1_quant) & sampling_location == 'uganda, queen elizabeth national park',2.5,birth_pulse_1_quant)) %>%
  # for the hipposiderod spp in uganda going to use the h. caffer estimate from uganda
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
  mutate(birth_pulse_1_quant = ifelse(species == 'pteropus scapulatus' & is.na(birth_pulse_1_quant) & sampling_location =='australia, elsey national park ',3.5,birth_pulse_1_quant)) %>%
  #our location for 'pteropus scapulatus' is from northern australia, which google thinks is a bar in eastern australia, so have to do this one by hand i guess
  mutate(birth_pulse_1_quant = ifelse(species == 'pteropus vampyrus' & is.na(birth_pulse_1_quant & sampling_location =='indonesia, sumatera'),3.0,birth_pulse_1_quant)) %>%
  #going to use the thailand, malaysia estimates for the missing indonesia estimates
  mutate(birth_pulse_1_quant = ifelse(species == 'rousettus leschenaulti' & is.na(birth_pulse_1_quant & sampling_location =='vietnam'),3.0,birth_pulse_1_quant))  %>%
  #we have estimate for vietnam, not sure why this isn't working
  mutate(birth_pulse_1_quant = ifelse(species == 'eidolon helvum' & sampling_location =='annobon island',8.5,birth_pulse_1_quant))
#the eidolon helvum results are funky because they are part of the same ecoregion as all the other islands that peel works on, so probably better to just do this by hand

#............................................remove extra columns and remove experimental studies ...............

seroprevalence_join2 <- seroprevalence_join %>%
  dplyr::select(-c(north_final, south_final, west_final, east_final, coordinate_box)) %>%
  unique() %>%
  filter(study_type !="Experimental ")

seroprevalence_join2 <- seroprevalence_join2 %>%
  mutate(birth_pulse_2_quant = ifelse(birth_pulse_2_quant=='NA', NA, as.numeric(birth_pulse_2_quant)))

#..........................collapse studies that fall across multiple ecoregions to obtain one birth pulse estimate per study....

seroprevalence_join3 <- seroprevalence_join2 %>%  
  group_by_at(vars(-ECO_NAME, -birth_pulse_1_quant, -birth_pulse_2_quant)) %>%
  mutate(birth_pulse_1_quant_new= mean(birth_pulse_1_quant, na.rm=TRUE), birth_pulse_2_quant_new = mean(birth_pulse_2_quant, na.rm=TRUE)) %>%
  #mutate(birth_pulse_1_quant= mean(birth_pulse_1_quant, na.rm=TRUE),birth_pulse_2_quant = mean(birth_pulse_2_quant, na.rm=TRUE), birth_pulse_1_quant_new_std = sd(birth_pulse_1_quant, na.rm=TRUE)) %>%
  summarise(birth_pulse_1_quant_new = mean(birth_pulse_1_quant, na.rm=TRUE),birth_pulse_2_quant_new = mean(birth_pulse_2_quant, na.rm=TRUE), birth_pulse_1_quant_new_std = sd(birth_pulse_1_quant, na.rm=TRUE)) %>%
  unique() %>%
  ungroup() %>%
  mutate(birth_pulse_two_cat = ifelse(is.na(birth_pulse_2_quant_new), 'Annual Birth Event', 'Biannual Birth Event'))#%>% 
#.........moving this part to an earlier script
#mutate(seroprevalence_percentage = ifelse(is.na(seroprevalence_percentage), successes/sample_size, seroprevalence_percentage)) %>%
#filter(!is.na(seroprevalence_percentage)) 

#.............................making the adjustments to our time series data for those studies that reported a sampling range

seroprevalence_join4 <- seroprevalence_join3 %>%
  mutate(month = round_date(sampling_date_single_time_point, unit= 'month'))%>%
  mutate(month = as.numeric(format(month, "%m"))) %>%
  mutate(difftime = as.numeric(difftime(end_of_sampling, start_of_sampling))) %>%
  filter(difftime < 365/12 | is.na(difftime)) %>% #total arbitrary filter, but we need to come up with some time point
  mutate(end_of_sampling_month = round_date(end_of_sampling, unit= 'months')) %>%
  mutate(start_of_sampling_month = round_date(start_of_sampling, unit = 'months')) %>%
  mutate(end_of_sampling_month = as.numeric(format(end_of_sampling_month, "%m"))) %>%
  mutate(year = ifelse(!is.na(start_of_sampling), as.numeric(format(start_of_sampling, "%Y")), as.numeric(format(sampling_date_single_time_point, "%Y")))) %>%
  mutate(start_of_sampling_month = as.numeric(format(start_of_sampling_month, "%m"))) %>%
  mutate(month = ifelse(!is.na(end_of_sampling_month), (end_of_sampling_month + start_of_sampling_month)/2, month))

seroprevalence_join5 <- seroprevalence_join4 %>%
  filter(!(is.na(birth_pulse_1_quant_new))) %>%
  mutate(month.dc.birthpulse_original = month - birth_pulse_1_quant_new) %>%
  mutate(month.dc.birthpulse = month - birth_pulse_1_quant_new) %>%
  mutate(month.dc.birthpulse = ifelse(month.dc.birthpulse > 6, month.dc.birthpulse - 12,
                                      ifelse(month.dc.birthpulse < -6, month.dc.birthpulse + 12, month.dc.birthpulse))) %>%
  mutate(seroprevalence_percentage = ifelse(seroprevalence_percentage == 0, exp(0), seroprevalence_percentage)) %>% #need to decide what we want to set 0 seroprevalence to
  mutate(seroprevalence_percentage = seroprevalence_percentage/100) %>%
  #mutate(odds = seroprevalence_percentage / (1-seroprevalence_percentage)) %>%
  #mutate(logit = log(odds)) %>%
  mutate(seroprevalence_percentage_lt = log(seroprevalence_percentage)) %>%
  mutate(substudy = paste(title, species, sex, methodology, age_class, sampling_location, year, sep = ', ')) %>%
  #mutate(birth_pulse_2_quant = ifelse(birth_pulse_2_quant=='NA', NA, birth_pulse_2_quant)) %>%
  mutate(birth_pulse_two_cat = as.factor(as.character(birth_pulse_two_cat))) %>%
  mutate(substudy = as.factor(as.character(substudy))) %>%
  mutate(title = as.factor(as.character(title))) %>%
  mutate(virus = as.factor(as.character(virus))) %>%
  mutate(species = as.factor(as.character(species))) %>%
  mutate(methodology = as.factor(as.character(methodology))) 

# ..................................... random scatter of dates...............................

# runif_month<-runif(n=nrow(seroprevalence_join5), -6,6)
# seroprevalence_join5 <- cbind(seroprevalence_join5,runif_month) %>%
#   mutate(runif_month1= round(as.numeric(runif_month)))
# x <-as.data.frame(unique(seroprevalence_join5$title))
# 
# x<-seroprevalence_join5%>%
#   filter(title=='Prevalence of Henipavirus and Rubulavirus Antibodies in Pteropid Bats, Papau New Guinea')

# ..................................... graphing ...............................
plot<-seroprevalence_join5 %>%
  #filter(last_name_of_first_author!="Wacharapluesadee") %>%
  mutate(month.dc.birthpulse_original = ifelse(month.dc.birthpulse_original < 0, 12 - abs(month.dc.birthpulse_original), month.dc.birthpulse_original)) %>%
  ggplot(aes(x= month.dc.birthpulse_original, y= seroprevalence_percentage, colour = substudy, group = substudy, text = paste('group', substudy))) +
  geom_point(aes(size=sample_size)) +
  geom_line(size=1.5) +
  #geom_ribbon(alpha = .2, aes(x= month.dc.birthpulse,
  #                            ymin=seroprevalence_per_lower_bound,
  #                            ymax=seroprevalence_per_upper_bound,
  #                            fill= group)) +
  theme(legend.position="none") +
  ylab('change in seroprevalence from mean, adjusted by study') +
  xlab('month since (or until) birth pulse') +
  ggtitle(paste(c("Filovirus & Henipavirus Seroprevalence; Data from Meta Analysis"))) +
  facet_wrap(~methodology, scales="free")
#facet_grid(~methodology)
plot
#ggplotly(x)

##................................trying alex's idea...................................

seroprevalence_join6c <- seroprevalence_join5 %>%
  mutate(month.dc.birthpulse_original = ifelse(month.dc.birthpulse_original < 0, 12 - abs(month.dc.birthpulse_original), month.dc.birthpulse_original)) %>%
  mutate(methodology = ifelse(methodology == 'PCR based method', 'PCR based method', "Antibody Based Method")) %>%
  mutate(successes = sample_size * seroprevalence_percentage) %>%
  mutate(failures = sample_size-successes) %>%
  data.frame() 

a <- as.data.frame(table(seroprevalence_join6c$substudy)) %>%
  mutate(substudy = Var1) %>%
  select(-c(Var1))

seroprevalence_join7c <- full_join(seroprevalence_join6c, a) %>%
  filter(Freq > 1)

df.successes <- seroprevalence_join7c[rep(row.names(seroprevalence_join7c), seroprevalence_join7c$successes), 1:34] %>%
  mutate(successes = 1)
df.failures <- seroprevalence_join7c[rep(row.names(seroprevalence_join7c), seroprevalence_join7c$failures), 1:34] %>%
  mutate(successes = 0)
seroprevalence_join8c <- rbind(df.successes,df.failures)

seroprevalence_join9c <- seroprevalence_join8c %>%
  mutate(parabola = month.dc.birthpulse_original * (12-month.dc.birthpulse_original))

#test one
#glmer(successes ~ parabola + (1|substudy), family=binomial(link='logit'), data= seroprevalence_join9) %>% summary()

seroprevalence_join9c %>%
  ggplot(aes(x = parabola, y= seroprevalence_percentage, colour = substudy, group = substudy)) + 
  geom_point(aes(size=sample_size)) +
  geom_line(size=1.5) + 
  theme(legend.position="none") + 
  facet_wrap(~methodology, scales = 'free')

nore <- glm(successes ~ parabola + virus + methodology + parabola*methodology , family=binomial(link='logit'), data= seroprevalence_join9c)  
re<- glmer(successes ~ parabola + virus + methodology + parabola*methodology + (1|substudy), family=binomial(link='logit'), data= seroprevalence_join9c)  
anova( re, nore)


