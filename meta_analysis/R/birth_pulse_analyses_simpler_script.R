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

load("~/BZDEL/meta_analysis/data/Bat_Birth_Pulse_Data_final_alternative.Rdata")
load("~/Desktop/BDEL/BZDEL/meta_analysis/data/Bat_Birth_Pulse_Data_final_alternative.Rdata")
#load("~/BZDEL/Data/MetaAnalysis/seroprevalence_ecoregions.Rdata")
load("~/BZDEL/meta_analysis/data/seroprevalence_ecoregions_alternative.Rdata")
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

x<-Bat_Birth_Pulse_Data_final %>%
  dplyr::select(c(species)) %>%
  mutate(p = 'birth') %>%
  unique()

y<-seroprevalence_x_final %>%
  dplyr::select(c(species)) %>%
  mutate(p1='seroprev') %>%
  unique()

xy <- full_join(x,y)
xy <- xy[order(xy$species),]

x<-seroprevalence_x_final %>% #find out whats going with dates
  filter(is.na(single_sampling_point) & is.na(start_of_sampling))
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

seroprevalence_join2 <- seroprevalence_join %>%
  select(-c(north_final, south_final, west_final, east_final, coordinate_box)) %>%
  unique() %>%
  filter(study_type !="Experimental ")

seroprevalence_join2 <- seroprevalence_join2 %>%
  mutate(birth_pulse_2_quant = ifelse(birth_pulse_2_quant=='NA', NA, as.numeric(birth_pulse_2_quant)))

seroprevalence_join3 <- seroprevalence_join2 %>%  
  group_by_at(vars(-ECO_NAME, -birth_pulse_1_quant, -birth_pulse_2_quant)) %>%
  mutate(birth_pulse_1_quant_new= mean(birth_pulse_1_quant, na.rm=TRUE), birth_pulse_2_quant_new = mean(birth_pulse_2_quant, na.rm=TRUE)) %>%
  #mutate(birth_pulse_1_quant= mean(birth_pulse_1_quant, na.rm=TRUE),birth_pulse_2_quant = mean(birth_pulse_2_quant, na.rm=TRUE), birth_pulse_1_quant_new_std = sd(birth_pulse_1_quant, na.rm=TRUE)) %>%
  summarise(birth_pulse_1_quant_new = mean(birth_pulse_1_quant, na.rm=TRUE),birth_pulse_2_quant_new = mean(birth_pulse_2_quant, na.rm=TRUE), birth_pulse_1_quant_new_std = sd(birth_pulse_1_quant, na.rm=TRUE)) %>%
  unique() %>%
  ungroup() %>%
  mutate(birth_pulse_two_cat = ifelse(is.na(birth_pulse_2_quant_new), 'Annual Birth Event', 'Biannual Birth Event'))%>% 
  mutate(seroprevalence_percentage = ifelse(is.na(seroprevalence_percentage), successes/sample_size, seroprevalence_percentage)) %>%
  filter(!is.na(seroprevalence_percentage)) 

#.............................making the adjustments to our time series data

# also join

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

# ..................................... graphing ...............................
x<-seroprevalence_join5 %>%
  filter(methodology == 'PCR based method') %>%
  #filter(birth_pulse_two_cat==TRUE) %>%
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
  facet_grid(methodology~birth_pulse_two_cat)

xy<-x %>%
  filter(sampling_location=='indonesia, sumatera')

ggplotly(x)
#.....................binomial stuffss for pcr...............................................

seroprevalence_join7 <- seroprevalence_join5 %>%
  mutate(month.dc.birthpulse_original = ifelse(month.dc.birthpulse_original < 0, 12 - abs(month.dc.birthpulse_original), month.dc.birthpulse_original)) %>%
  filter(methodology == 'PCR based method') %>%
  mutate(failures = sample_size-successes) %>%
  data.frame()

df.successes <- seroprevalence_join7[rep(row.names(seroprevalence_join7), seroprevalence_join7$successes), 1:34] %>%
  mutate(successes = 1)
df.failures <- seroprevalence_join7[rep(row.names(seroprevalence_join7), seroprevalence_join7$failures), 1:34] %>%
  mutate(successes = 0)
seroprevalence_join8 <- rbind(df.successes,df.failures)

#gam8_alt<-gam(successes ~ s(month.dc.birthpulse_original, by=birth_pulse_two_cat)  + virus + s(substudy, bs='re') + s(title, bs='re'), data = seroprevalence_join8, method="REML", family=binomial(link='logit'))


gam8<-gam(successes   ~ s(month.dc.birthpulse_original)  + s(substudy, bs='re') , data = seroprevalence_join8, method="REML", family=binomial(link='logit'))
#gam8.title<-gam(successes   ~ s(month.dc.birthpulse_original)  + s(substudy, bs='re') +s(title, bs='re'), data = seroprevalence_join8, method="REML", family=binomial(link='logit'))
#gam8_alt<-gam(successes ~ s(month.dc.birthpulse_original, substudy, bs='fs')  + s(substudy, bs='re') , data = seroprevalence_join8, method="REML", family=binomial(link='logit'))
gam8.1<-gam(successes ~                                  + s(substudy, bs='re') , data = seroprevalence_join8, method="REML", family=binomial(link='logit'))
#gam8.2<-gam(successes ~ s(month.dc.birthpulse_original)  + s(substudy, bs='re')                       , data = seroprevalence_join8, method="REML", family=binomial(link='logit'))
gam8.3<-gam(successes ~ s(month.dc.birthpulse_original)                         , data = seroprevalence_join8, method="REML", family=binomial(link='logit'))
gam8.null<-gam(successes ~ 1, data = seroprevalence_join8, method="REML", family=binomial(link='logit'))

summary(gam8_alt)
#basing this from one of simon woods suggestions: 
#http://grokbase.com/t/r/r-help/11ba8hjbn5/r-sum-of-the-deviance-explained-by-each-term-in-a-gam-model-does-not-equal-to-the-deviance-explained-by-the-full-model

dev.birth_pulse <- abs(deviance(gam8.1)-deviance(gam8))/deviance(gam8.null) ## prop explained by
#dev.title <- abs(deviance(gam8.2)-deviance(gam8))/deviance(gam8.null) ## prop explained by
dev.substudy <- abs(deviance(gam8.3)-deviance(gam8))/deviance(gam8.null) ## prop explained by

dev.birth_pulse_rde <- round(100*dev.birth_pulse/(dev.birth_pulse + dev.substudy),2)
dev.substudy_rde <-round(100*dev.substudy/(dev.birth_pulse + dev.substudy),2)
#dev.title_rde <- round(100*dev.title/(dev.title+dev.birth_pulse + dev.substudy),2)

summary_gam8 <- summary(gam8)
visreg(gam8, "month.dc.birthpulse_original", gg=TRUE)
gam.check(gam8)

# #.....................binomial stuffss for serology...............................................

seroprevalence_join7a <- seroprevalence_join5 %>%
  mutate(month.dc.birthpulse_original = ifelse(month.dc.birthpulse_original < 0, 12 - abs(month.dc.birthpulse_original), month.dc.birthpulse_original)) %>%
  filter(methodology == 'nAb based method' | methodology == "non nAb based method") %>%
  mutate(failures = sample_size-successes) %>%
  data.frame() 

df.successes <- seroprevalence_join7a[rep(row.names(seroprevalence_join7a), seroprevalence_join7a$successes), 1:34] %>%
  mutate(successes = 1)
df.failures <- seroprevalence_join7a[rep(row.names(seroprevalence_join7a), seroprevalence_join7a$failures), 1:34] %>%
  mutate(successes = 0)
seroprevalence_join8a <- rbind(df.successes,df.failures) %>%
  mutate(last_name_of_first_author = as.factor(last_name_of_first_author))

gam8a     <-gam(successes ~ s(month.dc.birthpulse_original) + s(substudy, bs='re')    , data = seroprevalence_join8a, method="REML", family=binomial(link='logit'))


gam8a     <-gam(successes ~ s(month.dc.birthpulse_original) + s(substudy, bs='re')    , data = seroprevalence_join8a, method="REML", family=binomial(link='logit'))
#gam8a_alt     <-gam(successes ~ s(month.dc.birthpulse_original, substudy, bs = 'fs')  + s(substudy, bs='re')   , data = seroprevalence_join8a, method="REML", family=binomial(link='logit'))
gam8a.1   <-gam(successes ~                                   s(substudy, bs='re')    , data = seroprevalence_join8a, method="REML", family=binomial(link='logit'))
gam8a.2   <-gam(successes ~ s(month.dc.birthpulse_original)                           , data = seroprevalence_join8a, method="REML", family=binomial(link='logit'))
#gam8a.3   <-gam(successes ~ s(month.dc.birthpulse_original) + s(substudy, bs = 're')                     , data = seroprevalence_join8a, method="REML", family=binomial(link='logit'))
gam8a.null<-gam(successes ~ 1                                                                            , data = seroprevalence_join8a, method="REML", family=binomial(link='logit'))

dev.birth_pulse_a <- abs(deviance(gam8a.1)-deviance(gam8a))/deviance(gam8a.null) ## prop explained by
dev.substudy_a    <- abs(deviance(gam8a.2)-deviance(gam8a))/deviance(gam8a.null) ## prop explained by
#dev.title_a       <- abs(deviance(gam8a.3)-deviance(gam8a))/deviance(gam8a.null) ## prop explained by

dev.birth_pulse_rde_a <- round(100*dev.birth_pulse_a/(dev.birth_pulse_a + dev.substudy_a),2)
dev.substudy_rde_a <-round(100*dev.substudy_a/(dev.birth_pulse_a + dev.substudy_a),2)
#dev.title_rde_a <- round(100*dev.title_a/(dev.title_a+dev.birth_pulse_a + dev.substudy_a),2)

summary_gam8_a <- summary(gam8a)

#.................................making regression table..............................................

# study_table <- matrix(c(
#   '','Term',"Value", "Z Statistic", "Chi-Sq Statistic", "P-value", "Effective DOF", "Total Deviance Explained", "Relative Deviance Explained",
#   #...........................pcr data.........................................................
#   paste("PCR Prevalence Data (n =", summary_gam8$n,")", sep=""),"","","","","","", paste(round(100*summary_gam8$dev.expl,2), "%") ,"",
#                       "","Intercept", round(summary_gam8$p.table[1,1],2), round(summary_gam8$p.table[1,3],2),"","","","","",
#   "","Months since birth event","","",round(summary_gam8$s.table[1,3],2), round(summary_gam8$s.table[1,4],2),round(summary_gam8$s.table[1,1],2),"", paste(dev.birth_pulse_rde,"%"),
#                  "", "Substudy","","",round(summary_gam8$s.table[2,3],2), round(summary_gam8$s.table[2,4],2),round(summary_gam8$s.table[2,1],2),"", paste(dev.substudy_rde,"%"),
#                  "",  "Title","","",  round(summary_gam8$s.table[3,3],2), round(summary_gam8$s.table[3,4],2),round(summary_gam8$s.table[3,1],2),"", paste(dev.title_rde,"%"),
#   #..........................serology data.....................................................
#   paste("Ab seroprevalence Data (n =", summary_gam8_a$n,")",sep=""),"","","","","","", paste(round(100*summary_gam8_a$dev.expl,2), "%") ,"",
#                       "","Intercept", round(summary_gam8$p.table[1,1],2), round(summary_gam8_a$p.table[1,3],2),"","","","","",
#   "","Months since birth event","","",round(summary_gam8_a$s.table[1,3],2),round(summary_gam8_a$s.table[1,4],2),round(summary_gam8_a$s.table[1,1],2),"", paste(dev.birth_pulse_rde_a,"%"),
#                  "", "Substudy","","",round(summary_gam8_a$s.table[2,3],2),round(summary_gam8_a$s.table[2,4],2),round(summary_gam8_a$s.table[2,1],2),"", paste(dev.substudy_rde_a,"%"),
#                  "",  "Title","","",  round(summary_gam8_a$s.table[3,3],2),round(summary_gam8_a$s.table[3,4],2),round(summary_gam8_a$s.table[3,1],2),"", paste(dev.title_rde_a,"%")), ncol = 9, byrow = TRUE)
# 
# study_table


study_table <- matrix(c(
  '','Term',"Value", "Z Statistic", "Chi-Sq Statistic", "P-value", "Effective DOF", "Total Deviance Explained", "Relative Deviance Explained",
  #...........................pcr data.........................................................
  paste("PCR Prevalence Data (n =", summary_gam8$n,")", sep=""),"","","","","","", paste(round(100*summary_gam8$dev.expl,2), "%") ,"",
  "","Intercept", round(summary_gam8$p.table[1,1],2), round(summary_gam8$p.table[1,3],2),"",round(summary_gam8$p.table[1,4],2),"","","",
  "","Months since birth event","","",round(summary_gam8$s.table[1,3],2), round(summary_gam8$s.table[1,4],2),round(summary_gam8$s.table[1,1],2),"", paste(dev.birth_pulse_rde,"%"),
  "", "Substudy","","",round(summary_gam8$s.table[2,3],2), round(summary_gam8$s.table[2,4],2),round(summary_gam8$s.table[2,1],2),"", paste(dev.substudy_rde,"%"),
  #"",  "Title","","",  round(summary_gam8$s.table[3,3],2), round(summary_gam8$s.table[3,4],2),round(summary_gam8$s.table[3,1],2),"", paste(dev.title_rde,"%"),
  #..........................serology data.....................................................
  paste("Ab seroprevalence Data (n =", summary_gam8_a$n,")",sep=""),"","","","","","", paste(round(100*summary_gam8_a$dev.expl,2), "%") ,"",
  "","Intercept", round(summary_gam8_a$p.table[1,1],2), round(summary_gam8_a$p.table[1,3],2),"",round(summary_gam8_a$p.table[1,4],2),"","","",
  "","Months since birth event","","",round(summary_gam8_a$s.table[1,3],2),round(summary_gam8_a$s.table[1,4],2),round(summary_gam8_a$s.table[1,1],2),"", paste(dev.birth_pulse_rde_a,"%"),
  "", "Substudy","","",round(summary_gam8_a$s.table[2,3],2),round(summary_gam8_a$s.table[2,4],2),round(summary_gam8_a$s.table[2,1],2),"", paste(dev.substudy_rde_a,"%")
  #,"",  "Title","","",  round(summary_gam8_a$s.table[3,3],2),round(summary_gam8_a$s.table[3,4],2),round(summary_gam8_a$s.table[3,1],2),"", paste(dev.title_rde_a,"%")),
  ), ncol = 9, byrow = TRUE)

study_table

write.csv(study_table, file= '/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/GAM.csv')
# #...............................anti join...........................................
# seroprevalence_join_no_month <- seroprevalence_join %>%
#   filter(!(is.na(birth_pulse_1_quant)))
#          
# seroprevalence_anti_join <- anti_join(seroprevalence_x_final, Bat_Birth_Pulse_Data_final) %>%
#   filter(methodology != 'PCR based method')  %>%
#   mutate(difftime = as.numeric(difftime(end_of_sampling, start_of_sampling))) %>%
#   mutate(difftime = ifelse(difftime > 1000, 1000, difftime)) %>%
#   filter(difftime < 365/12 | is.na(difftime))  %>% #total arbitrary filter, but we need to come up with some time point
#   filter(!(row_number %in% seroprevalence_join_no_month$row_number))

##................................trying alex's idea...................................

# seroprevalence_join7 <- seroprevalence_join5 %>%
#   mutate(month.dc.birthpulse_original = ifelse(month.dc.birthpulse_original < 0, 12 - abs(month.dc.birthpulse_original), month.dc.birthpulse_original)) %>%
#   filter(methodology == 'PCR based method') %>%
#   mutate(failures = sample_size-successes) %>%
#   data.frame()
# #
# df.successes <- seroprevalence_join7[rep(row.names(seroprevalence_join7), seroprevalence_join7$successes), 1:34] %>%
#   mutate(successes = 1)
# df.failures <- seroprevalence_join7[rep(row.names(seroprevalence_join7), seroprevalence_join7$failures), 1:34] %>%
#   mutate(successes = 0)
# seroprevalence_join8 <- rbind(df.successes,df.failures)
# #
# seroprevalence_join9 <- seroprevalence_join8 %>%
#   mutate(parabola = month.dc.birthpulse_original * (12-month.dc.birthpulse_original))
# 
# #test one
# glmer(successes ~ parabola, family=binomial(link='logit'), data= seroprevalence_join9) %>% summary()

# x<-as.data.frame(t(table(seroprevalence_join9$successes, seroprevalence_join9$parabola)))
# x<-x%>%
#   spread(Var2, Freq) 
# x %>%
#   mutate(p = `1`/`0`) %>%
#   ggplot(aes(x=Var1, y= p, size=`1`+`0`)) +
#   geom_point()
# 
# x<-as.data.frame(t(table(seroprevalence_join9$successes, seroprevalence_join9$month.dc.birthpulse_original)))
# x<-x%>%
#   spread(Var2, Freq) 
# 
# x %>%
#   mutate(p = `1`/`0`) %>%
#   ggplot(aes(x=Var1, y= p, size=`1`+`0`)) +
#   geom_point()
# 
# y <-seroprevalence_join7 %>%
#   filter(month.dc.birthpulse_original==0)
# 
# #..........................trying for all  data.........................................
# 

# seroprevalence_join9 <- seroprevalence_join8 %>%
#   mutate(parabola = month.dc.birthpulse_original * (12-month.dc.birthpulse_original))
# 
# seroprevalence_join9a <- seroprevalence_join8a %>%
#   mutate(parabola = month.dc.birthpulse_original * (12-month.dc.birthpulse_original))
# # # 
# library(lme4)
# # #first test
# glmer(successes ~ parabola + (1|substudy) + (1|title), family=binomial(link='logit'), data= seroprevalence_join9) %>% summary()

# # x<-as.data.frame(t(table(seroprevalence_join9$successes, seroprevalence_join9$parabola)))
# # x<-x%>%
# #   spread(Var2, Freq)
# # x%>%
# #   mutate(p = `1`/`0`) %>%
# #   ggplot(aes(x=Var1, y= p, size=`1`+`0`)) +
# #   geom_point()
#
# x<-as.data.frame(t(table(seroprevalence_join9$successes, seroprevalence_join9$month.dc.birthpulse_original)))
# x<-x%>%
#   spread(Var2, Freq) 
# x %>%
#   mutate(p = `1`/`0`) %>%
#   ggplot(aes(x=Var1, y= p, size=`1`+`0`)) +
#   geom_point()
# 
# x<- seroprevalence_join3 %>%
#   mutate(seroprevalence_percentage = round(seroprevalence_percentage)) %>%
#   group_by(methodology) %>%
#   mutate(sum = n()) %>%
#   ungroup() %>%
#   group_by(methodology, seroprevalence_percentage, sum) %>%
#   summarise(n=n()) %>%
#   mutate(n=n/sum) %>%
#   filter(seroprevalence_percentage == 0)
# 
# seroprevalence_join3 %>%
#   ggplot(aes(seroprevalence_percentage)) +
#   geom_histogram() + 
#   facet_grid(~methodology_kw) 


