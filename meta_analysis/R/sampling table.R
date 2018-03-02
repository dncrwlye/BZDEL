
#.............ehhhh gonna make a table describing the amount of longitudinal studies.........
library(tidyverse)
library(stringi)
library(readxl)
library(tidyverse)

setwd('/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis')
load(file='data/seroprevalence.Rdata')


seroprevalence.hnv <- seroprevalence %>% filter(virus == 'Henipavirus')
#seroprevalence.filo <- seroprevalence %>% filter(virus == 'Filovirus')

#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................
#................................henipavirus.........................................................

seroprevalence.hnv %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 1 | date_diff < 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling, end_of_sampling)%>%
  unique() %>%
  group_by(substudy_non_annual) %>%
  summarise(counts= n()) %>%
  nrow()

explicit_longitudinal.hnv.4<-seroprevalence.hnv %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 1 | date_diff < 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling, end_of_sampling)%>%
  unique() %>%
  group_by(substudy_non_annual) %>%
  summarise(counts= n()) %>%
  filter(counts >= 4) %>%
  ungroup() %>%
  unique()

explicit_longitudinal.hnv.3<-seroprevalence.hnv %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 1 | date_diff < 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling, end_of_sampling)%>%
  unique() %>%
  group_by(substudy_non_annual) %>%
  summarise(counts= n()) %>%
  filter(counts == 3) %>%
  ungroup() %>%
  unique()

explicit_longitudinal.hnv.2<-seroprevalence.hnv %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 1 | date_diff < 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling, end_of_sampling)%>%
  unique() %>%
  group_by(substudy_non_annual) %>%
  summarise(counts= n()) %>%
  filter(counts == 2) %>%
  ungroup() %>%
  unique()

single_time_points_but_decent_range.hnv<-seroprevalence.hnv %>%
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

nrow(single_time_points_but_decent_range.hnv) + nrow(explicit_longitudinal.hnv.2)+ nrow(explicit_longitudinal.hnv.3)+ nrow(explicit_longitudinal.hnv.4)

pooled_estimates_just_horrible.hnv<-seroprevalence.hnv %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 0 & date_diff >= 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling, end_of_sampling, date_diff) %>%
  unique()
  # group_by(substudy_non_annual, sampling_date_single_time_point, start_of_sampling, end_of_sampling, date_diff) %>%
  # summarise(counts= n()) 
  
nrow(single_time_points_but_decent_range.hnv) + nrow(explicit_longitudinal.hnv.2) + nrow(explicit_longitudinal.hnv.3)+ nrow(explicit_longitudinal.hnv.4) + nrow(pooled_estimates_just_horrible.hnv)
length(unique(seroprevalence.hnv$substudy_non_annual))

mean(pooled_estimates_just_horrible$date_diff)
sd(pooled_estimates_just_horrible$date_diff)


#..........................filovirus....................................
#..........................filovirus....................................
#..........................filovirus....................................
#..........................filovirus....................................
#..........................filovirus....................................
#..........................filovirus....................................
#..........................filovirus....................................
#..........................filovirus....................................
#..........................filovirus....................................
#..........................filovirus....................................


seroprevalence.filo %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 1 | date_diff < 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling, end_of_sampling)%>%
  unique() %>%
  group_by(substudy_non_annual) %>%
  summarise(counts= n()) %>%
  nrow()

explicit_longitudinal.filo.4<-seroprevalence.filo %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 1 | date_diff < 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling, end_of_sampling)%>%
  unique() %>%
  group_by(substudy_non_annual) %>%
  summarise(counts= n()) %>%
  filter(counts >= 4) %>%
  ungroup() %>%
  unique()

explicit_longitudinal.filo.3<-seroprevalence.filo %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 1 | date_diff < 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling, end_of_sampling)%>%
  unique() %>%
  group_by(substudy_non_annual) %>%
  summarise(counts= n()) %>%
  filter(counts == 3) %>%
  ungroup() %>%
  unique()

explicit_longitudinal.filo.2<-seroprevalence.filo %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 1 | date_diff < 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling, end_of_sampling)%>%
  unique() %>%
  group_by(substudy_non_annual) %>%
  summarise(counts= n()) %>%
  filter(counts == 2) %>%
  ungroup() %>%
  unique()

single_time_points_but_decent_range.filo<-seroprevalence.filo %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 1 | date_diff < 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling) %>%
  unique() %>%
  group_by(substudy_non_annual) %>%
  summarise(counts= n()) %>%
  filter(counts == 1) %>%
  ungroup() %>%
  unique()

nrow(single_time_points_but_decent_range.filo) + nrow(explicit_longitudinal.filo.2) + nrow(explicit_longitudinal.filo.3) + nrow(explicit_longitudinal.filo.4)

pooled_estimates_just_horrible.filo<-seroprevalence.filo %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 0 & date_diff >= 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling, end_of_sampling, date_diff) %>%
  unique() %>%
  group_by(substudy_non_annual, date_diff) %>%
  summarise(counts= n()) %>%
  ungroup() %>%
  unique()




single_time_points_but_decent_range.filo$substudy_non_annual %in% pooled_estimates_just_horrible.filo$substudy_non_annual
explicit_longitudinal.filo$substudy_non_annual %in% pooled_estimates_just_horrible.filo$substudy_non_annual

x<- explicit_longitudinal.filo %>%
  filter(substudy_non_annual %in% pooled_estimates_just_horrible.filo$substudy_non_annual)
  
y <- seroprevalence.filo %>%
  filter(substudy_non_annual %in% x$substudy_non_annual)
  length(unique(seroprevalence.filo$substudy_non_annual))


paper_breakdown_table <- matrix(c(
  "","N","%","Mean","SD",
  "Total Sampling Sub Units", length(unique(seroprevalence$substudy_non_annual)), "","","",
  "Explicitly Longitudinal", nrow(explicit_longitudinal), round(nrow(explicit_longitudinal)/length(unique(seroprevalence$substudy_non_annual)),3)*100,"","",
  "Single Time Point", nrow(single_time_points_but_decent_range), round(nrow(single_time_points_but_decent_range)/length(unique(seroprevalence$substudy_non_annual)),3)*100,"","",
  "Study Longitudinal, Reported Single Time Point", nrow(pooled_estimates_just_horrible), round(nrow(pooled_estimates_just_horrible)/length(unique(seroprevalence$substudy_non_annual)),3)*100,round(mean(pooled_estimates_just_horrible$date_diff),1),round(sd(pooled_estimates_just_horrible$date_diff),1)
  ), byrow=TRUE,ncol=5)
paper_breakdown_table                                  

kable(paper_breakdown_table, align = 'l', col.names = c('', 'Count', '%'))
x
cat(x)






