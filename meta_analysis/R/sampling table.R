
#.............ehhhh gonna make a table describing the amount of longitudinal studies.........
library(tidyverse)
library(stringi)
library(readxl)
library(tidyverse)

length(unique(seroprevalence$title))
length(unique(seroprevalence$substudy_non_annual))

explicit_longitudinal<-seroprevalence %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 1 | date_diff < 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling) %>%
  unique() %>%
  group_by(substudy_non_annual) %>%
  summarise(counts= n()) %>%
  filter(counts >= 2) %>%
  ungroup() %>%
  #dplyr::select(title) %>%
  unique()

single_time_points_but_decent_range<-seroprevalence %>%
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

pooled_estimates_just_horrible<-seroprevalence %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 0 & date_diff >= 365/12) %>%
  dplyr::select(substudy_non_annual, sampling_date_single_time_point, start_of_sampling, date_diff) %>%
  unique() %>%
  group_by(substudy_non_annual, date_diff) %>%
  summarise(counts= n()) %>%
  ungroup() %>%
  unique()

mean(pooled_estimates_just_horrible$date_diff)
sd(pooled_estimates_just_horrible$date_diff)

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






