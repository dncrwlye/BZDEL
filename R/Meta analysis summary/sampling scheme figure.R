#.................. sampling scheme................
library(ggplot2)
library(magrittr)
library(knitr)
library(dplyr)


x<-seq(0,365,by=1)
amp <- rnorm(n=length(x), mean=.2)
amp <- runif(n=length(x))
y<-amp*sin(.04*x)+amp
y<-y*sd(y)
plot(y, type='l')

smooth_y <- vector()
start = 20
exp_value =1.1
for (i in seq(start,length(y)))
{
l <- i-(start-1)
vector_exp <- 1/(exp_value^(2:(start+1)))
smooth_y[i] <- sum(vector_exp * y[i:l]) * (1/sum(vector_exp))
}
plot(y, type='l')
lines(smooth_y, col='red', lwd=5)

smooth_yt <- as.data.frame((as.matrix(smooth_y))) %>%
  mutate(x = row_number())

rect_sampling_point_unique <- c(30, 30*3,  30*5,  30*7,  30*9, 30*11)
rect_sampling_point_unique <- data.frame(
  xmin = rect_sampling_point_unique,
  xmax = rect_sampling_point_unique + 8,
  ymin = 0,
  ymax = 1
)
rect_sampling_point_long <- c(40)
rect_sampling_point_long <- data.frame(
  xmin = rect_sampling_point_long,
  xmax = rect_sampling_point_long + 200,
  ymin = .2,
  ymax = .6
)

singe_time_point <- c(300)
singe_time_point <- data.frame(
  xmin = singe_time_point,
  xmax = singe_time_point + 10,
  ymin = 0,
  ymax = 1
)

ggplot() +
  geom_rect(data=rect_sampling_point_unique, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
            fill='blue', alpha=0.2) +
  geom_rect(data=rect_sampling_point_long, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
            fill='green', alpha=0.4) +
  geom_rect(data=singe_time_point, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
            fill='black', alpha=0.4) +
  geom_line(data=smooth_yt, aes(y=V1, x=x, col='seroprevalence')) +
  ylab('seroprevalence') +
  xlab('julian days') +
  #ggtitle('sampling schemes of henipavirus/filovirus seroprevalence studies present in published literature') +
  labs(caption = "Sampling schemes present in published studies of henipavirus/filoviru sereprevalence in bats.
                  Black represents a single sampling event.
                  Blue represents a longitudinal sampling at bimonthly intervals. 
                  Green represents a potentially longitudinal sampling, but the study reports only a single mean estimate of seroprevalence for a given sample period.", parse = TRUE)


#.............ehhhh gonna make a table describing the amount of longitudinal studies.........
library(tidyverse)
library(stringi)
library(readxl)
MetaAnalysis_Data_New_Version <- read_excel("~/Dropbox_gmail/Dropbox/bat virus meta-analysis/MetaAnalysis Data New Version.xlsx", 
                                            col_types = c("text", "numeric", "text", 
                                                          "text", "text", "text", "text", "text", 
                                                          "text", "text", "text", "text", "text", 
                                                          "text", "text", "text", "text", "numeric", "numeric", 
                                                          "numeric", "numeric", "numeric", 
                                                          "numeric", "numeric", "numeric", 
                                                          "numeric", "text", "numeric", "numeric", 
                                                          "numeric", "numeric", "date", "date", 
                                                          "date", "text", "text", "text", "text", 
                                                          "numeric"))

seroprevalence <- MetaAnalysis_Data_New_Version %>%
  filter(outcome == 'Seroprevalence') %>%
  dplyr::select(title, last_name_of_first_author, virus, study_type, study_design, methodology, species, sex, age_class, sampling_location, sampling_location_two, sample_size, seroprevalence_percentage, single_sampling_point, sampling_date_single_time_point, start_of_sampling, end_of_sampling) %>%
  mutate(virus = ifelse((virus == "Ebola" | 
                           virus == "Marburg" | 
                           virus == "Zaire Ebola"|
                           virus == "Sudan virus" | 
                           virus == "Zaire Ebolavirus" | 
                           virus == "Zaire Ebola "|
                           virus == "Reston Ebola"), "Filovirus", 
                        ifelse(virus  == "Nipah" | virus == "Hendra" | virus  == "Henipavirus", "Henipavirus", virus))) %>%
  filter(virus != 'Tioman') %>%
  mutate(sampling_location = tolower(sampling_location)) %>%
  mutate(country = stri_extract_first_regex(sampling_location, '[a-z]+')) %>%
  mutate(species = tolower(species)) %>%
  mutate(species = trimws(species)) %>%
  mutate(species = stri_extract_first_regex(species, '[a-z]+ [a-z]+')) %>%
  mutate(methodology = ifelse(methodology == 'PCR'| methodology == 'RT-PCR'| methodology == "RT-PCR  (urine)"| methodology == "RT-PCR (Oro-pharangyeal swab)", 'PCR based method', 
                              ifelse(methodology == "ELISA" | methodology == "ELISA + WB" | methodology == "Luminex", 'non nAb based method',
                                     ifelse(methodology == "VNT"| methodology == "SNT"| methodology == "Unclear, Presumably Neutralizing Antibodies", "nAb based method", methodology)))) %>%
  mutate(successes = round((1/100)*seroprevalence_percentage * sample_size,0)) 

seroprevalence_total <- seroprevalence %>%
  dplyr::select(title) %>% 
  unique()

explicit_longitudinal<-seroprevalence %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 1 | date_diff < 365/6) %>%
  select(title, sampling_date_single_time_point, start_of_sampling, sampling_location, species) %>%
  unique() %>%
  group_by(title, sampling_location, species) %>%
  summarise(counts= n()) %>%
  filter(counts >= 2) %>%
  ungroup() %>%
  select(title) %>%
  unique()
 
single_time_points_but_decent_range<-seroprevalence %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 1 | date_diff < 365/6) %>%
  select(title, sampling_date_single_time_point, start_of_sampling, sampling_location, species) %>%
  unique() %>%
  group_by(title, sampling_location, species) %>%
  summarise(counts= n()) %>%
  filter(counts == 1) %>%
  ungroup() %>%
  select(title) %>%
  unique()

pooled_estimates_just_horrible<-seroprevalence %>%
  mutate(date_diff = end_of_sampling - start_of_sampling) %>%
  filter(single_sampling_point == 0 & date_diff >= 365/6) %>%
  select(title, sampling_date_single_time_point, date_diff, start_of_sampling, sampling_location, species) %>%
  unique() %>%
  group_by(title, sampling_location, species, date_diff) %>%
  summarise(counts= n()) %>%
  ungroup() %>%
  select(title) %>%
  unique()

11+16+15

unique(seroprevalence$title)


paper_breakdown_table <- matrix(c(
                                  "Total Studies", nrow(seroprevalence_total), "-",
                                  "Explicitly Longitudinal", nrow(explicit_longitudinal), round(nrow(explicit_longitudinal)/nrow(seroprevalence_total),3)*100,
                                  "Single Time Point", nrow(single_time_points_but_decent_range), round(nrow(single_time_points_but_decent_range)/nrow(seroprevalence_total),3)*100,
                                  "Study Longitudinal, Reported Single Time Point", nrow(pooled_estimates_just_horrible), round(nrow(pooled_estimates_just_horrible)/nrow(seroprevalence_total),3)*100)
                                  , byrow=TRUE,ncol=3)
paper_breakdown_table                                  
                                  
                                  
kable(paper_breakdown_table, align = 'l', col.names = c('', 'Count', '%'))
x
cat(x)











