#.................. sampling scheme................
library(ggplot2)
library(magrittr)
library(knitr)
library(dplyr)


x<-seq(0,365,by=1)
x <- c(rep(-2, 60), x[1:100], rep(-2, 200), x[150:250], rep(-2, 60))

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

longitudinal_sampling_point <-  smooth_yt[seq(1, nrow(smooth_yt)/2, 30),]
single_point <-  smooth_yt[275,]
longitudinal_pooled_sampling_point <-  smooth_yt[sample((nrow(smooth_yt)/2):nrow(smooth_yt), size=6),]

longitudinal_pooled_sampling_point_box <- data.frame(
  xmin = min(longitudinal_pooled_sampling_point$x),
  xmax = max(longitudinal_pooled_sampling_point$x),
  ymin = mean(longitudinal_pooled_sampling_point$V1)-var(longitudinal_pooled_sampling_point$V1),
  ymax = mean(longitudinal_pooled_sampling_point$V1+var(longitudinal_pooled_sampling_point$V1))
)
ggplot() +
  # geom_rect(data=rect_sampling_point_long, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
  #           fill='green', alpha=0.4) +
  geom_line(aes(longitudinal_sampling_point$x,longitudinal_sampling_point$V1), col= 'blue', alpha=0.5, size = 1) + 
  geom_line(data=smooth_yt, aes(y=V1, x=x), col='orange') +
  ylab('') +
  xlab('') +
  geom_point(aes(longitudinal_sampling_point$x,longitudinal_sampling_point$V1), col= 'blue', alpha=0.8, size = 5) + 
  geom_point(aes(longitudinal_pooled_sampling_point$x,longitudinal_pooled_sampling_point$V1), col= 'green', alpha=0.8, size = 5) + 
  geom_point(aes(single_point$x,single_point$V1), col= 'black', alpha=0.8, size = 5) + 
  geom_rect(data=longitudinal_pooled_sampling_point_box, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
  fill='green', alpha=0.4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

  
  
  
  #+
  # labs(caption = "Sampling schemes present in published studies of henipavirus/filoviru sereprevalence in bats.
  #                 Black represents a single sampling event.
  #                 Blue represents a longitudinal sampling at bimonthly intervals. 
  #                 Green represents a potentially longitudinal sampling, but the study reports only a single mean estimate of seroprevalence for a given sample period.", parse = TRUE)

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











