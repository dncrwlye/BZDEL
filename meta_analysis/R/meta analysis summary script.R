library(tidyverse)
library(stringi)
library(readxl)
library(binom)
library(lubridate)
MetaAnalysis_Data_New_Version_ <- read_excel("~/Dropbox_gmail/Dropbox/bat virus meta-analysis/MetaAnalysis Data New Version .xlsx", 
                                           col_types = c("text", "numeric", "text", 
                                                                "text", "text", "text", "text", "text", 
                                                                "text", "text", "text", "text", "text", 
                                                                "text", "numeric", "text", "text", 
                                                                "text", "text", "text", "text", "text", 
                                                                "numeric", "text", "text", "text", 
                                                                "text", "text", "text", "date", "date", 
                                                                "date", "text", "text", "text", "text", 
                                                                "numeric", "text"))



seroprevalence <- MetaAnalysis_Data_New_Version_ %>%
  filter(Outcome == 'Seroprevalence') %>%
  select(Title, First_Author_last_name, Virus, Study, study_design, Methodology, Species, Sex, `Age Class`, `Sampling Location`, `Sample Size`, `Seroprevalence (%)`, `Single Time Point Sampling  [True/False]`, `Sampling Date (Single Time Point)`, `Start of Sampling`, `End of Sampling`) %>%
  mutate(Virus = ifelse((Virus == "Ebola" | 
                        Virus == "Marburg" | 
                        Virus == "Sudan Virus" | 
                        Virus == "Zaire Ebolavirus" | 
                        Virus == "Reston Ebola"), "Filovirus", Virus)) %>%
  mutate(global_location = ifelse(grepl('Thailand|Malaysia|China|Cambodia', `Sampling Location`), 'East Asia', 
                           ifelse(grepl('Brazil', `Sampling Location`), 'South America', 
                           ifelse(grepl('Bangladesh|India',`Sampling Location`), 'Central Asia', 
                           ifelse(grepl('Australia|Indonesia|Phillippines', `Sampling Location`), "Austranesia",
                           ifelse(grepl('Ghana|DRC|Gabon|Annobon|Uganda|Congo|Zambia|Madagascar',`Sampling Location`), "sub_saharan_africa", NA))))))    
    

serology_graph_function <- function(Virus_of_Interest)
{
  seroprevalence_henipavirus<- seroprevalence %>%
    filter(Virus == Virus_of_Interest) %>%
    mutate(seroprevalence_per = as.numeric(`Seroprevalence (%)`)) %>%
    filter(!is.na(seroprevalence_per)) %>%
    mutate(seroprevalence_per = round(seroprevalence_per, 2)) %>%
    mutate(sample_size = as.numeric(`Sample Size`)) %>%
    mutate(successes = round((seroprevalence_per/100) * sample_size,0)) %>%
    filter(!is.na(successes)) 
  
  seroprevalence_henipavirus_single_time_point <- seroprevalence_henipavirus %>%
    filter(`Single Time Point Sampling  [True/False]` == 1) %>%
    mutate(month = as.numeric(format(`Sampling Date (Single Time Point)`, "%m"))) %>%
    group_by(month, Title, First_Author_last_name, global_location) %>%
    summarise(sample_size = sum(sample_size), successes = sum (successes))  %>%
    mutate(seroprevalence_per = 100 * (successes/sample_size)) %>%
    mutate(seroprevalence_per_upper_bound = 100*binom.exact(successes, sample_size, conf.level = .95)$upper) %>%
    mutate(seroprevalence_per_lower_bound = 100*binom.exact(successes, sample_size,conf.level = .95)$lower) %>%
    ungroup() 
  
  seroprevalence_henipavirus_unclear_time_point <- seroprevalence_henipavirus %>%
    filter(`Single Time Point Sampling  [True/False]` == 0) %>%
    filter(!is.na(`Start of Sampling`)) %>%
    mutate(difftime = as.numeric(difftime(`End of Sampling`, `Start of Sampling`))) %>%
    mutate(difftime = ifelse(difftime > 1000, 1000, difftime)) %>%
    mutate(difftime_sqrt = sqrt(difftime)) %>% #going to divide any long sampling period by the square root of the sampling period
    mutate(month_1 = ifelse(difftime < 364, as.numeric(format(`Start of Sampling`, "%m")), 1)) %>%
    mutate(month_2 = ifelse(difftime < 364, as.numeric(format(`End of Sampling`, "%m")), 12)) %>%
    group_by(month_1, month_2, difftime, difftime_sqrt, Title, First_Author_last_name, global_location) %>%
    summarise(sample_size = sum(sample_size), successes = sum (successes))%>%
    mutate(seroprevalence_per = 100 * (successes/sample_size)) %>%
    mutate(seroprevalence_per_upper_bound = 100*binom.exact(successes, sample_size, conf.level = .95)$upper) %>%
    mutate(seroprevalence_per_lower_bound = 100*binom.exact(successes, sample_size,conf.level = .95)$lower) %>%
    mutate(seroprevalence_per = ifelse(difftime > 32, seroprevalence_per/difftime_sqrt, seroprevalence_per)) %>%
    mutate(seroprevalence_per_upper_bound = ifelse(difftime > 32, seroprevalence_per_upper_bound/difftime_sqrt, seroprevalence_per_upper_bound)) %>%
    mutate(seroprevalence_per_lower_bound = ifelse(difftime > 32, seroprevalence_per_lower_bound/difftime_sqrt, seroprevalence_per_lower_bound)) %>%
    ungroup() %>%
    select(-c(difftime, difftime_sqrt)) %>%
    mutate(Title = paste(Title, "xc"))
  
  x <- seroprevalence_henipavirus_unclear_time_point[,c(2:10)]
  y <-seroprevalence_henipavirus_unclear_time_point[,c(1,3:10)]
  colnames(y) <- colnames(x)
  seroprevalence_henipavirus_unclear_time_point <- rbind(y,x)
  colnames(seroprevalence_henipavirus_unclear_time_point) <- colnames(seroprevalence_henipavirus_single_time_point)
  
  seroprevalence_graph <- rbind(seroprevalence_henipavirus_unclear_time_point, seroprevalence_henipavirus_single_time_point)
  
  graph <- seroprevalence_graph %>%
    ungroup() %>%
    ggplot() +
    geom_point(aes(x= month, y= seroprevalence_per, colour = Title)) +
    geom_line(lwd=1, aes(x= month, y=seroprevalence_per, colour = Title)) +
    geom_ribbon(alpha = .2, aes(x= month, 
                                ymin=seroprevalence_per_lower_bound, 
                                ymax=seroprevalence_per_upper_bound,
                                fill= Title)) +
    theme(legend.position="none") +
    ylab('seroprevalence') +
    ggtitle(paste(c(Virus_of_Interest, "Seroprevalence; Data from Meta Analysis"))) +
    facet_wrap(~global_location)
    
  
    return(graph)
}


Nipah <- serology_graph_function('Nipah')
Hendra <- serology_graph_function('Hendra')
Filovirus <- serology_graph_function('Filovirus')

