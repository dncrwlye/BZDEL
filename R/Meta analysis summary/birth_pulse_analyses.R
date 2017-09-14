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


filovirus.seroprevalence<-seroprevalence %>%
  filter(Virus=='Filovirus') %>%
  mutate(`Sampling Location` = tolower(`Sampling Location`)) %>%
  mutate(location.birth.pulse = stri_extract_first_regex(`Sampling Location`, '[a-z]+'))
