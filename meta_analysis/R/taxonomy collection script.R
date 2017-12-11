#........................round 1 of taxonomyfying...............................................

library(stringi)
library(tidyverse)
library(taxize)
#.......... clean bat phylogeny

batphy_unique_species <- batphy$species

batphy <- batphy %>%
  mutate(filo_surv = ifelse(is.na(filo_surv), 0, filo_surv)) %>%
  mutate(species = as.character(species)) %>%
  mutate(species = ifelse(species == "triaenops auritus", "paratriaenops auritus",
                          ifelse(species == "artibeus incomitatus", "dermanura incomitata", species)))



#classifications <- classification((batphy_unique_species), db = 'itis', rows=1)

m <- matrix(0,1105,2) %>% data.frame()
for(i in 1:1105)
{
  m[i,1] <- paste(as.data.frame(classifications[i])[1], collapse = "; ")
  print(i)
}

colnames(m) <- c('tax', 'species')
m <- m %>%
  mutate(tax = gsub("c\\(", "", tax)) %>%
  mutate(tax = gsub("\\\"", "", tax)) %>%
  mutate(tax = gsub(",", ";", tax)) %>%
  mutate(species = stri_extract_last_regex(tax,'[A-Z][a-z]+ [a-z]+')) %>%
  mutate(species = tolower(species))

batphy.x <- batphy %>%
  filter(!(species %in% m$species))

#........................round 2 of taxonomyfying...............................................

#classifications_round2 <- classification((batphy.x$species), db = 'itis')

m.r2 <- matrix(0,15,2) %>% data.frame()
for(i in 1:15)
{
  m.r2[i,1] <- paste(as.data.frame(classifications_round2[i])[1], collapse = "; ")
  print(i)
}

colnames(m.r2) <- c('tax', 'species')

m.r2 <- m.r2 %>%
  mutate(tax = gsub("c\\(", "", tax)) %>%
  mutate(tax = gsub("\\\"", "", tax)) %>%
  mutate(tax = gsub(",", ";", tax)) %>%
  mutate(species = stri_extract_last_regex(tax,'[A-Z][a-z]+ [a-z]+')) %>%
  mutate(species = tolower(species))

batphy.x <- batphy %>%
  filter(!(species %in% m$species)) %>%
  filter(!(species %in% m.r2$species)) %>%
  mutate(species = as.character(species))

#........................round 3 of taxonomyfying...............................................

#classifications_round3 <- classification((batphy.x$species), db = "wiki", wiki_site = "pedia")

m.r3 <- matrix(0,11,2) %>% data.frame()
for(i in 1:11)
{
  m.r3[i,1] <- paste(as.data.frame(classifications_round3[i])[1], collapse = "; ")
  print(i)
}

colnames(m.r3) <- c('tax', 'species')

m.r3 <- m.r3 %>%
  mutate(tax = gsub("c\\(", "", tax)) %>%
  mutate(tax = gsub("\\\"", "", tax)) %>%
  mutate(tax = gsub(",", ";", tax)) %>%
  mutate(species = stri_extract_last_regex(tax,'[A-Z][a-z]+ [a-z]+')) %>%
  mutate(species = tolower(species))

batphy.x <- batphy %>%
  filter(!(species %in% m$species)) %>%
  filter(!(species %in% m.r2$species)) %>%
  filter(!(species %in% m.r3$species)) %>%
  mutate(species = as.character(species))

#........................round 4 of taxonomyfying...............................................

#classifications_round4 <- classification((batphy.x$species), db = "ncbi")

m.r4 <- matrix(0,8,2) %>% data.frame()
for(i in 1:8)
{
  m.r4[i,1] <- paste(as.data.frame(classifications_round4[i])[1], collapse = "; ")
  print(i)
}

colnames(m.r4) <- c('tax', 'species')

m.r4 <- m.r4 %>%
  mutate(tax = gsub("c\\(", "", tax)) %>%
  mutate(tax = gsub("\\\"", "", tax)) %>%
  mutate(tax = gsub(",", ";", tax)) %>%
  mutate(species = stri_extract_last_regex(tax,'[A-Z][a-z]+ [a-z]+')) %>%
  mutate(species = tolower(species))

batphy.x <- batphy %>%
  filter(!(species %in% m$species)) %>%
  filter(!(species %in% m.r2$species)) %>%
  filter(!(species %in% m.r3$species)) %>%
  filter(!(species %in% m.r4$species)) %>%
  mutate(species = as.character(species))

#........................round 5 of taxonomyfying...............................................

#classifications_round5 <- classification((batphy.x$species), db = "ncbi")

m.r5 <- matrix(0,3,2) %>% data.frame()
for(i in 1:3)
{
  m.r5[i,1] <- paste(as.data.frame(classifications_round5[i])[1], collapse = "; ")
  print(i)
}

colnames(m.r5) <- c('tax', 'species')

m.r5 <- m.r5 %>%
  mutate(tax = gsub("c\\(", "", tax)) %>%
  mutate(tax = gsub("\\\"", "", tax)) %>%
  mutate(tax = gsub(",", ";", tax)) %>%
  mutate(species = stri_extract_last_regex(tax,'[A-Z][a-z]+ [a-z]+')) %>%
  mutate(species = tolower(species))

batphy.x <- batphy %>%
  filter(!(species %in% m$species)) %>%
  filter(!(species %in% m.r2$species)) %>%
  filter(!(species %in% m.r3$species)) %>%
  filter(!(species %in% m.r4$species)) %>%
  filter(!(species %in% m.r5$species)) %>%
  mutate(species = as.character(species))

batphy <- batphy %>%
  mutate(filo_surv = ifelse(is.na(filo_surv), 0, filo_surv)) %>%
  mutate(species = as.character(species)) %>%
  mutate(species = ifelse(species == "triaenops auritus", "paratriaenops auritus",
                          ifelse(species == "artibeus incomitatus", "dermanura incomitata", species)))

#........................round 6 of taxonomyfying...............................................

#classifications_round6 <- classification((batphy.x$species), db = "wiki")

m.r6 <- matrix(0,1,2) %>% data.frame()
for(i in 1:1)
{
  m.r6[i,1] <- paste(as.data.frame(classifications_round6[i])[1], collapse = "; ")
  print(i)
}

colnames(m.r6) <- c('tax', 'species')

m.r6[1,2] <- "glauconycteris superba"

batphy.x <- batphy %>%
  filter(!(species %in% m$species)) %>%
  filter(!(species %in% m.r2$species)) %>%
  filter(!(species %in% m.r3$species)) %>%
  filter(!(species %in% m.r4$species)) %>%
  filter(!(species %in% m.r5$species)) %>%
  filter(!(species %in% m.r6$species)) %>%
  mutate(species = as.character(species))

batphy <- batphy %>%
  mutate(filo_surv = ifelse(is.na(filo_surv), 0, filo_surv)) %>%
  mutate(species = as.character(species)) %>%
  mutate(species = ifelse(species == "triaenops auritus", "paratriaenops auritus",
                          ifelse(species == "artibeus incomitatus", "dermanura incomitata", species)))

m.rx <- rbind(m,  m.r2, m.r3, m.r4, m.r5, m.r6) %>% unique()
remove <- c('19','1132','1140','951','1112','1130','1138')
m.rx1 <- m.rx[!rownames(m.rx) %in% remove, ]

batphy1 <- left_join(batphy, m.rx1) %>%
  rename(index = X)

batphy1 <- batphy1 %>%
  mutate(tax = gsub("c\\(", "", tax)) %>%
  mutate(tax = gsub("\"", "", tax)) %>%
  mutate(tax = gsub(",", ";", tax))

classifications_Neoromicia_somalica  <- classification('Neoromicia somalica', db = "ncbi")

batphy1[339,'tax'] <- paste(as.data.frame(classifications_Neoromicia_somalica[1])[1], collapse = "; ")

setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")

save(batphy1, file = 'data/bat_taxonomy_data.Rdata')
# 
# 
# inspect <- batphy1 %>%
#   filter(!grepl('Chiroptera', tax))
