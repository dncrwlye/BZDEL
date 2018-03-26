#............................script to prepare data for phylofactor...........................
library(stringr)
load('data/seroprevalence.Rdata')

seroprevalence.multi.species <- seroprevalence %>%
  filter(str_detect(species, '([a-z]+ [a-z]+),|([a-z]+ [a-z]+) &')==TRUE) %>%
  mutate(species = stri_split_regex(species,',|&')) %>%
  unnest(species) %>%
  mutate(species = trimws(species)) %>% 
  mutate(original.multiple.species = TRUE)

seroprevalence.single.species <- seroprevalence %>%
  filter(str_detect(species, '([a-z]+ [a-z]+),|([a-z]+ [a-z]+) &')==FALSE) %>% 
  mutate(original.multiple.species = FALSE)

seroprevalence.phylo.analysis <- rbind(seroprevalence.multi.species, seroprevalence.single.species) #%>%
  #filter(str_detect(species, " species")==FALSE) 
  
unique(seroprevalence.phylofactor$species) %>% as.vector() %>% sort()

save(seroprevalence.phylo.analysis, file = 'data/seroprevalence.phylo.analysis.Rdata')


load(file = 'data/seroprevalence.phylo.analysis.Rdata')
