#this is going to tack on the number of bats positive for NIPAH ONLY onto becker's batphy file 
#this for that images of indians bats, boy this is a lot of work! 

tryCatch(setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis"),
         error=function(e) setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/"))

library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
require(pscl)
require(boot)
library(pscl)

load(file = 'data/batphy.Rdata')
load('data/seroprevalence.Rdata')

x <- filter(seroprevalence) %>%
  filter(species == "miniopterus schreibersii")
data <- seroprevalence %>%
  filter(sample_size < 9999)
rm(seroprevalence)

## make seroprevalence accessible
data$serop=data$seroprevalence_percentage
data$serop=as.numeric(data$serop)/100

## clean
data$seroprevalence_percentage=NULL

## sample size accessible
data$sample=data$sample_size
data$sample=as.numeric(data$sample)

## make detection
data$detection=ifelse(data$serop>0,1,0)

## if methodology==isolation
data$detection=ifelse(data$methodology=="isolation",1,data$detection)

## trim data down to just species, detection, and virus
data=data[c("species","virusold","detection","sample", "serop")]

data$species=plyr::revalue(data$species,c("rhinolophus refulgens"="rhinolophus lepidus",
                                          "hipposideros cf caffer"="hipposideros caffer",
                                          "hipposideros cf caffer/ruber"="hipposideros caffer",
                                          "hipposideros cf ruber"="hipposideros ruber",
                                          "natalus lanatus"="natalus mexicanus",
                                          "myonycteris leptodon"="myonycteris torquata",
                                          "nanonycteris veldkampiii"="nanonycteris veldkampii",
                                          "tadarida plicata"="chaerephon plicatus",
                                          "chaerephon plicata"="chaerephon plicatus"))

## remove unapplicable names
data=data[-which(data$species=="bat"),]
data=data[-which(data$species=="vespertilioniformes species"),]
data=data[-which(data$species=="vespertilionidae species1"),]
data=data[-which(data$species=="tadarida species"),]
data=data[-which(data$species=="synconycterus species"),]
data=data[-which(data$species=="scotorepens species"),]
data=data[-which(data$species=="rhinolophus species"),]
data=data[-which(data$species=="pteropus species"),]
data=data[-which(data$species=="nyctophilus species"),]
data=data[-which(data$species=="myotis species"),]
data=data[-which(data$species=="mops species"),]
data=data[-which(data$species=="pipistrellus species"),]
data=data[-which(data$species=="molossideo species1"),]
data=data[-which(data$species=="miniopterus species"),]
data=data[-which(data$species=="micropteropus/nanonycteris"),]
data=data[-which(data$species=="microchiroptere species"),]
data=data[-which(data$species=="hipposideros species"),]
data=data[-which(data$species=="glossophaginae species3"),]
data=data[-which(data$species=="glossophaginae species1"),]
data=data[-which(data$species=="epomophorus species"),]
data=data[-which(data$species=="eonycteris spelaea & rousettus species"),]
data=data[-which(data$species=="cynopterus species"),]
data=data[-which(data$species=="chalinobus species"),]
data=data[-which(data$species=="carollia species"),]
data=data[-which(data$species=="bats*"),]
data=data[-which(data$species=="nycteris species"),]
data=data[-which(data$species=="hypsignathus monstrosus, epomops franqueti, myonycteris torquata"),]
data=data[-which(data$species=="pteropus alecto, pteropus poliocephalus"),]
data=data[-which(data$species=="pteropus alecto, pteropus poliocephalus, pteropus scapulatus"),]
data=data[-which(data$species=="pteropus alecto, pteropus scapulatus"),]
data=data[-which(data$species=="pteropus conspicillatus, pteropus scapulatus"),]
data=data[-which(data$species=="pteropus conspicillatus & dobsonia magna"),]
data=data[-which(data$species=="unknown species"),]

## fix names to match our phylogeny
data$species=plyr::revalue(data$species,c("lissonycteris angolensis"="myonycteris angolensis",
                                          "nyctalus plancyi"="Nnctalus velutinus",
                                          "rhinolophus stheno"="rhinolophus microglobosus",
                                          "natalus stramineus"="natalus stramineus mexicanus",
                                          "neoromicia somalicus"="neoromicia malagasyensis",
                                          "scotophilus heathi"="scotophilus heathii",
                                          "pipistrellus cf nanus/nanulus" = "neoromicia nanus",
                                          "artibeus planirostris" = "artibeus jamaicensis",
                                          "enchisthenes hartii" = "artibeus hartii",
                                          "myotis capaccini" = "myotis capaccinii",
                                          "peroptery kappleri" = "peropteryx kappleri",
                                          "rhinolphus mehelyi" =  "rhinolophus mehelyi",
                                          "casinycteris ophiodon" = "scotonycteris ophiodon",
                                          "hipposideros pygmeus" = "hipposideros pygmaeus",
                                          "miniopterus magnate" = "miniopterus magnater",
                                          "natalus mexicanus" = "natalus stramineus mexicanus",
                                          "taphozous saccolaimus" = "saccolaimus saccolaimus"
))

x <- data[!(data$species %in% batphy$species),]
x <- x$species %>% unique() %>% as.data.frame()

data1 <- data %>% 
  dplyr:: mutate(number_pos = (serop * sample)) %>%
  group_by(species, virusold) %>%
  #dplyr::summarise(total_pos = sum(number_pos, na.rm=TRUE), total_sampled = sum(sample, na.rm=TRUE)) %>%
  dplyr::summarise(total_pos = sum(number_pos, na.rm=TRUE)) %>%
  spread(virusold, total_pos)

data1 <- data1 %>%
  select(species, Nipah) %>%
  mutate(Nipah = ifelse(Nipah > 0, 1, 0))

batphy <- left_join(batphy, data1)

save(batphy, file = 'data/batphy.Rdata')
rm(list=ls())


