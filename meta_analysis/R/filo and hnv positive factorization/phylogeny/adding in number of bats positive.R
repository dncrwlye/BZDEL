#this is going to tack on the number of bats positive per virus onto becker's batphy file 

tryCatch(setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis"),
         error=function(e) setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/"))

library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
require(pscl)
require(boot)
library(pscl)

load(file = "data/batphy")
load('data/seroprevalence.Rdata')

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
data=data[c("species","virus","detection","sample", "serop")]

data$species=revalue(data$species,c("rhinolophus refulgens"="rhinolophus lepidus",
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
data$species=revalue(data$species,c("lissonycteris angolensis"="myonycteris angolensis",
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
  dplyr:: mutate(x  = (serop * sample))%>%
  group_by(species) %>%
  summarise(sum = sum(x))
  
hist(data1$sum, breaks = seq(0,400,by=5))


x <- data$serop * data$sample
mean(x, na.rm=TRUE)
var(x, na.rm=TRUE)
hist(x)

