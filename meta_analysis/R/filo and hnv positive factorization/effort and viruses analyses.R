setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")

# Libraries ---------------------------------------------------------------

setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
library(Cairo)
par(family="Times")   
CairoFonts(regular = 'Times-12')


load(file='data/phylofactor work spaces/hnv_workspace')
load("data/bat_taxonomy_data.Rdata")

Data.hnv <- Data %>%
  mutate(virus = "hnv")

load(file='data/phylofactor work spaces/filo_workspace')
load("data/bat_taxonomy_data.Rdata")

Data.filo <- Data %>%
  mutate(virus = "filo")

data.filo.hnv.bind <- rbind(Data.hnv, Data.filo)

rm(list=(ls()[ls()!="data.filo.hnv.bind"]))

data.filo.hnv.bind$virus

glm(Z ~virus,family=binomial(link='logit'),data=data.filo.hnv.bind) %>% summary()
glm(Z ~virus + log_effort,family=binomial(link='logit'),data=data.filo.hnv.bind) %>% summary()


# Geography Analysis ----
# ... hnv ----
rm(list=ls())

load(file='data/phylofactor work spaces/hnv_workspace_sample_no_sample_all_bat_dataset')
load("data/bat_taxonomy_data.Rdata")
batphy_for_rotl_update <- read_csv("~/Desktop/BDEL/BZDEL/meta_analysis/data/batphy_for_rotl update.csv")


country.factor.data.frame <- matrix(0,1,8) %>% as.data.frame()
colnames(country.factor.data.frame) <- unique(batphy1$region)

par.reminder = pf$groups[[10]][[2]]

for (i in 1:10)
{
  indexes = pf$groups[[i]][[1]]
  print(indexes %in% par.reminder)
  species.group <- tree$tip.label[indexes]
x <- batphy_for_rotl_update %>% 
    filter(tree_species %in% species.group) %>%
    select(region) %>%
    group_by(region) %>%
    summarise(n = as.numeric(n())) %>%
    ungroup() 
mm <- as.character(unlist(x[,1]))
x = as.data.frame(x[,-1])
x <- x %>%
    t() %>%
    as.data.frame() 
colnames(x)
colnames(x) <- as.character(unlist(mm))
country.factor.data.frame <- dplyr::bind_rows(country.factor.data.frame, x)
}

country.factor.data.frame <- country.factor.data.frame[-1,]
country.factor.data.frame <- country.factor.data.frame[,-5]
country.factor.data.frame <- country.factor.data.frame[,-c(8,9)]
country.factor.data.frame[is.na(country.factor.data.frame)] <- 0
chisq.test(country.factor.data.frame) %>% summary()
str(country.factor.data.frame)


indexes.hi.sampling <- (c(pf$bins[[2]],
  pf$bins[[4]],
  pf$bins[[6]],
  pf$bins[[7]],
  pf$bins[[8]],
  pf$bins[[9]],
  pf$bins[[10]],
  pf$bins[[11]]
  ))

indexes.low.sampling <- (c(pf$bins[[3]],
                          pf$bins[[5]],
                          pf$bins[[1]]))

pf$groups[[4]][[1]] %in% pf$groups[[2]][[1]] 

length(indexes.low.sampling)
length(indexes.hi.sampling)

sum(indexes.low.sampling %in% indexes.hi.sampling)

length(indexes.low.sampling) + length(indexes.hi.sampling)
length(tree$tip.label)

species.hi.sampling <- tree$tip.label[indexes.hi.sampling]
species.low.sampling <- tree$tip.label[indexes.low.sampling]

batphy_for_rotl_update <- batphy_for_rotl_update %>%
  mutate(hi.low.sampling = ifelse(tree_species %in% species.hi.sampling, 1, 
                                  ifelse(tree_species %in% species.low.sampling, 0, NA))) %>%
  filter(!is.na(hi.low.sampling))

glm(hi.low.sampling ~  region,family=binomial(link='logit'),data=batphy_for_rotl_update) %>% summary()
glm(hnv_surv ~  region,family=binomial(link='logit'),data=batphy1) %>% summary()

# ... filo ----
rm(list=ls())

load(file='data/phylofactor work spaces/filo_workspace_sample_no_sample_all_bat_dataset')
load(file="data/bat_taxonomy_data.Rdata")
source('R/taxonomy group name function.R')
batphy_for_rotl_update <- read_csv("~/Desktop/BDEL/BZDEL/meta_analysis/data/batphy_for_rotl update.csv")
taxonomy <- batphy1 %>%
  select(c(species, tax)) 
country.factor.data.frame <- matrix(0,1,8) %>% as.data.frame()
colnames(country.factor.data.frame) <- unique(batphy1$region)

par.reminder = pf$groups[[10]][[2]]

for (i in 1:10)
{
  indexes = pf$groups[[i]][[1]]
  print(indexes %in% par.reminder)
  species.group <- tree$tip.label[indexes]
  x <- batphy_for_rotl_update %>% 
    filter(tree_species %in% species.group) %>%
    select(region) %>%
    group_by(region) %>%
    summarise(n = as.numeric(n())) %>%
    ungroup() 
  mm <- as.character(unlist(x[,1]))
  x = as.data.frame(x[,-1])
  x <- x %>%
    t() %>%
    as.data.frame() 
  colnames(x)
  colnames(x) <- as.character(unlist(mm))
  country.factor.data.frame <- dplyr::bind_rows(country.factor.data.frame, x)
}

country.factor.data.frame <- country.factor.data.frame[-1,]
country.factor.data.frame <- country.factor.data.frame[,-5]
country.factor.data.frame <- country.factor.data.frame[,-c(8,9)]
country.factor.data.frame[is.na(country.factor.data.frame)] <- 0
chisq.test(country.factor.data.frame) %>% summary()
str(country.factor.data.frame)


indexes.hi.sampling <- (c(pf$bins[[3]],
                          pf$bins[[4]],
                          pf$bins[[5]],
                          pf$bins[[7]],
                          pf$bins[[8]],
                          pf$bins[[9]],
                          pf$bins[[10]]
))

indexes.low.sampling <- (c(pf$bins[[2]],
                           pf$bins[[6]],
                           pf$bins[[1]]))

length(indexes.low.sampling)
length(indexes.hi.sampling)

sum(indexes.low.sampling %in% indexes.hi.sampling)

length(indexes.low.sampling) + length(indexes.hi.sampling)
length(tree$tip.label)

species.hi.sampling <- tree$tip.label[indexes.hi.sampling]
species.low.sampling <- tree$tip.label[indexes.low.sampling]

batphy_for_rotl_update <- batphy_for_rotl_update %>%
  mutate(hi.low.sampling = ifelse(tree_species %in% species.hi.sampling, 1, 
                                  ifelse(tree_species %in% species.low.sampling, 0, NA))) %>%
  filter(!is.na(hi.low.sampling))

glm(hi.low.sampling ~  region,family=binomial(link='logit'),data=batphy_for_rotl_update) %>% summary()
glm(filo_surv ~  region,family=binomial(link='logit'),data=batphy1) %>% summary()


#...mas preguntas ----
#there is a consistency between the species that were actually not sampled, and the clades for the filoviruses.
#this wasn't the case for the henipaviruses. 
#what i think this is saying is that the clades for the henipaviruses are picking out clades
#of bats with a diverse geographies, while for the filoviruses there is a strong correlation
#between the geographies at the clades

pf.tree <- pf.tree(pf, lwd=1, factors = 1:10, branch.length = "none", bg.color = NA)

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue') +
  ggtree::geom_label2(aes(label=node), size=2, color="darkred", alpha=0.5) +
  ggtree::geom_hilight(node=1086, fill = 'steelblue')

species.list.1086 <- ggtree::get.offspring.tip(pf$tree, node=1086)
  
batphy_for_rotl_update <- batphy_for_rotl_update %>%
  mutate(node.1086 = ifelse(tree_species %in% species.list.1086, 1, 0))
 
glm(node.1086 ~  region,family=binomial(link='logit'),data=batphy_for_rotl_update) %>% summary()

#yup, looks like 1086 has a lot of latin american bats :P

species.list.1086.tf <- gsub("_", " ", tolower(species.list.1086))


group_taxonomy_list <- as.data.frame(taxonomy[match(species.list.1086.tf,taxonomy[,1]),2])
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            x = group_taxonomy_list$`taxonomy[match(species.list.1086.tf, taxonomy[, 1]), 2]`) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list
group_taxonomy_list %>% table()
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))



