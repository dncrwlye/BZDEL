library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
library(ggtree)
par(family="Times")   

#setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis")
setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
#bat_spp_India <- read_csv("C:/Users/r83c996/Desktop/bat spp India.csv")
bat_spp_India <- read_csv('data/phylofactor work spaces/bat spp India.csv')
bat_spp_Kerala <- read_csv('data/phylofactor work spaces/Bats_Kerala.csv')

load('data/phylofactor work spaces/bat_tree')
load("data/bat_taxonomy_data.Rdata")

bat_spp_India %>%
  filter(!(MSW05_Binomial %in% batphy1$unique_name))

bat_spp_India$MSW05_Binomial =plyr::revalue(bat_spp_India$MSW05_Binomial,c("Myotis_blythii"="Myotis_oxygnathus",
                                                                           "Myotis_peytoni"="Myotis_montivagus"))

bat_spp_India %>%
  filter(!(MSW05_Binomial %in% batphy1$unique_name))

#only Eptesicus_bottae doesn't appear in the database, I can't find another name for this 
#species that is the database either. Also, according to the internets it doesn't actually 
#live in in India so maybe thats not a problem? funky face

Data <- batphy1 %>%
  filter(unique_name %in% bat_spp_India$MSW05_Binomial) %>% #only use bats that are in India, obvi gonna loose a couple 
  mutate(hnv_samps = ifelse(is.na(hnv_samps), 0, hnv_samps)) %>% #remove NAs
  mutate(hnv_positive = ifelse(is.na(hnv_positive), 0, hnv_positive)) %>% #remove NAs
  mutate(log_hnv_samps = log(hnv_samps)) %>% #log transform sampling effort
  select(c(hnv_samps, hnv_positive, species, log_hnv_samps)) %>%
  dplyr::rename(Species = species) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  mutate(Species = stri_trans_totitle(Species)) %>%
  mutate(Sample = 1)

# remove tips from tree that are in our Data file
Data$Species %in% bat_tree$tip.label
tree <- ape::drop.tip(bat_tree,bat_tree$tip.label[!(bat_tree$tip.label %in% Data$Species)])
rm(batphy1, bat_tree)

names(Data) <- c('effort','Z', 'Species','log_effort', 'Sample')

# ncores = 4
# tot.reps=500
# reps.per.worker=round(tot.reps/ncores)
# sampling_effort = TRUE
# 
# source('R/filo and hnv positive factorization/null simulations script.R')
# save(list=ls(),file='data/phylofactor work spaces/indian bats only niv_workspace')

# pf <- gpf(Data,tree,frmla.phylo=Z~phylo,nfactors=10,family=binomial,algorithm='phylo')
# pf$factors
# ncores = 4
# tot.reps=500
# reps.per.worker=round(tot.reps/ncores)
# sampling_effort = FALSE
# source('R/filo and hnv positive factorization/null simulations script.R')
# save(list=ls(),file='data/phylofactor work spaces/indian bats no sampling effort only niv_workspace')

#rm(list=ls())

#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................

#load(file='data/phylofactor work spaces/indian bats no sampling effort only niv_workspace')
bat_spp_India <- bat_spp_India %>%
  mutate(MSW05_Binomial = gsub(" ", "_", MSW05_Binomial))

bat_spp_Kerala <- bat_spp_Kerala %>%
  mutate(MSW05_binomial = gsub(" ", "_", MSW05_binomial))

Data <- Data %>%
  mutate(IB = ifelse(Species %in% bat_spp_India$MSW05_Binomial, 1, 
                     ifelse(!(Species %in% bat_spp_India$MSW05_Binomial), 0, NA))) %>%
  mutate(Kerala = ifelse(Species %in% bat_spp_Kerala$MSW05_binomial, 1,
                         ifelse(!(Species %in% bat_spp_Kerala$MSW05_binomial), 0, NA)))

load('data/seroprevalence.Rdata')

# pcr.pos.hnv <- seroprevalence %>%
#   filter(methodology == 'PCR based method') %>%
#   filter(seroprevalence_percentage > 0) %>%
#   filter(virusold == 'Nipah') %>%
#   select(species) %>%
#   unique()
# 
# neg.hnv <- seroprevalence %>%
#   filter(seroprevalence_percentage == 0) %>%
#   filter(virus == 'Henipavirus') %>%
#   select(species) %>%
#   unique()

Data <- Data %>%
  mutate(species.mutate = tolower(gsub("_", " ", Species))) %>%
  #mutate(all.neg = ifelse(species.mutate %in% neg.hnv$species & effort > 0, 1, 0))%>%
  mutate(all.neg = ifelse(effort > 0 & Z == 0, 1, 0)) %>%
  mutate(effort_cat = ifelse(effort > 0, 1,0))

Data <- Data %>%
  mutate(color = ifelse(Kerala == 1, "black", "grey"))

source('r/stupid_pointless_internal_node_color_function.R')
library('phylobase')
internal.edges<- (internal.edge.color(tree, Data))

row.names(internal.edges) <- internal.edges[,1]
#internal.edges <- internal.edges[,-c(1)] %>% as.data.frame()
g1 = as(tree, 'phylo4')
length(nodeId(g1, "internal")) == length(row.names(internal.edges))

g2 = phylo4d(g1, Data)
nodeData(g2) <- internal.edges

load("data/bat_taxonomy_data.Rdata")

taxonomy <- batphy1 %>%
  select(c(species, tax)) 
source('R/taxonomy group name function.R')

#...............................................................................

ggtree.object <- ggtree::ggtree(g2, layout= 'circular',branch.length = "none",aes(color=I(color))) 

jj <- nrow(Data)

d <- data.frame(x=ggtree.object$data[1:jj,'x'],
                xend=ggtree.object$data[1:jj,'x']  + 1 + Data$log_effort,
                y=ggtree.object$data[1:jj,'y'],
                yend=ggtree.object$data[1:jj,'y'] )

d

d[is.finite(d$xend)==FALSE,] <- 0

#g2 <- ggtree::ggtree(tree,layout='circular') 

ggtree.object$data$label <- gsub("_"," ", ggtree.object$data$label)

i <- .5
ii <- 5
iii <- 2

ggtree.object$data$Z[ggtree.object$data$Z==0]  <- NA
ggtree.object$data$effort_cat[ggtree.object$data$effort_cat==0]  <- NA

ggtree.object + 
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_cladelabel(node=118, label="Vespertilionidae", 
                          color="black", angle=12.16-180, offset=ii, offset.text = iii) +
  ggtree::geom_cladelabel(node=200, label="Pteropodidae", 
                          color="black", angle=-68.91, offset=ii, offset.text = iii) +
  ggtree::geom_tippoint(aes(size=effort_cat), pch =1) +
  ggtree::geom_tippoint(aes(size=Z), color = 'blue', alpha = .8) +
  geom_segment(data= d,aes(x=x + i ,y=y,xend=xend +i ,yend=yend, size= Data$effort_cat,
                           colour = 'orange', alpha = (.5 + Data$Kerala))) +
  ggtree::geom_tiplab(offset=i, aes(angle=angle)) 

ggsave("figures/indian bats only info sampling effort no tip lab circular.png", bg = "transparent", height = 18, width = 18)

Data %>% filter(effort_cat == 1)

#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................
#...................Adding in Nipah Only Data....................................................

setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
#bat_spp_India <- read_csv("C:/Users/r83c996/Desktop/bat spp India.csv")
bat_spp_India <- read_csv('data/phylofactor work spaces/bat spp India.csv')
bat_spp_Kerala <- read_csv('data/phylofactor work spaces/Bats_Kerala.csv')

load('data/phylofactor work spaces/bat_tree')
#load("data/bat_taxonomy_data.Rdata")
load("data/batphy.Rdata")

bat_spp_India %>%
  filter(!(MSW05_Binomial %in% batphy$unique_name))

bat_spp_India$MSW05_Binomial =plyr::revalue(bat_spp_India$MSW05_Binomial,c("Myotis_blythii"="Myotis_oxygnathus",
                                                                           "Myotis_peytoni"="Myotis_montivagus"))

bat_spp_India %>%
  filter(!(MSW05_Binomial %in% batphy$unique_name))

#only Eptesicus_bottae doesn't appear in the database, I can't find another name for this 
#species that is the database either. Also, according to the internets it doesn't actually 
#live in in India so maybe thats not a problem? funky face

Data <- batphy %>%
  filter(unique_name %in% bat_spp_India$MSW05_Binomial) %>% #only use bats that are in India, obvi gonna loose a couple 
  mutate(Nipah.samps = ifelse(is.na(Nipah.samps), 0, Nipah.samps)) %>% #remove NAs
  mutate(Nipah.Pos = ifelse(is.na(Nipah.Pos), 0, Nipah.Pos)) %>% #remove NAs
  mutate(log_Nipah.samps = log(Nipah.samps)) %>% #log transform sampling effort
  select(c(Nipah.samps, Nipah.Pos, species, log_Nipah.samps)) %>%
  dplyr::rename(Species = species) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  mutate(Species = stri_trans_totitle(Species)) %>%
  mutate(Sample = 1)

# remove tips from tree that are in our Data file
Data$Species %in% bat_tree$tip.label
tree <- ape::drop.tip(bat_tree,bat_tree$tip.label[!(bat_tree$tip.label %in% Data$Species)])
rm(batphy, bat_tree)

names(Data) <- c('effort','Z', 'Species','log_effort', 'Sample')

#load(file='data/phylofactor work spaces/indian bats no sampling effort only niv_workspace')
bat_spp_India <- bat_spp_India %>%
  mutate(MSW05_Binomial = gsub(" ", "_", MSW05_Binomial))

bat_spp_Kerala <- bat_spp_Kerala %>%
  mutate(MSW05_binomial = gsub(" ", "_", MSW05_binomial))

Data <- Data %>%
  mutate(IB = ifelse(Species %in% bat_spp_India$MSW05_Binomial, 1, 
                     ifelse(!(Species %in% bat_spp_India$MSW05_Binomial), 0, NA))) %>%
  mutate(Kerala = ifelse(Species %in% bat_spp_Kerala$MSW05_binomial, 1,
                         ifelse(!(Species %in% bat_spp_Kerala$MSW05_binomial), 0, NA)))

load('data/seroprevalence.Rdata')

Data <- Data %>%
  mutate(species.mutate = tolower(gsub("_", " ", Species))) %>%
  #mutate(all.neg = ifelse(species.mutate %in% neg.hnv$species & effort > 0, 1, 0))%>%
  mutate(all.neg = ifelse(effort > 0 & Z == 0, 1, 0)) %>%
  mutate(effort_cat = ifelse(effort > 0, 1,0))

Data <- Data %>%
  mutate(color = ifelse(Kerala == 1, "black", "grey"))

source('r/stupid_pointless_internal_node_color_function.R')
library('phylobase')
internal.edges<- (internal.edge.color(tree, Data))

row.names(internal.edges) <- internal.edges[,1]
#internal.edges <- internal.edges[,-c(1)] %>% as.data.frame()
g1 = as(tree, 'phylo4')
length(nodeId(g1, "internal")) == length(row.names(internal.edges))

g2 = phylo4d(g1, Data)
nodeData(g2) <- internal.edges

#...............................................................................

ggtree.object <- ggtree::ggtree(g2, layout= 'circular',branch.length = "none",aes(color=I(color))) 

jj <- nrow(Data)

d <- data.frame(x=ggtree.object$data[1:jj,'x'],
                xend=ggtree.object$data[1:jj,'x']  + 1 + Data$log_effort,
                y=ggtree.object$data[1:jj,'y'],
                yend=ggtree.object$data[1:jj,'y'] )

d

d[is.finite(d$xend)==FALSE,] <- 0

#g2 <- ggtree::ggtree(tree,layout='circular') 

ggtree.object$data$label <- gsub("_"," ", ggtree.object$data$label)

i <- .5
ii <- 5
iii <- 2

ggtree.object$data$Z[ggtree.object$data$Z==0]  <- NA
ggtree.object$data$effort_cat[ggtree.object$data$effort_cat==0]  <- NA

ggtree.object + 
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_cladelabel(node=118, label="Vespertilionidae", 
                          color="black", angle=12.16-180, offset=ii, offset.text = iii) +
  ggtree::geom_cladelabel(node=200, label="Pteropodidae", 
                          color="black", angle=-68.91, offset=ii, offset.text = iii) +
  ggtree::geom_tippoint(aes(size=effort_cat), pch =1) +
  ggtree::geom_tippoint(aes(size=Z), color = 'blue', alpha = .8) +
  geom_segment(data= d,aes(x=x + i ,y=y,xend=xend +i ,yend=yend, size= Data$effort_cat,
                           colour = 'orange', alpha = (.5 + Data$Kerala))) +
  ggtree::geom_tiplab(offset=i, aes(angle=angle)) 

ggsave("figures/indian bats only info sampling effort no tip lab circular NIPAH.png", bg = "transparent", height = 18, width = 18)

Data %>% filter(effort_cat == 1)
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............
# Making those silly circles...............

# .................................pteropodidae bats......................................

species.list <- ggtree::get.offspring.tip(tree, node=200)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ",
                            replacement = "",
                            group_taxonomy_list[,1])
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(tree, node=200)
mean(ggtree.object$data[ggtree.object$data$label %in% species.list,]$angle-90)

#.................................vester bats......................................

species.list <- ggtree::get.offspring.tip(tree, node=119)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ",
                            replacement = "",
                            group_taxonomy_list[,1])
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(tree, node=119)
mean(ggtree.object$data[ggtree.object$data$label %in% species.list,]$angle-90)

# 
# 
# 
# 
# 
# 
# 
