## code for bat filovirus and HNV sampling effort phylofactor
## daniel.becker3@montana.edu

## clean environment & plots
#rm(list=ls()) 
graphics.off()

## packages
library(rotl)
library(caper)
library(ape)

## read in batphy
setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
batphy=read.csv("data/batphy_for_rotl update.csv",header=T)


## get tree with otts
bat_tree=tol_induced_subtree(ott_ids=batphy$ott_id)

## remove ott information from the tips
bat_tree$tip.label=strip_ott_ids(bat_tree$tip.label)

## assign branch lengths
bat_tree=compute.brlen(bat_tree,method="Grafen")

## sort batphy into same order as bat_tree
bat_tree=makeLabel(bat_tree)
batphy=batphy[match(bat_tree$tip.label,batphy$tree_species),]

## give row names
rownames(batphy)=batphy$tree_species

## use caper to join (takes a while)
cdata=comparative.data(phy=bat_tree,data=batphy,names.col=tree_species,vcv=T,na.omit=F,warn.dropped=T)

## assign colors for sampled or not
cdata$data$survey_col=ifelse(cdata$data$surveyed==1,"blueviolet",NA)

## assign colors for viral detection
cdata$data$survey_col=ifelse(cdata$data$surveyed==0,NA,ifelse(cdata$data$both_surv==1,"blueviolet",
                                                              ifelse(cdata$data$filo_surv==1,
                                                                     "cornflowerblue","firebrick3")))

## trim to africa only
afdata=cdata[which(cdata$data$region=="Africa"),]

## plot trees
par(mar=c(0,0,0,0),oma=c(0,0,0,0),mfrow=c(1,2))

## full phylogeny
plot(cdata$phy,type="fan",show.tip.label=F,edge.width=1)

## add labels based on sampling or not
tiplabels(pch=21,col=cdata$data$survey_col,bg=cdata$data$survey_col,cex=1)

## africa only
plot(afdata$phy,type="fan",show.tip.label=F,edge.width=1)

africa_tree <- afdata$phy

## add labels based on sampling or not
tiplabels(pch=21,col=afdata$data$survey_col,bg=afdata$data$survey_col,cex=1)



#.............add in filovirus positivity....................
# 
# library(tidyverse)
# setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
# load(file='data/seroprevalence.Rdata')
# seroprevalence_species <- seroprevalence %>%
#   mutate(seroprevalence_percentage_cat = ifelse(seroprevalence_percentage == 0, 0, 1)) %>%
#   filter(!(is.na(species))) %>%
#   dplyr::group_by(species, virus, seroprevalence_percentage_cat) %>%
#   summarise(count = n()) %>%
#   ungroup() %>%
#   filter(virus == 'Filovirus') %>%
#   spread(seroprevalence_percentage_cat, count) %>%
#   select(c(species, `1`, `0`)) %>%
#   mutate(filovirus_positive_cat  = ifelse(is.na(`1`) & !is.na(`0`), 0, 
#                                           ifelse(!is.na(`1`), 1, NA)))
# 
# batphy <- full_join(batphy, seroprevalence_species) %>%
#   mutate(filovirus_positive_cat = ifelse(is.na(filovirus_positive_cat), 0, filovirus_positive_cat))
