setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis")
library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
load('data/phylofactor work spaces/bat_tree')
load("data/bat_taxonomy_data.Rdata")
ncores = 4
tot.reps=200
reps.per.worker=round(tot.reps/ncores)

Data <- batphy1 %>%
  mutate(hnv_surv_bin = ifelse(is.na(hnv_surv), 0, hnv_surv)) %>%
  mutate(hnv_surv_bin = ifelse(hnv_surv_bin == 1, 1, 0)) %>%
  mutate(log_hnv_samps = log(hnv_samps)) %>%
  select(c(hnv_surv_bin, species, log_hnv_samps))%>%
  rename(Species = species) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  mutate(Species = stri_trans_totitle(Species)) %>%
  mutate(Sample = 1) %>%
  unique()

Data[duplicated(Data$Species),]

tree<-bat_tree
tree <- ape::drop.tip(tree,tree$tip.label[!(tree$tip.label %in% Data$Species)])
Data <- Data[Data$Species %in% tree$tip.label,]

rm(batphy1, bat_tree)

names(Data) <- c('Z', 'Species', 'log_effort', 'Sample')

n_factors = 10

pf <- phylofactor::twoSampleFactor(Data$Z, tree, method='Fisher', n_factors ,ncores = 3)


source('R/filo and hnv positive factorization/null simulations script all bats try 2.R')

save(list=ls(),file='data/phylofactor work spaces/hnv_workspace_sample_no_sample_all_bat_dataset')








rm(list=ls())

load('data/phylofactor work spaces/bat_tree')
load("data/bat_taxonomy_data.Rdata")
ncores = 4
tot.reps=200
reps.per.worker=round(tot.reps/ncores)
n_factors = 10

Data <- batphy1 %>%
  mutate(filo_surv_bin = ifelse(is.na(filo_surv), 0, filo_surv)) %>%
  mutate(filo_surv_bin = ifelse(filo_surv_bin == 1, 1, 0)) %>%
  mutate(log_filo_samps = log(filo_samps)) %>%
  select(c(filo_surv_bin, species, log_filo_samps))%>%
  rename(Species = species) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  mutate(Species = stri_trans_totitle(Species)) %>%
  mutate(Sample = 1) %>%
  unique()

Data[duplicated(Data$Species),]

tree<-bat_tree
tree <- ape::drop.tip(tree,tree$tip.label[!(tree$tip.label %in% Data$Species)])
Data <- Data[Data$Species %in% tree$tip.label,]

rm(batphy1, bat_tree)

names(Data) <- c('Z', 'Species', 'log_effort', 'Sample')

#pf <- gpf(Data,tree,frmla.phylo=Z~phylo,nfactors=10,family=binomial,algorithm='phylo')
pf <- phylofactor::twoSampleFactor(Data$Z, tree, method='Fisher', n_factors ,ncores = 3)

pf$factors
# 
source('R/filo and hnv positive factorization/null simulations script all bats try 2.R')

save(list=ls(),file='data/phylofactor work spaces/filo_workspace_sample_no_sample_all_bat_dataset')

rm(list=ls())


#............................................ visualization
