library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
load('data/phylofactor work spaces/bat_tree')
load("data/bat_taxonomy_data.Rdata")

Data <- batphy1 %>%
  filter(hnv_samps > 0) %>%
  mutate(log_hnv_samps = log(hnv_samps)) %>%
  select(c(hnv_samps, hnv_positive, species, log_hnv_samps)) %>%
  rename(Species = species) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  mutate(Species = stri_trans_totitle(Species)) %>%
  mutate(Sample = 1)

tree <- ape::drop.tip(bat_tree,bat_tree$tip.label[!(bat_tree$tip.label %in% Data$Species)])

rm(batphy1, bat_tree)

names(Data) <- c('effort','Z', 'Species','log_effort', 'Sample')

model.sampling.effort <- glm(Z~log_effort,family=binomial,data=Data)

Data <- Data %>%
  mutate(effort.fit = coef(model.sampling.effort)['log_effort']*Data$log_effort)

rm(model.sampling.effort)

pf <- gpf(Data,tree,frmla=Z~offset(effort.fit),
          frmla.phylo=Z~offset(effort.fit)+phylo,
          family=binomial,
          nfactors=10,
          algorithm='phylo')
pf$factors
# 
source('R/filo and hnv positive factorization/null simulations script.R')

save(list=ls(),file='data/phylofactor work spaces/hnv_workspace')

rm(list=ls())


