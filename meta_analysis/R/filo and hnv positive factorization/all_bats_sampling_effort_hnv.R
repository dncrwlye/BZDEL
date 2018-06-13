#setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis")
setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
library(pscl)
library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
require(pscl)
require(boot)
load('data/phylofactor work spaces/bat_tree')
load("data/bat_taxonomy_data.Rdata")
ncores = 4
tot.reps=200
reps.per.worker=round(tot.reps/ncores)

Data <- batphy1 %>%
  mutate(hnv_surv_bin = ifelse(is.na(hnv_surv), 0, hnv_surv)) %>%
  mutate(hnv_surv_bin = ifelse(hnv_surv_bin == 1, 1, 0)) %>%
  mutate(log_hnv_samps = log(hnv_samps)) %>%
  select(c(hnv_samps, hnv_surv_bin, species, log_hnv_samps))%>%
  dplyr::rename(Species = species) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  mutate(Species = stri_trans_totitle(Species)) %>%
  mutate(Sample = 1) %>%
  unique()

Data[duplicated(Data$Species),]

tree<-bat_tree
tree <- ape::drop.tip(tree,tree$tip.label[!(tree$tip.label %in% Data$Species)])
Data <- Data[Data$Species %in% tree$tip.label,]

#rm(batphy1, bat_tree)

names(Data) <- c('Z.poisson',"Z.binom", 'Species', 'log_effort', 'Sample')

Data <- Data[Data$Species %in% tree$tip.label,]
n_factors = 10
Data[is.na(Data)] <- 0

# obj.fcn <- function(fit,grp,tree,PartitioningVariables,model.fcn,phyloData,...){
#   return(abs(summary(fit)$coefficients$zero['phyloS','z value']))
# }

model.fcn <- function(formula,data,...){
   fit <- tryCatch(zeroinfl(formula,data,...),
                   error=function(e) NA)
  #fit <- do.call
  return(fit)
  }

obj.fcn <- function(fit,grp,tree,PartitioningVariables,model.fcn,phyloData,...){
if (!'zeroinfl' %in% class(fit))
  {
    return(0)
  }
else 
  {
    return(abs(summary(fit)$coefficients$zero['phyloS','z value']))
  }
}

pf <- gpf(Data,tree,Z.poisson~phylo,nfactors=10,algorithm = 'phylo',
          model.fcn = model.fcn,objective.fcn = obj.fcn,
          cluster.depends = {library(pscl)},
          dist = "negbin")

# pf <- gpf(Data,tree,Z.poisson~phylo,nfactors=2,
#           algorithm = 'phylo',
#           ncores=1,
#           family=poisson
#           )
# 
# hist(Data$Z.poisson, breaks=seq(0,65, by=1))
#   
# )# pf <- gpf(Data,tree,
# #           frmla.phylo=Z~phylo,
# #           family=binomial,
# #           nfactors=10,
# #           algorithm='phylo')
# 
# 
# pf <- phylofactor::twoSampleFactor(Data$Z, tree, method='Fisher', n_factors ,ncores = 3)
# 
# source('R/filo and hnv positive factorization/null simulations script all bats try 2.R')
# 
# save(list=ls(),file='data/phylofactor work spaces/hnv_workspace_sample_no_sample_all_bat_dataset')
# 
# 
# rm(list=ls())
# 
# load('data/phylofactor work spaces/bat_tree')
# load("data/bat_taxonomy_data.Rdata")
# ncores = 4
# tot.reps=200
# reps.per.worker=round(tot.reps/ncores)
# n_factors = 10
# 
# Data <- batphy1 %>%
#   mutate(filo_surv_bin = ifelse(is.na(filo_surv), 0, filo_surv)) %>%
#   mutate(filo_surv_bin = ifelse(filo_surv_bin == 1, 1, 0)) %>%
#   mutate(log_filo_samps = log(filo_samps)) %>%
#   select(c(filo_surv_bin, species, log_filo_samps))%>%
#   rename(Species = species) %>%
#   mutate(Species = gsub(" ", "_", Species)) %>%
#   mutate(Species = stri_trans_totitle(Species)) %>%
#   mutate(Sample = 1) %>%
#   unique()
# 
# Data[duplicated(Data$Species),]
# 
# tree<-bat_tree
# tree <- ape::drop.tip(tree,tree$tip.label[!(tree$tip.label %in% Data$Species)])
# Data <- Data[Data$Species %in% tree$tip.label,]
# 
# rm(batphy1, bat_tree)
# 
# names(Data) <- c('Z', 'Species', 'log_effort', 'Sample')
# 
# #pf <- gpf(Data,tree,frmla.phylo=Z~phylo,nfactors=10,family=binomial,algorithm='phylo')
# pf <- phylofactor::twoSampleFactor(Data$Z, tree, method='Fisher', n_factors ,ncores = 3)
# 
# pf$factors
# # 
# source('R/filo and hnv positive factorization/null simulations script all bats try 2.R')

save(list=ls(),file='data/phylofactor work spaces/filo_workspace_sample_no_sample_all_bat_dataset')

rm(list=ls())


#............................................ visualization
