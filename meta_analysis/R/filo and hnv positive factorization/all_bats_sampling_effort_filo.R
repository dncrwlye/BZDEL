
tryCatch(setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis"),
                error=function(e) setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/"))

library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
require(pscl)
require(boot)
library(pscl)
load('data/phylofactor work spaces/bat_tree')
load("data/bat_taxonomy_data.Rdata")

Data <- batphy1 %>%
  mutate(filo_surv_bin = ifelse(is.na(filo_surv), 0, filo_surv)) %>%
  mutate(filo_surv_bin = ifelse(filo_surv_bin == 1, 1, 0)) %>%
  mutate(log_filo_samps = log(filo_samps)) %>%
  select(c(filo_samps, filo_surv_bin, species, log_filo_samps))%>%
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

source('R/filo and hnv positive factorization/negative binomial no zero inflation tests for all bats sampling effort.R')
source('R/filo and hnv positive factorization/null simulations script all bats try 4.R')

save(list=ls(),file='data/phylofactor work spaces/neg binom no zeroinfl filo_workspace_try2')

rm(list=ls())
