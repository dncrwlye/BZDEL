####### Replicating Olival's GAMs with and without removal of non-zoonotic clades.
### Below, we've copied & pasted their file "04-fit-models.R"
library(mgcv)
library(dplyr)
library(stringi)
library(parallel)
library(purrr)
library(ggplot2)
library(viridis)
library(knitr)
library(svglite)
library(xlsx)
require(MuMIn)
source("Olival/R/model_reduction.R")
source("Olival/R/fit_gam.R")
source("Olival/R/relative_contributions.R")
source("Olival/R/cross_validation.R")
source("Olival/R/logp.R")
set.seed(0)

for (original in c(T,F)){
  db <- readRDS("Olival/intermediates/postprocessed_database.rds")
  
  # our modification
  if (!original){
    V <- read.xlsx('Datasets/Olival_w_phylobin.xlsx',sheetIndex = 1)
    db$viruses <- cbind(db$viruses,V[,c('pfIsZ','pfIsZS','nonZ')])
    db$viruses <- db$viruses[!as.logical(db$viruses$nonZ),]
  }
  
  hosts <- db$hosts
  viruses <- db$viruses
  associations <- db$associations
  rm(db)
  
  #---- All zoonoses ----
  
  data_set = hosts %>%
    filter(hMarOTerr == "Terrestrial",
           hWildDomFAO == "wild",
           !is.na(PdHoSa.cbCst_order))
  
  outcome_variable = "NSharedWithHoSa"
  
  model_family = poisson
  
  #  Create dummy variables for orders to use as random effects
  dummys = as.data.frame(with(data_set, model.matrix(~hOrder))[,-1])
  data_set = cbind(data_set, dummys)
  dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
  names(dummy_terms) <- names(dummys)
  
  ## Create data.frame of all possible models
  terms = list(
    mass = "s(hMassGramsPVR, bs = 'tp', k=7)",
    interaction = c(
      "s(HabAreaCropLn, bs = 'tp', k=7)   + s(HabAreaCropChgLn, bs = 'tp', k=7)",
      "s(HabAreaGrassLn, bs = 'tp', k=7)  + s(HabAreaGrassChgLn, bs = 'tp', k=7)",
      "s(HabAreaUrbanLn, bs = 'tp', k=7)  + s(HabAreaUrbanChgLn, bs = 'tp', k=7)",
      "s(HabInhabitedLn, bs = 'tp', k=7)  + s(HabInhabitedChgLn, bs = 'tp', k=7)",
      "s(TotHumPopLn, bs = 'tp', k=7) + s(TotHumPopChgLn, bs = 'tp', k=7) + s(UrbRurPopRatioLn, bs = 'tp', k=7) + s(UrbRurPopRatioChg, bs = 'tp', k=7)",
      "s(HumPopDensLn, bs = 'tp', k=7) + s(HumPopDensLnChg, bs = 'tp', k=7) + s(UrbRurPopRatioLn, bs = 'tp', k=7) + s(UrbRurPopRatioChg, bs = 'tp', k=7)"),
    interaction2 = "s(hHuntedIUCN, bs='re')",
    interaction3 = "s(hArtfclHbttUsrIUCN, bs='re')",
    phylo_distance = c("s(PdHoSa.cbCst, bs = 'tp', k=7)", "s(PdHoSaSTPD, bs = 'tp', k=7)"),
    bias = c("s(hDiseaseZACitesLn, bs = 'tp', k=7)"),
    offset = "offset(LnTotNumVirus)",
    stringsAsFactors=FALSE)
  
  terms = c(dummy_terms, terms)
  
  
  if (!original){
    all_zoonoses_new=fit_all_gams(data_set,
                                  outcome_variable,
                                  poisson,
                                  terms)
  } else {
    all_zoonoses = fit_all_gams(data_set,
                                outcome_variable,
                                poisson,
                                terms)
  }
}

# ---- Viral Traits GAM - All Associations ----
for (original in c(T,F)){
  db <- readRDS("Olival/intermediates/postprocessed_database.rds")
  
  # our modification
  if (!original){
    V <- read.xlsx('Datasets/Olival_w_phylobin.xlsx',sheetIndex = 1)
    db$viruses <- cbind(db$viruses,V[,c('pfIsZ','pfIsZS','nonZ')])
    db$viruses <- db$viruses[!as.logical(db$viruses$nonZ),]
  }
  
  hosts <- db$hosts
  viruses <- db$viruses
  associations <- db$associations
  rm(db)
  PD_centers = names(viruses)[stri_detect_regex(names(viruses), "^(cb|st)_.*(?<!stringent)$")]
  
  model_family = binomial
  
  
  data_set = viruses %>%
    mutate_each_(funs("logp"), vars=PD_centers)
  names(data_set)[names(data_set) %in% PD_centers] <- paste0(PD_centers, "Ln")
  PD_centers <- paste0(PD_centers, "Ln")
  
  dummys = as.data.frame(with(viruses, model.matrix(~vFamily))[,-1])
  data_set = cbind(data_set, dummys)
  dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
  names(dummy_terms) <- names(dummys)
  
  
  terms = list(
    PD     = paste0("s(", PD_centers, ", bs='tp', k=7)"),
    bias   = c("s(vPubMedCitesLn, bs='tp', k=7)", "s(vWOKcitesLn, bs='tp', k=7)"),
    strand = c("s(RNA, bs='re')", "s(SS, bs='re')", "s(vCytoReplicTF, bs='re')"),
    vector = "s(Vector, bs='re')",
    env = "s(Envelope, bs='re')",
    genome = "s(vGenomeAveLengthLn, bs='tp', k=7)",
    stringsAsFactors=FALSE)
  
  if (original){
    vtraits = fit_all_gams(data_set,
                         outcome_variable = "IsZoonotic",
                         model_family = binomial,
                         terms)
    data_set_old <- data_set
  } else {
    data_set_new <- data_set
    vtraits_new = fit_all_gams(data_set,
                           outcome_variable = "IsZoonotic",
                           model_family = binomial,
                           terms)
  }
}


save(list=ls(),file='Workspaces/phylofactor_comparisons')
