# null simulations script -------
#setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis")
#setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
require(pscl)
require(boot)
library(tictoc)
#pf <- readRDS("data/phylofactor work spaces/hnv_phylofactor_object_negbin")

Data.random <- pf2$Data
tree <- pf2$tree
ncores = 4
nfactors = 10
#reps.per.worker=round(tot.reps/ncores)

randomPF <- function(pf){
  Data.random$Z.poisson <- sample(Data.random$Z.poisson, replace=TRUE)
  #pf.random <- gpf(Data.random,tree,Z.poisson~phylo,nfactors=10,algorithm = 'phylo',
  #               model.fcn = model.fcn2,objective.fcn = obj.fcn2,
  #               cluster.depends='library(MASS)', 
  #               ncores = 4)
  
  pf.random <- tryCatch(gpf(Data.random,tree,Z.poisson~phylo,nfactors=nfactors,algorithm = 'phylo',
    model.fcn = model.fcn2,objective.fcn = obj.fcn2,
    cluster.depends='library(MASS)', 
    ncores = 4), error = function(e) NULL)
  
  
  #summaries.null <- (lapply(pf.random$models, "[[", "model"))
  #loglik.null <- lapply(summaries.null, 
  #                      function(data){
  #                        MASS::glm.nb(Z.poisson~1, data)
  #                      })
  
  #logLik.models <- unlist(lapply(pf.random$models,logLik))
  #loglik.null <- unlist(sapply(loglik.null, logLik))
  #return(logLik.models-loglik.null)
  
  return(unlist(lapply(pf.random$models, "[[", "null.deviance")) - unlist(lapply(pf.random$models, "[[", "deviance")))
}

randomPF2 <- function(pf){
  Data.random$Z.poisson <- rnbinom(n =length(Data$Z.poisson),
                                   mu = mean(Data$Z.poisson), 
                                   size = mean(Data$Z.poisson)^2 / (var(Data$Z.poisson) - mean(Data$Z.poisson)))
  #pf.random <- gpf(Data.random,tree,Z.poisson~phylo,nfactors=10,algorithm = 'phylo',
  #                 model.fcn = model.fcn2,objective.fcn = obj.fcn2,
  #                 cluster.depends='library(MASS)', 
  #                 ncores = 4)
  
  pf.random <- tryCatch(gpf(Data.random,tree,Z.poisson~phylo,nfactors=nfactors, algorithm = 'phylo',
                            model.fcn = model.fcn2,objective.fcn = obj.fcn2,
                            cluster.depends='library(MASS)', 
                            ncores = 4), error = function(e) NULL)
  
  #summaries.null <- (lapply(pf.random$models, "[[", "model"))
  #loglik.null <- lapply(summaries.null, 
  #                      function(data){
  #                        MASS::glm.nb(Z.poisson~1, data)
  #                      })
  
  #logLik.models <- unlist(lapply(pf.random$models,logLik))
  #loglik.null <- unlist(sapply(loglik.null, logLik))
  
  #return(logLik.models-loglik.null)
  return(unlist(lapply(pf.random$models, "[[", "null.deviance")) - unlist(lapply(pf.random$models, "[[", "deviance")))
  
}

randomPFs <- function(reps,pf){
  obj <- matrix(NA,nrow=reps,ncol=nfactors)
  for (rr in 1:reps){
    obj[rr,] <- randomPF(pf)
  }
  return(obj)
}

randomPFs2 <- function(reps,pf){
  obj <- matrix(NA,nrow=reps,ncol=nfactors)
  for (rr in 1:reps){
    obj[rr,] <- randomPF2(pf)
  }
  return(obj)
}

#.....................tests...............................
tot.reps = 100
OBJ1 <- lapply(tot.reps, randomPFs, pf2)
 
tot.reps = 10
OBJ2.1 <- lapply(tot.reps, randomPFs2, pf2)
OBJ2.2 <- lapply(tot.reps, randomPFs2, pf2)
OBJ2.3 <- lapply(tot.reps, randomPFs2, pf2)
OBJ2.4 <- lapply(tot.reps, randomPFs2, pf2)
OBJ2.5 <- lapply(tot.reps, randomPFs2, pf2)
OBJ2.6 <- lapply(tot.reps, randomPFs2, pf2)
OBJ2.7 <- lapply(tot.reps, randomPFs2, pf2)
OBJ2.8 <- lapply(tot.reps, randomPFs2, pf2)
OBJ2.9 <- lapply(tot.reps, randomPFs2, pf2)
OBJ2.10 <-lapply(tot.reps, randomPFs2, pf2)

dfs <- mget(ls(pattern = "OBJ2"))
dfs <- lapply(dfs, '[[', 1)
OBJ.final.2 <- do.call(rbind, dfs)
                       
rm(list = ls(pattern = "OBJ2"))

#x <- randomPF(pf)
#y <- randomPFs(2,pf )
#...............................................................................


