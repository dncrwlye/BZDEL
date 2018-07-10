# null simulations script -------
setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis")
#setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
require(pscl)
require(boot)
library(tictoc)
#pf <- readRDS("data/phylofactor work spaces/hnv_phylofactor_object_negbin")

model.fcn <- function(formula,data,...){
  fit <- tryCatch(pscl::zeroinfl(formula,data,...),
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
    fit2 <- zeroinfl(Z.poisson~1,data = fit$model,dist='negbin')
    # return(abs(summary(fit)$coefficients$zero['phyloS','z value']))
    fit$loglik-fit2$loglik %>% return()
  }
}
Data <- pf$Data
tree <- pf$tree
ncores = 
#reps.per.worker=round(tot.reps/ncores)

randomPF <- function(pf){
  Data$Z.poisson <- sample(Data$Z.poisson)
  pf.random <- gpf(Data,tree,Z.poisson~phylo,nfactors=10,algorithm = 'phylo',
                   model.fcn = model.fcn,objective.fcn = obj.fcn,
                   cluster.depends='library(pscl)', ncores = 4,
                   dist = "negbin")
  summaries <- lapply(pf.random$models,summary)
  loglik <- unlist(sapply(summaries, "[", "loglik"))
  
  summaries.null <- (lapply(pf.random$models, "[[", "model"))
  loglik.null <- lapply(summaries.null, 
                        function(data){
                          pscl::zeroinfl(Z.poisson~1, data)
                        })
  loglik.null <- unlist(sapply(loglik.null, "[", "loglik"))
  return(loglik-loglik.null)
}

randomPFs <- function(reps,pf){
  obj <- matrix(NA,nrow=reps,ncol=pf$nfactors)
  for (rr in 1:reps){
    obj[rr,] <- randomPF(pf)
  }
  return(obj)
}
#.....................tests...............................

OBJ1 <- lapply(tot.reps, randomPFs, pf)
OBJ2 <- lapply(tot.reps, randomPFs, pf)

#x <- randomPF(pf)
#y <- randomPFs(2,pf )
#...............................................................................


