# null simulations script -------
#setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis")
#setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
require(pscl)
require(boot)

pf <- readRDS("meta_analysis/data/phylofactor work spaces/hnv_phylofactor_object_negbin")

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
    return(abs(summary(fit)$coefficients$zero['phyloS','z value']))
  }
}


Data <- pf$Data
tree <- pf$tree
ncores=1
tot.reps=2
reps.per.worker=round(tot.reps/ncores)

randomPF <- function(pf){
  Data$Z.poisson <- sample(Data$Z.poisson)
  pf.random <- gpf(Data,tree,Z.poisson~phylo,nfactors=10,algorithm = 'phylo',
            model.fcn = model.fcn,objective.fcn = obj.fcn,
            cluster.depends='library(pscl)', #ncores = ncores,
            dist = "negbin")
  summaries <- lapply(pf.random$models,summary)
  loglik <- unlist(sapply(summaries, "[", "loglik"))
  
  return(loglik)
}

randomPFs <- function(reps,pf){
  obj <- matrix(NA,nrow=reps,ncol=pf$nfactors)
  for (rr in 1:reps){
    obj[rr,] <- randomPF(pf)
  }
  return(obj)
}

#y <- randomPFs(2,pf )

reps <- as.list(rep(reps.per.worker,ncores))
cl <- phyloFcluster(ncores=ncores)


clusterExport(cl,
              varlist = 
                   c('pf',
                   'randomPF',
                   'randomPFs',
                   'model.fcn',
                   'obj.fcn',
                   'Data',
                   'tree',
                   'ncores'
                   ))

#OBJ <- lapply(X=reps,FUN= randomPFs,pf=pf)

OBJ<-parLapply(cl, reps, randomPFs, pf)


#...............................................................................


