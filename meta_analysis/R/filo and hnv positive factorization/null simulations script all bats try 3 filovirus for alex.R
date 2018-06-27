# null simulations script -------
#setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis")
setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
require(pscl)
require(boot)
library(tictoc)
pf <- readRDS("data/phylofactor work spaces/hnv_phylofactor_object_negbin")

model.fcn <- function(formula,data,...){
  fit <- tryCatch(pscl::zeroinfl(formula,data,...),
                  error=function(e) NA)
  #fit <- do.call
  return(fit)
}

model.fcn2 <- function(formula,data,...){
  fit <- tryCatch(zeroinfl(formula,data,...),
                  error=function(e) NA)
  #fit <- do.call
  return(fit)
}


randomPF <- function(pf){
  Data$Z.poisson <- sample(Data$Z.poisson)
  pf.random <- gpf(Data,tree,Z.poisson~phylo,nfactors=10,algorithm = 'phylo',
                   model.fcn = model.fcn,objective.fcn = obj.fcn,
                   cluster.depends='library(pscl)', #ncores = ncores,
                   dist = "negbin")
  summaries <- lapply(pf.random$models,summary)
  loglik <- unlist(sapply(summaries, "[", "loglik"))
  
  summaries.null <- (lapply(pf$models, "[[", "model"))
  loglik.null <- lapply(summaries.null, 
                        function(data){
                          pscl::zeroinfl(Z.poisson~1, data)
                        })
  loglik.null <- unlist(sapply(loglik.null, "[", "loglik"))
  loglik-loglik.null
  
  return(loglik-loglik.null)
}

obj.fcn2 <- function(fit,grp,tree,PartitioningVariables,model.fcn,phyloData,...){
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
#cl <- phyloFcluster(ncores=ncores)
stopCluster(cl)

#clusterExport(cl,
#               varlist = 
#                    c('pf',
#                    'randomPF',
#                    'randomPFs',
#                    'model.fcn',
#                    'obj.fcn',
#                    'Data',
#                    'tree',
#                    'ncores'
#                    ))

OBJ <- lapply(X=reps,FUN= randomPFs,pf=pf)

#OBJ<-parLapply(cl, reps, randomPFs, pf)

tic()
dum=gpf(pf$Data,pf$tree,frmla.phylo = Z.poisson~phylo, nfactors=1,
    model.fcn=model.fcn2,objective.fcn=obj.fcn2,
    algorithm='phylo',dist='negbin')
toc()
#...............................................................................
tic()
dum2=gpf(pf$Data,pf$tree,frmla.phylo = Z.binom~phylo, nfactors=1,
        family=binomial,
        algorithm='phylo')
toc()

