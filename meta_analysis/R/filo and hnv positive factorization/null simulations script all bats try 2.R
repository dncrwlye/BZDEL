# null simulations script try 2!  -------

pf$pvals

fit.effort <- glm(Z~1,data=Data,family=binomial)
Z.probs <- predict(fit.effort,type='response')


randomPF <- function(pf){
  Data$Z <- rbinom(nrow(pf$Data),size=1,prob=Z.probs)
  pf <- phylofactor::twoSampleFactor(Data$Z, tree, method='Fisher', n_factors)
  return(pf$pvals)
}

randomPFs <- function(reps,pf){
  obj <- matrix(NA,nrow=reps,ncol=pf$nfactors)
  for (rr in 1:reps){
    obj[rr,] <- randomPF(pf)
  }
  return(obj)
}

reps <- as.list(rep(reps.per.worker,ncores))
cl <- phyloFcluster(ncores=ncores)

clusterExport(cl,c('randomPF','randomPFs','Z.probs','Data','tree', 'n_factors'))

OBJ <- parLapply(cl,reps,randomPFs,pf=pf)






