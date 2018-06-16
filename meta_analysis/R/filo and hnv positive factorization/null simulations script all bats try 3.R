# null simulations script -------

#getDV <- function(ss) ss['phylo','Deviance']
summaries <- lapply(pf$models,summary)
deviances <- sapply(summaries, "[", "loglik")

plot(deviances,type='o',lwd=2,cex=2)

randomPF <- function(pf){
  Data$Z.poisson <- sample(Data$Z.poisson)
  pf.random <- gpf(Data,tree,Z.poisson~phylo,nfactors=10,algorithm = 'phylo',
            model.fcn = model.fcn,objective.fcn = obj.fcn,
            cluster.depends='library(pscl)', ncores = ncores,
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

reps <- as.list(rep(reps.per.worker,ncores))
cl <- phyloFcluster(ncores=ncores)

clusterExport(cl,c('randomPF',
                   'randomPFs',
                   'phylobin',
                   'model.fcn',
                   'Data',
                   'tree',
                   'ncores',
                   'obj.fcn'
                   
                   ))

OBJ <- parLapply(cl,reps,randomPFs,pf=pf)

ltoMat <- function(OBJ){
  nr <- sum(sapply(OBJ,nrow))
  nc <- ncol(OBJ[[1]])
  Mat <- matrix(NA,nrow=nr,ncol=nc)
  for (l in 1:length(OBJ)){
    Mat[reps.per.worker*(l-1)+1:reps.per.worker,] <- OBJ[[l]]
  }
  return(Mat)
}

S <- ltoMat(OBJ) 

mm <- min(c(S))
mx <- max(c(c(S),Obj))

ecdf(S[,1])(Obj[1])



save(list=ls(), file='data/phylofactor work spaces/poisson hnv workspace' )

