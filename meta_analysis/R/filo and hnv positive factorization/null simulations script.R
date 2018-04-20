# null simulations script -------

getDV <- function(ss) ss['phylo','Deviance']
summaries <- lapply(pf$models,anova)
deviances <- sapply(summaries,getDV)

plot(deviances,type='o',lwd=2,cex=2)

pb <- phylofactor:::phylobin

globalObj <- function(factor,pf){
  DF <- cbind(pf$Data,'phylo'=factor(pb(bins(pf$basis[,1:factor,drop=F]))))
  fit <- pf$model.fcn(pf$frmla.phylo,data=DF,family=binomial)
  return(anova(fit)['phylo','Deviance'])
}

Obj <- sapply(1:pf$nfactors,globalObj,pf=pf)

plot(Obj,type='o',lwd=2,cex=2)

if (sampling_effort == TRUE)
{
fit.effort <- glm(Z~log_effort,data=Data,family=binomial)
Z.probs <- predict(fit.effort,type='response')
}

if (sampling_effort == FALSE)
{
  fit.effort <- glm(Z~1,data=Data,family=binomial)
  Z.probs <- predict(fit.effort,type='response')
}

randomPF <- function(pf){
  pf$Data$Z <- rbinom(nrow(pf$Data),size=1,prob=Z.probs)
  pf <- gpf(pf$Data,pf$tree,frmla.phyl=pf$frmla.phylo,nfactors=pf$nfactors,family=binomial,algorithm='phylo')
  return(sapply(1:pf$nfactors,globalObj,pf=pf))
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

clusterExport(cl,c('randomPF','randomPFs','pb','globalObj','Z.probs'))

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


