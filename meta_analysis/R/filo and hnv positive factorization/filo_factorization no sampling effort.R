library(phylofactor)
library(parallel)
setwd('/Users/buckcrowley/Desktop/')

load('Z.filo_pos.filo1')
load('bat_tree.filo.pos')

tree <- bat_tree.filo.pos
Data <- Z.filo_pos.filo1


names(Data)[c(1,2)] <- c('effort','Z')
Data$effort <- log(Data$effort)

pf <- gpf(Data,tree,frmla.phylo=Z~phylo,nfactors=10,family=binomial,algorithm='phylo')
pf$factors
# 
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

fit.effort <- glm(Z~effort,data=Data,family=binomial)
Z.probs <- predict(fit.effort,type='response')

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


reps <- as.list(rep(20,7))
cl <- phyloFcluster(ncores=1)

clusterExport(cl,c('randomPF','randomPFs','pb','globalObj','Z.probs'))


OBJ <- parLapply(cl,reps,randomPFs,pf=pf)

ltoMat <- function(OBJ){
  nr <- sum(sapply(OBJ,nrow))
  nc <- ncol(OBJ[[1]])
  Mat <- matrix(NA,nrow=nr,ncol=nc)
  for (l in 1:length(OBJ)){
    Mat[20*(l-1)+1:20,] <- OBJ[[l]]
  }
  return(Mat)
}

S <- ltoMat(OBJ) 


mm <- min(c(S))
mx <- max(c(c(S),Obj))


plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,15))
for (rr in 1:140){
  lines(c(S[rr,1],diff(S[rr,])),col=rgb(0,0,0,0.2))
}

points(deviances,col='red',pch=16,cex=2)


save(list=ls(),file='filo_workspace')
