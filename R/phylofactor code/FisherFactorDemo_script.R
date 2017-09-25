library(phylofactor)
library(ape)
library(phangorn)
library(phytools)
###### Customized fisher-test phylofactorization



# Producing tree & dataset ----------------------------------------------------------

set.seed(1)
tree <- ape::rtree(100)
#nodes 173 high, 154 low, 108 high probabilities of spillover
nds <- c(181,130,112)
probs <- c(.9,.1,.8,.5)
grps <- lapply(as.list(nds),FUN = function(n,t) Descendants(t,n,'tips')[[1]],t=tree)
grps[[4]] <- setdiff(1:100,unlist(grps))

Z <- numeric(100)
names(Z) <- tree$tip.label
set.seed(1)
for (i in 1:4){
  nms <- tree$tip.label[grps[[i]]]
  Z[nms] <- rbinom(length(grps[[i]]),1,probs[i])
}

plot.phylo(tree,use.edge.length = F,show.tip.label = F)
nodelabels(probs[1:3],nds)





# Functions - may later be incorporated into package ----------------------
########### Our test:
FisherTest <- function(grps,tree,Z){
  s1 <- Z[tree$tip.label[grps[[1]]]]
  s2 <- Z[tree$tip.label[grps[[2]]]]
  n1 <- sum(s1)
  n2 <- sum(s2)
  ft <- fisher.test(matrix(c(n1,length(s1)-n1,n2,length(s2)-n2),ncol=2),alternative='two.sided')
  return(ft$p.value)
}

p.fcn <- function(pvals,lambda){
  return(pvals^lambda)
}

fisherFactor <- function(Z,tree,nfactors,random=F,prob.fcn=p.fcn,lambda=-1,null=F){
  treeList <- list(tree)
  binList <- list(1:ape::Ntip(tree))
  Grps=getGroups(tree)
  output <- NULL
  
  if(null){
    p <- sum(Z)/length(Z)
    Z <- rbinom(length(Z),1,prob = p)
    names(Z) <- tree$tip.label
  }
  
  pfs=0
  while (pfs < min(length(Z)-1,nfactors)){

    if (pfs>=1){
      treeList <- updateTreeList(treeList,binList,grp,tree,skip.check=T)
      binList <- updateBinList(binList,grp)
      Grps <- getNewGroups(tree,treeList,binList)
    }
    
    pvals <- sapply(Grps,FUN=FisherTest,tree=tree,Z=Z)
    
    if (!random){
      best.grp <- Grps[[which.min(pvals)]]
      P <- min(pvals)
    } else {
      probs <- prob.fcn(pvals,lambda)
      ix <- sample(length(Grps),1,prob = probs)
      best.grp <- Grps[[ix]]
      P <- pvals[ix]
    }
    
    output$pvals <- c(output$pvals,P)
    
    grp <- getLabelledGrp(tree=tree,Groups=best.grp)
    output$groups <- c(output$groups,list(best.grp))
    
    grpInfo <- matrix(c(names(grp)),nrow=2)
    output$factors <- cbind(output$factors,grpInfo)
    
    pfs=pfs+1
  }
  
  return(output)
}




# Demo - deterministically pull three factors from Z ----------------------

Y <- fisherFactor(Z,tree,3,random=F)

#get binned phylogenetic units (BPUs)
V <- NULL
for (i in 1:length(Y$groups)){
 V <- cbind(V,ilrvec(Y$groups[[i]],Ntip(tree)))
}

BPUs <- bins(V)


# parallelization of stochastic phylofactorization ------------------------

library(parallel)
library(tictoc)

reps=600
tic()
cl <- phyloFcluster(6)
clusterExport(cl=cl,varlist=c('p.fcn','FisherTest','fisherFactor'))
FACTORS <- parLapply(cl,1:reps,fun=function(r,Z,tree,nfactors) fisherFactor(Z,tree,nfactors,random=T,lambda=-1),Z=Z,tree=tree,nfactors=10)

NullFactors <- parLapply(cl,1:reps,fun=function(r,Z,tree,nfactors) fisherFactor(Z,tree,nfactors,random=T,null=T,lambda=-1),Z=Z,tree=tree,nfactors=10)

stopCluster(cl)

toc()


save(list=ls(),file='C:/Users/Alculus/Documents/Bozeman/Spillover_Phylofactorization/FisherFactorDemo_workspace')


## extract P values
Ps <- lapply(FACTORS,FUN=function(x) x$pvals)
Pnull <- lapply(NullFactors,FUN=function(x) x$pvals)


## extract a measure of how well the whole phylofactorization fit: the residual deviance from regression of binned phylogenetic unit (BPU) on Z
binResidDev <- function(Y,Z){
  V <- NULL
  bResidDev <- numeric(length(Y$groups))
  for (i in 1:length(Y$groups)){
    V <- V %>% cbind(.,ilrvec(Y$groups[[i]],length(Z)))
    B <- bins(V)
    X <- numeric(length(Z))
    for (j in 1:length(B)){
      X[B[[j]]] <- j
    }
    X <- as.factor(X)
    bglm <- glm(Z~X,family = binomial(link='logit'))
    aa <- anova(bglm,test = 'Chisq')
    bResidDev[i] <- aa$`Resid. Dev`[2]
  }
  return(bResidDev)
}

## Another measure of 'total fit' - F-statistic from regression of BPU on Z
binF <- function(Y,Z){
  V <- NULL
  bF <- numeric(length(Y$groups))
  for (i in 1:length(Y$groups)){
    V <- V %>% cbind(.,ilrvec(Y$groups[[i]],length(Z)))
    B <- bins(V)
    X <- numeric(length(Z))
    for (j in 1:length(B)){
      X[B[[j]]] <- j
    }
    X <- as.factor(X)
    bglm <- glm(Z~X,family = binomial(link='logit'))
    s <- aov(bglm) %>% summary
    bF[i] <- s[[1]]['X','F value']
  }
  return(bF)
}


## Calculate Aggregate objective functions
BResidDev <- lapply(FACTORS,binResidDev,Z=Z)
BnullResidDev<- lapply(NullFactors,binResidDev,Z=Z)

BF <- lapply(FACTORS,binF,Z=Z)
BnullF<- lapply(NullFactors,binF,Z=Z)


## Plot stochastic phylofactorizations vs. null
par(mfrow=c(2,1))
ymin <- min(c(unlist(BResidDev),unlist(BnullResidDev)))
ymax <- max(c(unlist(BResidDev),unlist(BnullResidDev)))
for (r in 1:reps){
  if (r==1){
    plot(BResidDev[[1]],log='y',ylim=c(ymin,ymax))
  } else {
    lines(BResidDev[[r]])
  }
}
for (r in 1:reps){
  lines(BnullResidDev[[r]],col='grey',lwd=0.01)
}

M <- unlist(BResidDev) %>% matrix(.,nrow=reps,byrow=T)
Mnull <- unlist(BnullResidDev) %>% matrix(.,nrow=reps,byrow=T)

QRD <- apply(M,MARGIN=2,FUN=function(x) quantile(x,probs = c(0.025,0.975)))
QnullRD <- apply(Mnull,MARGIN=2,FUN=function(x) quantile(x,probs = c(0.025,0.975)))
lines(QnullRD[1,],col='green',lwd=4)
lines(QnullRD[2,],col='green',lwd=4)
lines(QRD[1,],col='red',lwd=4)
lines(QRD[2,],col='red',lwd=4)


ymin <- min(c(unlist(BF),unlist(BnullF)))
ymax <- max(c(unlist(BF),unlist(BnullF)))
for (r in 1:reps){
  if (r==1){
    plot(BF[[1]],log='y',ylim=c(ymin,ymax))
  } else {
    lines(BF[[r]])
  }
}
for (r in 1:reps){
  lines(BnullF[[r]],col='grey')
}


M <- unlist(BF) %>% matrix(.,nrow=reps,byrow=T)
Mnull <- unlist(BnullF) %>% matrix(.,nrow=reps,byrow=T)

QF <- apply(M,MARGIN=2,FUN=function(x) quantile(x,probs = c(0.025,0.975)))
QnullF <- apply(Mnull,MARGIN=2,FUN=function(x) quantile(x,probs = c(0.025,0.975)))
lines(QnullF[1,],col='green',lwd=4)
lines(QnullF[2,],col='green',lwd=4)
lines(QF[1,],col='red',lwd=4)
lines(QF[2,],col='red',lwd=4)

save(list=ls(),file='FisherFactorDemo_workspace')
