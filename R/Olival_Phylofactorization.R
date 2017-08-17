### Downloading Olival et al. dataset and scanning NCBI FTP genomes - who do we have genomes for?
library(RCurl)
library(ape)
library(phylofactor)
library(parallel)
source('R/fisherFactor.R')
parpf <- function(reps,pf,phat){
  Pvals <- matrix(NA,nrow=reps,ncol=pf$nfactors)
  for (nn in 1:reps){
    Z <- rbinom(Ntip(pf$tree),size=1,prob=phat)
    names(Z) <- pf$tree$tip.label
    pf <- fisherFactor(Z,pf$tree,nfactors=10)
    Pvals[nn,] <- pf$pvals
  }
  return(Pvals)
}

load('Workspaces/el_jefe.rdata')
VOlival <- read.csv('Data/HP3/data/viruses.csv')
tree <- el_jefe_grande  %>% di2multi ## resolve polytomies
all.equal(tree$tip.label,sapply(as.list(VOlival$vVirusNameCorrected),toString))




# Phylofactorization of IsZoonotic ----------------------------------------
Z <- VOlival$IsZoonotic
names(Z) <- tree$tip.label
phat <- mean(Z)     #fraction of viruses which are zoonotic in Olival data

pf.IsZoonotic <- fisherFactor(Z,tree,nfactors=10)
pf.IsZoonotic$factors

nms <- lapply(pf.IsZoonotic$bins,FUN=function(a,tree) tree$tip.label[a],tree=tree)


################ Null Simulation ##############
ncores=7
reps=50
cl <- phyloFcluster(ncores)
clusterExport(cl,varlist = c('fisherFactor','parpf'))
reps=as.list(rep(reps,ncores))
Ps <- parLapply(cl,reps,fun=parpf,pf=pf.IsZoonotic,phat=mean(VOlival$IsZoonotic))
stopCluster(cl)
for (i in 1:length(Ps)){
  if (i==1){Pvals <- Ps[[i]]} else{Pvals <- rbind(Pvals,Ps[[i]])}
}
rm('Ps')

tiff('Figures/Default_FisherFactors_of_IsZoonotic.tiff',height=800,width=1600)
  par(mfrow=c(1,2))
  binPhyloPlot(pf.IsZoonotic,factor=8,type='radial',use.edge.length=F,edge.width=2,main='Phylogenetic Factors of Spillover',cex.main=3,cex.lab=2)
  plot(pf.IsZoonotic$pvals,log='y',type='l',lwd=2,col='red',ylim=c(1e-8,1),xlab='Factor',ylab='P-value',main='Fisher-test Phylofactorization',cex.main=3,cex.lab=2,cex.axis=2)
  for (nn in 1:nrow(Pvals)){
    lines(Pvals[nn,],col=adjustcolor('black',alpha.f=0.1))
  }
  # matplot(t(Pvals),col='grey',type='l')
  lines(pf.IsZoonotic$pvals,lwd=2,col='red')
  legend(4,1e-4,legend=c('Olival et al. Data','Null Simulations'),col=c('black','grey'),lwd=c(2,1))
dev.off()

Clades <- lapply(pf.IsZoonotic$groups,FUN=function(g) g[[1]])
Clades <- c(Clades,pf.IsZoonotic$groups[[10]][2])
CladeNames <- lapply(Clades,FUN=function(g,tree) tree$tip.label[g],tree=tree)

Z <- VOlival$IsZoonotic
probs <- lapply(Clades,FUN=function(b,Z) mean(Z[b]),Z)
## all very low probabilities
## factoring out the first monophyletic clades
## increases the fraction of viruses that are zoonotic
## by 50%
## the remainder has P=0.43. 
##Then, a bunch of Unavirus, equine encaphalitis etc., have P=0.64,
##finally, Bovine leukemia and three primate T=lymphotropic viruses have P=1

# library(xlsx)
# 
# for (i in 1:11){
#   if (i<11){
#     nm <- paste('Factor ',toString(i),', phat=',toString(probs[[i]]),sep='')
#   } else {
#     nm <- paste('Paraphyletic Remaidner, phat=',toString(probs[[i]]),sep='')
#   }
#   write.xlsx(VOlival[VOlival$vVirusNameCorrected %in% CladeNames[[i]],],file='Data/Olival_Taxonomic_phylofactors.xlsx',
#              append=(!(i==1)),sheetName = nm)
# }

########## creating phylobin variable ##############
DF <- VOlival
DF$phyloBin <- numeric(nrow(DF))
for (i in 1:11){
  DF$phyloBin[DF$vVirusNameCorrected %in% CladeNames[[i]]] <- i
}
DF$phyloBin <- as.factor(DF$phyloBin)

write.xlsx(DF,file='Data/Olival_Dataset_with_phyloBin_Variable.xlsx')
save(list=ls(),file='Virus_Phylofactorization_workspace')












# Phylofactorization by Human-Association ---------------------------------

Z <- VOlival$IsHoSa
names(Z) <- VOlival$vVirusNameCorrected
phat <- mean(Z)
pf.IsHoSa <- fisherFactor(Z,tree,nfactors=10)
class(pf.IsHoSa) <- 'phylofactor'
pf.IsHoSa$tree <- tree
HumanBins <- getBins(pf.IsHoSa) 
HumanBinNms <- lapply(HumanBins,FUN=function(a,b) b[a],b=names(Z))


########### Null factors
ncores=7
reps=100
cl <- phyloFcluster(ncores)
clusterExport(cl,varlist = c('fisherFactor','parpf'))
reps=as.list(rep(reps,ncores))
Ps <- parLapply(cl,reps,fun=parpf,pf=pf.IsHoSa,phat=mean(VOlival$IsHoSa))
stopCluster(cl)
for (i in 1:length(Ps)){
  if (i==1){Pvals <- Ps[[i]]} else{Pvals <- rbind(Pvals,Ps[[i]])}
}
rm('Ps')


tiff('Figures/Default_FisherFactors_of_IsHoSa.tiff',height=800,width=1600)
  par(mfrow=c(1,2))
  binPhyloPlot(pf.IsHoSa,factor=8,type='radial',use.edge.length=F,edge.width=2,main='Phylogenetic Factors of Human Viruses',cex.main=3,cex.lab=2)
  plot(pf.IsHoSa$pvals,log='y',type='l',lwd=2,col='red',ylim=c(1e-8,1),xlab='Factor',ylab='P-value',main='Fisher-test Phylofactorization',cex.main=3,cex.lab=2,cex.axis=2)
  for (nn in 1:reps){
    lines(HPvals[nn,],col='grey',col=adjustcolor('black',alpha.f=0.1))
  }
  lines(pf.IsHoSa$pvals,lwd=2,col='red')
  legend(4,1e-4,legend=c('Olival et al. Data','Null Simulations'),col=c('black','grey'),lwd=c(2,1))
dev.off()


hosts <- read.csv('Data/HP3/data/hosts.csv')
cytree <- read.tree('Data/HP3/data/cytb_supertree.tree')    #not sure what this is
stree <- read.tree('Data/HP3/data/supertree_mammals.tree')  #mammalian supertree
htree <- drop.tip(stree,setdiff(stree$tip.label,hosts$hHostNameFinal))


























# Genomes -----------------------------------------------------------------


genomes.url <- 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/'
filenames = getURL(genomes.url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
filenames <- strsplit(filenames, "\r\n")
filenames = unlist(filenames)
filenames

as.list(VOlival$vVirusNameCorrected) %>% 
  sapply(.,FUN=function(v) any(grepl(v,filenames))) %>% sum

329/586   #percent of viruses for which we have genomes - about 56%. 

