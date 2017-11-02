### Downloading Olival et al. dataset and scanning NCBI FTP genomes - who do we have genomes for?
library(RCurl)
library(ape)
library(phylofactor)
library(parallel)
library(dplyr)
library(ggpubr)
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

parf <- function(reps,pf,phat){
  Fs <- matrix(NA,nrow=reps,ncol=pf$nfactors)
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




# Visualization -----------------------------------------------------------

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

GG <- pf.tree(pf,group.map=data.frame(1:8,rep(1,8)))
ddf <- data.frame('Pvals'=c(t(Pvals)),'factor'=rep(1:10,nrow(Pvals)),'replicate'=rep(1:nrow(Pvals),each=10))
obs <- data.frame('Pvals'=pf.IsZoonotic$pvals,'factor'=1:10,'replicate'=NA)
ddf <- rbind(ddf,obs)

cols <- c('Null Simulations'='black','Observed Phylofactorization'='red')

PLOT <- ggplot(ddf,aes(factor,Pvals))+
          stat_summary(data = ddf[!is.na(ddf$replicate),],
                       geom='ribbon',
                       fun.ymin=function(x) quantile(x,0.025),
                       fun.ymax=function(x) quantile(x,0.975),
                       alpha=0.2,lwd=2)+
          stat_summary(data = ddf[!is.na(ddf$replicate),],
                       aes(x=factor,y=Pvals,color='Null Simulations'),geom='line',fun.y = mean,lwd=2)+
          geom_line(data=ddf[is.na(ddf$replicate),],
                    aes(y=Pvals,color='Observed Phylofactorization'),lwd=2)+
          scale_y_continuous(trans='log',breaks = 10^(c(-8:0)))+
          scale_x_continuous(breaks=1:10)+
          scale_fill_continuous(cols)+
          theme(legend.position = c(.5,.5),title = element_text(size=18))+
          ggtitle('Phylofactorization - observed vs. null')

GG <- pf.tree(pf.IsZoonotic,
              group.map = data.frame(1:8,rep(1,8)),
              bg.color='grey',
              bg.alpha=0.5)$ggplot
  ggarrange(GG,PLOT,ncol = 2,nrow=1)
ggsave(filename='Figures/FisherFactors_of_IsZoonotic.tiff',height=8,width=16)


Clades <- lapply(pf.IsZoonotic$groups[1:8],FUN=function(g) g[[1]])
Clades <- c(Clades,pf.IsZoonotic$groups[[8]][2])
CladeNames <- lapply(Clades,FUN=function(g,tree) tree$tip.label[g],tree=tree)

Z <- VOlival$IsZoonotic
alphas=rep(0.7,pf$nfactors+1)
phat <- mean(Z)
col <-colorRamp(c('blue','red'))(phat)/255
col <- rgb(col[1],col[2],col[3])
ggZ <- ggtree(tree,layout='circular')+
  geom_hilight(Ntip(tree)+1,fill=col,alpha=0.7)+
  ggtitle('0 Factors')+
  theme(title=element_text(size=18))

probs <- sapply(Clades,FUN=function(b,Z) mean(Z[b]),Z)
mycolors <- colorRamp(c('blue','red'))(probs)/255
mycolors <- apply(mycolors,1,FUN=function(x) rgb(x[1],x[2],x[3]))
colfcn <- function(n){
  return(mycolors)
}
gg <- pf.tree(pf.IsZoonotic,Grps=Clades,color.fcn = colfcn,rootnode=T,alphas=alphas)$ggplot
gg
gg <- gg+ggtitle('8-Factors')+
  theme(title = element_text(size=18))+
  scale_fill_gradient(low='blue',high='red')

probs2 <- sapply(pf.IsZoonotic$bins,FUN=function(a,Z) mean(Z[a]),Z=Z)
mycolors <- colorRamp(c('blue','red'))(probs2)/255
mycolors <- apply(mycolors,1,FUN=function(x) rgb(x[1],x[2],x[3]))
colfcn <- function(n){
  return(mycolors)
}
gg2 <- pf.tree(pf.IsZoonotic,Grps=pf.IsZoonotic$bins,color.fcn=colfcn,rootnode=T,alphas=alphas)$ggplot
gg2 <- gg2+ggtitle('10-factors')+
  theme(title = element_text(size=18))+
  scale_fill_gradient(low='blue',high='red')

ggarrange(ggZ,gg,gg2,ncol=3,nrow=1)
ggsave(filename = 'Figures/Fraction_zoonotic_fisherFactorization.tiff',height=6,width=16)


# Appending & Writing phyloBin results ------------------------------------

########## creating phylobin variable ##############
DF <- VOlival
DF$phyloBin <- numeric(nrow(DF))
for (i in 1:11){
  DF$phyloBin[DF$vVirusNameCorrected %in% CladeNames[[i]]] <- i
}
DF$phyloBin <- as.factor(DF$phyloBin)

write.xlsx(DF,file='Data/Olival_Dataset_with_phyloBin_Variable.xlsx')
save(list=ls(),file='Virus_Phylofactorization_workspace')




# Phylofactorization by stringent Zoon ------------------------------------
library(phylofactor)
library(parallel)
library(dplyr)
library(ggpubr)
source('R/fisherFactor.R')
load('Virus_Phylofactorization_workspace')
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

Z <- VOlival$IsZoonotic.stringent
names(Z) <- tree$tip.label
phat <- mean(Z)     #fraction of viruses which are zoonotic in Olival data

pf.Stringent <- fisherFactor(Z,tree,nfactors=10)

ncores=7
reps=50
cl <- phyloFcluster(ncores)
clusterExport(cl,varlist = c('fisherFactor','parpf'))
reps=as.list(rep(reps,ncores))
Ps <- parLapply(cl,reps,fun=parpf,pf=pf.Stringent,phat=phat)
stopCluster(cl)
for (i in 1:length(Ps)){
  if (i==1){Pvals.Stringent <- Ps[[i]]} else{Pvals.Stringent <- rbind(Pvals.Stringent,Ps[[i]])}
}
rm('Ps')



GG <- pf.tree(pf.Stringent)
ddf <- data.frame('Pvals'=c(t(Pvals.Stringent)),'factor'=rep(1:10,nrow(Pvals.Stringent)),'replicate'=rep(1:nrow(Pvals.Stringent),each=10))
obs <- data.frame('Pvals'=pf.Stringent$pvals,'factor'=1:10,'replicate'=NA)
ddf <- rbind(ddf,obs)

cols <- c('Null Simulations'='black','Observed Phylofactorization'='red')

PLOT <- ggplot(ddf,aes(factor,Pvals))+
  stat_summary(data = ddf[!is.na(ddf$replicate),],
               geom='ribbon',
               fun.ymin=function(x) quantile(x,0.05),
               fun.ymax=function(x) quantile(x,1),
               alpha=0.2,lwd=2)+
  stat_summary(data = ddf[!is.na(ddf$replicate),],
               aes(x=factor,y=Pvals,color='Null Simulations'),geom='line',fun.y = mean,lwd=2)+
  geom_line(data=ddf[is.na(ddf$replicate),],
            aes(y=Pvals,color='Observed Phylofactorization'),lwd=2)+
  scale_y_continuous(trans='log',breaks = 10^(c(-8:0)))+
  scale_x_continuous(breaks=1:10)+
  scale_fill_continuous(cols)+
  theme(legend.position = c(.5,.5),title = element_text(size=18))+
  ggtitle('Phylofactorization - observed vs. null')

GG <- pf.tree(pf.Stringent,
              bg.color='grey',
              bg.alpha=0.5)$ggplot
ggarrange(GG,PLOT,ncol = 2,nrow=1)
ggsave(filename='Figures/FisherFactors_of_IsZoonotic_Stringent_small.tiff',height=4,width=8)


Clades <- lapply(pf.Stringent$groups[1:10],FUN=function(g) g[[1]])
Clades <- c(Clades,pf.Stringent$groups[[9]][2])
CladeNames <- lapply(Clades,FUN=function(g,tree) tree$tip.label[g],tree=tree)

probs <- sapply(Clades,FUN=function(b,Z) mean(Z[b]),Z)
cbind(1:11,probs,sapply(Clades,length),sapply(Clades,length)*probs)
## the final, paraphyletic bin is only 22% zoonotic.
## We have partitioned a great deal of non-zoonotic clades, much like before.

getName <- function(Clade,VOlival){
  V <- VOlival[Clade,2:5,drop=F]
  ss=NULL
  if (nrow(V)==1){
    tx <- names(V)
    tx <- sapply(tx,FUN=function(x) gsub('v','',x))
    V <- sapply(V,toString)
    V <- mapply(paste,tx,V,sep='_')
    ss <- paste(V,collapse='; ')
  } else {
    for (i in ncol(V):1){
      if (length(unique(V[,i]))==1){
        if (!is.na(V[1,i])){
          tx <- gsub('v','',names(V[1,])[i])
          ss <- paste(tx,toString(V[1,i]),sep='_')
          break
        }
      }
    }
  }
  if (is.null(ss)){
    ss <- 'No Common Name'
  }
  return(ss)
}

sCladeNames <- sapply(Clades,getName,VOlival=VOlival)
output <- cbind(sCladeNames,probs,sapply(Clades,length),sapply(Clades,length)*probs)
colnames(output) <- c('Name','Fraction_Zoonotic','no_viruses','no_zoonotic')
write.csv(output,file='Olival_Stringent_phylogenetic_bins.csv')

save(list=ls(),file='Virus_Phylofactorization_workspace')









# Appending VOlival with phylobins ----------------------------------------
library(phylofactor)
library(xlsx)
load('Virus_Phylofactorization_workspace')


DF <- VOlival

for (pp in 1:2){
  if (pp==1){
    DF$pfIsZ <- numeric(nrow(DF))
    nm <- 'pfIsZ'
    B <- pf.IsZoonotic$bins
  } else {
    DF$pfIsZS <- numeric(nrow(DF))
    nm <- 'pfIsZS'
    B <- bins(pf.Stringent$basis[,1:9])
  }
  
  for (i in 1:length(B)){
    DF[,nm][1:nrow(DF) %in% B[[i]]] <- i
  }
  DF[,nm] <- as.factor(DF[,nm])
}

### now we append the nonZ clades, i.e. those in bins. In this case, we remove all those in bins 2-6
pf <- pf.IsZoonotic
sapply(pf$groups[1:5],FUN=function(g,Z) mean(Z[g[[1]]]), Z=VOlival$IsZoonotic)

nonZ <- sapply(pf$groups[1:5],FUN=function(g) g[[1]]) %>% unlist

DF$nonZ <- numeric(nrow(DF))
DF$nonZ[nonZ] <- 1

write.xlsx(DF,file='Data/Olival_w_phylobin.xlsx')

































# Phylofactorization by Human-Association ---------------------------------
# 
# Z <- VOlival$IsHoSa
# names(Z) <- VOlival$vVirusNameCorrected
# phat <- mean(Z)
# pf.IsHoSa <- fisherFactor(Z,tree,nfactors=10)
# class(pf.IsHoSa) <- 'phylofactor'
# pf.IsHoSa$tree <- tree
# HumanBins <- getBins(pf.IsHoSa) 
# HumanBinNms <- lapply(HumanBins,FUN=function(a,b) b[a],b=names(Z))
# 
# 
# ########### Null factors
# ncores=7
# reps=100
# cl <- phyloFcluster(ncores)
# clusterExport(cl,varlist = c('fisherFactor','parpf'))
# reps=as.list(rep(reps,ncores))
# Ps <- parLapply(cl,reps,fun=parpf,pf=pf.IsHoSa,phat=mean(VOlival$IsHoSa))
# stopCluster(cl)
# for (i in 1:length(Ps)){
#   if (i==1){Pvals <- Ps[[i]]} else{Pvals <- rbind(Pvals,Ps[[i]])}
# }
# rm('Ps')
# 
# 
# tiff('Figures/Default_FisherFactors_of_IsHoSa.tiff',height=800,width=1600)
#   par(mfrow=c(1,2))
#   binPhyloPlot(pf.IsHoSa,factor=8,type='radial',use.edge.length=F,edge.width=2,main='Phylogenetic Factors of Human Viruses',cex.main=3,cex.lab=2)
#   plot(pf.IsHoSa$pvals,log='y',type='l',lwd=2,col='red',ylim=c(1e-8,1),xlab='Factor',ylab='P-value',main='Fisher-test Phylofactorization',cex.main=3,cex.lab=2,cex.axis=2)
#   for (nn in 1:reps){
#     lines(HPvals[nn,],col='grey',col=adjustcolor('black',alpha.f=0.1))
#   }
#   lines(pf.IsHoSa$pvals,lwd=2,col='red')
#   legend(4,1e-4,legend=c('Olival et al. Data','Null Simulations'),col=c('black','grey'),lwd=c(2,1))
# dev.off()
# 
# 
# hosts <- read.csv('Data/HP3/data/hosts.csv')
# cytree <- read.tree('Data/HP3/data/cytb_supertree.tree')    #not sure what this is
# stree <- read.tree('Data/HP3/data/supertree_mammals.tree')  #mammalian supertree
# htree <- drop.tip(stree,setdiff(stree$tip.label,hosts$hHostNameFinal))
# 
# 
# 
# 
# 




















# 
# # Genomes -----------------------------------------------------------------
# 
# # 
# # genomes.url <- 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/'
# # filenames = getURL(genomes.url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
# # filenames <- strsplit(filenames, "\r\n")
# # filenames = unlist(filenames)
# # filenames
# 
# ncbi_viruses <- read.table('Data/NCBI_viral_genomes.tbl',fill = T,header=T)
# 
# as.list(VOlival$vVirusNameCorrected) %>% 
#   sapply(.,FUN=function(v) any(grepl(v,filenames))) %>% sum
# 
# 329/586   #percent of viruses for which we have genomes - about 56%. 
# 
# ix <- as.list(VOlival$vVirusNameCorrected) %>% 
#   sapply(.,FUN=function(v) any(grepl(v,filenames))) %>% which
# 
# genomes <- numeric(Ntip(tree))
# genomes[ix]=2
# cols <- rep('white',length(genomes))
# cols[ix] <- 'blue'
# names(genomes) <- VOlival$vVirusNameCorrected
# gng <- ggtree(tree,layout='circular')+
#   geom_tippoint(size=genomes,col=cols)+
#   ggtitle('Genomes available on NCBI')+
#   theme(title=element_text(size=18))
#   
# 
# gng
# ggsave(filename='Figures/Genomes_available_on_NCBI.tiff',height=8,width=8)
