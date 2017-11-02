############ Olival cb_dist_noHoSa_maxLn factorization
setwd('~/Bozeman/Spillover_Phylofactorization/')
library(phylofactor)
library(xlsx)
library(stringi)
library(dplyr)
library(parallel)
library(ggplot2)
library(ggpubr)
library(viridis)
source('R/fisherFactor.R')
source('R/visualization_fcns.R')
source("~/Bozeman/BZDEL/Data/Olival/R/model_reduction.R")
source("~/Bozeman/BZDEL/Data/Olival/R/fit_gam.R")
source("~/Bozeman/BZDEL/Data/Olival/R/relative_contributions.R")
source("~/Bozeman/BZDEL/Data/Olival/R/cross_validation.R")
source("~/Bozeman/BZDEL/Data/Olival/R/logp.R")
parpf <- function(reps,pf){
  Pvals <- matrix(NA,nrow=reps,ncol=pf$nfactors)
  for (nn in 1:reps){
    Z <- sample(pf$Data,size=length(pf$Data),replace=T)
    names(Z) <- pf$tree$tip.label
    pf <- fisherFactor(Z,pf$tree,nfactors=9,TestFunction = WilcoxTest)
    Pvals[nn,] <- pf$pvals
  }
  return(Pvals)
}


load('Workspaces/el_jefe.rdata')
db <- readRDS("~/Bozeman/BZDEL/Data/Olival/intermediates/postprocessed_database.rds")
tree <- el_jefe_grande  %>% di2multi ## resolve polytomies
all.equal(tree$tip.label,sapply(as.list(db$viruses$vVirusNameCorrected),toString))

hosts <- db$hosts
viruses <- db$viruses
associations <- db$associations
rm(db)
PD_centers = names(viruses)[stri_detect_regex(names(viruses), "^(cb|st)_.*(?<!stringent)$")]

data_set = viruses %>%
  mutate_each_(funs("logp"), vars=PD_centers)
names(data_set)[names(data_set) %in% PD_centers] <- paste0(PD_centers, "Ln")

Z <- data_set$cb_dist_noHoSa_maxLn
names(Z) <- data_set$vVirusNameCorrected


WilcoxTest <- function(grps,tree,Z){
  s1 <- Z[tree$tip.label[grps[[1]]]]
  s2 <- Z[tree$tip.label[grps[[2]]]]
  pval <- tryCatch(wilcox.test(s1,s2)$p.value,
                   error=function(e) 1)
  return(pval)
}

pf.cb_dist_noHoSa_maxLn <- fisherFactor(Z,tree,nfactors=9,TestFunction = WilcoxTest)



################ Null Simulation ##############
ncores=7
reps=50
cl <- phyloFcluster(ncores)
clusterExport(cl,varlist = c('fisherFactor','parpf','WilcoxTest'))
reps=as.list(rep(reps,ncores))
Ps <- parLapply(cl,reps,fun=parpf,pf=pf.cb_dist_noHoSa_maxLn)
stopCluster(cl)
rm('cl')
for (i in 1:length(Ps)){
  if (i==1){Pvals <- Ps[[i]]} else{Pvals <- rbind(Pvals,Ps[[i]])}
}
rm('Ps')



################ Visualization #################

GG <- pf.tree(pf.cb_dist_noHoSa_maxLn)
ddf <- data.frame('Pvals'=c(t(Pvals)),'factor'=rep(1:pf.cb_dist_noHoSa_maxLn$nfactors,nrow(Pvals)),'replicate'=rep(1:nrow(Pvals),each=pf.cb_dist_noHoSa_maxLn$nfactors))
obs <- data.frame('Pvals'=pf.cb_dist_noHoSa_maxLn$pvals,'factor'=1:9,'replicate'=NA)
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

GG <- pf.tree(pf.cb_dist_noHoSa_maxLn,
              bg.color='grey',
              bg.alpha=0.5)$ggplot
ggarrange(GG,PLOT,ncol = 2,nrow=1)
ggsave(filename='Figures/Wilcox_factorization_cb_dist.tiff',height=8,width=16)




# Bin-differences of cb_dist ----------------------------------------------

B <- bins(pf.cb_dist_noHoSa_maxLn$basis)
B_means <- sapply(B,FUN=function(b,z) mean(z[b],na.rm=T),z=Z)
nms <- sapply(B,getName,VOlival=data_set)
nms[5] <- "Order_Mononegavirales"
cbind(B_means,nms)
nm.order <- nms[order(B_means,decreasing = T)]

colfcn <- function(n)  rainbow(10)
GG <- pf.tree(pf.cb_dist_noHoSa_maxLn,Grps = B,color.fcn = colfcn,top.layer = 5,top.alpha = 0.6,bg.alpha = 1)
GTREE <- GG$ggplot+
  ggtitle('Phylogenetic Breadth Factorization')

phylobin <- function(B){
  bn <- numeric(sum(sapply(B,length)))
  for (i in 1:length(B)){
    bn[B[[i]]] <- i
  }
  return(bn)
}

DF <- data.frame('Phylogenetic_Breadth'=Z,'name'=factor(nms[phylobin(B)],levels = rev(nm.order)))

cols <- rainbow(11)[order(B_means,decreasing = F)]
alphas <- rep()
PLOT <- ggplot(DF,aes(x=name,y=Phylogenetic_Breadth,fill=name))+
          geom_boxplot()+
          scale_fill_manual(values=cols,guide=guide_legend(reverse=T))+
          coord_flip()+
          ggtitle('Phylogenetic Breadth')

ggarrange(GTREE,PLOT,ncol=2)
ggsave('~/Bozeman/BZDEL/Figures/Olival phylofactorization/cb_dist_Phylofactorization.tiff',height=4,width=12)
save(list=ls(),file='cb_dist_phylofactorization')
save(list=ls(),file='~/Bozeman/BZDEL/Data/Olival/cb_dist_phylofactorization')  ##let's also put a copy in BZDEL

