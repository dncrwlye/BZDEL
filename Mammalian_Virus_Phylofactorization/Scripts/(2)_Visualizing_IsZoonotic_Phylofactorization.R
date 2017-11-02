library(phylofactor)
library(ggtree)
source('R/visualization_fcns.R')

load('Virus_Phylofactorization_workspace')

pf <- pf.IsZoonotic
Z <- VOlival$IsZoonotic
# pf <- pf.Stringent
# Z <- VOlival$IsZoonotic.stringent

pp=pf.tree(pf)

Taxonomy <- cbind(VOlival$vVirusNameCorrected,VOlival[,c(2:5,1)])
Taxonomy[Taxonomy=='Unassigned']=NA
Taxonomy <- apply(Taxonomy,2,FUN=function(x) sapply(x,toString))
IDs <- Taxonomy[,1]
length(unique(IDs))==length(IDs)
colnames(Taxonomy)[c(1,6)] <- 'ID'

colnames(Taxonomy) <- gsub('v','',colnames(Taxonomy))



# nms <- sutax(CladeNames,Taxonomy,trim=T)
B <- bins(pf$basis[,1:10])
B <- B[2:11]
nms <- sapply(B,getName,VOlival=VOlival)
probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
# nms[1] <- 'Remainder'
names(nms) <- probs
if (all.equal(Z,VOlival$IsZoonotic.stringent)){
  nms[10] <- 'Genus_Lyphocryptovirus'
}

Legend <- pp$legend
Legend$names <- nms
P <- sapply(probs,FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend$names <- mapply(paste,Legend$names,probs)
plot.new()
legend('topleft',legend=Legend$names,fill=Legend$colors,cex=1)





# The Plot ----------------------------------------------------------------


library(grid)
library(gridBase)
library(ggplot2)
gtr <- pp$ggplot+geom_tippoint(size=3*Z,col='black')

plot.new()
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,2)))

pushViewport(viewport(layout.pos.col = 1))
print(gtr, newpage = FALSE)
popViewport()

#Draw bsae plot
pushViewport(viewport(layout.pos.col = 2))
par(fig = gridFIG(), new = TRUE)
plot.new()
legend('topleft',legend=Legend$names,fill=Legend$colors,cex=1.8)
popViewport()
# ggsave('Figures/IsZoonotic_fig.png')
# ggsave('Figures/IsZoonotic_Stringent_fig.png')








# Mammal body mass ---------------------------------------------------------
rm(list=ls())
source('R/visualization_fcns.R')
source('R/fisherFactor.R')
##### Can we replicate this for mammalian body mass?

library(phylofactor)
library(ape)
library(ggplot2)

htree <- read.tree('Data/HP3/data/supertree_mammals.tree')
hmass <- read.csv('Data/HP3/data/intermediate/PVR_cytb_hostmass.csv')
hosts <- read.csv('Data/HP3/data/hosts.csv')

htree <- drop.tip(htree,setdiff(htree$tip.label,hmass$hHostNameFinal))
# Tax <- data.frame('Species IDs'=hosts$hHostNameFinal,
#                   'Taxonomy'= apply(hosts[,2:5],1,
#                                     FUN=function(x)
#                                       paste('o_',x[1],'; ',
#                                             'f_',x[2],'; ',
#                                             'g_',x[3],'; ',
#                                             's_',x[4],'; ',
#                                             sep='')))
# Tax <- Tax[Tax$Species.IDs %in% htree$tip.label,]
Tax <- as.data.frame(cbind(hosts[,1],hosts[,c(2:5,1)]))
colnames(Tax) <- c('ID',gsub('h','',colnames(Tax[2:6])))

Z <- log(hmass$hMassGrams)
names(Z) <- hmass$hHostNameFinal
Z <- Z[htree$tip.label]

TestFunction <- function(grps,tree,Z){
  v <- ilrvec(grps,length(Z))
  return(1/abs(v %*% as.matrix(Z,ncol=1)))
}

pf.mass <- fisherFactor(Z,htree,nfactors=5,TestFunction = TestFunction)

Clades <- bins(pf.mass$basis)
CladeNames <- lapply(Clades,FUN=function(ix,tree) tree$tip.label[ix],tree=pf.mass$tree)
nms <- sutax(CladeNames,Tax,trim=T)

pp=pf.tree(pf.mass,layout='fan',alphas=1,method = 'factors')
# pp$ggplot

sapply(Clades,length)
### Clade--Factor map: Clades[1,2,3,4]=factor[4,1]
#by clades:
# nms <- c('Remainder','Fereuungulata','Primates','Boreoeutheria')
# by factor:
nms <- c('Fereuungulata','Primates','Boreoeutheria')
M <- sapply(pf.mass$groups,FUN=function(grp,Z) median(exp(Z[grp[[1]]]))/1e3,Z=Z)

Legend <- pp$legend
Legend$names <- sapply(nms,FUN=function(x) if(length(x)>1){paste(x,collapse='\n')}else{x})
# M <- lapply(Clades,FUN=function(x,Z) mean(exp(Z[x]))/1e3,Z=Z)
M <- sapply(M[1:3],FUN=function(x,ls) paste('mass=',toString(signif(x,digits = 2)),'Kg',sep=''),ls=ls)
fn <- function(y,x){
  ls <- nchar(x)+nchar(y)
  paste(y,paste(rep(' ',30-ls),collapse=''),x,sep=' ')
}
Legend$names <- mapply(fn,Legend$names,M)


library(grid)
library(gridBase)
library(ggplot2)
gtr <- pp$ggplot+geom_tippoint(size=6*(Z/max(Z))^2,col='black')

plot.new()
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,2)))

pushViewport(viewport(layout.pos.col = 1))
print(gtr, newpage = FALSE)
popViewport()

#Draw bsae plot
pushViewport(viewport(layout.pos.col = 2))
par(fig = gridFIG(), new = TRUE)
plot.new()
legend('left',legend=Legend$names,fill=Legend$colors,cex=1.8)
popViewport()



# higher factors ----------------------------------------------------------

pf.mass <-  fisherFactor(Z,htree,nfactors=5,TestFunction = TestFunction)
pf.mass$factors
B <- bins(pf.mass$basis)
sapply(B,length)
nms <- lapply(B,FUN=function(a,tr) tr$tip.label[a],tr=htree)
Bnms <- sutax(nms,Tax)
