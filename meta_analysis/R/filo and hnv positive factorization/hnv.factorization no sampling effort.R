library(phylofactor)
library(parallel)
setwd('/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis')

load('Z.hnv.pos')
load('bat_tree.hnv.pos')

tree <- bat_tree.hnv.pos
Data <- Z.hnv_pos.hnv1


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


plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,25))
for (rr in 1:140){
  lines(c(S[rr,1],diff(S[rr,])),col=rgb(0,0,0,0.2))
}

points(deviances,col='red',pch=16,cex=2)

ecdf(S[,1])(Obj[1])

save(list=ls(),file='hnv_workspace_no_sampling_effort')

load(file='R/filo and hnv positive factorization/hnv_workspace_no_sampling_effort')



#............................define our colors....................................................
#.................................................................................................
#............................1:  Yangochiroptera....."#FF0000FF"..................................
#............................2:  Myotis.............."#FF9900FF"..................................
#............................3:  Miniopterus........."#CCFF00FF"..................................
#............................4:  Pteropodidae........"#33FF00FF"..................................
#............................5:  Emballonuroidea....."#00FF66FF"..................................
#............................6:  Rhinolophus........."#00FFFFFF"..................................
#............................7:  Pipistrellini......."#0066FFFF"..................................
#............................8:  Pteropus............"#3300FFFF"..................................
#............................9:  Pteropodinae........"#CC00FFFF".................................
#............................10: Scotophilus kuhlii.."#FF0099FF"..................................
#............................11: Eidolon............."#FF0099FF"..................................
#............................12: Hipposideros........"#00FF66FF"..................................
#............................13: C. perspicillata...."#00FFFFFF"..................................
#............................14: D. rotundus........."#00FFFFFF"..................................
#............................14: Carollia............"#FF9900FF"..................................

names.storage <- list()

for (i in 2:(11))
{
  indexes = pf$bins[[i]]
  species <- gsub("_", " ", tolower(tree$tip.label[indexes]))
  group_taxonomy_list <- as.data.frame(taxonomy[match(species,taxonomy[,1]),2])
  names.storage[i-1] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
  print(i)
}


B <- bins(pf$basis[,1:10])
B <- B[2:11]
Z <- Data$Z
probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

colfcn <- function(n) return(c("#CC00FFFF"))

pf.tree <- pf.tree(pf, lwd=1, factors = 1, color.fcn=colfcn, branch.length = "none")
pf.tree$ggplot +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue')  

Legend <- pf.tree$legend
Legend$names <- names.storage[1]
P <- sapply(probs[1],FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend$names <- mapply(paste,Legend$names,P)
plot.new()
plot.new()
legend('topleft',legend=Legend$names,fill=Legend$colors,cex=1)






