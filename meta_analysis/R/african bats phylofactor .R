setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")


library(phylofactor)
library(grid)
library(gridBase)
library(stringi)
#source('R/other phylofactor scripts/fisherFactor.R')
source('R/taxonomy group name function.R')
#source('R/other phylofactor scripts/visualization_fcns.R')
library(tidyverse)
library(taxize)

# WilcoxTest <- function(grps,tree,Z){
#   s1 <- na.omit(Z[grps[[1]]])
#   s2 <- na.omit(Z[grps[[2]]])
#   wt <- tryCatch(wilcox.test(s1,s2)$p.value,error=function(e) 1)
#   return(wt)
# }


#..................which african bats have been sampled for filoviruses
#load('data/bat_taxonomy_data.Rdata') 
load('data/bat_taxonomy_data.Rdata') #oops typo

#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.......................okay going to filter out the results for which we dont have taxonomies
#.......................but really we shouldn't do this and need to figure out the missing taxonomies 

batphy1.africa <- dplyr::filter(batphy1, region == 'Africa') 
batphy1.africa<-batphy1.africa %>%
  mutate(index = row.names(batphy1.africa))

Z <- batphy1.africa$filo_surv

n_factors = 7
#pf <- fisherFactor(Z, africa_tree,10,TestFunction=WilcoxTest,ncores = 1)

pf <- twoSampleFactor(Z, africa_tree, method='Fisher', n_factors ,ncores = 1)

taxonomy_africa<- batphy1.africa[,c('species','tax')]

#names(smry)
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#......................START: alex's slack suggestion.............................................

nms <- list()
  
for (i in 2:length(pf$bins))
{
  indexes = pf$bins[[i]]
  species <- gsub("_", " ", tolower(africa_tree$tip.label[indexes]))
  group_taxonomy_list <- as.data.frame(taxonomy_africa[match(species,taxonomy_africa[,1]),2])
  
  nms[i] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))

  print(i)
}
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#......................END: alex's slack suggestion...............................................

B <- bins(pf$basis[,1:n_factors])
B <- B[2:n_factors] 
nms1 <- nms[2:(n_factors+1)]

## remove paraphyletic bin

probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
names(nms1) <- probs


#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#...............Null Simulation..........................................................

null_simulations <- pf.nullsim(pf, reps = 50, nfactors = n_factors, output ='pvals')

#nfactors <- pf$nfactors
obs <- data.frame('Pvals'=pf$pvals,
                  'factor'=1:n_factors
                  )
null_simulations_matrix<- matrix(0,50,n_factors)
for (i in 1:50)
{
  null_simulations_matrix[i,] <- null_simulations[[i]]$pvals
}
null_simulations_matrix <- as.data.frame(t(null_simulations_matrix))
ddf <- cbind(obs,null_simulations_matrix)

ddf %>%
  tidyr::gather(group, pvalue, c(1,3:22)) %>%
  mutate(color = ifelse(group =='Pvals', "observed", 'simulated')) %>%
  ggplot() +
  geom_line(aes(x = factor, y= log(pvalue), group = group, color = color))

#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
# The Plot ---------------------------------------------------------------

pf.final <- twoSampleFactor(Z, africa_tree, 1 ,ncores = 1) 

pp=pf.tree(pf.final)

gtr <- pp$ggplot+
  ggtree::geom_tippoint(size=3*Z,col='blue') #+ 
  #ggtree::geom_tiplab(aes(angle=angle))

plot.new()
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,2)))
pushViewport(viewport(layout.pos.col = 1))
print(gtr, newpage = FALSE)
popViewport()

#Draw bsae plot
pushViewport(viewport(layout.pos.col = 2))
#par(fig = gridFIG(), new = TRUE)
plot.new()
legend('topleft',legend=Legend[1,'names'],fill=Legend[1,'colors'],cex=1.8)
popViewport()
ggsave('Figures/Sampled_fig.png', gtr, width=15, height=15)

#....................................................................................
#....................................................................................
#....................................................................................
#....................................................................................
#....................................................................................
#....................................................................................
#....................................................................................
#....................................................................................
#....................................................................................
#....................................................................................
#....................................................................................
#....................................................................................
#.................comparison of henipa and filo......................................

batphy1 <- batphy1 %>%
  mutate(hnv_surv = ifelse(is.na(hnv_surv), 0, hnv_surv)) 

taxonomy<- batphy1[,c('species','tax')]

n_factors = 10

#........................henipaviruses...............................

n_factors.hvn = 7

Z.hnv <- batphy1$hnv_surv
pf.hvn <- phylofactor::twoSampleFactor(Z.hnv, bat_tree, method='Fisher', n_factors.hvn ,ncores = 1)

nms.hnv <- list()
for (i in 2:(n_factors.hvn+1))
{
  indexes = pf.hvn$bins[[i]]
  species <- gsub("_", " ", tolower(bat_tree$tip.label[indexes]))
  group_taxonomy_list <- as.data.frame(taxonomy[match(species,taxonomy[,1]),2])
  nms.hnv[i-1] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
  print(i)
}

B.hnv <- bins(pf.hvn$basis[,1:n_factors.hvn])
B.hnv <- B.hnv[2:(n_factors.hvn+1)] 
nms1.hvn <- nms.hnv[1:(n_factors.hvn)]

## remove paraphyletic bin

probs.hnv <- sapply(B.hnv,FUN=function(ix,Z.hnv) mean(Z.hnv[ix]),Z.hnv=Z.hnv) %>% signif(.,digits=2)
names(nms1.hvn) <- probs.hnv

pf.tree.hnv <- pf.tree(pf.hvn, lwd=0.5, , branch.length = "none")
pf.tree.hnv$ggplot +
  ggtree::geom_tippoint(size=3*Z.hnv,col='blue')  

Legend.hnv <- pf.tree.hnv$legend
Legend.hnv$names <- nms1.hvn
P.hnv <- sapply(probs.hnv,FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend.hnv$names <- mapply(paste,Legend.hnv$names,P.hnv)
plot.new()
legend('topleft',legend=Legend.hnv$names,fill=Legend.hnv$colors,cex=1)

#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#...............Null Simulation..........................................................

#null_simulations <- pf.nullsim(pf.hvn, reps = 50, nfactors = n_factors, output ='pvals')
null_simulations.hvn <- readRDS('data/my_phylofactor_object_hvn_nullsim')

#nfactors <- pf$nfactors
obs.hvn <- data.frame('Pvals'=pf.hvn$pvals,
                  'factor'=1:n_factors
)
null_simulations_matrix.hvn<- matrix(0,1000,n_factors)
for (i in 1:1000)
{
  null_simulations_matrix.hvn[i,] <- null_simulations.hvn[[i]]$pvals
}
null_simulations_matrix.hvn <- as.data.frame(t(null_simulations_matrix.hvn))
ddf.hvn <- cbind(obs.hvn,null_simulations_matrix.hvn)

ddf.hvn %>%
  tidyr::gather(group, pvalue, c(1,3:1002)) %>%
  mutate(color = ifelse(group =='Pvals', "observed", 'simulated')) %>%
  ggplot() +
  geom_line(aes(x = factor, y= log(pvalue), group = group, color = color)) +
  ggtitle('henipavirus virus simulations') + 
  scale_x_continuous(limits=c(1, 10), breaks = seq(1,10, by= 1))

#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#...................................filo viruses......................................

Z.filo <- batphy1$filo_surv

pf.filo <- phylofactor::twoSampleFactor(Z.filo, method='Fisher', bat_tree, n_factors ,ncores = 1)

#nms.filo <- matrix(0,n_factors+1,2)
nms.filo <- list()
for (i in 2:(n_factors+1))
{
  indexes = pf.filo$bins[[i]]
  species_x <- gsub("_", " ", tolower(bat_tree$tip.label[indexes]))
  group_taxonomy_list <- as.data.frame(taxonomy[match(species_x,taxonomy[,1]),2])
  #nms.filo[i-1, 1] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
  nms.filo[i-1] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
  
  print(i)
  
  p<-batphy1 %>%
  filter(tolower(species) %in% species_x) %>%
  group_by(filo_surv) %>%
  summarise(n = n()) 
  
  print(as.numeric(p[2,'n'] / (p[2,'n'] + p[1,'n'])))
}

B.filo <- bins(pf.filo$basis[,1:n_factors])
B.filo <- B.filo[2:(n_factors+1)] 
nms1.filo <- nms.filo[1:(n_factors)]

#nms1.hvn <- nms.hnv[2:(n_factors+1)]
## remove paraphyletic bin

probs.filo <- sapply(B.filo,FUN=function(ix,Z.filo) mean(Z.filo[ix]),Z=Z.filo) %>% signif(.,digits=2)

#probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

names(nms1.filo) <- probs.filo

#pf.tree.filo$legend$colors[1] <- pf.tree.filo$legend$colors[7]
pf.tree.filo <- pf.tree(pf.filo, lwd=0.0025, branch.length = "none") 

pf.tree.filo$ggplot
pf.tree.filo$ggplot +
  #ggtree::geom_treescale(linesize=.02, col='green') +  
  ggtree::geom_tippoint(size=5*Z.filo,col='blue')  

Legend.filo <- pf.tree.filo$legend
Legend.filo$names <- nms1.filo
P.filo <- sapply(probs.filo,FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend.filo$names <- mapply(paste,Legend.filo$names,P.filo)
plot.new()
legend('topleft',legend=Legend.filo$names,fill=Legend.filo$colors,cex=1)


#legend('topleft',legend=Legend[1,'names'],fill=Legend[1,'colors'],cex=1.8)

#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#...............Null Simulation..........................................................

#null_simulations <- pf.nullsim(pf.filo, reps = 50, nfactors = n_factors, output ='pvals')
#null_simulations.filo <- readRDS('data/my_phylofactor_object_nullsim')


#nfactors <- pf$nfactors
obs <- data.frame('Pvals'=pf.filo$pvals,
                  'factor'=1:n_factors
                  )
null_simulations_matrix<- matrix(0,50,n_factors)
for (i in 1:50)
{
  null_simulations_matrix[i,] <- null_simulations.filo[[i]]$pvals
}
null_simulations_matrix <- as.data.frame(t(null_simulations_matrix))
ddf.filo <- cbind(obs,null_simulations_matrix)

ddf.filo %>%
  tidyr::gather(group, pvalue, c(1,3:22)) %>%
  mutate(color = ifelse(group =='Pvals', "observed", 'simulated')) %>%
  ggplot() +
  geom_line(aes(x = factor, y= log(pvalue), group = group, color = color)) +
  ggtitle('filo virus simulations') +
  scale_x_continuous(limits=c(0, 10), breaks = seq(1,10, by= 1))
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................

#Yinpterochiroptera is the group that is not Yangochiroptera 
