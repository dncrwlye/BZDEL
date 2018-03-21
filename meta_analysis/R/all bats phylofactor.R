setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")

library(phylofactor)
library(grid)
library(gridBase)
library(stringi)
#source('R/other phylofactor scripts/fisherFactor.R')
source('R/taxonomy group name function.R')
source('R/bat filo and hnv sampling phylogeny update.R')
#source('R/other phylofactor scripts/visualization_fcns.R')
library(tidyverse)
library(taxize)

#..................obtain a couple more outcome variables..................................
load('data/bat_taxonomy_data.Rdata') 

load(file='data/seroprevalence.Rdata')
seroprevalence_filo_binary <- seroprevalence %>%
  mutate(seroprevalence_percentage_cat = ifelse(seroprevalence_percentage == 0, 0, 1)) %>%
  filter(!(is.na(species))) %>%
  dplyr::group_by(species, virus, seroprevalence_percentage_cat) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  filter(virus == 'Filovirus') %>%
  spread(seroprevalence_percentage_cat, count) %>%
  mutate(filovirus_positive_cat  = ifelse(is.na(`1`) & !is.na(`0`), 0,
                                          ifelse(!is.na(`1`), 1, NA))) %>%
  dplyr::select(c(species, filovirus_positive_cat)) %>%
  unique()

batphy1 <- left_join(batphy1, seroprevalence_filo_binary) 

seroprevalence_hnv_binary <- seroprevalence %>%
  mutate(seroprevalence_percentage_cat = ifelse(seroprevalence_percentage == 0, 0, 1)) %>%
  filter(!(is.na(species))) %>%
  dplyr::group_by(species, virus, seroprevalence_percentage_cat) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  filter(virus == 'Henipavirus') %>%
  spread(seroprevalence_percentage_cat, count) %>%
  mutate(henipavirus_positive_cat  = ifelse(is.na(`1`) & !is.na(`0`), 0,
                                            ifelse(!is.na(`1`), 1, NA))) %>%
  dplyr::select(c(species, henipavirus_positive_cat)) %>%
  unique()

batphy1 <- left_join(batphy1, seroprevalence_hnv_binary) 

batphy1 <- batphy1 %>%
  mutate(hnv_surv = ifelse(is.na(hnv_surv), 0, hnv_surv)) %>%
  mutate(filo_samps = ifelse(is.na(filo_samps), 0, filo_samps)) %>%
  mutate(hnv_samps =  ifelse(is.na(hnv_samps), 0, hnv_samps)) 

taxonomy<- batphy1[,c('species','tax')]

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

#........................henipaviruses...............................

n_factors.hnv = 10

Z.hnv <- batphy1$hnv_surv

pf.hnv <- phylofactor::twoSampleFactor(Z.hnv, bat_tree, method='Fisher', n_factors.hnv ,ncores = 1)

pf.hnv.pval.storage <- list()
nms.hnv <- list()
for (i in 2:(n_factors.hnv+1))
{
  indexes = pf.hnv$bins[[i]]
  species_x <- gsub("_", " ", tolower(bat_tree$tip.label[indexes]))
  group_taxonomy_list <- as.data.frame(taxonomy[match(species_x,taxonomy[,1]),2])
  nms.hnv[i-1] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
  print(i)
  
  p<-batphy1 %>%
    filter(tolower(species) %in% species_x) %>%
    group_by(hnv_surv) %>%
    summarise(n = n())
  
  pf.hnv.pval.storage[i-1] <-  (as.numeric(p[2,'n'] / (p[2,'n'] + p[1,'n'])))
}

B.hnv <- bins(pf.hnv$basis[,1:n_factors.hnv])
B.hnv <- B.hnv[2:(n_factors.hnv+1)] 
nms1.hnv <- nms.hnv[1:(n_factors.hnv)]

## remove paraphyletic bin

probs.hnv <- sapply(B.hnv,FUN=function(ix,Z.hnv) mean(Z.hnv[ix]),Z.hnv=Z.hnv) %>% signif(.,digits=2)
names(nms1.hnv) <- probs.hnv

#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#........................................................................................
#...............Null Simulation..........................................................
# 
null_simulations.hnv <- readRDS('data/my_phylofactor_object_hvn_nullsim')

#nfactors <- pf$nfactors
obs.hnv <- data.frame('Pvals'=pf.hnv$pvals,
                      'factor'=1:n_factors.hnv
)
null_simulations_matrix.hnv<- matrix(0,1000,n_factors.hnv)
for (i in 1:1000)
{
  null_simulations_matrix.hnv[i,] <- null_simulations.hnv[[i]]$pvals
}
null_simulations_matrix.hnv <- as.data.frame(t(null_simulations_matrix.hnv))
ddf.hnv <- cbind(obs.hnv,null_simulations_matrix.hnv)

ddf.hnv %>%
  tidyr::gather(group, pvalue, c(1,3:1002)) %>%
  mutate(color = ifelse(group =='Pvals', "observed", 'simulated')) %>%
  ggplot() +
  geom_line(aes(x = factor, y= log(pvalue), group = group, color = color)) +
  ggtitle('henipavirus virus simulations') +
  scale_x_continuous(limits=c(1, 10), breaks = seq(1,10, by= 1))

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

colfcn.hnv <- function(n) return(c("#FF0000FF", "#FF9900FF","#3300FFFF","#33FF00FF","#00FFFFFF","#FF0099FF","#00FF66FF"))

pf.tree.hnv <- pf.tree(pf.hnv, lwd=1, factors = 1:7, color.fcn=colfcn.hnv, branch.length = "none")
pf.tree.hnv$ggplot +
  ggtree::geom_tippoint(size=10*Z.hnv,col='blue')  

Legend.hnv <- pf.tree.hnv$legend
Legend.hnv$names <- nms1.hnv[1:7]
P.hnv <- sapply(probs.hnv[1:7],FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend.hnv$names <- mapply(paste,Legend.hnv$names,P.hnv)
plot.new()
plot.new()
legend('topleft',legend=Legend.hnv$names,fill=Legend.hnv$colors,cex=1)


#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................
#...................................filo viruses......................................

Z.filo <- batphy1$filo_surv
n_factors.filo = 10
pf.filo <- phylofactor::twoSampleFactor(Z.filo, method='Fisher', bat_tree, n_factors.filo ,ncores = 1)

#nms.filo <- matrix(0,n_factors+1,2)
nms.filo <- list()
for (i in 2:(n_factors.filo+1))
{
  indexes = pf.filo$bins[[i]]
  species_x <- gsub("_", " ", tolower(bat_tree$tip.label[indexes]))
  group_taxonomy_list <- as.data.frame(taxonomy[match(species_x,taxonomy[,1]),2])
  nms.filo[i-1] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
  
  print(i)
  
  p<-batphy1 %>%
    filter(tolower(species) %in% species_x) %>%
    group_by(filo_surv) %>%
    summarise(n = n()) 
  
  print(as.numeric(p[2,'n'] / (p[2,'n'] + p[1,'n'])))
}

B.filo <- bins(pf.filo$basis[,1:n_factors.filo])
B.filo <- B.filo[2:(n_factors.filo+1)] 
nms1.filo <- nms.filo[1:(n_factors.filo)]

#nms1.hnv <- nms.hnv[2:(n_factors+1)]
## remove paraphyletic bin

probs.filo <- sapply(B.filo,FUN=function(ix,Z.filo) mean(Z.filo[ix]),Z=Z.filo) %>% signif(.,digits=2)
names(nms1.filo) <- probs.filo

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
null_simulations.filo <- readRDS('data/my_phylofactor_object_filo_nullsim')

#nfactors <- pf$nfactors
obs <- data.frame('Pvals'=pf.filo$pvals,
                  'factor'=1:n_factors.filo
)
null_simulations_matrix<- matrix(0,1000,n_factors.filo)
for (i in 1:1000)
{
  null_simulations_matrix[i,] <- null_simulations.filo[[i]]$pvals
}
null_simulations_matrix <- as.data.frame(t(null_simulations_matrix))
ddf.filo <- cbind(obs,null_simulations_matrix)

ddf.filo %>%
  tidyr::gather(group, pvalue, c(1,3:1002)) %>%
  mutate(color = ifelse(group =='Pvals', "observed", 'simulated')) %>%
  ggplot() +
  geom_line(aes(x = factor, y= log(pvalue), group = group, color = color)) +
  ggtitle('filo virus simulations') +
  scale_x_continuous(limits=c(1, 10), breaks = seq(1,10, by= 1))

#hard to make any sense out of it, i suppse choose 7 or 8, 7 to be conservative?
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................


#............................define our colors..................................................
#...............................................................................................
#............................1:  Yangochiroptera.....#FF0000FF..................................
#............................2:  Myotis..............#FF9900FF..................................
#............................3:  Miniopterus.........#CCFF00FF..................................
#............................4:  Pteropodidae........#33FF00FF..................................
#............................5:  Emballonuroidea.....#FF004DFF..................................
#............................6:  Rhinolophus.........#00FFFFFF..................................
#............................7:  Pipistrellini.......#0066FFFF..................................
#............................8:  Pteropus............#3300FFFF..................................
#............................9:  Pteropodinae........#CC00FFFF".................................
#............................10: Scotophilus kuhlii..#FF0099FF..................................

colfcn.filo <- function(n) return(c("#FF0000FF", "#FF9900FF" ,"#CCFF00FF" ,"#33FF00FF", "#FF00E6FF", "#00FFFFFF" ,"#0066FFFF")) #, "#3300FFFF")) ,"#CC00FFFF" ,"#FF0099FF"))

pf.tree.filo <- pf.tree(pf.filo, lwd=1, color.fcn = colfcn.filo, factors = 1:7, branch.length = "none")

pf.tree.filo$ggplot +
  ggtree::geom_tippoint(size=10*Z.filo,col='blue')  

Legend.filo <- pf.tree.filo$legend
Legend.filo$names <- nms1.filo[1:7]
P.filo <- sapply(probs.filo,FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend.filo$names <- mapply(paste,Legend.filo$names,P.filo[1:7])
plot.new()
plot.new()
legend('topleft',legend=Legend.filo$names,fill=Legend.filo$colors,cex=1)

#legend('topleft',legend=Legend[1,'names'],fill=Legend[1,'colors'],cex=1.8)



#Yinpterochiroptera is the group that is not Yangochiroptera 
\