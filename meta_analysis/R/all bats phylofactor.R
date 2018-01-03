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

#..............obtain a couple more outcome variables..................................
#................load bat phylogeny data...................................

load('data/bat_taxonomy_data.Rdata') 
load(file='data/seroprevalence.Rdata')

rm(batphy)

seroprevalence_filo_binary <- seroprevalence %>%
  #mutate(seroprevalence_percentage_cat = ifelse(seroprevalence_percentage == 0, 0, 1)) %>%
  filter(!(is.na(species))) %>%
  #dplyr::group_by(species, virus, seroprevalence_percentage_cat) %>%
  #summarise(count = n()) %>%
  dplyr::group_by(species, virus) %>%
  summarise(filovirus_positive_cat=mean(seroprevalence_percentage)) %>% #keep the same name as before just to make things easier
  ungroup() %>%
  filter(virus == 'Filovirus') %>%
  select(-virus)#%>%
  #spread(seroprevalence_percentage_cat, count) %>%
  #mutate(filovirus_positive_cat  = ifelse(is.na(`1`) & !is.na(`0`), 0, ifelse(!is.na(`1`), 1, NA))) %>%
  #dplyr::select(c(species, filovirus_positive_cat)) %>%
  #unique()

batphy1 <- left_join(batphy1, seroprevalence_filo_binary) 

# seroprevalence_hnv_binary <- seroprevalence %>%
#   mutate(seroprevalence_percentage_cat = ifelse(seroprevalence_percentage == 0, 0, 1)) %>%
#   filter(!(is.na(species))) %>%
#   dplyr::group_by(species, virus, seroprevalence_percentage_cat) %>%
#   summarise(count = n()) %>%
#   ungroup() %>%
#   filter(virus == 'Henipavirus') %>%
#   spread(seroprevalence_percentage_cat, count) %>%
#   mutate(henipavirus_positive_cat  = ifelse(is.na(`1`) & !is.na(`0`), 0,
#                                             ifelse(!is.na(`1`), 1, NA))) %>%
#   dplyr::select(c(species, henipavirus_positive_cat)) %>%
#   unique()

seroprevalence_hnv_binary <- seroprevalence %>%
  #mutate(seroprevalence_percentage_cat = ifelse(seroprevalence_percentage == 0, 0, 1)) %>%
  filter(!(is.na(species))) %>%
  #dplyr::group_by(species, virus, seroprevalence_percentage_cat) %>%
  #summarise(count = n()) %>%
  dplyr::group_by(species, virus) %>%
  summarise(henipavirus_positive_cat=mean(seroprevalence_percentage)) %>%
  ungroup() %>%
  filter(virus == 'Henipavirus') %>%
  select(-virus)
  #%>%
#spread(seroprevalence_percentage_cat, count) %>%
#mutate(filovirus_positive_cat  = ifelse(is.na(`1`) & !is.na(`0`), 0, ifelse(!is.na(`1`), 1, NA))) %>%
#dplyr::select(c(species, filovirus_positive_cat)) %>%
#unique()

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

n_factors.hnv = 7

Z.hnv <- batphy1$hnv_surv

pf.hnv <- phylofactor::twoSampleFactor(Z.hnv, bat_tree, method='Fisher', n_factors.hnv ,ncores = 1)

nms.hnv <- list()
for (i in 2:(n_factors.hnv+1))
{
  indexes = pf.hnv$bins[[i]]
  species <- gsub("_", " ", tolower(bat_tree$tip.label[indexes]))
  group_taxonomy_list <- as.data.frame(taxonomy[match(species,taxonomy[,1]),2])
  nms.hnv[i-1] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
  print(i)
}

B.hnv <- bins(pf.hnv$basis[,1:n_factors.hnv])
B.hnv <- B.hnv[2:(n_factors.hnv+1)] 
nms1.hnv <- nms.hnv[1:(n_factors.hnv)]

## remove paraphyletic bin

probs.hnv <- sapply(B.hnv,FUN=function(ix,Z.hnv) mean(Z.hnv[ix]),Z.hnv=Z.hnv) %>% signif(.,digits=2)
names(nms1.hnv) <- probs.hnv

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
#............................9:  Pteropodidae........"#CC00FFFF".................................
#............................10: Scotophilus kuhlii.."#FF0099FF"..................................
#............................11: Eidolon............."#FF0099FF"..................................
#............................12: Hipposideros........"#00FF66FF"..................................
#............................13: C. perspicillata...."#00FFFFFF"..................................
#............................14: D. rotundus........."#00FFFFFF"..................................
#............................14: Carollia............"#00FFFFFF"..................................

colfcn.hnv <- function(n) return(c("#FF0000FF", "#00FFFFFF","#3300FFFF","#33FF00FF","#00FFFFFF","#FF0099FF","#00FF66FF"))

pf.tree.hnv <- pf.tree(pf.hnv, lwd=1, color.fcn=colfcn.hnv, branch.length = "none")
pf.tree.hnv$ggplot +
  ggtree::geom_tippoint(size=10*Z.hnv,col='blue')  

Legend.hnv <- pf.tree.hnv$legend
Legend.hnv$names <- nms1.hnv
P.hnv <- sapply(probs.hnv,FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend.hnv$names <- mapply(paste,Legend.hnv$names,P.hnv)
plot.new()
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
# 
# null_simulations.hnv <- readRDS('data/my_phylofactor_object_hvn_nullsim')
# 
# #nfactors <- pf$nfactors
# obs.hnv <- data.frame('Pvals'=pf.hnv$pvals,
#                   'factor'=1:n_factors
# )
# null_simulations_matrix.hnv<- matrix(0,1000,n_factors)
# for (i in 1:1000)
# {
#   null_simulations_matrix.hnv[i,] <- null_simulations.hnv[[i]]$pvals
# }
# null_simulations_matrix.hnv <- as.data.frame(t(null_simulations_matrix.hnv))
# ddf.hnv <- cbind(obs.hnv,null_simulations_matrix.hnv)
# 
# ddf.hnv %>%
#   tidyr::gather(group, pvalue, c(1,3:1002)) %>%
#   mutate(color = ifelse(group =='Pvals', "observed", 'simulated')) %>%
#   ggplot() +
#   geom_line(aes(x = factor, y= log(pvalue), group = group, color = color)) +
#   ggtitle('henipavirus virus simulations') + 
#   scale_x_continuous(limits=c(1, 10), breaks = seq(1,10, by= 1))

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
#............................9:  Pteropodidae........#CC00FFFF".................................
#............................10: Scotophilus kuhlii..#FF0099FF..................................

colfcn.filo <- function(n) return(c("#FF0000FF", "#FF9900FF" ,"#CCFF00FF" ,"#33FF00FF", "#FF00E6FF", "#00FFFFFF" ,"#0066FFFF", "#3300FFFF")) #,"#CC00FFFF" ,"#FF0099FF"))

pf.tree.filo <- pf.tree(pf.filo, lwd=1, color.fcn = colfcn.filo, factors = 1:8, branch.length = "none")

pf.tree.filo$ggplot +
  ggtree::geom_tippoint(size=10*Z.filo,col='blue')  

Legend.filo <- pf.tree.filo$legend
Legend.filo$names <- nms1.filo[1:8]
P.filo <- sapply(probs.filo,FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend.filo$names <- mapply(paste,Legend.filo$names,P.filo[1:8])
plot.new()
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
#.................................................................................................
#.................................................................................................
#.................................................................................................
#.................................................................................................

#Yinpterochiroptera is the group that is not Yangochiroptera 


#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................
#.............add in filovirus and henipavirus POSITIVITY......................................

#setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")

#.............henipavirus positivity....................

n_factors.hnv.pos = 10

#Z.hnv_pos.hnv <- batphy1$henipavirus_positive_cat
#pf.hnv_pos <- phylofactor::twoSampleFactor(Z.hnv_pos, bat_tree, method='Fisher', n_factors.hnv.pos ,ncores = 1)

Z.hnv_pos.hnv <- batphy1[,c('tree_species', 'hnv_samps', 'henipavirus_positive_cat')] 

Z.hnv_pos.hnv1  <- Z.hnv_pos.hnv %>%
  rename(Species = tree_species) %>%
  mutate(Species = as.character(Species)) %>%
  mutate(Sample =1) %>%
  filter(!(is.na(henipavirus_positive_cat)))
         
dropped.tips.hnv.drop.tips <- Z.hnv_pos.hnv %>%
  rename(Species = tree_species) %>%
  filter((is.na(henipavirus_positive_cat))) %>%
  mutate(Species = as.character(Species)) %>%
  select(Species) %>%
  as.vector()

bat_tree.hnv.pos <- ape::drop.tip(bat_tree,dropped.tips.hnv.drop.tips$Species)

pf.hnv_pos <- gpf(Z.hnv_pos.hnv1,bat_tree.hnv.pos,frmla=henipavirus_positive_cat~hnv_samps+phylo,nfactors=10,mStableAgg=F)

pf.tree(pf.hnv_pos)

pf.hnv_pos.pval.storage <- list()
nms.hnv.pos <- list()
for (i in 2:(n_factors.hnv.pos+1))
{
  indexes = pf.hnv_pos$bins[[i]]
  species_x <- gsub("_", " ", tolower(bat_tree.hnv.pos$tip.label[indexes]))
  group_taxonomy_list <- as.data.frame(taxonomy[match(species_x,taxonomy[,1]),2])
  nms.hnv.pos[i-1] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
  #print(i)
  print(gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list))))
  p<-batphy1 %>%
   filter(tolower(species) %in% species_x) %>%
   group_by(henipavirus_positive_cat) %>%
   summarise(n = n())
  
  pf.hnv_pos.pval.storage[i-1] <-  (as.numeric(p[2,'n'] / (p[2,'n'] + p[1,'n'])))
}

B.hnv.pos <- bins(pf.hnv_pos$basis[,1:n_factors.hnv.pos])
B.hnv.pos <- B.hnv.pos[2:(n_factors.hnv.pos+1)] 
nms1.hnv.pos <- nms.hnv.pos[1:(n_factors.hnv.pos)]

## remove paraphyletic bin

probs.hnv.pos <- sapply(B.hnv.pos,FUN=function(ix,Z.hnv_pos) mean(Z.hnv_pos[ix]),Z.hnv_pos=Z.hnv_pos) %>% signif(.,digits=2)
names(nms1.hnv.pos) <- probs.hnv.pos

#..........................NULL SIMULATIONS for + HENIPAVIRUSES.................................................
#..........................NULL SIMULATIONS for + HENIPAVIRUSES.................................................
#..........................NULL SIMULATIONS for + HENIPAVIRUSES.................................................
#..........................NULL SIMULATIONS for + HENIPAVIRUSES.................................................
#..........................NULL SIMULATIONS for + HENIPAVIRUSES.................................................
#..........................NULL SIMULATIONS for + HENIPAVIRUSES.................................................
#..........................NULL SIMULATIONS for + HENIPAVIRUSES.................................................
#..........................NULL SIMULATIONS for + HENIPAVIRUSES.................................................
#..........................NULL SIMULATIONS for + HENIPAVIRUSES.................................................
#..........................NULL SIMULATIONS for + HENIPAVIRUSES.................................................
#..........................NULL SIMULATIONS for + HENIPAVIRUSES.................................................
#..........................NULL SIMULATIONS for + HENIPAVIRUSES.................................................
#..........................NULL SIMULATIONS for + HENIPAVIRUSES.................................................
#..........................NULL SIMULATIONS for + HENIPAVIRUSES.................................................
#..........................NULL SIMULATIONS for + HENIPAVIRUSES.................................................

#only 3 bins are of interest, the rest are species single groups so for the null simulations lets set it to 3
nfactors.hnv.random = 3 

randomData <- function()
{
  rbinom(nrow(Z.hnv_pos.hnv1), size = 1, prob = sum(Z.hnv_pos.hnv1$henipavirus_positive_cat)/nrow(Z.hnv_pos.hnv1))
}
NoReplicates <- 100
p_storage <- matrix(0, NoReplicates, nfactors.hnv.random)

for (replicate in 1:NoReplicates){
  Z <- randomData()  ## define how you want to make this random data
  DF <- Z.hnv_pos.hnv1 %>%
    mutate(henipavirus_positive_cat = Z) %>%
    rename(henipavirus_positive_cat_random = henipavirus_positive_cat)
  pf <- gpf(DF,bat_tree.hnv.pos,frmla=henipavirus_positive_cat_random~hnv_samps+phylo,nfactors=nfactors.hnv.random,mStableAgg=F)
  
  for (i in 2:(nfactors.hnv.random+1))
  {
    indexes = pf$bins[[i]]
    p<-DF[c(indexes),] 
    p <-sum(p$henipavirus_positive_cat_random) / nrow(p)
    
    p_storage[replicate, i-1] <- p
  }
}

pf.hnv_pos.pval.storage <- matrix(1,3,1) #in this case all the group pvals were 1

obs <- data.frame('factor'=1:nfactors.hnv.random,
                  'Pvals'=pf.hnv_pos.pval.storage) %>% t()
                  
ddf.hnv.null.sim <- as.data.frame(t(rbind(obs, p_storage)))
  
ddf.hnv.null.sim %>%
  tidyr::gather(group, pvalue, c(2:100)) %>%
  mutate(color = ifelse(group =='Pvals', "observed", 'simulated')) %>%
  ggplot() +
  geom_line(aes(x = factor, y= log(pvalue), group = group, color = color))


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
#............................9:  Pteropodidae........"#CC00FFFF".................................
#............................10: Scotophilus kuhlii.."#FF0099FF"..................................
#............................11: Eidolon............."#FF0099FF"..................................
#............................12: Hipposideros........"#00FF66FF"..................................

colfcn.hnv.pos <- function(n) return(c("#33FF00FF"))

pf.tree.hnv.pos <- pf.tree(pf.hnv_pos, lwd=1, factors =1, color.fcn = colfcn.hnv.pos, branch.length = "none")
pf.tree.hnv.pos$ggplot +
  ggtree::geom_tippoint(size=10*Z.hnv_pos.hnv1$henipavirus_positive_cat,col='blue')  +
  ggtree::geom_tippoint(size=2*log(Z.hnv_pos.hnv1$hnv_samps),col='red')  
  
Legend.hnv.pos <- pf.tree.hnv.pos$legend
Legend.hnv.pos$names <- nms1.hnv.pos[1]
P.hnv.pos <- sapply(probs.hnv.pos,FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend.hnv.pos$names <- mapply(paste,Legend.hnv.pos$names,P.hnv.pos[1])
plot.new()
plot.new()
legend('topleft',legend=Legend.hnv.pos$names,fill=Legend.hnv.pos$colors,cex=1)

#.........................................add in filovirus positivity...............................
#.........................................add in filovirus positivity...............................
#.........................................add in filovirus positivity...............................
#.........................................add in filovirus positivity...............................
#.........................................add in filovirus positivity...............................
#.........................................add in filovirus positivity...............................
#.........................................add in filovirus positivity...............................
#.........................................add in filovirus positivity...............................
#.........................................add in filovirus positivity...............................
#.........................................add in filovirus positivity...............................
#.........................................add in filovirus positivity...............................
#.........................................add in filovirus positivity...............................
#.........................................add in filovirus positivity...............................
#.........................................add in filovirus positivity...............................
#.........................................add in filovirus positivity...............................
#.........................................add in filovirus positivity...............................

n_factors.filo.pos = 10

Z.joined.filo.pos <- batphy1[,c('tree_species', 'filo_samps', 'filovirus_positive_cat')] 

Z.joined.filo.pos.1  <- Z.joined.filo.pos %>%
  rename(Species = tree_species) %>%
  mutate(Species = as.character(Species)) %>%
  mutate(Sample =1) %>%
  filter(!(is.na(filovirus_positive_cat)))

Z.joined.filo.pos.drop.tips<- Z.joined.filo.pos %>%
  rename(Species = tree_species) %>%
  filter((is.na(filovirus_positive_cat))) %>%
  mutate(Species = as.character(Species)) %>%
  select(Species) %>%
  as.vector()

bat_tree.filo.pos<- ape::drop.tip(bat_tree,Z.joined.filo.pos.drop.tips$Species)

pf.filo_pos <- gpf(Z.joined.filo.pos.1,bat_tree.filo.pos,frmla=filovirus_positive_cat~filo_samps+phylo,nfactors=10,mStableAgg=F)

pf.tree(pf.hnv_pos)

#Z.filo_pos <- batphy1$filovirus_positive_cat
#pf.filo_pos <- phylofactor::twoSampleFactor(Z.filo_pos, bat_tree, method='Fisher', n_factors.filo.pos ,ncores = 1)

nms.filo.pos <- list()
pf.filo_pos.pval.storage <- list()

for (i in 2:(n_factors.filo.pos+1))
{
  indexes = pf.filo_pos$bins[[i]]
  species_x <- gsub("_", " ", tolower(bat_tree.filo.pos$tip.label[indexes]))
  group_taxonomy_list <- as.data.frame(taxonomy[match(species_x,taxonomy[,1]),2])
  nms.filo.pos[i-1] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
  print(i)
  
  print(gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list))))
  p<-batphy1 %>%
    filter(tolower(species) %in% species_x) %>%
    group_by(filovirus_positive_cat) %>%
    summarise(n = n()) 
  
  pf.hnv_pos.pval.storage[i] <-  (as.numeric(p[2,'n'] / (p[2,'n'] + p[1,'n'])))
  print(p)
}

B.filo.pos <- bins(pf.filo_pos$basis[,1:n_factors.filo.pos])
B.filo.pos <- B.filo.pos[2:(n_factors.filo.pos+1)] 
nms1.filo.pos <- nms.filo.pos[1:(n_factors.filo.pos)]

## remove paraphyletic bin

probs.filo.pos <- sapply(B.filo.pos,FUN=function(ix,Z.filo_pos) mean(Z.filo_pos[ix]),Z.filo_pos=Z.filo_pos) %>% signif(.,digits=2)
names(nms1.filo.pos) <- probs.filo.pos

#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................
#.............................NULL SIMULATIONS + FILOVIRUSES...........................................


randomData <- function()
{
  rbinom(nrow(Z.filo.pos1), size = 1, prob = sum(Z.filo.pos1$filovirus_positive_cat)/nrow(Z.filo.pos1))
}
NoReplicates <- 100
p_storage <- matrix(0, NoReplicates, 10)

for (replicate in 1:NoReplicates){
  Z <- randomData()  ## define how you want to make this random data
  DF <- Z.filo.pos1 %>%
    mutate(filovirus_positive_cat = Z) %>%
    rename(filovirus_positive_cat_random = filovirus_positive_cat)
  pf <- gpf(DF,bat_tree.filo.pos,frmla=filovirus_positive_cat_random~filo_samps+phylo,nfactors=10,mStableAgg=F)
  
  for (i in 2:11)
  {
    indexes = pf$bins[[i]]
    p<-DF[c(indexes),] 
    p <-sum(p$filovirus_positive_cat_random) / nrow(p)
    
    p_storage[replicate, i-1] <- p
  }
}

pf.filo_pos.pval.storage <- matrix(1,10,1) #in this case all the group pvals were 1

obs <- data.frame('factor'=1:n_factors.filo.pos,
                  'Pvals'=pf.filo_pos.pval.storage) %>% t()

ddf.filo.null.sim <- as.data.frame(t(rbind(obs, p_storage)))

ddf.filo.null.sim %>%
  tidyr::gather(group, pvalue, c(2:100)) %>%
  mutate(color = ifelse(group =='Pvals', "observed", 'simulated')) %>%
  ggplot() +
  geom_line(aes(x = factor, y= log(pvalue), group = group, color = color))


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
#............................9:  Pteropodidae........"#CC00FFFF".................................
#............................10: Scotophilus kuhlii.."#FF0099FF"..................................
#............................11: Eidolon............."#FF0099FF"..................................
#............................12: Hipposideros........"#00FF66FF"..................................
#............................13: C. perspicillata...."#00FFFFFF"..................................
#............................14: D. rotundus........."#0066FFFF"..................................

colfcn.filo.pos <- function(n) return(c("#33FF00FF"))

pf.tree.filo.pos <- pf.tree(pf.filo_pos, lwd=1, factors =1, color.fcn = colfcn.filo.pos, branch.length = "none")
pf.tree.filo.pos$ggplot +
  ggtree::geom_tippoint(size=10*Z.joined.filo.pos.1$filovirus_positive_cat,col='blue')  +
  ggtree::geom_tippoint(size=log(Z.joined.filo.pos.1$filo_samps),col='red')  
  
Legend.filo.pos <- pf.tree.filo.pos$legend
Legend.filo.pos$names <- nms1.filo.pos[1]
P.filo.pos <- sapply(probs.filo.pos,FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend.filo.pos$names <- mapply(paste,Legend.filo.pos$names,P.filo.pos[1])
plot.new()
plot.new()
legend('topleft',legend=Legend.filo.pos$names,fill=Legend.filo.pos$colors,cex=1)
