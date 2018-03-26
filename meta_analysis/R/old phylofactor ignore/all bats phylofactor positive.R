setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")

library(phylofactor)
library(grid)
library(gridBase)
library(stringi)
#source('R/other phylofactor scripts/fisherFactor.R')
source('R/taxonomy group name function.R')
source('R/bat filo and hnv sampling phylogeny.R')
#source('R/other phylofactor scripts/visualization_fcns.R')
library(tidyverse)
library(taxize)

#..................obtain a couple more outcome variables..................................
load('data/bat_taxonomy_data.Rdata') 

load(file='data/seroprevalence.Rdata')
load(file='data/seroprevalence.phylo.analysis.Rdata')

seroprevalence_filo_binary <- seroprevalence.phylo.analysis %>%
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

seroprevalence_hnv_binary <- seroprevalence.phylo.analysis %>%
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

#.............add in filovirus and henipavirus positivity....................
#.............add in filovirus and henipavirus positivity....................
#.............add in filovirus and henipavirus positivity....................
#.............add in filovirus and henipavirus positivity....................
#.............add in filovirus and henipavirus positivity....................
#.............add in filovirus and henipavirus positivity....................
#.............add in filovirus and henipavirus positivity....................
#.............add in filovirus and henipavirus positivity....................
#.............add in filovirus and henipavirus positivity....................
#.............add in filovirus and henipavirus positivity....................
#.............add in filovirus and henipavirus positivity....................
#.............add in filovirus and henipavirus positivity....................
#.............add in filovirus and henipavirus positivity....................
#.............add in filovirus and henipavirus positivity....................

#setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")

#.............henipavirus positivity....................

n_factors.hnv.pos = 10

#Z.hnv_pos <- batphy1$henipavirus_positive_cat
#pf.hnv_pos <- phylofactor::twoSampleFactor(Z.hnv_pos, bat_tree, method='Fisher', n_factors.hnv.pos ,ncores = 1)

Z.hnv_pos <- batphy1[,c('tree_species', 'hnv_samps', 'henipavirus_positive_cat')] 
taxonomy<- batphy1[,c('species','tax')]


Z.hnv_pos.1  <- Z.hnv_pos %>%
  mutate(Species = as.character(tree_species)) %>%
  select(-c(tree_species)) %>%
  mutate(Sample =1) %>%
  mutate(log_hnv_samps = log(hnv_samps)) %>%
  filter(!(is.na(henipavirus_positive_cat)))

dropped.tips.hnv.drop.tips <- Z.hnv_pos %>%
  mutate(Species = as.character(tree_species)) %>%
  select(-c(tree_species)) %>%
  filter((is.na(henipavirus_positive_cat))) %>%
  mutate(Species = as.character(Species)) %>%
  select(Species) %>%
  as.vector()

hist(Z.hnv_pos.1$hnv_samps)

bat_tree.hnv.pos <- ape::drop.tip(bat_tree,dropped.tips.hnv.drop.tips$Species)
#############################################HNV POSITIVITY NO SAMPLING EFFORT######################################
#############################################HNV POSITIVITY NO SAMPLING EFFORT######################################
#############################################HNV POSITIVITY NO SAMPLING EFFORT######################################
#############################################HNV POSITIVITY NO SAMPLING EFFORT######################################
#############################################HNV POSITIVITY NO SAMPLING EFFORT######################################
#############################################HNV POSITIVITY NO SAMPLING EFFORT######################################
#############################################HNV POSITIVITY NO SAMPLING EFFORT######################################
#############################################HNV POSITIVITY NO SAMPLING EFFORT######################################
#############################################HNV POSITIVITY NO SAMPLING EFFORT######################################
#############################################HNV POSITIVITY NO SAMPLING EFFORT######################################
#############################################HNV POSITIVITY NO SAMPLING EFFORT######################################
#############################################HNV POSITIVITY NO SAMPLING EFFORT######################################


pf.hnv_pos.no.effort <- gpf(Z.hnv_pos.1,bat_tree.hnv.pos,frmla.phylo =henipavirus_positive_cat~phylo,nfactors=10,algorithm='phylo',family=binomial)

#save(Z.hnv_pos.1, file = '/Users/buckcrowley/Desktop/Z.hnv_pos.1')
#save(bat_tree.hnv.pos, file = '/Users/buckcrowley/Desktop/bat_tree.hnv.pos')

pf.hnv_pos.pval.storage.no.effort <- list()
nms.hnv.pos.no.effort <- list()
for (i in 2:(n_factors.hnv.pos+1))
{
  indexes = pf.hnv_pos.no.effort$bins[[i]]
  species_x <- gsub("_", " ", tolower(bat_tree.hnv.pos$tip.label[indexes]))
  group_taxonomy_list <- as.data.frame(taxonomy[match(species_x,taxonomy[,1]),2])
  nms.hnv.pos.no.effort[i-1] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
  #print(i)
  print(gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list))))
}

#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################
#############################################HNV POSITIVITY WITH SAMPLING EFFORT######################################

pf.hnv_pos <- gpf(Z.hnv_pos.1,bat_tree.hnv.pos,frmla.phylo =henipavirus_positive_cat~log_hnv_samps+phylo,nfactors=10,algorithm='phylo',family=binomial)

#save(Z.hnv_pos.1, file = '/Users/buckcrowley/Desktop/Z.hnv_pos.1')
#save(bat_tree.hnv.pos, file = '/Users/buckcrowley/Desktop/bat_tree.hnv.pos')

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
    filter(!(tolower(species) %in% species_x))
  
  p<-Z.hnv_pos.1[c(-indexes),] 
  fisher.matrix <- matrix(c(sum(p$henipavirus_positive_cat == 1),
                            sum(p$henipavirus_positive_cat == 0),  
                            sum(Z.hnv_pos.1$henipavirus_positive_cat == 1),
                            sum(Z.hnv_pos.1$henipavirus_positive_cat == 0)), nrow=2)
  p.value <-fisher.test(fisher.matrix)$p.value 
  
  pf.hnv_pos.pval.storage[i-1] <- p.value
  
}

# #an alternative way of getting p values...i guess??
# pf.hnv_pos.pval.storage.alt <- list()
# pf.hnv_pos.pval.storage.alt.aov <- list()
# for (i in 1:10)
# {
#   # print(i)
#   pf.hnv_pos.pval.storage.alt[i] <- pchisq(pf.hnv_pos$models[[i]]$null.deviance- pf.hnv_pos$models[[i]]$deviance, df=(pf.hnv_pos$models[[1]]$df.null-pf.hnv_pos$models[[1]]$df.residual))
#   x<-summary(aov(pf.hnv_pos$models[[i]])) 
#   pf.hnv_pos.pval.storage.alt.aov[i] <- x[[1]]$`Pr(>F)`[2]
# }
# 
# #okay now that is over
# 
# B.hnv.pos <- bins(pf.hnv_pos$basis[,1:n_factors.hnv.pos])
# B.hnv.pos <- B.hnv.pos[2:(n_factors.hnv.pos+1)] 
# nms1.hnv.pos <- nms.hnv.pos[1:(n_factors.hnv.pos)]

## remove paraphyletic bin

# probs.hnv.pos <- sapply(B.hnv.pos,FUN=function(ix,Z.hnv_pos) mean(Z.hnv_pos[ix]),Z.hnv_pos=Z.hnv_pos) %>% signif(.,digits=2)
# names(nms1.hnv.pos) <- probs.hnv.pos

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


#null simulation are run in another script now. the results are....
#no factors are significant for Henipaviruses after accounting for sampling effort :0 !!!


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

pf.tree.hnv.pos <- pf.tree(pf.hnv_pos, lwd=1, branch.length = "none") #again, its easier to make the tree and just assign a factor that is really hard to see, since none are signigant, than actually make a tree wtih no factors
pf.tree.hnv.pos$ggplot +
  ggtree::geom_tippoint(size=10*Z.hnv_pos.1$henipavirus_positive_cat,col='blue')   #+
#ggtree::geom_tippoint(size=3*log(Z.hnv_pos.1$hnv_samps),col='red')  

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

Z.filo.pos <- batphy1[,c('tree_species', 'filo_samps', 'filovirus_positive_cat')] 

Z.filo.pos.1  <- Z.filo.pos %>%
  mutate(Species = as.character(tree_species)) %>%
  select(-c(tree_species)) %>%
  #dplyr::rename(Species = tree_species) %>%
  mutate(Sample =1) %>%
  mutate(log_filo_samps = log(filo_samps)) %>%
  filter(!(is.na(filovirus_positive_cat)))

Z.filo.pos.drop.tips<- Z.filo.pos %>%
  mutate(Species = tree_species) %>%
  select(-c(tree_species)) %>%
  filter((is.na(filovirus_positive_cat))) %>%
  mutate(Species = as.character(Species)) %>%
  select(Species) %>%
  as.vector()

bat_tree.filo.pos<- ape::drop.tip(bat_tree,Z.filo.pos.drop.tips$Species)

save(Z.filo.pos.1, file = '/Users/buckcrowley/Desktop/Z.filo.pos.1')
save(bat_tree.filo.pos, file = '/Users/buckcrowley/Desktop/bat_tree.filo.pos')

pf.filo_pos <- gpf(Z.filo.pos.1,bat_tree.filo.pos,frmla.phylo =filovirus_positive_cat~log_filo_samps+phylo,nfactors=10, algorithm='phylo',family=binomial)

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
    filter(!(tolower(species) %in% species_x))
  
  p<-Z.filo.pos.1[c(-indexes),] 
  fisher.matrix <- matrix(c(sum(p$filovirus_positive_cat == 1),
                            sum(p$filovirus_positive_cat == 0),  
                            sum(Z.filo.pos.1$filovirus_positive_cat == 1),
                            sum(Z.filo.pos.1$filovirus_positive_cat == 0)), nrow=2)
  p.value <-fisher.test(fisher.matrix)$p.value 
  
  pf.filo_pos.pval.storage[i-1] <- p.value
}

#an alternative way of getting p values...i guess??
pf.filo_pos.pval.storage.alt <- list()
pf.filo_pos.pval.storage.alt.aov <- list()
for (i in 1:10)
{
  x<-summary(aov(pf.filo_pos$models[[i]])) 
  pf.filo_pos.pval.storage.alt.aov[i] <- x[[1]]$`Pr(>F)`[2]
  pf.filo_pos.pval.storage.alt[i] <- pchisq(pf.filo_pos$models[[i]]$null.deviance- pf.filo_pos$models[[i]]$deviance, df=(pf.filo_pos$models[[1]]$df.null-pf.filo_pos$models[[1]]$df.residual))
}

B.filo.pos <- bins(pf.filo_pos$basis[,1:n_factors.filo.pos])
B.filo.pos <- B.filo.pos[2:(n_factors.filo.pos+1)] 
nms1.filo.pos <- nms.filo.pos[1:(n_factors.filo.pos)]

## remove paraphyletic bin

# probs.filo.pos <- sapply(B.filo.pos,FUN=function(ix,Z.filo_pos) mean(Z.filo_pos[ix]),Z.filo_pos=Z.filo_pos) %>% signif(.,digits=2)
# names(nms1.filo.pos) <- probs.filo.pos

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

#null simulations now run in another script, only first factor is signifigant 

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

colfcn.filo.pos <- function(n) return(c("#CC00FFFF"))

pf.tree.filo.pos <- pf.tree(pf.filo_pos, lwd=1, factors =1, color.fcn = colfcn.filo.pos, branch.length = "none")


pf.tree.filo.pos$ggplot +
  ggtree::geom_tippoint(size=10*Z.filo.pos.1$filovirus_positive_cat,col='blue')  +
  ggtree::geom_tippoint(size=3*log(Z.filo.pos.1$filo_samps),col='red')  


Legend.filo.pos <- pf.tree.filo.pos$legend
Legend.filo.pos$names <- nms1.filo.pos[1]
P.filo.pos <- sapply(probs.filo.pos,FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend.filo.pos$names <- mapply(paste,Legend.filo.pos$names,P.filo.pos[1])
plot.new()
plot.new()
legend('topleft',legend=Legend.filo.pos$names,fill=Legend.filo.pos$colors,cex=1)
