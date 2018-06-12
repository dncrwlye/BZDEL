#setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis")
setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")

library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
load('data/phylofactor work spaces/bat_tree')
load("data/bat_taxonomy_data.Rdata")
species_NiV_data_stratified_by_detection <- read_csv("C:/Users/r83c996/Dropbox/Nipah paper PLOS Outbreaks/Resources/species NiV data_stratified by detection.csv")
ncores = 4
tot.reps=200
reps.per.worker=round(tot.reps/ncores)

x<-species_NiV_data_stratified_by_detection[,c('species',"positive_all")]

x$species =plyr::revalue(x$species,c("scotophilus heathi"="scotophilus heathii"))
                                                                           
batphy1 <- left_join(batphy1, x)

xx <- batphy1 %>% 
  filter(!is.na(positive_all))

jj <- x %>%
  filter(!(species %in% batphy1$species))



# Data <- batphy1 %>%
#   filter(!is.na(positive_all)) %>%
#   mutate(log_hnv_samps = log(hnv_samps)) %>%
#   select(c(hnv_samps, positive_all, species, log_hnv_samps))%>%
#   rename(Species = species) %>%
#   mutate(Species = gsub(" ", "_", Species)) %>%
#   mutate(Species = stri_trans_totitle(Species)) %>%
#   mutate(Sample = 1) %>%
#   unique()

Data <- batphy1 %>%
  mutate(log_hnv_samps = log(hnv_samps)) %>%
  filter(virus =='Henipavirus') %>%
  filter(grepl('thailand|malaysia|china|indonesia|bangladesh|papua new guinea|vietnam|india|timor', sampling_location)) %>%
  select(species) %>%
  unique() %>%
  select(c(hnv_samps, positive_all, species, log_hnv_samps))%>%
  rename(Species = species) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  mutate(Species = stri_trans_totitle(Species)) %>%
  mutate(Sample = 1) %>%
  unique()

Data[duplicated(Data$Species),]

tree<-bat_tree
tree <- ape::drop.tip(tree,tree$tip.label[!(tree$tip.label %in% Data$Species)])
Data <- Data[Data$Species %in% tree$tip.label,]

rm(batphy1, bat_tree)

names(Data) <- c('effort', 'Z', 'Species', 'log_effort', 'Sample')


model.sampling.effort <- glm(Z~log_effort,family=binomial,data=Data)

Data <- Data %>%
  mutate(effort.fit = coef(model.sampling.effort)['log_effort']*Data$log_effort)

rm(model.sampling.effort)

pf <- gpf(Data,tree,frmla=Z~offset(effort.fit),
          frmla.phylo=Z~offset(effort.fit)+phylo,
          family=binomial,
          nfactors=10,
          algorithm='phylo')
pf$factors
# 
ncores = 4
tot.reps=500
reps.per.worker=round(tot.reps/ncores)
sampling_effort = TRUE
source('R/filo and hnv positive factorization/null simulations script.R')



plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,15))
for (rr in 1:497){
  lines(c(S[rr,1],diff(S[rr,])),col=rgb(0,0,0,0.2))
}
points(deviances,col='red',pch=16,cex=2)

bb <- (t(diff(t(S))))
bb <- cbind(S[,1],bb)
Ojb.bb <- c(Obj[1],diff(Obj))
for (i in 1:10)
{
  print(ecdf(bb[,i])(Ojb.bb[i]))
}
names.storage
pf.tree <- pf.tree(pf, lwd=1, factors = 8, branch.length = "none", bg.color = NA)

pf.tree$ggplot+
ggtree::geom_tippoint(size=10*Data$Z,col='blue')
  

save(list=ls(),file='data/phylofactor work spaces/niv_only_workspace')
load('data/phylofactor work spaces/niv_only_workspace')
1-ecdf(X)(anova(pf$models[[1]])['phylo','Deviance'])



pf$call
x <- pf$models[[1]] %>% summary()
1/exp(x$coefficients[2,1])

exp(-2.157)

taxonomy <- batphy1 %>%
  select(c(species, tax)) 

names.storage <- list()

for (i in 1:10)
{
  indexes = pf$groups[[i]][[1]]
  species <- gsub("_", " ", tolower(tree$tip.label[indexes]))
  #print(species)
  #print(i)
  group_taxonomy_list <- as.data.frame(taxonomy[match(species,taxonomy[,1]),2])
  names.storage[i] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
}

factor.1 <- pf$groups[[1]][[1]]
paraphyletic.remainder <-     pf$groups[[1]][[2]]


