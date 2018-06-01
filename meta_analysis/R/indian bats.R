library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
par(family="Times")   
CairoFonts(regular = 'Times-12')

setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis")

bat_spp_India <- read_csv("C:/Users/r83c996/Desktop/bat spp India.csv")
load('data/phylofactor work spaces/bat_tree')
load("data/bat_taxonomy_data.Rdata")

bat_spp_India %>%
  filter(!(MSW05_Binomial %in% batphy1$unique_name))

bat_spp_India$MSW05_Binomial =plyr::revalue(bat_spp_India$MSW05_Binomial,c("Myotis_blythii"="Myotis_oxygnathus"))

bat_spp_India %>%
  filter(!(MSW05_Binomial %in% batphy1$unique_name))

#only Eptesicus_bottae doesn't appear in the database, I can't find another name for this 
#species that is the database either. Also, according to the internets it doesn't actually 
#live in in India so maybe thats not a problem?

revised_NiV_positive_bats_uncleaned <- read_csv("~/BZDEL/meta_analysis/data/phylofactor work spaces/revised NiV positive bats_uncleaned.csv")

revised_NiV_positive_bats_uncleaned %>%
  filter(!(species %in% batphy1$species)) 

revised_NiV_positive_bats_uncleaned$species =
  plyr::revalue(revised_NiV_positive_bats_uncleaned$species,c("ballionycterus maculata"="balionycteris maculata",
                                                            "scotophilus heathi" = "scotophilus heathii"
                                                              ))

revised_NiV_positive_bats_uncleaned %>%
  filter(!(species %in% batphy1$species)) 

#i'm okay with this for now

batphy1 <- left_join(batphy1, revised_NiV_positive_bats_uncleaned)
nrow(revised_NiV_positive_bats_uncleaned) ==  nrow(filter(batphy1, !is.na(hnvpos))) + nrow(revised_NiV_positive_bats_uncleaned %>% filter(!(species %in% batphy1$species)) )

Data <- batphy1 %>%
  rename(niv_pos = hnvpos) %>%
  select(-X1) %>% 
  mutate(log_hnv_samps = log(hnv_samps)) %>%
  select(c(hnv_samps, niv_pos, species, log_hnv_samps)) %>%
  rename(Species = species) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  mutate(Species = stri_trans_totitle(Species)) %>%
  mutate(Sample = 1)

Data <- Data %>%
  filter(!is.na(niv_pos))

tree <- ape::drop.tip(bat_tree,bat_tree$tip.label[!(bat_tree$tip.label %in% Data$Species)])

rm(batphy1, bat_tree)

names(Data) <- c('effort','Z', 'Species','log_effort', 'Sample')

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

save(list=ls(),file='data/phylofactor work spaces/indian bats only niv_workspace')

pf <- gpf(Data,tree,frmla.phylo=Z~phylo,nfactors=10,family=binomial,algorithm='phylo')
pf$factors
# 
ncores = 4
tot.reps=500
reps.per.worker=round(tot.reps/ncores)
sampling_effort = FALSE
source('R/filo and hnv positive factorization/null simulations script.R')

save(list=ls(),file='data/phylofactor work spaces/indian bats no sampling effort only niv_workspace')


#............phylofactor visualization script...................................

# Libraries ---------------------------------------------------------------

setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
library(Cairo)
par(family="Times")   
CairoFonts(regular = 'Times-12')

# hnv all bats factors-----------------------------------------------

rm(list=ls())

load(file='data/phylofactor work spaces/indian bats no sampling effort only niv_workspace')

plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,25))
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

rm(list=ls())

load(file='data/phylofactor work spaces/indian bats only niv_workspace')

plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,25))
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


pf.tree <- pf.tree(pf, lwd=1, branch.length = "none", bg.color = NA)

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue') 
 

taxonomy <- batphy1 %>%
  select(c(species, tax)) 

names.storage <- list()

for (i in 1:10)
{
  indexes = pf$groups[[i]][[1]]
  species <- gsub("_", " ", tolower(tree$tip.label[indexes]))
  group_taxonomy_list <- as.data.frame(taxonomy[match(species,taxonomy[,1]),2])
  names.storage[i] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
  print(i)
}









