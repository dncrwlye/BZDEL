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


Data <- batphy1 %>%
  filter(unique_name %in% bat_spp_India$MSW05_Binomial) %>%
  filter(hnv_samps > 0) %>%
  mutate(log_hnv_samps = log(hnv_samps)) %>%
  select(c(hnv_samps, hnv_positive, species, log_hnv_samps)) %>%
  rename(Species = species) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  mutate(Species = stri_trans_totitle(Species)) %>%
  mutate(Sample = 1)

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

ncores = 4
tot.reps=500
reps.per.worker=round(tot.reps/ncores)
sampling_effort = TRUE

source('R/filo and hnv positive factorization/null simulations script.R')
save(list=ls(),file='data/phylofactor work spaces/indian bats only niv_workspace')


pf <- gpf(Data,tree,frmla.phylo=Z~phylo,nfactors=10,family=binomial,algorithm='phylo')
pf$factors
ncores = 4
tot.reps=500
reps.per.worker=round(tot.reps/ncores)
sampling_effort = FALSE
source('R/filo and hnv positive factorization/null simulations script.R')
save(list=ls(),file='data/phylofactor work spaces/indian bats no sampling effort only niv_workspace')

rm(list=ls())

#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................
#....................signifigant factors .............................................

load(file='data/phylofactor work spaces/indian bats no sampling effort only niv_workspace')

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

load(file='data/phylofactor work spaces/indian bats only niv_workspace')

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


load("data/bat_taxonomy_data.Rdata")
source('R/taxonomy group name function.R')

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

names.storage
factor.1 <- pf$groups[[1]][[1]]
paraphyletic.remainder <-     pf$groups[[1]][[2]]

Z <- Data$Z
factor.1.p <- sapply(factor.1,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
paraphyletic.remainder.p <- sapply(paraphyletic.remainder,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

length(factor.1)
length(paraphyletic.remainder.p)

factor.1.prob <- sum(factor.1.p)/length(factor.1.p)
paraphyletic.remainder.prob <- sum(paraphyletic.remainder.p)/length(paraphyletic.remainder.p)

# 0.547619 for Pteropodidae

# B <- lapply(pf$groups, function(l) l[[1]])
# # B <- bins(pf$basis[,1:10])
# # B <- B[2:11]
# Z <- Data$Z
# probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
# 


names.storage[[1]]

#factors from image 1 

# factor group    colors
# 1       Pteropodidae          1 #440154FF
# 2       Vespertilionoidea     1 #482878FF
# 3       Pipistrellini         1 #3E4A89FF
# 4       Myotis                1 #31688EFF
# 5       Pteropus              1 #26828EFF
# 6       Rhinolophus           1 #1F9E89FF
# 7       Pteropus 2????        1 #35B779FF
# 8       Phyllostomini         1 #6DCD59FF
# 9       Molossus              1 #B4DE2CFF
# 10      Carollia              1 #FDE725FF
# 11      Molossinae            1 #B4DE2CFF
# 12      Miniopterus           1 #FDE725FF
# 13      Vespertilionini       1 #FF9900FF
# 14      Macroglossus          1 #33FF00FF

colfcn <- function(n) return(c("#440154FF"))

pf.tree <- pf.tree(pf, lwd=1, factors = 1, color.fcn=colfcn, branch.length = "none", bg.color = NA)

d <- data.frame(x=pf.tree$ggplot$data[1:143,'x'] + .5,
                xend=pf.tree$ggplot$data[1:143,'x'] + 1+ Data$log_effort,
                y=pf.tree$ggplot$data[1:143,'y'],
                yend=pf.tree$ggplot$data[1:143,'y'] )

Data <- Data %>%
  mutate(ib = ifelse(Species %in% text2$V3, 1, 0))

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue')  +
  #ggtree::geom_tippoint(size=5*Data$ib,col='red') +
  geom_segment(data= d,aes(x=x,y=y,xend=xend,yend=yend, size= Data$log_effort, colour = 'blue'))

ggsave("figures/hnv sampling effort tree w indian bats.png", bg = "transparent", height = 18, width = 18)

#something strange is going on, we have a factor 2 that is only one tip, but 
#in the bins it contains many many species. This is odd....
Legend <- pf.tree$legend
Legend$names <- names.storage[1]
P <- sapply(probs[1],FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend$names <- mapply(paste,Legend$names,P)
Cairo(file='figures/hnv sampling effort tree legend.jpg', 
      type="png",
      units="in", 
      width=5, 
      height=4, 
      pointsize=12, 
      dpi=200)
plot.new()
plot.new()
legend('topleft', legend=Legend$names,fill=Legend$colors,cex=1, bg="transparent", bty="n")
dev.off()







#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
#######bunch of more bat bs................
load('data/seroprevalence.Rdata')

data.set.for.barbara <- seroprevalence %>%
  filter(virus == 'Henipavirus') %>%
  mutate(seroprevalence_percentage_cat = ifelse(seroprevalence_percentage > 0, 1, 0)) %>%
  group_by(species, virus, virusold, seroprevalence_percentage_cat, methodology) %>%
  summarise(n=n()) %>%
  filter(seroprevalence_percentage_cat == 1)

data.set.for.barbara.2 <- seroprevalence %>%
  filter(virus == 'Henipavirus') %>%
  mutate(seroprevalence_percentage_cat = ifelse(seroprevalence_percentage > 0, 1, 0)) %>%
  group_by(species, virus, virusold, seroprevalence_percentage_cat, methodology) %>%
  summarise(n=n()) %>%
  filter(seroprevalence_percentage_cat == 0) %>%
  filter(!(species %in% data.set.for.barbara$species))

data.set.for.barbara <- rbind(data.set.for.barbara, data.set.for.barbara.2)

rm(data.set.for.barbara.2)

hnv_pos_for_barabara <- read_csv("~/Desktop/hnv.pos.for.barabara.csv")

x <- data.set.for.barbara %>%
  filter(!(species %in% hnv_pos_for_barabara$species)) %>%
  filter(spe)

data.set.for.barbara <- data.set.for.barbara %>%
  filter(species %in% hnv_pos_for_barabara$species)

write_csv(data.set.for.barbara, path = "~/Desktop/data.set.for.barbara.update.csv")

x<- seroprevalence %>%
  filter(virus == 'Henipavirus') %>%
  filter(species == 'rousettus leschenaultii')

x <- seroprevalence %>%
  filter(grepl('diadema',species))

bat_spp_India <- read_csv("~/Desktop/BDEL/BZDEL/meta_analysis/data/phylofactor work spaces/bat spp India.csv")
load("data/bat_taxonomy_data.Rdata")
taxonomy <- batphy1 %>%
  select(c(species, tax)) 

unique(seroprevalence$sampling_location)

india.outbreaks <- seroprevalence %>%
  filter(grepl('india', sampling_location)) %>%
  filter(virus == 'Henipavirus')

india.outbreaks %>%
  select(species) %>%
  unique()

india.outbreaks %>%
  select(sampling_location, sampling_location_two) %>%
  unique()

bat_spp_India <- bat_spp_India %>%
  mutate(species.tnsfrm = tolower(MSW05_Binomial)) %>%
  mutate(species.tnsfrm = gsub("_", " ", species.tnsfrm)) 

bat_spp_India.join <- left_join(bat_spp_India, seroprevalence, by = c("species.tnsfrm" = "species")) 
bat_spp_India.join <- left_join(bat_spp_India.join, taxonomy, by =c('species.tnsfrm' = 'species'))

bat_spp_India.join <- bat_spp_India.join %>%
  filter(virus != 'Filovirus'| is.na(virus))

bat_spp_India.join <- bat_spp_India.join %>%
  mutate(Pteropodidae = ifelse(grepl('Pteropodidae', tax,), TRUE, FALSE))

x <- bat_spp_India.join %>%
  filter(Pteropodidae == TRUE) %>%
  #filter(seroprevalence_percentage >0) %>%
  select(species.tnsfrm, methodology, seroprevalence_percentage) %>%
  unique()



bat_spp_India.join %>%
  select(sampling_location, sampling_location_two, sampling_location_three) %>%
  unique()

bat_spp_India.join.pcr <- bat_spp_India.join %>%
  filter(methodology == "PCR based method")

bat_spp_India.join <- bat_spp_India.join %>%
  mutate(Vespertilionoidea = ifelse(grepl('Vespertilionoidea', tax,), TRUE, FALSE))

bat_spp_India.join %>%
  filter(Vespertilionoidea == TRUE) %>%
  filter(is.na(virus)) %>%
  select(species.tnsfrm) %>%
  unique()

#...........evenmore!!!!.......
#...........evenmore!!!!.......
#...........evenmore!!!!.......
#...........evenmore!!!!.......
#...........evenmore!!!!.......
#...........evenmore!!!!.......
#...........evenmore!!!!.......
#...........evenmore!!!!.......

library(tidyverse)
load('data/seroprevalence.Rdata')
batphy_for_rotl_update <- read_csv("~/Desktop/BDEL/BZDEL/meta_analysis/data/batphy_for_rotl update.csv")

nipah.question <- seroprevalence %>%
  filter(virusold == 'Nipah'|virusold == "Henipavirus") 

nipah.question.missing <- nipah.question%>%
  filter(!(species %in% batphy_for_rotl_update$species))

unique(nipah.question.missing$species)

nipah.question$species=plyr::revalue(nipah.question$species,c("rhinolophus refulgens"="rhinolophus lepidus",
                                                              "hipposideros cf caffer"="hipposideros caffer",
                                                              "hipposideros cf caffer/ruber"="hipposideros caffer",
                                                              "hipposideros cf ruber"="hipposideros ruber",
                                                              "natalus lanatus"="natalus mexicanus",
                                                              "myonycteris leptodon"="myonycteris torquata",
                                                              "nanonycteris veldkampiii"="nanonycteris veldkampii",
                                                              "tadarida plicata"="chaerephon plicatus"))

nipah.question$species=plyr::revalue(nipah.question$species,c("Lissonycteris angolensis"="Myonycteris angolensis",
                                                              "Nyctalus plancyi"="Nyctalus velutinus",
                                                              "Rhinolophus stheno"="Rhinolophus microglobosus",
                                                              "Natalus stramineus"="Natalus stramineus mexicanus",
                                                              "Neoromicia somalicus"="Neoromicia malagasyensis"))

nipah.question.missing <- nipah.question%>%
  filter(!(species %in% batphy_for_rotl_update$species))

unique(nipah.question.missing$species)

nipah.question <- nipah.question %>%
  mutate(seroprevalence_prevalence_percentage.cat  = ifelse(seroprevalence_percentage >0, 1, 0)) 

nipah.question <- left_join(nipah.question, batphy_for_rotl_update[,c('species', 'region')])

nipah.question.table <-  nipah.question %>%       
  select(seroprevalence_prevalence_percentage.cat, species, region) %>%
  group_by(species, region) %>%
  summarise(nipah.pos.binary=sum(seroprevalence_prevalence_percentage.cat, na.rm=TRUE))  %>%
  mutate(nipah.pos.binary = ifelse(nipah.pos.binary > 0, 1, 0))


revised_HNV_positive_bats_uncleaned <- read_csv("~/Dropbox_gmail/Dropbox/bat virus meta-analysis/revised HNV positive bats_uncleaned.csv")
revised_NiV_positive_bats_uncleaned <- read_csv("~/Desktop/revised NiV positive bats_uncleaned.csv")

revised_NiV_positive_bats_uncleaned <-revised_NiV_positive_bats_uncleaned %>%
  rename(niv.pos = hnvpos)

x <- full_join(revised_HNV_positive_bats_uncleaned, revised_NiV_positive_bats_uncleaned)

y <- left_join(x, batphy_for_rotl_update[,c('species', 'region')])

y <- y %>%
  rename(wilson.reeder.region = region)



write.csv(y, file = "/Users/buckcrowley/Desktop/hnv.pos.for.barabara.csv")


















