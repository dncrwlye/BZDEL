library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
library(ggtree)
par(family="Times")   

#setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis")
setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
#bat_spp_India <- read_csv("C:/Users/r83c996/Desktop/bat spp India.csv")
bat_spp_India <- read_csv('data/phylofactor work spaces/bat spp India.csv')
bat_spp_Kerala <- read_csv('data/phylofactor work spaces/Bats_Kerala.csv')

load('data/phylofactor work spaces/bat_tree')
load("data/bat_taxonomy_data.Rdata")

bat_spp_India %>%
  filter(!(MSW05_Binomial %in% batphy1$unique_name))

bat_spp_India$MSW05_Binomial =plyr::revalue(bat_spp_India$MSW05_Binomial,c("Myotis_blythii"="Myotis_oxygnathus",
                                                                           "Myotis_peytoni"="Myotis_montivagus"))

bat_spp_India %>%
  filter(!(MSW05_Binomial %in% batphy1$unique_name))

#only Eptesicus_bottae doesn't appear in the database, I can't find another name for this 
#species that is the database either. Also, according to the internets it doesn't actually 
#live in in India so maybe thats not a problem? funky face

Data <- batphy1 %>%
  filter(unique_name %in% bat_spp_India$MSW05_Binomial) %>% #only use bats that are in India, obvi gonna loose a couple 
  mutate(hnv_samps = ifelse(is.na(hnv_samps), 0, hnv_samps)) %>% #remove NAs
  mutate(hnv_positive = ifelse(is.na(hnv_positive), 0, hnv_positive)) %>% #remove NAs
  mutate(log_hnv_samps = log(hnv_samps)) %>% #log transform sampling effort
  select(c(hnv_samps, hnv_positive, species, log_hnv_samps)) %>%
  dplyr::rename(Species = species) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  mutate(Species = stri_trans_totitle(Species)) %>%
  mutate(Sample = 1)

# remove tips from tree that are in our Data file
Data$Species %in% bat_tree$tip.label
tree <- ape::drop.tip(bat_tree,bat_tree$tip.label[!(bat_tree$tip.label %in% Data$Species)])
rm(batphy1, bat_tree)

names(Data) <- c('effort','Z', 'Species','log_effort', 'Sample')

# ncores = 4
# tot.reps=500
# reps.per.worker=round(tot.reps/ncores)
# sampling_effort = TRUE
# 
# source('R/filo and hnv positive factorization/null simulations script.R')
# save(list=ls(),file='data/phylofactor work spaces/indian bats only niv_workspace')

# pf <- gpf(Data,tree,frmla.phylo=Z~phylo,nfactors=10,family=binomial,algorithm='phylo')
# pf$factors
# ncores = 4
# tot.reps=500
# reps.per.worker=round(tot.reps/ncores)
# sampling_effort = FALSE
# source('R/filo and hnv positive factorization/null simulations script.R')
# save(list=ls(),file='data/phylofactor work spaces/indian bats no sampling effort only niv_workspace')

#rm(list=ls())

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

#load(file='data/phylofactor work spaces/indian bats no sampling effort only niv_workspace')
bat_spp_India <- bat_spp_India %>%
  mutate(MSW05_Binomial = gsub(" ", "_", MSW05_Binomial))

bat_spp_Kerala <- bat_spp_Kerala %>%
  mutate(MSW05_binomial = gsub(" ", "_", MSW05_binomial))

Data <- Data %>%
  mutate(IB = ifelse(Species %in% bat_spp_India$MSW05_Binomial, 1, 
                     ifelse(!(Species %in% bat_spp_India$MSW05_Binomial), 0, NA))) %>%
  mutate(Kerala = ifelse(Species %in% bat_spp_Kerala$MSW05_binomial, 1,
                         ifelse(!(Species %in% bat_spp_Kerala$MSW05_binomial), 0, NA)))

load('data/seroprevalence.Rdata')

# pcr.pos.hnv <- seroprevalence %>%
#   filter(methodology == 'PCR based method') %>%
#   filter(seroprevalence_percentage > 0) %>%
#   filter(virusold == 'Nipah') %>%
#   select(species) %>%
#   unique()
# 
# neg.hnv <- seroprevalence %>%
#   filter(seroprevalence_percentage == 0) %>%
#   filter(virus == 'Henipavirus') %>%
#   select(species) %>%
#   unique()

Data <- Data %>%
  mutate(species.mutate = tolower(gsub("_", " ", Species))) %>%
  #mutate(all.neg = ifelse(species.mutate %in% neg.hnv$species & effort > 0, 1, 0))%>%
  mutate(all.neg = ifelse(effort > 0 & Z == 0, 1, 0)) %>%
  mutate(effort_cat = ifelse(effort > 0, 1,0))

Data <- Data %>%
  mutate(color = ifelse(Kerala == 1, "black", "grey"))

source('r/stupid_pointless_internal_node_color_function.R')
library('phylobase')
internal.edges<- (internal.edge.color(tree, Data))

row.names(internal.edges) <- internal.edges[,1]
#internal.edges <- internal.edges[,-c(1)] %>% as.data.frame()
g1 = as(tree, 'phylo4')
length(nodeId(g1, "internal")) == length(row.names(internal.edges))

g2 = phylo4d(g1, Data)
nodeData(g2) <- internal.edges

load("data/bat_taxonomy_data.Rdata")

taxonomy <- batphy1 %>%
  select(c(species, tax)) 
source('R/taxonomy group name function.R')

#...............................................................................


jj <- nrow(Data)

d <- data.frame(x=ggtree.object$data[1:jj,'x'],
                xend=ggtree.object$data[1:jj,'x']  + 1 + Data$log_effort,
                y=ggtree.object$data[1:jj,'y'],
                yend=ggtree.object$data[1:jj,'y'] )

d

d[is.finite(d$xend)==FALSE,] <- 0

#g2 <- ggtree::ggtree(tree,layout='circular') 

ggtree.object <- ggtree::ggtree(g2, layout= 'circular',branch.length = "none",aes(color=I(color))) 

ggtree.object$data$label <- gsub("_"," ", ggtree.object$data$label)

i <- .5
ii <- 5
iii <- 2

ggtree.object$data$Z[ggtree.object$data$Z==0]  <- NA
ggtree.object$data$effort_cat[ggtree.object$data$effort_cat==0]  <- NA

ggtree.object + 
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_cladelabel(node=118, label="Vespertilionidae", 
                          color="black", angle=12.16-180, offset=ii, offset.text = iii) +
  ggtree::geom_cladelabel(node=200, label="Pteropodidae", 
                          color="black", angle=-68.91, offset=ii, offset.text = iii) +
  ggtree::geom_tippoint(aes(size=effort_cat), pch =1) +
  ggtree::geom_tippoint(aes(size=Z), color = 'blue', alpha = .8) +
  geom_segment(data= d,aes(x=x + i ,y=y,xend=xend +i ,yend=yend, size= Data$effort_cat,
                           colour = 'orange', alpha = (.5 + Data$Kerala))) +
  ggtree::geom_tiplab(offset=i, aes(angle=angle)) 

ggsave("figures/indian bats only info sampling effort no tip lab circular.png", bg = "transparent", height = 18, width = 18)

Data %>% filter(effort_cat == 1)

#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
#...................Alex's binomial check....................................................
# 
# rm(list=ls())
# bat_spp_India <- read_csv("C:/Users/r83c996/Desktop/bat spp India.csv")
# bat_spp_India$MSW05_Binomial =plyr::revalue(bat_spp_India$MSW05_Binomial,c("Myotis_blythii"="Myotis_oxygnathus"))
# 
# load("data/bat_taxonomy_data.Rdata")
# source('R/taxonomy group name function.R')
# 
# taxonomy <- batphy1 %>%
#   select(c(unique_name, species, tax)) 
# 
# bat_spp_India <- left_join(bat_spp_India, taxonomy, by = c("MSW05_Binomial" = "unique_name"))
# 
# 
# filter(bat_spp_India, grepl("Pteropodidae", tax)) %>% nrow()
# 
# #12 out of 112
# 
# rm(list=ls())
# 
# load(file='data/phylofactor work spaces/hnv_workspace_no_sampling_effort')
# load("data/bat_taxonomy_data.Rdata")
# 
# taxonomy <- batphy1 %>%
#   select(c(unique_name, species, tax)) 
# 
# Data <- left_join(Data, taxonomy, by = c("Species" = "unique_name"))
# 
# filter(Data, grepl("Pteropodidae", tax)) %>% nrow()
# 
# #41 out of 143
# 
# x <- matrix(c(12, 112-12, 41, 143-41), nrow=2, byrow=TRUE)
# x
# 
# fisher.test(x)
# chisq.test(x) 
# 
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# #..........figure........................................................................
# 
# 
# rm(list=ls())
# 
# load(file='data/phylofactor work spaces/hnv_workspace_no_sampling_effort')
# load("data/bat_taxonomy_data.Rdata")
# source('R/taxonomy group name function.R')
# 
# 
# load('data/seroprevalence.Rdata')
# 
# pcr.pos.hnv <- seroprevalence %>%
#   filter(methodology == 'PCR based method') %>%
#   filter(seroprevalence_percentage > 0) %>%
#   filter(virus == 'Henipavirus')
# 
# Data <- Data %>%
#   mutate(species.mutate = tolower(gsub("_", " ", Species))) %>%
#   mutate(pcr.pos = ifelse(species.mutate %in% pcr.pos.hnv$species, 1,0))
# 
# colfcn <- function(n) return(c("#33FF00FF", "#440154FF", "#FF9900FF"))
# pf.tree <- pf.tree(pf, lwd=1, factors = 1:3, color.fcn=colfcn, branch.length = "none", bg.color = NA)
# 
# pf.tree$ggplot +
#   ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
#   ggtree::geom_tippoint(size=10*Data$Z,col='blue') +
#   ggtree::geom_tippoint(size=5*Data$pcr.pos,col='red') 
# 
# 
# ggtree::geom_cladelabel(node=179, label="Noctilionoidea", 
#                         color="black", angle=222, offset=2, offset.text = 2) 
# 
# 
# 
# #######bunch of more bat bs................
# #######bunch of more bat bs................
# #######bunch of more bat bs................
# #######bunch of more bat bs................
# #######bunch of more bat bs................
# #######bunch of more bat bs................
# #######bunch of more bat bs................
# #######bunch of more bat bs................
# #######bunch of more bat bs................
# #######bunch of more bat bs................
# #######bunch of more bat bs................
# load('data/seroprevalence.Rdata')
# 
# data.set.for.barbara <- seroprevalence %>%
#   filter(virus == 'Henipavirus') %>%
#   mutate(seroprevalence_percentage_cat = ifelse(seroprevalence_percentage > 0, 1, 0)) %>%
#   group_by(species, virus, virusold, seroprevalence_percentage_cat, methodology) %>%
#   summarise(n=n()) %>%
#   filter(seroprevalence_percentage_cat == 1)
# 
# data.set.for.barbara.2 <- seroprevalence %>%
#   filter(virus == 'Henipavirus') %>%
#   mutate(seroprevalence_percentage_cat = ifelse(seroprevalence_percentage > 0, 1, 0)) %>%
#   group_by(species, virus, virusold, seroprevalence_percentage_cat, methodology) %>%
#   summarise(n=n()) %>%
#   filter(seroprevalence_percentage_cat == 0) %>%
#   filter(!(species %in% data.set.for.barbara$species))
# 
# data.set.for.barbara <- rbind(data.set.for.barbara, data.set.for.barbara.2)
# 
# rm(data.set.for.barbara.2)
# 
# hnv_pos_for_barabara <- read_csv("~/Desktop/hnv.pos.for.barabara.csv")
# 
# x <- data.set.for.barbara %>%
#   filter(!(species %in% hnv_pos_for_barabara$species)) %>%
#   filter(spe)
# 
# data.set.for.barbara <- data.set.for.barbara %>%
#   filter(species %in% hnv_pos_for_barabara$species)
# 
# write_csv(data.set.for.barbara, path = "~/Desktop/data.set.for.barbara.update.csv")
# 
# x<- seroprevalence %>%
#   filter(virus == 'Henipavirus') %>%
#   filter(species == 'rousettus leschenaultii')
# 
# x <- seroprevalence %>%
#   filter(grepl('diadema',species))
# 
# bat_spp_India <- read_csv("~/Desktop/BDEL/BZDEL/meta_analysis/data/phylofactor work spaces/bat spp India.csv")
# load("data/bat_taxonomy_data.Rdata")
# taxonomy <- batphy1 %>%
#   select(c(species, tax)) 
# 
# unique(seroprevalence$sampling_location)
# 
# india.outbreaks <- seroprevalence %>%
#   filter(grepl('india', sampling_location)) %>%
#   filter(virus == 'Henipavirus')
# 
# india.outbreaks %>%
#   select(species) %>%
#   unique()
# 
# india.outbreaks %>%
#   select(sampling_location, sampling_location_two) %>%
#   unique()
# 
# bat_spp_India <- bat_spp_India %>%
#   mutate(species.tnsfrm = tolower(MSW05_Binomial)) %>%
#   mutate(species.tnsfrm = gsub("_", " ", species.tnsfrm)) 
# 
# bat_spp_India.join <- left_join(bat_spp_India, seroprevalence, by = c("species.tnsfrm" = "species")) 
# bat_spp_India.join <- left_join(bat_spp_India.join, taxonomy, by =c('species.tnsfrm' = 'species'))
# 
# bat_spp_India.join <- bat_spp_India.join %>%
#   filter(virus != 'Filovirus'| is.na(virus))
# 
# bat_spp_India.join <- bat_spp_India.join %>%
#   mutate(Pteropodidae = ifelse(grepl('Pteropodidae', tax,), TRUE, FALSE))
# 
# x <- bat_spp_India.join %>%
#   filter(Pteropodidae == TRUE) %>%
#   #filter(seroprevalence_percentage >0) %>%
#   select(species.tnsfrm, methodology, seroprevalence_percentage) %>%
#   unique()
# 
# 
# 
# bat_spp_India.join %>%
#   select(sampling_location, sampling_location_two, sampling_location_three) %>%
#   unique()
# 
# bat_spp_India.join.pcr <- bat_spp_India.join %>%
#   filter(methodology == "PCR based method")
# 
# bat_spp_India.join <- bat_spp_India.join %>%
#   mutate(Vespertilionoidea = ifelse(grepl('Vespertilionoidea', tax,), TRUE, FALSE))
# 
# bat_spp_India.join %>%
#   filter(Vespertilionoidea == TRUE) %>%
#   filter(is.na(virus)) %>%
#   select(species.tnsfrm) %>%
#   unique()
# 
# #...........evenmore!!!!.......
# #...........evenmore!!!!.......
# #...........evenmore!!!!.......
# #...........evenmore!!!!.......
# #...........evenmore!!!!.......
# #...........evenmore!!!!.......
# #...........evenmore!!!!.......
# #...........evenmore!!!!.......
# 
#library(tidyverse)
load('data/seroprevalence.Rdata')
batphy_for_rotl_update <- read_csv("~/Desktop/BDEL/BZDEL/meta_analysis/data/batphy_for_rotl update.csv")

nipah.question <- seroprevalence %>%
  filter(virusold == 'Nipah'|virusold == "Henipavirus")

nipah.question.missing <- nipah.question%>%
  filter(!(species %in% batphy_for_rotl_update$species))

unique(nipah.question.missing$species)
# 
# nipah.question$species=plyr::revalue(nipah.question$species,c("rhinolophus refulgens"="rhinolophus lepidus",
#                                                               "hipposideros cf caffer"="hipposideros caffer",
#                                                               "hipposideros cf caffer/ruber"="hipposideros caffer",
#                                                               "hipposideros cf ruber"="hipposideros ruber",
#                                                               "natalus lanatus"="natalus mexicanus",
#                                                               "myonycteris leptodon"="myonycteris torquata",
#                                                               "nanonycteris veldkampiii"="nanonycteris veldkampii",
#                                                               "tadarida plicata"="chaerephon plicatus"))
# 
# nipah.question$species=plyr::revalue(nipah.question$species,c("Lissonycteris angolensis"="Myonycteris angolensis",
#                                                               "Nyctalus plancyi"="Nyctalus velutinus",
#                                                               "Rhinolophus stheno"="Rhinolophus microglobosus",
#                                                               "Natalus stramineus"="Natalus stramineus mexicanus",
#                                                               "Neoromicia somalicus"="Neoromicia malagasyensis"))
# 
# nipah.question.missing <- nipah.question%>%
#   filter(!(species %in% batphy_for_rotl_update$species))
# 
# unique(nipah.question.missing$species)
# 
# nipah.question <- nipah.question %>%
#   mutate(seroprevalence_prevalence_percentage.cat  = ifelse(seroprevalence_percentage >0, 1, 0)) 
# 
# nipah.question <- left_join(nipah.question, batphy_for_rotl_update[,c('species', 'region')])
# 
# nipah.question.table <-  nipah.question %>%       
#   select(seroprevalence_prevalence_percentage.cat, species, region) %>%
#   group_by(species, region) %>%
#   summarise(nipah.pos.binary=sum(seroprevalence_prevalence_percentage.cat, na.rm=TRUE))  %>%
#   mutate(nipah.pos.binary = ifelse(nipah.pos.binary > 0, 1, 0))
# 
# 
# revised_HNV_positive_bats_uncleaned <- read_csv("~/Dropbox_gmail/Dropbox/bat virus meta-analysis/revised HNV positive bats_uncleaned.csv")
# revised_NiV_positive_bats_uncleaned <- read_csv("~/Desktop/revised NiV positive bats_uncleaned.csv")
# 
# revised_NiV_positive_bats_uncleaned <-revised_NiV_positive_bats_uncleaned %>%
#   rename(niv.pos = hnvpos)
# 
# x <- full_join(revised_HNV_positive_bats_uncleaned, revised_NiV_positive_bats_uncleaned)
# 
# y <- left_join(x, batphy_for_rotl_update[,c('species', 'region')])
# 
# y <- y %>%
#   rename(wilson.reeder.region = region)
# 
# 
# 
# write.csv(y, file = "/Users/buckcrowley/Desktop/hnv.pos.for.barabara.csv")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# .................................pteropodidae bats......................................

species.list <- ggtree::get.offspring.tip(tree, node=200)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ",
                            replacement = "",
                            group_taxonomy_list[,1])
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(tree, node=200)
mean(ggtree.object$data[ggtree.object$data$label %in% species.list,]$angle-90)

#.................................vester bats......................................

species.list <- ggtree::get.offspring.tip(tree, node=119)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ",
                            replacement = "",
                            group_taxonomy_list[,1])
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(tree, node=119)
mean(ggtree.object$data[ggtree.object$data$label %in% species.list,]$angle-90)

# 
# 
# 
# 
# 
# 
# 
