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

revised_NiV_positive_bats_uncleaned <- read_csv("~/BZDEL/meta_analysis/data/phylofactor work spaces/revised HNV positive bats_uncleaned.csv")

revised_NiV_positive_bats_uncleaned %>%
  filter(!(species %in% batphy1$species)) 

revised_NiV_positive_bats_uncleaned$species =
  plyr::revalue(revised_NiV_positive_bats_uncleaned$species,c("ballionycterus maculata"="balionycteris maculata",
                                                              "scotophilus heathi" = "scotophilus heathii"
  ))

x <- revised_NiV_positive_bats_uncleaned %>%
  filter(!(species %in% batphy1$species)) 

#i'm okay with this for now

batphy1 <- left_join(batphy1, revised_NiV_positive_bats_uncleaned)
nrow(revised_NiV_positive_bats_uncleaned) ==  nrow(filter(batphy1, !is.na(hnvpos))) + nrow(revised_NiV_positive_bats_uncleaned %>% filter(!(species %in% batphy1$species)) )

Data <- batphy1 %>%
  filter(unique_name %in% bat_spp_India$MSW05_Binomial) %>%
  rename(niv_pos = hnvpos) %>%
  mutate(niv_sample = ifelse(is.na(niv_pos), 0, 1)) %>%
  select(-X1) %>% 
  mutate(log_hnv_samps = log(hnv_samps)) %>%
  select(c(hnv_samps, niv_pos, species, niv_sample)) %>%
  rename(Species = species) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  mutate(Species = stri_trans_totitle(Species)) %>%
  mutate(Sample = 1)

names(Data) <- c('effort','niv.pos', 'Species','Z', 'Sample')

tree <- ape::drop.tip(bat_tree,bat_tree$tip.label[!(bat_tree$tip.label %in% Data$Species)])

pf <- gpf(Data,tree,frmla.phylo=Z~phylo,nfactors=10,family=binomial,algorithm='phylo')

ncores = 4
tot.reps=1000
reps.per.worker=round(tot.reps/ncores)
sampling_effort = FALSE
source('R/filo and hnv positive factorization/null simulations script.R')

load(file='data/phylofactor work spaces/indian bats sample no sample niv_workspace')


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


load("data/bat_taxonomy_data.Rdata")
source('R/taxonomy group name function.R')

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

group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            x = group_taxonomy_list$`taxonomy[match(species, taxonomy[, 1]), 2]`) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list %>% table()


factor.1 <- pf$groups[[1]][[1]]
factor.2 <- pf$groups[[2]][[1]]
paraphyletic.remainder <-     pf$groups[[2]][[2]]

length(factor.1)
length(factor.2)
length(paraphyletic.remainder)

Z <- Data$Z
factor.1.p <- sapply(factor.1,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.2.p <- sapply(factor.2,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
paraphyletic.remainder.p <- sapply(paraphyletic.remainder,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

factor.1.prob <- sum(factor.1.p)/length(factor.1.p)
factor.2.prob <- sum(factor.2.p)/length(factor.2.p)
paraphyletic.remainder.prob <- sum(paraphyletic.remainder.p)/length(paraphyletic.remainder.p)

factor.1.prob
factor.2.prob
paraphyletic.remainder.prob

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
# 15      Rhinolophoidea        1 #00FF66FF
# 16      Pteropodinae          1 #33FF00FF
# 17      Hipposideros          1 #FF9900FF
load('data/seroprevalence.Rdata')


Data <- Data %>%
  mutate(species.mutate = tolower(gsub("_", " ", Species))) %>%
  mutate(pcr.pos = ifelse(species.mutate %in% pcr.pos.filo$species, 1,0))
colfcn <- function(n) return(c("#00FF66FF", "#440154FF"))

pf.tree <- pf.tree(pf, lwd=1, factors = 1:2, color.fcn=colfcn, branch.length = "none", bg.color = NA)























pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue') +
  ggtree::geom_tippoint(size=5*Data$niv.pos,col='RED') 
  

ggsave("figures/filo no sampling effort tree.png", bg = "transparent", height = 18, width = 18)

Legend <- pf.tree$legend
Legend$names <- c(names.storage[1], "Pteropodinae") #megaloglossus woermanni isn't contained within the Pteropodinae according to our taxonomy search, but it clearly is according to the tree so going to override that 
P <- sapply(probs[1:2],FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend$names <- mapply(paste,Legend$names,P)



Cairo(file='figures/filo no sampling effort tree legend.jpg', 
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

