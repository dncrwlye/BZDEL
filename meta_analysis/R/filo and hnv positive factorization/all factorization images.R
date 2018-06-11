
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

load(file='data/phylofactor work spaces/hnv_workspace_sample_no_sample_all_bat_dataset')

x<- cbind(as.matrix(OBJ[[1]]),as.matrix(OBJ[[2]]))
x<-rbind(as.data.frame(OBJ[[1]]),as.data.frame(OBJ[[2]]),as.data.frame(OBJ[[3]]),as.data.frame(OBJ[[4]]))
x <- t(x)

jpeg(filename = 'figures/hnv sample no sample all bats factors.jpg')
plot(pf$pvals, type='l',ylim = c(0,.05), col='red')
for (i in 1:200)
{
  lines(x[,i])
}
points(deviances,col='red',pch=16,cex=2)
dev.off()

plot(pf$pvals, type='l',ylim = c(0,.05), col='red')
for (i in 1:200)
{
  lines(x[,i])
}
#10 factors!!!! ????? :o

plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,15))
for (rr in 1:497){
  lines(c(S[rr,1],diff(S[rr,])),col=rgb(0,0,0,0.2))
}


# filo all bats factors-----------------------------------------------

rm(list=ls())

load(file='data/phylofactor work spaces/filo_workspace_sample_no_sample_all_bat_dataset')

x<- cbind(as.matrix(OBJ[[1]]),as.matrix(OBJ[[2]]))
x<-rbind(as.data.frame(OBJ[[1]]),as.data.frame(OBJ[[2]]),as.data.frame(OBJ[[3]]),as.data.frame(OBJ[[4]]))
x <- t(x)
jpeg(filename = 'figures/filo sample no sample all bats factors.jpg')
plot(pf$pvals, type='l',ylim = c(0,.05), col='red')
for (i in 1:200)
{
  lines(x[,i])
}
dev.off()

plot(pf$pvals, type='l',ylim = c(0,.05), col='red')
for (i in 1:200)
{
  lines(x[,i])
}
#10 factors again!!!! ????? :o

rm(list=ls())

# hnv SE+ factors-----------------------------------------------

load(file='data/phylofactor work spaces/hnv_workspace')

jpeg(filename = 'figures/hnv sampling effort factors.jpg')
plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,15))
for (rr in 1:497){
  lines(c(S[rr,1],diff(S[rr,])),col=rgb(0,0,0,0.2))
}
points(deviances,col='red',pch=16,cex=2)
dev.off()

bb <- (t(diff(t(S))))
bb <- cbind(S[,1],bb)
Ojb.bb <- c(Obj[1],diff(Obj))
for (i in 1:10)
{
  print(ecdf(bb[,i])(Ojb.bb[i]))
}
#1 factor 

plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,15))
for (rr in 1:497){
  lines(c(S[rr,1],diff(S[rr,])),col=rgb(0,0,0,0.2))
}


#going to go with 1 factor! 

rm(list=ls())

# hnv SE- factors-----------------------------------------------

load(file='data/phylofactor work spaces/hnv_workspace_no_sampling_effort')
jpeg(filename = 'figures/hnv no sampling effort factors.jpg')

plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,40))
for (rr in 1:497){
  lines(c(S[rr,1],diff(S[rr,])),col=rgb(0,0,0,0.2))
}
points(deviances,col='red',pch=16,cex=2)
dev.off()

plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,40))
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

#going to say....3 factors!!!

# filo SE+ factors-----------------------------------------------

load(file='data/phylofactor work spaces/filo_workspace')
jpeg(filename = 'figures/filo with sampling effort factors.jpg')

plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,15))
for (rr in 1:497){
  lines(c(S[rr,1],diff(S[rr,])),col=rgb(0,0,0,0.2))
}
points(deviances,col='red',pch=16,cex=2)
dev.off()

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

rm(list=ls())

#1 factor?

# filo SE- factors-----------------------------------------------

rm(list=ls())

load(file='data/phylofactor work spaces/filo_workspace_no_sampling_effort')
jpeg(filename = 'figures/filo no sampling effort factors.jpg')

plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,25))
for (rr in 1:497){
  lines(c(S[rr,1],diff(S[rr,])),col=rgb(0,0,0,0.2))
}
points(deviances,col='red',pch=16,cex=2)
dev.off()

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

#looks like....2 factors?

# FIGURE MAKING ----

# hnv all bats+ figure 1 ----

# all the factors....

rm(list=ls())
load(file='data/phylofactor work spaces/hnv_workspace_sample_no_sample_all_bat_dataset')
load("data/bat_taxonomy_data.Rdata")
source('R/taxonomy group name function.R')

taxonomy <- batphy1 %>%
  select(c(species, tax)) 

names.storage <- list()

for (i in 1:10)
{
  indexes = pf$groups[[i]][[1]]
  species <- gsub("_", " ", tolower(tree$tip.label[indexes]))
  if(is.na(match(species,taxonomy[,1])))
  {
    break
  }
  group_taxonomy_list <- as.data.frame(taxonomy[match(species,taxonomy[,1]),2])
  names.storage[i] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
  
  group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                              replacement = "", 
                              x = group_taxonomy_list$`taxonomy[match(species, taxonomy[, 1]), 2]`) 
  group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
  
  print(group_taxonomy_list %>% table())
  print("NEW GROUP
        NEW GROUP
        NEW GROUP
        NEW GROUP
        NEW GROUP
        NEW GROUP
        NEW GROUP
        NEW GROUP
        NEW GROUP
        NEW GROUP")
}

#usefull little thing
#group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; Yinpterochiroptera; Pteropodidae; ", 

factor.1 <- pf$groups[[1]][[1]]
factor.2 <- pf$groups[[2]][[1]]
factor.3 <- pf$groups[[3]][[1]]
factor.4 <- pf$groups[[4]][[1]]
factor.5 <- pf$groups[[5]][[1]]
factor.6 <- pf$groups[[6]][[1]]
factor.7 <- pf$groups[[7]][[1]]
factor.8 <- pf$groups[[8]][[1]]
factor.9 <- pf$groups[[9]][[1]]
factor.10 <- pf$groups[[10]][[1]]
paraphyletic.remainder <- pf$groups[[10]][[2]]

length(factor.1) +
length(factor.2) +
length(factor.3) +
length(factor.4) +
length(factor.5)+
length(factor.6)+
length(factor.7)+
length(factor.8)+
length(factor.9)+
length(factor.10)+
length(paraphyletic.remainder)

Z <- Data$Z
factor.1.p <- sapply(factor.1,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.2.p <- sapply(factor.2,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.3.p <- sapply(factor.3,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.4.p <- sapply(factor.4,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.5.p <- sapply(factor.5,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.6.p <- sapply(factor.6,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.7.p <- sapply(factor.7,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.8.p <- sapply(factor.8,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.9.p <- sapply(factor.9,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.10.p <- sapply(factor.10,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
paraphyletic.remainder.p <- sapply(paraphyletic.remainder,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

factor.1.prob <- sum(factor.1.p)/length(factor.1.p)
factor.2.prob <- sum(factor.2.p)/length(factor.2.p)
factor.3.prob <- sum(factor.3.p)/length(factor.3.p)
factor.4.prob <- sum(factor.4.p)/length(factor.4.p)
factor.5.prob <- sum(factor.5.p)/length(factor.5.p)
factor.6.prob <- sum(factor.6.p)/length(factor.6.p)
factor.7.prob <- sum(factor.7.p)/length(factor.7.p)
factor.8.prob <- sum(factor.8.p)/length(factor.8.p)
factor.9.prob <- sum(factor.9.p)/length(factor.9.p)
factor.10.prob <- sum(factor.10.p)/length(factor.10.p)
paraphyletic.remainder.prob <- sum(paraphyletic.remainder.p)/length(paraphyletic.remainder.p)

factor.1.prob
factor.2.prob
factor.3.prob
factor.4.prob
factor.5.prob
factor.6.prob
factor.7.prob
factor.8.prob
factor.9.prob
factor.10.prob
paraphyletic.remainder.prob


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

colfcn <- function(n) return(c("#440154FF",  
                               "#FF0066FF", 
                               "#3E4A89FF", 
                               "#31688EFF",
                               "#26828EFF", 
                               "#1F9E89FF", 
                               "#35B779FF",
                               "#6DCD59FF",
                               "#B4DE2CFF", 
                               "#FDE725FF"  
))


pf.tree <- pf.tree(pf, lwd=1, color.fcn = colfcn, factors = 1:10, branch.length = "none", bg.color = NA)
Data <- Data %>%
  mutate(ib = ifelse(Species %in% text2$V3, 1, 0))

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=20*Data$ib,col='red') + 
  ggtree::geom_tippoint(size=10*Data$Z,col='blue')  

ggsave("figures/hnv global sampling effort tree india bats.png", bg = "transparent", height = 18, width = 18)

#something strange is going on, we have a factor 2 that is only one tip, but 
#in the bins it contains many many species. This is odd....
Legend <- pf.tree$legend

Legend$names <- c('Pteropodidae: Rousettus, Epomophorus, ...', names.storage[2:10])

P <- sapply(probs[1:10],FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend$names <- mapply(paste,Legend$names,P)

Cairo(file='figures/hnv global sampling effort legend.jpg', 
      type="png",
      units="in", 
      width=6, 
      height=4, 
      pointsize=12, 
      dpi=200)
plot.new()
plot.new()
legend('topleft', legend=Legend$names,fill=Legend$colors,cex=1, bg="transparent", bty="n")
dev.off()
# filo all bats+ figure 1 ----

# all the factors....

rm(list=ls())
load(file='data/phylofactor work spaces/filo_workspace_sample_no_sample_all_bat_dataset')
load("data/bat_taxonomy_data.Rdata")
source('R/taxonomy group name function.R')

taxonomy <- batphy1 %>%
  select(c(species, tax)) 

names.storage <- list()

for (i in 1:10)
{
  indexes = pf$groups[[i]][[1]]
  species <- gsub("_", " ", tolower(tree$tip.label[indexes]))
  if(is.na(match(species,taxonomy[,1])))
  {
    break
  }
  group_taxonomy_list <- as.data.frame(taxonomy[match(species,taxonomy[,1]),2])
  names.storage[i] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
}

#usefull little thing
#group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; Yinpterochiroptera; Pteropodidae; ", 
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            x = group_taxonomy_list$`taxonomy[match(species, taxonomy[, 1]), 2]`) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list
group_taxonomy_list %>% table()

factor.1 <- pf$groups[[1]][[1]]
factor.2 <- pf$groups[[2]][[1]]
factor.3 <- pf$groups[[3]][[1]]
factor.4 <- pf$groups[[4]][[1]]
factor.5 <- pf$groups[[5]][[1]]
factor.6 <- pf$groups[[6]][[1]]
factor.7 <- pf$groups[[7]][[1]]
factor.8 <- pf$groups[[8]][[1]]
factor.9 <- pf$groups[[9]][[1]]
factor.10 <- pf$groups[[10]][[1]]
paraphyletic.remainder <- pf$groups[[10]][[2]]

length(factor.1) +
  length(factor.2) +
  length(factor.3) +
  length(factor.4) +
  length(factor.5)+
  length(factor.6)+
  length(factor.7)+
  length(factor.8)+
  length(factor.9)+
  length(factor.10)+
  length(paraphyletic.remainder)

Z <- Data$Z
factor.1.p <- sapply(factor.1,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.2.p <- sapply(factor.2,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.3.p <- sapply(factor.3,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.4.p <- sapply(factor.4,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.5.p <- sapply(factor.5,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.6.p <- sapply(factor.6,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.7.p <- sapply(factor.7,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.8.p <- sapply(factor.8,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.9.p <- sapply(factor.9,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.10.p <- sapply(factor.10,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
paraphyletic.remainder.p <- sapply(paraphyletic.remainder,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

factor.1.prob <- sum(factor.1.p)/length(factor.1.p)
factor.2.prob <- sum(factor.2.p)/length(factor.2.p)
factor.3.prob <- sum(factor.3.p)/length(factor.3.p)
factor.4.prob <- sum(factor.4.p)/length(factor.4.p)
factor.5.prob <- sum(factor.5.p)/length(factor.5.p)
factor.6.prob <- sum(factor.6.p)/length(factor.6.p)
factor.7.prob <- sum(factor.7.p)/length(factor.7.p)
factor.8.prob <- sum(factor.8.p)/length(factor.8.p)
factor.9.prob <- sum(factor.9.p)/length(factor.9.p)
factor.10.prob <- sum(factor.10.p)/length(factor.10.p)
paraphyletic.remainder.prob <- sum(paraphyletic.remainder.p)/length(paraphyletic.remainder.p)

factor.1.prob
factor.2.prob
factor.3.prob
factor.4.prob
factor.5.prob
factor.6.prob
factor.7.prob
factor.8.prob
factor.9.prob
factor.10.prob
paraphyletic.remainder.prob

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
#@colfcn <- function(n) return(c("#00FF66FF", "#440154FF"))

colfcn <- function(n) return(c("#FF0066FF", #Yangochiroptera ~ similar to Vespertilionoidea 
                               "#B4DE2CFF", #Molossinae
                               "#FDE725FF", #Miniopterus
                               "#31688EFF", #Myotis ~ Myotis
                               "#6DCD59FF", #Pteropodidae: Nyctimene, Dobsonia, Pteropus, ...
                               "#CCFF00FF", #Taphozoinae
                               "#440154FF", #Pteropodidae: Rousettus, Epomophorus, ...
                               "#FF9900FF", #Vespertilionini
                               "#FF0000FF", #Pteropus 
                               "#33FF00FF"  #Macroglossus 
))

pf.tree <- pf.tree(pf, lwd=1, factors = 1:10,  color.fcn=colfcn, branch.length = "none", bg.color = NA)

#species.list <- ggtree::get.offspring.tip(pf$tree, node=1086)

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue') +
  ggtree::geom_cladelabel(node=1086, label="Noctilionoidea", 
                          color="black", angle="auto", offset=2, offset.text = 2, align=TRUE) 

#geom_segment(data= d,aes(x=x,y=y,xend=xend,yend=yend, size= Data$log_effort, colour = 'blue'))


#ggtree::geom_label2(aes(label=node), size=2, color="darkred", alpha=0.5) +
#ggtree::geom_hilight(node=1086, fill = 'steelblue')

ggsave("figures/filo global sampling effort.png", bg = "transparent", height = 18, width = 18)

#something strange is going on, we have a factor 2 that is only one tip, but 
#in the bins it contains many many species. This is odd....
Legend <- pf.tree$legend
Legend$names <- names.storage[1:10]

Legend$names <- c(names.storage[1:4], 'Pteropodidae: Nyctimene, Dobsonia, Pteropus, ...',names.storage[6], 'Pteropodidae: Rousettus, Epomophorus, ...', names.storage[8:10])

P <- sapply(probs[1:10],FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend$names <- mapply(paste,Legend$names,P)

Cairo(file='figures/filo global sampling effort legend.jpg', 
      type="png",
      units="in", 
      width=6, 
      height=4, 
      pointsize=12, 
      dpi=200)
plot.new()
plot.new()
legend('topleft', legend=Legend$names,fill=Legend$colors,cex=1, bg="transparent", bty="n")
dev.off()

# Africa only filo figure 1 ----

# all the factors....

rm(list=ls())
load(file='data/phylofactor work spaces/filo_workspace_sample_no_sample_all_african_bat_dataset')
load("data/bat_taxonomy_data.Rdata")
source('R/taxonomy group name function.R')

taxonomy <- batphy1 %>%
  select(c(species, tax)) 

names.storage <- list()

for (i in 1:10)
{
  indexes = pf$groups[[i]][[1]]
  species <- gsub("_", " ", tolower(tree$tip.label[indexes]))
  if(is.na(match(species,taxonomy[,1])))
  {
    break
  }
  group_taxonomy_list <- as.data.frame(taxonomy[match(species,taxonomy[,1]),2])
  names.storage[i] <- gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
}

#usefull little thing
#group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; Yinpterochiroptera; Pteropodidae; ", 
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            x = group_taxonomy_list$`taxonomy[match(species, taxonomy[, 1]), 2]`) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list %>% table()

#5
#Nyctimene, Dobsonia, Pteropus, ...

#7
#Rousettus, Epomophorus, ...

B <- lapply(pf$groups, function(l) l[[1]])

# B <- bins(pf$basis[,1:10])
# B <- B[2:11]
Z <- Data$Z
probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

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
#@colfcn <- function(n) return(c("#00FF66FF", "#440154FF"))

colfcn <- function(n) return(c("#FF0066FF", #Yangochiroptera ~ similar to Vespertilionoidea 
                               "#B4DE2CFF", #Molossinae
                               "#FDE725FF", #Miniopterus
                               "#31688EFF", #Myotis ~ Myotis
                               "#440154FF", #Pteropodidae: Nyctimene, Dobsonia, Pteropus, ...
                               "#FF0000FF", #Taphozoinae
                               "#6DCD59FF", #Pteropodidae: Rousettus, Epomophorus, ...
                               "#FF9900FF", #Vespertilionini
                               "#CCFF00FF", #Pteropus 
                               "#33FF00FF"  #Macroglossus 
))

pf.tree <- pf.tree(pf, lwd=1, factors = 1:10,  color.fcn=colfcn, branch.length = "none", bg.color = NA)

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue')

ggsave("figures/filo global sampling effort.png", bg = "transparent", height = 18, width = 18)

#something strange is going on, we have a factor 2 that is only one tip, but 
#in the bins it contains many many species. This is odd....
Legend <- pf.tree$legend
Legend$names <- names.storage[1:10]

Legend$names <- c(names.storage[1:4], 'Pteropodidae: Nyctimene, Dobsonia, Pteropus, ...',names.storage[6], 'Pteropodidae: Rousettus, Epomophorus, ...', names.storage[8:10])

P <- sapply(probs[1:10],FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend$names <- mapply(paste,Legend$names,P)

Cairo(file='figures/filo global sampling effort legend.jpg', 
      type="png",
      units="in", 
      width=6, 
      height=4, 
      pointsize=12, 
      dpi=200)
plot.new()
plot.new()
legend('topleft', legend=Legend$names,fill=Legend$colors,cex=1, bg="transparent", bty="n")
dev.off()


# hnv SE+ figure 1 ----

# 1 factor 

rm(list=ls())
load(file='data/phylofactor work spaces/hnv_workspace')
load("data/bat_taxonomy_data.Rdata")
source('R/taxonomy group name function.R')

pf$call
x <- pf$models[[1]] %>% summary()
1/exp(x$coefficients[2,1])

exp(-2.157)

taxonomy <- batphy1 %>%
  select(c(species, tax)) 

names.storage <- list()

for (i in 1)
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
# hnv SE- figure 1 --------------------------------------------------------------------------------------

rm(list=ls())

load(file='data/phylofactor work spaces/hnv_workspace_no_sampling_effort')
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
factor.3 <- pf$groups[[3]][[1]]
paraphyletic.remainder <-     pf$groups[[3]][[2]]

length(factor.1)
length(factor.2)
length(factor.3)
length(paraphyletic.remainder)

sum(15 + 27 + 2 + 99)

Z <- Data$Z
factor.1.p <- sapply(factor.1,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.2.p <- sapply(factor.2,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.3.p <- sapply(factor.3,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
paraphyletic.remainder.p <- sapply(paraphyletic.remainder,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

factor.1.prob <- sum(factor.1.p)/length(factor.1.p)
factor.2.prob <- sum(factor.2.p)/length(factor.2.p)
factor.3.prob <- sum(factor.3.p)/length(factor.3.p)
paraphyletic.remainder.prob <- sum(paraphyletic.remainder.p)/length(paraphyletic.remainder.p)

factor.1.prob
factor.2.prob
factor.3.prob
paraphyletic.remainder.prob

# B <- lapply(pf$groups, function(l) l[[1]])
# # B <- bins(pf$basis[,1:10])
# # B <- B[2:11]
# Z <- Data$Z
# probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

#Pteropodinae 0.8666667

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
# 15      Pteropodinae          1 #33FF00FF
# 15      Hipposideros          1 #FF9900FF
load('data/seroprevalence.Rdata')

pcr.pos.hnv <- seroprevalence %>%
  filter(methodology == 'PCR based method') %>%
  filter(seroprevalence_percentage > 0) %>%
  filter(virus == 'Henipavirus')

Data <- Data %>%
  mutate(species.mutate = tolower(gsub("_", " ", Species))) %>%
  mutate(pcr.pos = ifelse(species.mutate %in% pcr.pos.hnv$species, 1,0))

colfcn <- function(n) return(c("#33FF00FF", "#440154FF", "#FF9900FF"))
pf.tree <- pf.tree(pf, lwd=1, factors = 1:3, color.fcn=colfcn, branch.length = "none", bg.color = NA)

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue') +
  ggtree::geom_tippoint(size=5*Data$pcr.pos,col='red') 

ggsave("figures/hnv no sampling effort tree.png", bg = "transparent", height = 18, width = 18)

Legend <- pf.tree$legend
Legend$names <- c('Pteropodinae: Pteropus, Eidolon...' ,names.storage[2:3])
P <- sapply(probs[1:3],FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend$names <- mapply(paste,Legend$names,P)
Cairo(file='figures/hnv no sampling effort tree legend.jpg', 
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
#.............................indian bats aside......................................
#.............................indian bats aside......................................
#.............................indian bats aside......................................
#.............................indian bats aside......................................
#.............................indian bats aside......................................
#.............................indian bats aside......................................
#.............................indian bats aside......................................
#.............................indian bats aside......................................
#.............................indian bats aside......................................
#.............................indian bats aside......................................
#.............................indian bats aside......................................
#.............................indian bats aside......................................
#.............................indian bats aside......................................
#.............................indian bats aside......................................
#.............................indian bats aside......................................
#.............................indian bats aside......................................
# indian bats aside ----

colfcn <- function(n) return(c("#33FF00FF"))

pf.tree <- pf.tree(pf, lwd=1, factors = 2, color.fcn=colfcn, branch.length = "none", bg.color = NA)

bat_spp_India <- read_csv("C:/Users/r83c996/Desktop/Bats_Kerala.csv")
#bat_spp_India <- read_csv("~/Desktop/BDEL/BZDEL/meta_analysis/data/phylofactor work spaces/bat spp India.csv")

bat_spp_India %>%
  filter(!(MSW05_binomial %in% batphy1$unique_name))
bat_spp_India$MSW05_binomial =plyr::revalue(bat_spp_India$MSW05_binomial,c("Myotis_blythii"="Myotis_oxygnathus"))
bat_spp_India %>%
  filter(!(MSW05_binomial %in% batphy1$unique_name))
bat_spp_India <- bat_spp_India %>%
  mutate(MSW05_binomial = gsub(" ", "_", MSW05_binomial))

Data <- Data %>%
  mutate(IB = ifelse(Species %in% bat_spp_India$MSW05_binomial, 1, 0))
pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_text2(aes(subset=!isTip,label=node))

# .................................vesper bats......................................
species.list <- ggtree::get.offspring.tip(pf$tree, node=179)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            group_taxonomy_list[,1]) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(pf$tree, node=179)
mean(pf.tree$ggplot$data[pf.tree$ggplot$data$label %in% species.list,]$angle-90)

# .................................pteropodidae bats......................................
species.list <- ggtree::get.offspring.tip(pf$tree, node=244)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            group_taxonomy_list[,1]) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(pf$tree, node=244)
mean(pf.tree$ggplot$data[pf.tree$ggplot$data$label %in% species.list,]$angle-90)

# .................................218 bats......................................
species.list <- ggtree::get.offspring.tip(pf$tree, node=218)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            group_taxonomy_list[,1]) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(pf$tree, node=244)
mean(pf.tree$ggplot$data[pf.tree$ggplot$data$label %in% species.list,]$angle-90)



pf.tree <- pf.tree(pf, lwd=1, factors = 10, color.fcn=colfcn, branch.length = "none", bg.color = NA)


d <- data.frame(x=pf.tree$ggplot$data[1:143,'x'] + 1,
                xend=pf.tree$ggplot$data[1:143,'x'] + 1 + Data$IB,
                y=pf.tree$ggplot$data[1:143,'y'],
                yend=pf.tree$ggplot$data[1:143,'y'] )

load('data/seroprevalence.Rdata')

pcr.pos.hnv <- seroprevalence %>%
  filter(methodology == 'PCR based method') %>%
  filter(seroprevalence_percentage > 0) %>%
  filter(virusold == 'Nipah') %>%
  select(species) %>%
  unique()

Data <- Data %>%
  mutate(species.mutate = tolower(gsub("_", " ", Species))) %>%
  mutate(pcr.pos = ifelse(species.mutate %in% pcr.pos.hnv$species, 1,0))

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue', alpha = .5) +
  ggtree::geom_cladelabel(node=179, label="Vespertilionidae", 
                        color="black", angle=222, offset=2, offset.text = 2) +
  ggtree::geom_cladelabel(node=244, label="Pteropodidae", 
                        color="black", angle=32, offset=2, offset.text = 2)  +
  geom_segment(data= d,aes(x=x,y=y,xend=xend,yend=yend, size= Data$IB, colour = 'purple')) +
  ggtree::geom_tippoint(size=5*Data$pcr.pos, shape = 17, col='red') +
  ggtree::geom_tiplab(aes(angle=angle))

ggsave("figures/indian bats info.png", bg = "transparent", height = 18, width = 18)


# filo SE + figure 1 --------------------------------------------------------------------------------------

#1 factor again 
rm(list=ls())

load(file='data/phylofactor work spaces/filo_workspace')
load("data/bat_taxonomy_data.Rdata")
source('R/taxonomy group name function.R')

taxonomy <- batphy1 %>%
  select(c(species, tax)) 

names.storage <- list()

for (i in 1)
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
group_taxonomy_list

factor.1 <- pf$groups[[1]][[1]]
paraphyletic.remainder <-     pf$groups[[1]][[2]]

length(factor.1)
length(paraphyletic.remainder)


Z <- Data$Z
factor.1.p <- sapply(factor.1,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
paraphyletic.remainder.p <- sapply(paraphyletic.remainder,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

factor.1.prob <- sum(factor.1.p)/length(factor.1.p)
paraphyletic.remainder.prob <- sum(paraphyletic.remainder.p)/length(paraphyletic.remainder.p)

factor.1.prob
paraphyletic.remainder.prob


B <- pf$groups[[1]][[1]]
# B <- bins(pf$basis[,1:10])
# B <- B[2:11]
Z <- Data$Z
probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
probs <- sum(probs)/length(probs)
#Rhinolophoidea at  0.125

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

#colors 3 "Rhinolophoidea (Rhinolophidae & Hipposideridae)"  #00FF66FF

load('data/seroprevalence.Rdata')

pcr.pos.filo <- seroprevalence %>%
  filter(methodology == 'PCR based method') %>%
  filter(seroprevalence_percentage > 0) %>%
  filter(virus == 'Filovirus')

Data <- Data %>%
  mutate(species.mutate = tolower(gsub("_", " ", Species))) %>%
  mutate(pcr.pos = ifelse(species.mutate %in% pcr.pos.filo$species, 1,0))

colfcn <- function(n) return(c("#00FF66FF"))

pf.tree <- pf.tree(pf, lwd=1, factors = 1, color.fcn=colfcn, branch.length = "none", bg.color = NA)

d <- data.frame(x=pf.tree$ggplot$data[1:104,'x'] + .5,
                xend=pf.tree$ggplot$data[1:104,'x'] + 1+ Data$log_effort,
                y=pf.tree$ggplot$data[1:104,'y'],
                yend=pf.tree$ggplot$data[1:104,'y'] )

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue') +
  ggtree::geom_tippoint(size=5*Data$pcr.pos, shape = 17, col='red') +
  geom_segment(data= d,aes(x=x,y=y,xend=xend,yend=yend, size= Data$log_effort, colour = 'blue'))

ggsave("figures/filo sampling effort tree.png", bg = "transparent", height = 18, width = 18)

Legend <- pf.tree$legend
Legend$names <- names.storage[1]
P <- sapply(probs[1],FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend$names <- mapply(paste,Legend$names,P)
Cairo(file='figures/filo sampling effort tree legend.jpg', 
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

# filo SE- figure 1 --------------------------------------------------------------------------------------

#likely 2 factors 

rm(list=ls())

load(file='data/phylofactor work spaces/filo_workspace_no_sampling_effort')
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

sum(32 + 12 + 60)
sum(32 + 72)

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

pcr.pos.filo <- seroprevalence %>%
  filter(methodology == 'PCR based method') %>%
  filter(seroprevalence_percentage > 0) %>%
  filter(virus == 'Filovirus')

Data <- Data %>%
  mutate(species.mutate = tolower(gsub("_", " ", Species))) %>%
  mutate(pcr.pos = ifelse(species.mutate %in% pcr.pos.filo$species, 1,0))
colfcn <- function(n) return(c("#00FF66FF", "#440154FF"))

pf.tree <- pf.tree(pf, lwd=1, factors = 1:2, color.fcn=colfcn, branch.length = "none", bg.color = NA)

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue') +
  ggtree::geom_tippoint(size=5*Data$pcr.pos, shape = 17, col='red') 

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

# supplemental anlayses --------------------------------------------------------------------------------------

rm(list=ls())

load("~/Desktop/BDEL/BZDEL/meta_analysis/data/seroprevalence.phylo.analysis.Rdata")
load(file='data/phylofactor work spaces/filo_workspace')
load("data/bat_taxonomy_data.Rdata")
source('R/taxonomy group name function.R')

taxonomy <- batphy1 %>%
  select(c(species, tax)) 

names.storage <- list()

i <- 1
indexes = pf$groups[[i]][[1]]
rhino.spp <- gsub("_", " ", tolower(tree$tip.label[indexes]))

rhino <- seroprevalence.phylo.analysis %>%
  filter(virus == 'Filovirus') %>%
  filter(species %in% rhino.spp) 

sum(rhino$sample_size) #how many of these bats have been sampled
sum(rhino$sample_size * rhino$seroprevalence_percentage/100) #how many have returned positive
unique(rhino$sampling_location)
unique(rhino$title)

rhino %>%
  group_by(sampling_location) %>%
  summarize(n=sum(sample_size))


rhino.pos <- rhino %>% 
  filter(seroprevalence_percentage > 0)

unique(rhino.pos$title)

rhino.2 <- batphy1 %>%
  filter(species %in% rhino.spp)%>%
  filter(filo_positive == 1)


#supp.compar. anlayses --------------------------------------------------------------------------------------

rm(list=ls())

load(file='data/phylofactor work spaces/hnv_workspace_sample_no_sample_all_bat_dataset')

plot(pf$pvals, type='l',ylim = c(0,.01), col='red')

load(file='data/phylofactor work spaces/filo_workspace_sample_no_sample_all_bat_dataset')

lines(pf$pvals, type='l', col='green')

#supp.positive compare anlayses --------------------------------------------------------------------------------------

#the p value of the null simulations 
load(file='data/phylofactor work spaces/hnv_workspace')

m1 <- nrow(Data)
E1 <- dim(pf$tree$edge)[1]
deviance.hnv <- anova(pf$models[[1]])$Deviance[2]

load(file='data/phylofactor work spaces/filo_workspace')

m2 <- nrow(Data)
E2 <- dim(pf$tree$edge)[1]
deviance.filo <- anova(pf$models[[1]])$Deviance[2]

m=1e5
F1 <- matrix(rchisq(E1*m,df = 1),ncol=E1)
F1 <- apply(F1,1,max)
F2 <- matrix(rchisq(E2*m,df = 1),ncol=E2)
F2 <- apply(F2,1,max)
hist(abs(F1-F2),breaks=100)
cdf <- ecdf(abs(F1-F2))
cdf(abs(deviance.hnv-deviance.filo))



