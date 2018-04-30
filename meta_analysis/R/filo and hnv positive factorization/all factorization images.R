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

rm(list=ls())

# filo all bats factors-----------------------------------------------

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
}

#usefull little thing
#group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; Yinpterochiroptera; Pteropodidae; ", 
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            x = group_taxonomy_list$`taxonomy[match(species, taxonomy[, 1]), 2]`) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list %>% table()


B <- lapply(pf$groups, function(l) l[[1]])

# B <- bins(pf$basis[,1:10])
# B <- B[2:11]
Z <- Data$Z
probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

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

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue')

ggsave("figures/hnv global sampling effort tree.png", bg = "transparent", height = 18, width = 18)

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

B <- pf$groups[[1]][[1]]
# B <- bins(pf$basis[,1:10])
# B <- B[2:11]
Z <- Data$Z
probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

probs <- sum(probs)/length(probs)
# 0.547619 for Pteropodidae

B <- lapply(pf$groups, function(l) l[[1]])
# B <- bins(pf$basis[,1:10])
# B <- B[2:11]
Z <- Data$Z
probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)


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

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue')  +
  geom_segment(data= d,aes(x=x,y=y,xend=xend,yend=yend, size= Data$log_effort, colour = 'blue'))

ggsave("figures/hnv sampling effort tree.png", bg = "transparent", height = 18, width = 18)

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


B <- lapply(pf$groups, function(l) l[[1]])
# B <- bins(pf$basis[,1:10])
# B <- B[2:11]
Z <- Data$Z
probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

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
colfcn <- function(n) return(c("#33FF00FF", "#440154FF", "#FF9900FF"))
pf.tree <- pf.tree(pf, lwd=1, factors = 1:3, color.fcn=colfcn, branch.length = "none", bg.color = NA)

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue') 
  
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

colfcn <- function(n) return(c("#00FF66FF"))

pf.tree <- pf.tree(pf, lwd=1, factors = 1, color.fcn=colfcn, branch.length = "none", bg.color = NA)

d <- data.frame(x=pf.tree$ggplot$data[1:104,'x'] + .5,
                xend=pf.tree$ggplot$data[1:104,'x'] + 1+ Data$log_effort,
                y=pf.tree$ggplot$data[1:104,'y'],
                yend=pf.tree$ggplot$data[1:104,'y'] )

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue') + 
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

B <- lapply(pf$groups, function(l) l[[1]])
# B <- bins(pf$basis[,1:10])
# B <- B[2:11]
Z <- Data$Z
probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)


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

colfcn <- function(n) return(c("#00FF66FF", "#440154FF"))

pf.tree <- pf.tree(pf, lwd=1, factors = 1:2, color.fcn=colfcn, branch.length = "none", bg.color = NA)

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue') 

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
