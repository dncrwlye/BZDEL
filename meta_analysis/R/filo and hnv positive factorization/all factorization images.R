#............phylofactor visualization script...................................

# Libraries ---------------------------------------------------------------

setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)

# hnv SE+ factors-----------------------------------------------

load(file='data/phylofactor work spaces/hnv_workspace')

jpeg(filename = 'figures/hnv sampling effort factors.jpg')
plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,15))
for (rr in 1:140){
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

plot((Obj), type = 'o')
for (rr in 1:140){
  lines((S[rr,]),col=rgb(0,0,0,0.2))
}

for (i in 1:10)
{
print(ecdf(S[,i])(Obj[i]))
}
#3 factors 

#going to go with 1 factor! 

rm(list=ls())

# hnv SE- factors-----------------------------------------------

load(file='data/phylofactor work spaces/hnv_workspace_no_sampling_effort')
jpeg(filename = 'figures/hnv no sampling effort factors.jpg')

plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,40))
for (rr in 1:140){
  lines(c(S[rr,1],diff(S[rr,])),col=rgb(0,0,0,0.2))
}
points(deviances,col='red',pch=16,cex=2)
dev.off()

plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,40))
for (rr in 1:140){
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

plot((Obj), type = 'o')
for (rr in 1:140){
  lines((S[rr,]),col=rgb(0,0,0,0.2))
}

for (i in 1:10)
{
  print(ecdf(S[,i])(Obj[i]))
}

rm(list=ls())

#going to say....1 factor!


# filo SE+ factors-----------------------------------------------

load(file='data/phylofactor work spaces/filo_workspace')
jpeg(filename = 'figures/filo with sampling effort factors.jpg')

plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,15))
for (rr in 1:140){
  lines(c(S[rr,1],diff(S[rr,])),col=rgb(0,0,0,0.2))
}
points(deviances,col='red',pch=16,cex=2)
dev.off()

plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,15))
for (rr in 1:140){
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

plot((Obj), type = 'o')
for (rr in 1:140){
  lines((S[rr,]),col=rgb(0,0,0,0.2))
}

for (i in 1:10)
{
  print(ecdf(S[,i])(Obj[i]))
}

rm(list=ls())

#1 factor?

# filo SE- factors-----------------------------------------------

rm(list=ls())

load(file='data/phylofactor work spaces/filo_workspace_no_sampling_effort')
jpeg(filename = 'figures/filo no sampling effort factors.jpg')

plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,25))
for (rr in 1:140){
  lines(c(S[rr,1],diff(S[rr,])),col=rgb(0,0,0,0.2))
}
points(deviances,col='red',pch=16,cex=2)
dev.off()

plot(c(Obj[1],diff(Obj)),type='o',lwd=2,cex=2,ylim=c(0,25))
for (rr in 1:140){
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

plot((Obj), type = 'o')
for (rr in 1:140){
  lines((S[rr,]),col=rgb(0,0,0,0.2))
}

for (i in 1:10)
{
  print(ecdf(S[,i])(Obj[i]))
}

rm(list=ls())

#looks like....0 factors?

# FIGURE MAKING ----

# hnv SE+ figure 1 ----

# 1 factor 

rm(list=ls())
load(file='data/phylofactor work spaces/hnv_workspace')
load("data/bat_taxonomy_data.Rdata")
source('R/taxonomy group name function.R')

taxonomy <- batphy1 %>%
  filter(hnv_samps > 0) %>%
  mutate(log_hnv_samps = log(hnv_samps)) %>%
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

B <- bins(pf$basis[,1:10])
B <- B[2:11]
Z <- Data$Z
probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

names.storage[[1]]

#colors 1 "Pteropodidae"  #CC00FFFF

colfcn <- function(n) return(c("#00EEEE"))

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
Legend$names <- names.storage[1:2]
P <- sapply(probs[1:2],FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend$names <- mapply(paste,Legend$names,P)
plot.new()
plot.new()
legend('topleft',legend=Legend$names,fill=Legend$colors,cex=1)

# hnv SE- figure 1 --------------------------------------------------------------------------------------

rm(list=ls())

load(file='data/phylofactor work spaces/hnv_workspace_no_sampling_effort')
load("data/bat_taxonomy_data.Rdata")
source('R/taxonomy group name function.R')

taxonomy <- batphy1 %>%
  filter(hnv_samps > 0) %>%
  mutate(log_hnv_samps = log(hnv_samps)) %>%
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

B <- bins(pf$basis[,1:10])
B <- B[2:11]
Z <- Data$Z
probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

#colors 1 "Pteropodidae"  #00FFFF
#colors 2 "Pteropodinae"  #7FFFD4

colfcn <- function(n) return(c("#7FFFD4"))
pf.tree <- pf.tree(pf, lwd=1, factors = 1, color.fcn=colfcn, branch.length = "none", bg.color = NA)

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue') 
  
ggsave("figures/hnv no sampling effort tree.png", bg = "transparent", height = 18, width = 18)

Legend <- pf.tree$legend
Legend$names <- names.storage[1]
P <- sapply(probs[1],FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend$names <- mapply(paste,Legend$names,P)
plot.new()
plot.new()
legend('topleft',legend=Legend$names,fill=Legend$colors,cex=1)

# filo SE + figure 1 --------------------------------------------------------------------------------------

#1 figure again 
rm(list=ls())

load(file='data/phylofactor work spaces/filo_workspace')
load("data/bat_taxonomy_data.Rdata")
source('R/taxonomy group name function.R')

taxonomy <- batphy1 %>%
  filter(filo_samps > 0) %>%
  mutate(log_filo_samps = log(filo_samps)) %>%
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

# B <- bins(pf$basis[,1:10])
# B <- B[2:11]
# Z <- Data$Z
# probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

#colors 1 "Pteropodidae"  #00FFFF
#colors 2 "Pteropodinae"  #7FFFD4
#colors 3 "Rhinolophoidea (Rhinolophidae & Hipposideridae)"  #FFB90F

colfcn <- function(n) return(c("#FFB90F"))

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
plot.new()
plot.new()
legend('topleft',legend=Legend$names,fill=Legend$colors,cex=1)

# filo SE- figure 1 --------------------------------------------------------------------------------------

#likely zero factors 

rm(list=ls())

load(file='data/phylofactor work spaces/filo_workspace_no_sampling_effort')
load("data/bat_taxonomy_data.Rdata")
source('R/taxonomy group name function.R')

taxonomy <- batphy1 %>%
  filter(filo_samps > 0) %>%
  mutate(log_hnv_samps = log(hnv_samps)) %>%
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

B <- bins(pf$basis[,1:10])
B <- B[2:11]
Z <- Data$Z
probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

#colors 1 "Pteropodidae"  #00FFFF
#colors 2 "Pteropodinae"  #7FFFD4

colfcn <- function(n) return(c("#7FFFD4"))
pf.tree <- pf.tree(pf, lwd=1, factors = 7, color.fcn=colfcn, branch.length = "none", bg.color = NA)

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_tippoint(size=10*Data$Z,col='blue') 

ggsave("figures/filo no sampling effort tree.png", bg = "transparent", height = 18, width = 18)

Legend <- pf.tree$legend
Legend$names <- names.storage[1]
P <- sapply(probs[1],FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend$names <- mapply(paste,Legend$names,P)
plot.new()
plot.new()
legend('topleft',legend=Legend$names,fill=Legend$colors,cex=1)
