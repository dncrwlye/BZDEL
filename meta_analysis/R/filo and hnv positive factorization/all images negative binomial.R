#negative binomial regression images 

tryCatch(setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis"),
         error=function(e) setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/"))

# FILO VIRUS IMAGES  ----

load(file='data/phylofactor work spaces/neg binom no zeroinfl filo_workspace_try2')

# null sim ----
# i now have two types of null simulations, the re shuffling method and 
# number 2, where i made new data from a distribution. lets see if it looks different!


# null sim 2 ----
# i now have two types of null simulations, the re shuffling method and 
# number 2, where i made new data from a distribution. lets see if it looks different!
x <- unlist(lapply(pf2$models,'[[', 'null.deviance')) - unlist(lapply(pf2$models,'[[', 'deviance'))

plot(x, type = 'l', col ='red')

OBJ1 <- OBJ1[[1]]
apply(unlist(OBJ1), 1, lines)

for (i in 1:10)
{
  print(ecdf(OBJ1[,i])(x[[i]]))
}

#.................Now we can try the method where we generated new data....

plot(x, 
     type ='l', 
     col = 'red',
     xlim=c(0, 11), 
     ylim=c(0, 150))

apply(unlist(OBJ.final.2), 1, lines)

for (i in 1:10)
{
  print(ecdf(OBJ.final.2[,i])(x[[i]]))
}

# factor names -----

taxonomy <- batphy1 %>%
  select(c(species, tax)) 
source('R/taxonomy group name function.R')

names.storage <- list()

for (i in 1:10)
{
  indexes = pf2$groups[[i]][[1]]
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
  
  #print(group_taxonomy_list %>% table())

}

#factor names 2-
#usefull little thing
#group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; Yinpterochiroptera; Pteropodidae; ", 

factor.1 <- pf2$groups[[1]][[1]]
factor.2 <- pf2$groups[[2]][[1]]
factor.3 <- pf2$groups[[3]][[1]]
paraphyletic.remainder <- pf2$groups[[1]][[2]][!pf2$groups[[1]][[2]] %in% factor.2 & !pf2$groups[[1]][[2]] %in% factor.3]

paraphyletic.remainder %in% factor.1
factor.1 %in% paraphyletic.remainder
paraphyletic.remainder %in% factor.2
factor.2 %in% paraphyletic.remainder
paraphyletic.remainder %in% factor.3
factor.3 %in% paraphyletic.remainder

length(factor.1) +
length(factor.2) +
length(factor.3) + 
length(paraphyletic.remainder)

Z <- Data$Z.poisson

factor.1.p <- sapply(factor.1,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.2.p <- sapply(factor.2,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.3.p <- sapply(factor.3,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
paraphyletic.remainder.p <- sapply(paraphyletic.remainder,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

f1.sps <- sum(factor.1.p) / length(factor.1.p)
f2.sps <- sum(factor.2.p) / length(factor.2.p)
f3.sps <- sum(factor.3.p) / length(factor.3.p)
pfr.sps <- sum(paraphyletic.remainder.p)/ length(paraphyletic.remainder.p)


f1.sps
f2.sps
f3.sps
pfr.sps

#so the three groups are the Pteropodidae, Noctilionoidea, and the Vespertilionoidea

pf.tree <- pf.tree(pf2, factors = 1:3, lwd=1, branch.length = "none", bg.color = NA)

# data for figure 

jj <- nrow(Data)
d <- data.frame(x=pf.tree$ggplot$data[1:jj,'x'],
                xend=pf.tree$ggplot$data[1:jj,'x'] + .2*(Data$Z.poisson),
                y=pf.tree$ggplot$data[1:jj,'y'],
                yend=pf.tree$ggplot$data[1:jj,'y'] )

# figure ----

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  ggtree::geom_cladelabel(node=1238, label="Emballonuroidea",
                           color="black", angle=35.98338, offset=4.5, offset.text = 2, align=TRUE) + 
  ggtree::geom_cladelabel(node=1085, label="Yangochiroptera",
                           color="black", angle=147.5069, offset=2.5, offset.text = 2, align=TRUE) + 
  geom_segment(data= d,aes(x=x,y=y,xend=xend,yend=yend, size= Data$Z.binom, colour = 'blue')) +
  ggtree::geom_cladelabel(node=1592, label="Rhinolophoidea",
                          color="black", angle=-62.40997, offset=3.5, offset.text = 2, align=TRUE) + 
  ggtree::geom_cladelabel(node=1700, label="Pteropodidae", 
                          color="black", angle=-5.069252, offset=3, offset.text = 2, align=TRUE)  +
  geom_segment(data= d,aes(x=x,y=y,xend=xend,yend=yend, size= Data$Z.binom, colour = 'blue'))

ggsave("figures/filo sampling effort tree neg binom.png", bg = "transparent", height = 18, width = 18)

rm(list = ls())

# HENIPAVIRUS ----

load(file='data/phylofactor work spaces/neg binom no zeroinfl hnv_workspace try 2')
# i now have two types of null simulations, the re shuffling method and 
# number 2, where i made new data from a distribution. lets see if it looks different!
#null sim  ----
x <- unlist(lapply(pf2$models,'[[', 'null.deviance')) - unlist(lapply(pf2$models,'[[', 'deviance'))
plot(x, type = 'l', col ='red')

OBJ1 <- OBJ1[[1]]
apply(unlist(OBJ1), 1, lines)
for (i in 1:10)
{
  print(ecdf(OBJ1[,i])(x[[i]]))
}

# null sim 2 ----

plot(x, 
     type ='l', 
     col = 'red',
     xlim=c(0, 11), 
     ylim=c(0, 150))

apply(unlist(OBJ.final.2), 1, lines)

for (i in 1:10)
{
  print(ecdf(OBJ.final.2[,i])(x[[i]]))
}


taxonomy <- batphy1 %>%
  select(c(species, tax)) 
source('R/taxonomy group name function.R')

# factor names  ---- 
names.storage <- list()

for (i in 1:10)
{
  indexes = pf2$groups[[i]][[1]]
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
  # print("NEW GROUP
  #       NEW GROUP
  #       NEW GROUP
  #       NEW GROUP
  #       NEW GROUP
  #       NEW GROUP
  #       NEW GROUP
  #       NEW GROUP
  #       NEW GROUP
  #       NEW GROUP")
}

names.storage

#factor names 2----
#usefull little thing
#group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; Yinpterochiroptera; Pteropodidae; ", 

factor.1 <- pf2$groups[[1]][[1]]
factor.2 <- pf2$groups[[2]][[1]]
paraphyletic.remainder <- pf2$groups[[2]][[2]][!pf2$groups[[2]][[2]] %in% factor.2]

paraphyletic.remainder %in% factor.1
factor.1 %in% paraphyletic.remainder
paraphyletic.remainder %in% factor.2
factor.2 %in% paraphyletic.remainder

length(factor.1) +
  length(factor.2) +
  length(paraphyletic.remainder)

Z <- Data$Z.poisson

factor.1.p <- sapply(factor.1,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
factor.2.p <- sapply(factor.2,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
paraphyletic.remainder.p <- sapply(paraphyletic.remainder,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)

f1.sps <- sum(factor.1.p) / length(factor.1.p)
f2.sps <- sum(factor.2.p) / length(factor.2.p)
pfr.sps <- sum(paraphyletic.remainder.p)/ length(paraphyletic.remainder.p)


f1.sps
f2.sps
f3.sps
pfr.sps

#print(ecdf(bb[,i])(Ojb.bb[i]))

#so the two groups are Pteropus and Yinpterochiroptera

# figure ----
pf.tree <- pf.tree(pf2, factors = 1:2, lwd=1, branch.length = "none", bg.color = NA)
# data for figure 
jj <- nrow(Data)
d <- data.frame(x=pf.tree$ggplot$data[1:jj,'x'],
                xend=pf.tree$ggplot$data[1:jj,'x'] + .2*(Data$Z.poisson),
                y=pf.tree$ggplot$data[1:jj,'y'],
                yend=pf.tree$ggplot$data[1:jj,'y'] )
pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  #ggtree::geom_tippoint(size=10*Data$Z.poisson,col='blue') +
  ggtree::geom_cladelabel(node=1700, label="Pteropodidae", 
                          color="black", angle=-5.069252, offset=2.5, offset.text = 2, align=TRUE) +
  ggtree::geom_cladelabel(node=1085, label="Yangochiroptera",
                          color="black", angle=147.5069, offset=1.2, offset.text = 2, align=TRUE) +
  ggtree::geom_cladelabel(node=1592, label="Rhinolophoidea",
                          color="black", angle=-62.40997, offset=3.5, offset.text = 2, align=TRUE) + 
  geom_segment(data= d,aes(x=x,y=y,xend=xend,yend=yend, size= Data$Z.binom, colour = 'blue')) 

ggsave("figures/hnv sampling effort tree neg binom.png", bg = "transparent", height = 18, width = 18)

rm(list = ls())

# node and angle data ----

# .......... Yinpterochiroptera (Rhinolophoidea & Pteropodidae) bats......................................

species.list <- ggtree::get.offspring.tip(pf2$tree, node=1591)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            group_taxonomy_list[,1]) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(pf2$tree, node=1591)
mean(pf.tree$ggplot$data[pf.tree$ggplot$data$label %in% species.list,]$angle-90)

# .......................... Vespertilionoidea bats......................................

species.list <- ggtree::get.offspring.tip(pf2$tree, node=1283)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            group_taxonomy_list[,1]) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(pf2$tree, node=1283)
mean(pf.tree$ggplot$data[pf.tree$ggplot$data$label %in% species.list,]$angle-90)

# .......................... Emballonuroidea bats......................................

species.list <- ggtree::get.offspring.tip(pf2$tree, node=1238)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            group_taxonomy_list[,1]) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(pf2$tree, node=1238)
mean(pf.tree$ggplot$data[pf.tree$ggplot$data$label %in% species.list,]$angle-90)

# .......... Pteropodidae bats......................................

species.list <- ggtree::get.offspring.tip(pf2$tree, node=1700)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            group_taxonomy_list[,1]) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(pf2$tree, node=1700)
mean(pf.tree$ggplot$data[pf.tree$ggplot$data$label %in% species.list,]$angle-90)

# .......................... Yangochiroptera bats......................................

species.list <- ggtree::get.offspring.tip(pf2$tree, node=1085)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            group_taxonomy_list[,1]) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(pf2$tree, node=1085)
mean(pf.tree$ggplot$data[pf.tree$ggplot$data$label %in% species.list,]$angle-90)

# .......................... Rhinolophoidea bats......................................

species.list <- ggtree::get.offspring.tip(pf2$tree, node=1592)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            group_taxonomy_list[,1]) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(pf2$tree, node=1592)
mean(pf.tree$ggplot$data[pf.tree$ggplot$data$label %in% species.list,]$angle-90)

