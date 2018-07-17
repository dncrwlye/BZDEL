
tryCatch(setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis"),
         error=function(e) setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/"))

library(phylofactor)
library(parallel)
library(tidyverse)
library(stringi)
require(pscl)
require(boot)
library(pscl)
load('data/phylofactor work spaces/bat_tree')
load("data/bat_taxonomy_data.Rdata")

Data <- batphy1 %>%
  mutate(hnv_surv_bin = ifelse(is.na(hnv_surv), 0, hnv_surv)) %>%
  mutate(hnv_surv_bin = ifelse(hnv_surv_bin == 1, 1, 0)) %>%
  mutate(log_hnv_samps = log(hnv_samps)) %>%
  select(c(hnv_samps, hnv_surv_bin, species, log_hnv_samps))%>%
  dplyr::rename(Species = species) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  mutate(Species = stri_trans_totitle(Species)) %>%
  mutate(Sample = 1) %>%
  unique()

Data[duplicated(Data$Species),]

tree<-bat_tree
tree <- ape::drop.tip(tree,tree$tip.label[!(tree$tip.label %in% Data$Species)])
Data <- Data[Data$Species %in% tree$tip.label,]

#rm(batphy1, bat_tree)

names(Data) <- c('Z.poisson',"Z.binom", 'Species', 'log_effort', 'Sample')

Data <- Data[Data$Species %in% tree$tip.label,]
n_factors = 10
Data[is.na(Data)] <- 0

source('R/filo and hnv positive factorization/negative binomial no zero inflation tests for all bats sampling effort.R')
source('R/filo and hnv positive factorization/null simulations script all bats try 4.R')

save(list=ls(),file='data/phylofactor work spaces/neg binom no zeroinfl hnv_workspace')

load(file='data/phylofactor work spaces/neg binom no zeroinfl hnv_workspace')
# plot null simulations against data ----


summaries.null <- (lapply(pf2$models, "[[", "model"))
loglik.null <- lapply(summaries.null, 
                      function(data){
                        MASS::glm.nb(Z.poisson~1, data)
                      })

logLik.models <- unlist(lapply(pf2$models,logLik))
loglik.null <- unlist(sapply(loglik.null, logLik))



plot(logLik.models-loglik.null, 
     type ='l', 
     col = 'red',
     xlim=c(0, 11), 
     ylim=c(0, 40))

apply(OBJ, 1, lines)

bb <- (t(diff(t(S))))
bb <- cbind(S[,1],bb)
Ojb.bb <- c(Obj[1],diff(Obj))

ll.model <- (logLik.models-loglik.null)
for (i in 1:10)
{
  print(ecdf(OBJ[,i])(ll.model[i]))
}

#print(ecdf(bb[,i])(Ojb.bb[i]))



pf.tree <- pf.tree(pf, lwd=1, branch.length = "none", bg.color = NA)
jj <- nrow(Data)
d <- data.frame(x=pf.tree$ggplot$data[1:jj,'x'],
                xend=pf.tree$ggplot$data[1:jj,'x'] + .2*(Data$Z.poisson),
                y=pf.tree$ggplot$data[1:jj,'y'],
                yend=pf.tree$ggplot$data[1:jj,'y'] )

pf.tree$ggplot +
  ggtree::theme_tree(bgcolor = NA, fgcolor = NA, plot.background = element_rect(fill = "transparent",colour = NA)) +
  #ggtree::geom_tippoint(size=10*Data$Z.poisson,col='blue') +
  geom_segment(data= d,aes(x=x,y=y,xend=xend,yend=yend, size= Data$Z.binom, colour = 'blue'))
