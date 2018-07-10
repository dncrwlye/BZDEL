
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
  mutate(filo_surv_bin = ifelse(is.na(filo_surv), 0, filo_surv)) %>%
  mutate(filo_surv_bin = ifelse(filo_surv_bin == 1, 1, 0)) %>%
  mutate(log_filo_samps = log(filo_samps)) %>%
  select(c(filo_samps, filo_surv_bin, species, log_filo_samps))%>%
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

model.fcn <- function(formula,data,...){
  fit <- tryCatch(pscl::zeroinfl(formula,data,...),
                  error=function(e) NA)
  #fit <- do.call
  return(fit)
}

obj.fcn <- function(fit,grp,tree,PartitioningVariables,model.fcn,phyloData,...){
  if (!'zeroinfl' %in% class(fit))
  {
    return(0)
  }
  else 
  {
    fit2 <- zeroinfl(Z.poisson~1,data = fit$model,dist='negbin')
    # return(abs(summary(fit)$coefficients$zero['phyloS','z value']))
    fit$loglik-fit2$loglik %>% return()
  }
}

tictoc::tic()
pf <- gpf(Data,tree,Z.poisson~phylo,nfactors=2,algorithm = 'phylo',
          model.fcn = model.fcn,objective.fcn = obj.fcn,
          ncores = 1,
          dist = "negbin", cluster.depends='library(pscl)')
tictoc::toc()

tot.reps=50

source('R/filo and hnv positive factorization/null simulations script all bats try 3.R')

save(list=ls(),file='data/phylofactor work spaces/neg binom filo_workspace')

load(file='data/phylofactor work spaces/neg binom filo_workspace')

# extract deviances from actual data ----
summaries <- lapply(pf$models,summary)
loglik <- unlist(sapply(summaries, "[", "loglik"))

summaries.null <- (lapply(pf$models, "[[", "model"))
loglik.null <- lapply(summaries.null, 
                      function(data){
                        pscl::zeroinfl(Z.poisson~1, data)
                      })
loglik.null <- unlist(sapply(loglik.null, "[", "loglik"))

# plot null simulations against data ----

OBJ <- rbind(OBJ1[[1]], OBJ2[[1]])


plot(loglik-loglik.null, type ='l', col = 'red')
apply(OBJ, 1, lines)

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
