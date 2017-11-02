setwd('C:/Users/Alculus/Documents/Bozeman/BZDEL/Data/Olival/')
library(mgcv)
library(stringi)
library(dplyr)
library(purrr)
library(ReporteRs)
library(visreg)
library(gridExtra)
library(ggpubr)
source("R/relative_contributions.R")
load('intermediates/phylofactor_comparisons')

Vmod <- vtraits$model[[1]]
Vmod_new <- vtraits_new$model[[1]]

summary(Vmod)
summary(Vmod_new)


############## Olival and Phylobin-removal datasets #############
D <- data_set_old
D$Envelope <- factor(D$Envelope)
D$Vector <- factor(D$Vector)
D$vCytoReplicTF <- factor(D$vCytoReplicTF)

D_new <-data_set_new
D_new$Envelope <- factor(D_new$Envelope)
D_new$Vector <- factor(D_new$Vector) 
D_new$vCytoReplicTF <- factor(D_new$vCytoReplicTF)

########### Olival gam ############
gam_old <- gam(formula=as.formula(Vmod$formula),data=D,family = binomial)
########### Olival formula, new dataset ##########
gam_compare <- gam(formula=as.formula(Vmod$formula),data=D_new,family = binomial)
tiff('~/Bozeman/BZDEL/Figures/Olival phylofactorization/Envelope_Vector_Comparison_w_and_without_phylobin_removal.tiff',height=1200,width=400)
  par(mfrow=c(3,2))
  visreg(gam_old,"Envelope",scale="response",main='Original',cex.main=2,cex.axis=2,cex.lab=2)
  visreg(gam_compare,"Envelope",scale='response',main='Post-Phylofactorization',cex.main=2,cex.axis=2,cex.lab=2,ylab=NA)
  visreg(gam_old,"Vector",scale='response',main='Original',cex.main=2,cex.axis=2,cex.lab=2)
  visreg(gam_compare,"Vector",scale='response',main='Post-Phylofactorization',cex.main=2,cex.axis=2,cex.lab=2,ylab=NA)
  visreg(gam_old,"cb_dist_noHoSa_maxLn",scale='response',main='Original',cex.main=2,cex.axis=2,cex.lab=2)
  visreg(gam_compare,"cb_dist_noHoSa_maxLn",scale='response',main='Post-Phylofactorization',cex.main=2,cex.axis=2,cex.lab=2,ylab=NA)
dev.off()

### GGplot will make a much more honest figure

clo <- function(x) x/sum(x)
ilogit <- function(x) 1/(1+exp(-x))

ix <- ! 1:nrow(D) %in% gam_old$na.action
preds <- ilogit(gam_old$linear.predictors-mean(gam_old$linear.predictors))
DF_old <- cbind(D[ix,c('Vector','Envelope','cb_dist_noHoSa_maxLn','IsZoonotic','IsZoonotic.stringent')],
                'version'=rep('Olival',sum(ix)),'Normalized_Prediction'=preds)

ix <-  ! 1:nrow(D_new) %in% gam_compare$na.action
preds <- ilogit(gam_compare$linear.predictors-mean(gam_compare$linear.predictors))
DF_new <- cbind(D_new[ix,c('Vector','Envelope','cb_dist_noHoSa_maxLn','IsZoonotic','IsZoonotic.stringent')],
              'version'=rep('Phylobin_removal',sum(ix)),'Normalized_Prediction'=preds)

DF <- rbind(DF_old,DF_new)
DF$version <- factor(DF$version)

ggV <- ggplot(DF,aes(x=Vector,y=Normalized_Prediction,fill=version))+
        geom_boxplot()+
        ggtitle(label='Vector-Borne')+
        facet_grid(.~version)

ggE <- ggplot(DF,aes(x=Envelope,y=Normalized_Prediction,fill=version))+
        geom_boxplot()+
        ggtitle(label='Enveloped')+
        facet_grid(.~version)

ggP <- ggplot(DF,aes(x=cb_dist_noHoSa_maxLn,y=Normalized_Prediction,fill=version,color=version))+
        geom_point()+
        scale_x_continuous(name='Non-Human Phylogenetic Breadth')+
        geom_smooth()+
        ggtitle(label='Phylogenetic Breadth')
        


# tiff('~/Bozeman/BZDEL/Figures/Olival phylofactorization/Envelope_Vector_Comparison_w_and_without_phylobin_removal.tiff',height=1200,width=400)       
# grid.arrange(ggV,ggE,ggP,layout.matrix=matrix(c(1,2,3),ncol=1))
ggarrange(ggV,ggE,ggP,ncol = 1,nrow=3)
ggsave('~/Bozeman/BZDEL/Figures/Olival phylofactorization/Envelope_Vector_Comparison_w_and_without_phylobin_removal.tiff',height=10,width=5)
ggarrange(ggV,ggE,ggP,ncol = 3,nrow=1)
ggsave('~/Bozeman/BZDEL/Figures/Olival phylofactorization/Envelope_Vector_Comparison_w_and_without_phylobin_removal_hamburger.tiff',height=4,width=14)
