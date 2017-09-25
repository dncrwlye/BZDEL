library(phytools)
library(tidyverse)
load("~/Desktop/BDEL/BZDEL/phylo_data_for_alex.rdata")
load("~/Desktop/BDEL/BZDEL/Olival_et_al_2017.Rdata")

#plotTree(Mononegavirales_hc_phylo,ftype="i")

viruses_herpes<- viruses %>%
  filter(vOrder == "Mononegavirales") 

row.names(viruses_herpes) <- viruses_herpes$vVirusNameCorrected 

viruses_herpes <- viruses_herpes %>%
  dplyr::select(c(IsZoonotic)) %>%
  mutate(IsZoonotic = as.factor(IsZoonotic))

svl<-as.matrix(viruses_herpes)[,1]

plotTree(Mononegavirales_hc_phylo,ftype="i")

cols<-setNames(palette()[1:length(unique(svl))],sort(unique(svl)))

tiplabels(pie=to.matrix(svl,sort(unique(svl))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(Herpesvirales_hc_phylo)),fsize=0.8)

########################



viruses_herpes<-viruses %>%
  filter(vOrder == "Nidovirales") 

row.names(viruses_herpes) <- viruses_herpes$vVirusNameCorrected 

viruses_herpes <- viruses_herpes %>%
  dplyr::select(c(IsZoonotic)) %>%
  mutate(IsZoonotic = as.factor(IsZoonotic))


svl<-as.matrix(viruses_herpes)[,1]

#plotTree(Mononegavirales_hc_phylo,ftype='i',type="fan",fsize=.3)

plotTree(Nidovirales_hc_phylo,fsize=.5)

#fitER<-ace(svl,Mononegavirales_hc_phylo,model="ER",type="discrete")
#fitER

dst<-Nidovirales_hc_phylo
#dst$edge.length[dst$edge.length==0]<-max(nodeHeights(Mononegavirales_hc_phylo))*1e-6
dst$edge.length[dst$edge.length==0]<-.01

fit3<-ace(svl,dst,type="discrete",model="ER")
fit3

cols<-setNames(palette()[1:length(unique(svl))],sort(unique(svl)))

tiplabels(pie=to.matrix(svl,sort(unique(svl))),piecol=cols,cex=0.25)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(Nidovirales_hc_phylo)),fsize=0.8)

nodelabels(pie=fit3$lik.anc, cex=.4)



######################################### 
######################################### 
######################################### 
######################################### 
######################################### 
######################################### 
######################################### 
######################################### 
######################################### 
######################################### 
######################################### 

virus_group='Mononegavirales'
phylo_group = Mononegavirales_hc_phylo

viruses_group_data<-olival_et_al_2017 %>%
    filter(vOrder == virus_group) 

row.names(viruses_group_data) <- viruses_group_data$vVirusNameCorrected 
  
viruses_group_data <- viruses_group_data %>%
    dplyr::select(c(IsZoonotic)) %>%
    mutate(IsZoonotic = as.factor(IsZoonotic))
  
svl<-as.matrix(viruses_group_data)[,1]
  
phylo_group_polytomies_resolved<- di2multi(phylo_group, 0.03)
  
#dst<-phylo_group_polytomies_resolved
#dst$edge.length[dst$edge.length==0]<-.01
 

fit3<-ace(svl,phylo_group_polytomies_resolved,type="discrete",model="ER")
fit3
fit4<-rerootingMethod(phylo_group_polytomies_resolved,svl,model="SYM")

tree<-plotTree(phylo_group_polytomies_resolved,fsize=.5, direction="leftwards")
  
cols<-setNames(palette()[1:length(unique(svl))],sort(unique(svl)))
  
tiplabels(pie=to.matrix(svl,sort(unique(svl))),piecol=cols,cex=0.25)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(phylo_group_polytomies_resolved)),fsize=0.8)
nodelabels(pie=fit3$lik.anc, cex=.4)
  
######################################### 
######################################### 
######################################### 
######################################### 
######################################### 
######################################### 
######################################### 
######################################### 
######################################### 
######################################### 
#########################################


