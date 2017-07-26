#devtools::install_github('reptalex/phylofactor')
#install.packages('caper')
library(phylofactor)
??phylofactor

load("~/Desktop/BDEL/BZDEL/Olival_et_al_2017.Rdata")
load("~/Desktop/BDEL/BZDEL/phylo_data_for_alex.rdata")

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
data <- as.matrix(phylo_group_polytomies_resolved$tip.label)



#phylofactor::getGroups(phylo_group_polytomies_resolved)

PhyloFactor(Data= data,tree=phylo_group_polytomies_resolved, X= as.factor(svl))
