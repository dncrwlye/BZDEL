library(ape)
library(tidyverse)

########################################
########################################
library(readr)
viruses <- read_csv("~/Desktop/BDEL/BZDEL/viruses.csv")

#first do Mononegavirales

taxonomy_tree <- function(order)
{
Mononegavirales<-viruses %>%
  filter(vOrder==order)
virus_names <- unique(Mononegavirales$vVirusNameCorrected)
matrix_mono <- matrix(100, length(virus_names), length(virus_names)) 
colnames(matrix_mono) <- virus_names
row.names(matrix_mono) <- virus_names
diag(matrix_mono) <- 0

for (i in 1:length(virus_names)) 
{
  indice1 <- which(Mononegavirales == virus_names[i])[[1]] 
  
  family <- as.character(Mononegavirales[i, 'vFamily'])
  family_vector <- which(Mononegavirales$vFamily == as.character(family))
  matrix_mono[i,c(family_vector)] = 75
  matrix_mono[c(family_vector),i] = 75
  
  subfamily <- as.character(Mononegavirales[i, 'vSubfamily'])
  #several viruses don't have subfamilies, need to exclude these
  if (!is.na(subfamily))
    {
    subfamily_vector <- which(Mononegavirales$vSubfamily == as.character(subfamily))
    matrix_mono[i,c(subfamily_vector)] = 50
    matrix_mono[c(subfamily_vector),i] = 50
    }
  #okay this should actually go last

  order <- as.character(viruses_no_order_subset[i, 'vGenus'])
  
  if(order != 'Unassigned')
  {
    order <- as.character(Mononegavirales[i, 'vGenus'])
  order_vector <- which(Mononegavirales$vGenus == as.character(order))
  matrix_mono[i,c(order_vector)] = 25
  matrix_mono[c(order_vector),i] = 25
  }
}

diag(matrix_mono) <- 0

x<-as.dist(matrix_mono)
hc <- hclust(x, "ave")
hc_phylo <- as.phylo(hc)
}


Mononegavirales_hc_phylo<-taxonomy_tree('Mononegavirales')
plot(Mononegavirales_hc_phylo)

Picornavirales_hc_phylo<-taxonomy_tree('Picornavirales')
plot(Picornavirales_hc_phylo)

Herpesvirales_hc_phylo<-taxonomy_tree('Herpesvirales')
plot(Picornavirales_hc_phylo)

Nidovirales_hc_phylo<-taxonomy_tree('Nidovirales')
plot(Picornavirales_hc_phylo)

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#############  ~~**ALTERNATIVE FUNCTION FOR THOSE VIRUSES W/O ORDER**~~  #################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

viruses_no_order <- viruses %>%
  filter(vOrder == 'Unassigned')

taxonomy_tree_alternative <- function(family)
{
  viruses_no_order_subset<-viruses_no_order %>%
    filter(vFamily==family)
  
  virus_names <- unique(viruses_no_order_subset$vVirusNameCorrected)
  matrix_mono <- matrix(100, length(virus_names), length(virus_names)) 
  colnames(matrix_mono) <- virus_names
  row.names(matrix_mono) <- virus_names
  diag(matrix_mono) <- 0
  
  for (i in 1:length(virus_names)) 
  {
    indice1 <- which(viruses_no_order_subset == virus_names[i])[[1]] 
    
    family <- as.character(viruses_no_order_subset[i, 'vFamily'])
    family_vector <- which(viruses_no_order_subset$vFamily == as.character(family))
    matrix_mono[i,c(family_vector)] = 75
    matrix_mono[c(family_vector),i] = 75
    
    subfamily <- as.character(viruses_no_order_subset[i, 'vSubfamily'])
    #several viruses don't have subfamilies, need to exclude these
    if (!is.na(subfamily))
    {
      subfamily_vector <- which(viruses_no_order_subset$vSubfamily == as.character(subfamily))
      matrix_mono[i,c(subfamily_vector)] = 50
      matrix_mono[c(subfamily_vector),i] = 50
    }
    
    #okay this should actually go last

    order <- as.character(viruses_no_order_subset[i, 'vGenus'])
    
    if(order != 'Unassigned')
    {
    order_vector <- which(viruses_no_order_subset$vGenus == as.character(order))
    matrix_mono[i,c(order_vector)] = 25
    matrix_mono[c(order_vector),i] = 25
    }
  }
  
  diag(matrix_mono) <- 0
  
  x<-as.dist(matrix_mono)
  hc <- hclust(x, "ave")
  hc_phylo <- as.phylo(hc)
}

Parvoviridae_phylo <- taxonomy_tree_alternative('Parvoviridae')
Polyomaviridae_phylo <- taxonomy_tree_alternative("Polyomaviridae")
Retroviridae_phylo <- taxonomy_tree_alternative('Retroviridae')
Reoviridae_phlyo <- taxonomy_tree_alternative('Reoviridae')
Asfarviridae_phylo <- taxonomy_tree_alternative('Asfarviridae')
Bunyaviridae_phylo <- taxonomy_tree_alternative('Bunyaviridae')
Arenaviridae_phylo <- taxonomy_tree_alternative('Arenaviridae')
Flaviviridae_phylo <- taxonomy_tree_alternative('Flaviviridae')
Hepadnaviridae_phylo <- taxonomy_tree_alternative('Hepadnaviridae')
Togaviridae_phylo <- taxonomy_tree_alternative('Togaviridae')
Astroviridae_phylo <- taxonomy_tree_alternative('Astroviridae')
Caliciviridae_phylo <- taxonomy_tree_alternative('Caliciviridae')
Papillomaviridae_phylo <- taxonomy_tree_alternative('Papillomaviridae')
Poxviridae_phylo <- taxonomy_tree_alternative('Poxviridae')
Orthomyxoviridae_phylo <- taxonomy_tree_alternative('Orthomyxoviridae')
Hepeviridae_phylo <- taxonomy_tree_alternative('Hepeviridae')
Picobirnaviridae_phylo <- taxonomy_tree_alternative('Picobirnaviridae')
Rhabdoviridae_phylo <- taxonomy_tree_alternative('Rhabdoviridae')
Circoviridae_phylo <- taxonomy_tree_alternative('Circoviridae')
Herpesviridae_phylo <- taxonomy_tree_alternative('Herpesviridae')
Anelloviridae_phylo <- taxonomy_tree_alternative('Anelloviridae')

#note: Asfarviridae_phylo, Hepeviridae_phylo, Picobirnaviridae_phylo, and Anelloviridae_phylo did not work 
#as they only had one virus in their family 
#also note there are herpes viruses in here with no order

x<-filter(viruses_no_order, vFamily == 'Herpesviridae')

save(Mononegavirales_hc_phylo,
     Picornavirales_hc_phylo,
     Herpesvirales_hc_phylo,
     Nidovirales_hc_phylo,
     Parvoviridae_phylo,
     Polyomaviridae_phylo,  
     Retroviridae_phylo, 
     Reoviridae_phlyo,  
     #Asfarviridae_phylo, 
     Bunyaviridae_phylo,  
     Arenaviridae_phylo,  
     Flaviviridae_phylo, 
     Hepadnaviridae_phylo,  
     Togaviridae_phylo,  
     Astroviridae_phylo,  
     Caliciviridae_phylo, 
     Papillomaviridae_phylo, 
     Poxviridae_phylo,  
     Orthomyxoviridae_phylo,  
     #Hepeviridae_phylo,  
     #Picobirnaviridae_phylo,  
     Rhabdoviridae_phylo, 
     Circoviridae_phylo,  
     Herpesviridae_phylo, 
     #Anelloviridae_phylo,  
     file = 'phylo_data_for_alex.rdata')
  

