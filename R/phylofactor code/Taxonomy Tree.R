library(ape)
library(tidyverse)
library(readr)

########################################
########################################
viruses <- read_csv("~/Desktop/BDEL/BZDEL/viruses.csv")

taxonomy_tree <- function(order)
{
Mononegavirales<-viruses %>% #I did Monogavirales first, then turned this into a function. Just ignore names meaning, its just a variable/dataset name now
  filter(vOrder==order)
virus_names <- unique(Mononegavirales$vVirusNameCorrected) #create a vector that contains all the unique viruses in the order you are curious about
matrix_mono <- matrix(100, length(virus_names), length(virus_names)) #create your distance matrix
colnames(matrix_mono) <- virus_names #label your distance matrix 
row.names(matrix_mono) <- virus_names
diag(matrix_mono) <- 0 #set the diaganol of the distance matrix to zero, this is required for Hclust package

for (i in 1:length(virus_names)) 
{
  indice1 <- which(Mononegavirales == virus_names[i])[[1]] # I'm not sure anymore if we need this
  
  family <- as.character(Mononegavirales[i, 'vFamily']) # pull the family name 
  if(family != 'Unassigned') # several families are unassigned, we need to not group those viruses together 
  {
  family_vector <- which(Mononegavirales$vFamily == as.character(family)) #create a vector of indices that indiciate which viruses are in the same family
  matrix_mono[i,c(family_vector)] = 75 #assign those locations on the matrix a 75
  matrix_mono[c(family_vector),i] = 75
  }
  subfamily <- as.character(Mononegavirales[i, 'vSubfamily'])
  #several viruses don't have subfamilies, need to exclude these
  if (!is.na(subfamily))
    {
    subfamily_vector <- which(Mononegavirales$vSubfamily == as.character(subfamily))
    matrix_mono[i,c(subfamily_vector)] = 50 #the viruses are more closely related, so assign them a 50
    matrix_mono[c(subfamily_vector),i] = 50
    }
  #okay this should actually go last

  order <- as.character(Mononegavirales[i, 'vGenus'])
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
plot(Herpesvirales_hc_phylo)

Nidovirales_hc_phylo<-taxonomy_tree('Nidovirales')
plot(Nidovirales_hc_phylo)

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
  

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#############  ~~**SUPER ALT. FUNCTION WHERE  WE COMBINE EVERYTHING**~~  #################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


taxonomy_tree_super_alt <- function()
{

  virus_names <- unique(viruses$vVirusNameCorrected)
  
  matrix_mono <- matrix(500, length(virus_names), length(virus_names)) 
  colnames(matrix_mono) <- virus_names
  row.names(matrix_mono) <- virus_names
  diag(matrix_mono) <- 0
  
  for (i in 1:length(virus_names)) 
  {
    indice1 <- which(viruses == virus_names[i])[[1]] 
    
    order <- as.character(viruses[i, 'vOrder'])
    if(order != 'Unassigned')    
    {
      order_vector <- which(viruses$vOrder == as.character(order))
      matrix_mono[i,c(order_vector)] = 250
      matrix_mono[c(order_vector),i] = 250
    }
    
    family <- as.character(viruses[i, 'vFamily'])
    if(family != 'Unassigned')    
    {
      family_vector <- which(viruses$vFamily == as.character(family))
      matrix_mono[i,c(family_vector)] = 125
      matrix_mono[c(family_vector),i] = 125
    }
    
    subfamily <- as.character(viruses[i, 'vSubfamily'])
    #several viruses don't have subfamilies, need to exclude these
    if (!is.na(subfamily))
    {
      subfamily_vector <- which(viruses$vSubfamily == as.character(subfamily))
      matrix_mono[i,c(subfamily_vector)] = 75
      matrix_mono[c(subfamily_vector),i] = 75
    }
    
    #okay this should actually go last
    genus <- as.character(viruses[i, 'vGenus'])
    if(genus != 'Unassigned')
    {
      genus_vector <- which(viruses$vGenus == as.character(genus))
      matrix_mono[i,c(genus_vector)] = 35
      matrix_mono[c(genus_vector),i] = 35
    }
  }
  
  diag(matrix_mono) <- 0
  
  x<-as.dist(matrix_mono)
  hc <- hclust(x, "ave")
  hc_phylo <- as.phylo(hc)
}

el_jefe_grande <-taxonomy_tree_super_alt()
plot(el_jefe_grande)

save(el_jefe_grande,  
     file = 'el_jefe.rdata')


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#############  ~~**SUPER ALT. FUNCTION WHERE  WE COMBINE EVERYTHING**~~  #################
####################   BUT THIS TIME FOR THE REDUCED DATASET  ############################
##########################################################################################
##########################################################################################
##########################################################################################


taxonomy_tree_super_alt <- function()
{
  #virus_names <- unique(olival_et_al_2017_reduced$virus_name_dc)
  virus_names <- (olival_et_al_2017_reduced$virus_name_dc)
  
  
  matrix_mono <- matrix(500, length(virus_names), length(virus_names)) 
  colnames(matrix_mono) <- virus_names
  row.names(matrix_mono) <- virus_names
  diag(matrix_mono) <- 0
  
  for (i in 1:length(virus_names)) 
  {
    indice1 <- which(olival_et_al_2017_reduced == virus_names[i])[[1]] 
    
    order <- as.character(olival_et_al_2017_reduced[i, 'vOrder'])
    if(order != 'Unassigned')    
    {
      order_vector <- which(olival_et_al_2017_reduced$vOrder == as.character(order))
      matrix_mono[i,c(order_vector)] = 250
      matrix_mono[c(order_vector),i] = 250
    }
    
    family <- as.character(olival_et_al_2017_reduced[i, 'vFamily'])
    if(family != 'Unassigned')    
    {
      family_vector <- which(olival_et_al_2017_reduced$vFamily == as.character(family))
      matrix_mono[i,c(family_vector)] = 125
      matrix_mono[c(family_vector),i] = 125
    }
    
    subfamily <- as.character(olival_et_al_2017_reduced[i, 'vSubfamily'])
    #several olival_et_al_2017_reduced don't have subfamilies, need to exclude these
    if (!is.na(subfamily))
    {
      subfamily_vector <- which(olival_et_al_2017_reduced$vSubfamily == as.character(subfamily))
      matrix_mono[i,c(subfamily_vector)] = 75
      matrix_mono[c(subfamily_vector),i] = 75
    }
    
    #okay this should actually go last
    genus <- as.character(olival_et_al_2017_reduced[i, 'vGenus'])
    if(genus != 'Unassigned')
    {
      genus_vector <- which(olival_et_al_2017_reduced$vGenus == as.character(genus))
      matrix_mono[i,c(genus_vector)] = 35
      matrix_mono[c(genus_vector),i] = 35
    }
  }
  
  diag(matrix_mono) <- 0
  
  x<-as.dist(matrix_mono)
  hc <- hclust(x, "ave")
  hc_phylo <- as.phylo(hc)
}

taxonomy_tree_olival_reduced_dataset <-taxonomy_tree_super_alt()
plot(taxonomy_tree_olival_reduced_dataset)

save(taxonomy_tree_olival_reduced_dataset,file = '/Users/buckcrowley/Desktop/BDEL/BZDEL/Data/Meta Analyses Papers/el_jefe_reduced.rdata')
