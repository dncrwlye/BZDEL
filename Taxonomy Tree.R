library(ape)

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
  indice1 <- which(Mononegavirales == virus_names[i])[[1]] 
  order <- as.character(Mononegavirales[i, 'vGenus'])
  order_vector <- which(Mononegavirales$vGenus == as.character(order))
  matrix_mono[i,c(order_vector)] = 25
  matrix_mono[c(order_vector),i] = 25
  
}

diag(matrix_mono) <- 0

x<-as.dist(matrix_mono)
hc <- hclust(x, "ave")
hc.phylo <- as.ph


plot(hc)
