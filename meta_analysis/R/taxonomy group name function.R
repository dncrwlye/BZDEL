taxonomy_group_name <- function(taxonomy_list)
{
  colnames(taxonomy_list) <- 'tax'
  taxonomy_list <- taxonomy_list %>%
    dplyr::filter(!(is.na(tax))) 
  
  if (nrow(taxonomy_list) !=1)
  {
  taxonomy_list <- taxonomy_list %>%
    mutate(tax = as.character(tax)) %>%
    dplyr::filter(nchar(tax)>145) %>%
    dplyr::filter(grepl('Animalia', tax))
  }
   
    if (nrow(taxonomy_list) ==1)
    { 
    colnames(taxonomy_list) <- 'tax'
    taxonomy_list <- taxonomy_list %>%
      mutate(tax = as.character(tax)) 
    str<-as.character(taxonomy_list[1,1])
    x <- unlist(strsplit(str, " "))
    integer <- length(x)
    return(paste(x[c(integer-1, integer)], sep = " ", collapse = " "))
    }
  
  str<-as.character(taxonomy_list[1,1])
  x <- unlist(strsplit(str, " "))
    for (zed in 1:nrow(taxonomy_list))
    #for (zed in 1:311)
      {
       stopwords = as.character(taxonomy_list[zed,1])
       stopwords <- unlist(strsplit(stopwords, " "))
       x <- x[x %in% stopwords]
       integer <- length(x)
       #print(zed)
       #print(x)
      }
  return(x[integer])
}
#taxonomy_list[311,1]
#taxonomy_list <- group_taxonomy_list
#taxonomy_group_name(taxonomy_africa_list)

