# 
# internal.edge.color = function(tree, data)
# {
#   nodes <-data.frame(unique(tree$edge[,1]))
#   nodes$color <- NA
# 
#   for (i in 1:nrow(nodes))
#   {
#   data.kerala <- data[data$Kerala==1,]
#     
#   tip.labels <- ape::extract.clade(tree, nodes[i,1])$tip.label
#   x <- tip.labels %in% data.kerala$Species 
#   if (TRUE %in% x)
#   {
#     nodes[i,2] <- 'black'
#   }
#   if (!(TRUE %in% x))
#   {
#     nodes[i,2] <- 'grey'
#   }
#   }
#   return(nodes)
# }

internal.edge.color = function(tree, data)
{
  ggtree.object <- ggtree::ggtree(tree) 
  nodes <-data.frame(unique(ggtree.object$data))
  nodes <- data.frame(nodes[c(nodes$isTip==FALSE),])
  nodes <- data.frame(nodes$node)
  nodes$color <- NA
  
  for (i in 1:nrow(nodes))
  {
    data.kerala <- data[data$Kerala==1,]
    
    tip.labels <-   species.list <- ggtree::get.offspring.tip(tree, node=nodes[i,1])
    x <- tip.labels %in% data.kerala$Species 
    if (TRUE %in% x)
    {
      nodes[i,2] <- 'black'
    }
    if (!(TRUE %in% x))
    {
      nodes[i,2] <- 'grey'
    }
  }
  
  
  return(nodes)
}





