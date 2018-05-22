species.list <- ggtree::get.offspring.tip(pf$tree, node=1591)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            group_taxonomy_list[,1]) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

#1591 Yinpterochiroptera 

species.list <- ggtree::get.offspring.tip(pf$tree, node=1700)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            group_taxonomy_list[,1]) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(pf$tree, node=1700)
abs(mean(pf.tree$ggplot$data[pf.tree$ggplot$data$label %in% species.list,]$angle+90))

# 1700 Pteropodidae

species.list <- ggtree::get.offspring.tip(pf$tree, node=1592)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            group_taxonomy_list[,1]) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(pf$tree, node=1592)
abs(mean(pf.tree$ggplot$data[pf.tree$ggplot$data$label %in% species.list,]$angle-90))


#1592 Rhinolophoidea

species.list <- ggtree::get.offspring.tip(pf$tree, node=1283)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            group_taxonomy_list[,1]) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(pf$tree, node=1283)
abs(mean(pf.tree$ggplot$data[pf.tree$ggplot$data$label %in% species.list,]$angle-90))


#1283 Vespertilionoidea

species.list <- ggtree::get.offspring.tip(pf$tree, node=1238)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            group_taxonomy_list[,1]) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(pf$tree, node=1238)

mean(pf.tree$ggplot$data[pf.tree$ggplot$data$label %in% species.list,]$angle-90)

#1238 Emballonuroidea

species.list <- ggtree::get.offspring.tip(pf$tree, node=1087)
species.list <- tolower(gsub("_", " ",species.list))
group_taxonomy_list <- as.data.frame(taxonomy[match(species.list,taxonomy[,1]),2])
gsub("\\)|;","", as.character(taxonomy_group_name(group_taxonomy_list)))
group_taxonomy_list <- gsub(pattern = "Animalia; Bilateria; Deuterostomia; Chordata; Vertebrata; Gnathostomata; Tetrapoda; Mammalia; Theria; Eutheria; Chiroptera; ", 
                            replacement = "", 
                            group_taxonomy_list[,1]) 
group_taxonomy_list <- gsub(pattern = "; [A-Z][a-z]+ [a-z]+)", "", group_taxonomy_list)
group_taxonomy_list

species.list <- ggtree::get.offspring.tip(pf$tree, node=1087)
mean(pf.tree$ggplot$data[pf.tree$ggplot$data$label %in% species.list,]$angle-90)

#1087 Noctilionoidea

