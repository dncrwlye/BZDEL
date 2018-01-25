pf.filo_pos <- gpf(Z.joined.filo.pos.1,bat_tree.filo.pos,frmla.phylo =filovirus_positive_cat~filo_samps+phylo,nfactors=10, algorithm='phylo',family=binomial)

indexes1 <- pf.filo_pos$bins[[1]]
indexes2 <-pf.filo_pos$bins[[2]]
species1<- (bat_tree.filo.pos$tip.label[indexes1]) %>% as.list()
species2 <- (bat_tree.filo.pos$tip.label[indexes2]) %>% as.list()
list <- c(list(species1, species2))

drop.tips <- Z.joined.filo.pos.1 %>%
  filter(!(Species %in% species1) & !(Species %in% species2)) %>%
  as.vector()

new.Data <- Z.joined.filo.pos.1 %>%
  filter(Species %in% species1 | Species %in% species2) %>%
  as.data.table()

setkey(new.Data, Species)

new.X <- Z.joined.filo.pos.1 %>%
  filter(Species %in% species1 | Species %in% species2) %>%
  select(c(Sample, filo_samps)) %>%
  as.data.table()

setkey(new.X, Sample)

new.tree <- ape::drop.tip(bat_tree.filo.pos,drop.tips$Species)

new.Data$Species %in% new.tree$tip.label
new.tree$tip.label %in% new.Data$Species

getObjective(grp=list,
             tree=new.tree,
             Data=new.Data,
             X=new.X,
             frmla=filovirus_positive_cat~filo_samps+phylo,
             mStableAgg=F)


  