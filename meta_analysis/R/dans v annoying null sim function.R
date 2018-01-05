dc.null.sim.function <- function (Z, outcome, NoReplicates, tree, predictor, n.factors, fisher.storage, deviance.storage, aov.storage)
{
  #Z.joined.filo.pos.1 = Z
  #filovirus_positive_cat = outcome
  #NoReplicates  100
  #bat_tree.filo.pos = tree
  #filo_samps = predictor
  #pf.filo_pos.pval.storage = fisher.storage
  #pf.filo_pos.pval.storage.alt = deviance.storage
  #pf.filo_pos.pval.storage.alt.aov = aov.storage
  
  randomData <- function()
  {
    rbinom(nrow(Z), size = 1, prob = sum(outcome)/nrow(Z))
  }
  p_storage <- matrix(0, NoReplicates, n.factors)
  p_storage.alt <- matrix(0, NoReplicates, n.factors)
  p_storage.alt.aov <- matrix(0, NoReplicates, n.factors)
  
  
  for (replicate in 1:NoReplicates){
    outcome <- randomData()  ## define how you want to make this random data
    DF <- Z %>%
      mutate(outcome = outcome) %>%
      rename(outcome_random = outcome)
    pf <- gpf(DF,tree,frmla=outcome_random~predictor+phylo,nfactors=n.factors,mStableAgg=F)
    
    for (i in 1:n.factors)
    {
      p_storage.alt[replicate, i] <- pchisq(pf$models[[i]]$null.deviance- pf$models[[i]]$deviance, df=(pf$models[[1]]$df.null-pf$models[[1]]$df.residual))
      
      x<-summary(aov(pf$models[[i]])) 
      p_storage.alt.aov[replicate, i] <- x[[1]]$`Pr(>F)`[2]
    }
    
    for (i in 2:11)
    {
      indexes = pf$bins[[i]]
      p<-DF[c(-indexes),] 
      fisher.matrix <- matrix(c(sum(p$outcome_random  == 1),
                                sum(p$outcome_random  == 0),  
                                sum(DF$outcome_random == 1),
                                sum(DF$outcome_random == 0)), nrow=2)
      p <-fisher.test(fisher.matrix)$p.value 
      
      p_storage[replicate, i-1] <- p
    }
  }
  
  obs <- rbind('factor'=1:n.factors, 'Pvals'= fisher.storage)
  
  ddf.filo.null.sim <- as.data.frame(t(rbind(obs, p_storage)))
  
 one <-  ddf.filo.null.sim %>%
    tidyr::gather(group, pvalue, c(2:(NoReplicates+2))) %>%
    mutate(color = ifelse(group =='Pvals', "observed", 'simulated')) %>%
    mutate(pvalue = as.numeric(pvalue)) %>%
    mutate(factor = as.numeric(factor)) %>%
    ggplot() +
    geom_line(aes(x = factor, y= log(pvalue), group = group, alpha = .5, color  = color))
  
obs.alt <- rbind('factor'=1:n.factors,
                   'Pvals'=deviance.storage)
  
  ddf.filo.null.sim.alt <- as.data.frame(t(rbind(obs.alt, p_storage.alt)))
  
  two <- ddf.filo.null.sim.alt %>%
    tidyr::gather(group, pvalue, c(2:(NoReplicates+2))) %>%
    mutate(color = ifelse(group =='Pvals', "observed", 'simulated')) %>%
    mutate(pvalue = as.numeric(pvalue)) %>%
    mutate(factor = as.numeric(factor)) %>%
    ggplot() +
    geom_line(aes(x = factor, y= (pvalue), group = group, color  = color))+
    ggtitle('positive filovirus simulations') +
    scale_x_continuous(limits=c(1, n.factors), breaks = seq(1,n.factors, by= 1))
  
  
  obs.alt.aov <- rbind('factor'=1:n.factors,
                       'Pvals'=aov.storage)
  
  ddf.filo.null.sim.alt.aov <- as.data.frame(t(rbind(obs.alt.aov, p_storage.alt.aov)))
  
 three <-  ddf.filo.null.sim.alt.aov %>%
    tidyr::gather(group, pvalue, c(2:(NoReplicates+2))) %>%
    mutate(color = ifelse(group =='Pvals', "observed", 'simulated')) %>%
    mutate(pvalue = as.numeric(pvalue)) %>%
    mutate(factor = as.numeric(factor)) %>%
    ggplot() +
    geom_line(aes(x = factor, y= (pvalue), group = group, color  = color))+
    ggtitle('positive filovirus simulations') +
    scale_x_continuous(limits=c(1, n.factors), breaks = seq(1,n.factors, by= 1))
 
 x <- list(one, two, three)
 return(x)
}
