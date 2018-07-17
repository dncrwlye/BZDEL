
model.fcn2 <- function(formula,data,...){
  fit <- tryCatch(MASS::glm.nb(formula,data,...),
                  error=function(e) NA)
  #fit <- do.call
  return(fit)
}

obj.fcn2 <- function(fit,grp,tree,PartitioningVariables,model.fcn,phyloData,...){
  #if (!'negbin' %in% class(fit) & !'glm' %in% class(fit) & !'lm' %in% class(fit))
  if (!'negbin' %in% class(fit))
  {
    return(0)
  }
  else 
  {
    fit2 <- MASS::glm.nb(Z.poisson~1,data = fit$model)
    fit$twologlik-fit2$twologlik %>% return()
    #fit$twologlik %>% return()
  }
}

tictoc::tic()
pf2 <- gpf(Data,tree,Z.poisson~phylo,nfactors=10,algorithm = 'phylo',
          model.fcn = model.fcn2,objective.fcn = obj.fcn2,
          cluster.depends='library(MASS)', 
          ncores = 4)
tictoc::toc()

# double check that the negative binomial is an improvement over the poisson ----
summaries.null <- (lapply(pf2$models, "[[", "model"))
poisson.models <- lapply(summaries.null, 
                         function(data){
                           glm(Z.poisson~phylo, data, family = "poisson")
                         })
logLik.poisson <- lapply(poisson.models,logLik)
logLik.negbinom <- lapply(pf2$models,logLik)

for (i in 1:length(logLik.negbinom))
{
pchisq(2 * (logLik.negbinom[[i]] - logLik.poisson[[i]]), df = 1, lower.tail = FALSE) %>% print()
}

#major improvement in negative binomial regression, I believe. At least for first 6 factors

