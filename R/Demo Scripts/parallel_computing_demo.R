library(magrittr)
library(parallel)
library(tictoc)

bar <- function(x) x

foo <- function(x){
  bar(x)^2
} 

x <- as.list(seq(1,1e7))

## serial

tic()
Sqr <- sapply(x,foo)
toc()

## parallel
ncores <- detectCores()-1
mycluster <- makeCluster(ncores)
clusterExport(mycluster,varlist='bar')
clusterEvalQ(cl = mycluster,
             expr = {
               library(nlme)
             })

tic()
Sqrpar <- parSapply(cl=mycluster,
                    X=x,
                    FUN=foo)
toc()



stopCluster(mycluster)


x <- 1:7
slp <- function(x) Sys.sleep(4)

tic()
sapply(x,Sys.sleep)
toc()

cl <- makeCluster(ncores)
clusterExport(cl,'slp')
tic()
parSapply(cl=cl,X=x,FUN=slp)
toc()

