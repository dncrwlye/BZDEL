require(abind)
library(tictoc)
### parameters
beta=matrix(2e-5)
K=matrix(1e4)
b=.4/365
sigma=1/6
g=1/7
tau=0.1
# C <- matrix(c(0,1,1,0),nrow=2)
C <- 0

mu=1/3650
Tmax=365
Y0=matrix(c(9999,0,1,0),ncol=1)
m=0

p=length(K)
if (!length(beta)==p){stop('beta and K are not same length')}

if (dim(C)[1]==dim(C)[2]){
  if (!dim(C)[1]==p){
    stop('dimension of C does not match beta')
  }
} else {stop('C is not square')}
### state matrix
if (is.null(Y0)){
  Y <- matrix(NA,nrow=4,ncol=p)
  for (i in 1:p){
    Y[,i] <- rmultinom(1,K[i],prob=c(.99,0,.01,0))
  }
} else {
  if (!all.equal(dim(Y0),c(4,p))){
    stop('Dimension of Y0 is not 4xlength(K)')
  }
  Y <- Y0
}

SEIR <-  c('S','E','I','R')
rownames(Y) <- SEIR

### indexes
births <- 1:p
SE <- (1:p)+p
EI <- (1:p)+2*p
IR <- (1:p)+3*p
death <- 4*p+1
migration <- 4*p+2

### propensities
lambdas <- numeric(4*p+2)
noLambdas <- length(lambdas)

### output
output <- NULL
output$time <- 0
output$SEIR <- Y


tt <- 0
  
tic()
while (tt<Tmax){
  
  N <- colSums(Y)
  Ntot <- sum(N)
  ################### update lambdas
  lambdas[births] <- max(b*(1-N/K),0)                      #births
  lambdas[SE] <- beta*Y['S',]*Y['I',]                      #S-->E
  lambdas[EI] <- sigma*Y['E',]                             #E-->I
  lambdas[IR] <- g*Y['I',]                                 #I-->R
  lambdas[death] <- mu*Ntot                                #deaths
  lambdas[migration]<- m*Ntot                              #migration
  
  S <- sum(lambdas)
  
  ################## draw timestep
  if (is.null(tau)){
    dt <- rexp(1,rate=S)
  } else {
    dt <- tau
  }
  ################## update time
  tt <- tt+dt
  output$time <- c(output$time,tt)
  
  ################## draw event
  if (is.null(tau)){
    events <- sample(noLambdas,size=1,prob = lambdas)
  } else {
    nevents <- rpois(1,lambda=tau*S)
    events <- sample(noLambdas,size=nevents,replace = T,prob=lambdas)
  }
  
  ################## update Y
  for (event in events){
    if (event %in% births){
      patch=event
      Y['S',patch] <- Y['S',patch]+1
    } else if (event %in% SE){
      patch=event-p
      Y[c('S','E'),patch]=Y[c('S','E'),patch]+c(-1,1)
    } else if (event %in% EI){
      patch=event-2*p
      Y[c('E','I'),patch] <- Y[c('E','I'),patch]+c(-1,1)
    } else if (event %in% IR){
      patch=event-3*p
      Y[c('I','R'),patch] <- Y[c('I','R'),patch]+c(-1,1)
    } else if (event %in% death){
      #sample patch
      patch <- sample(x=p,size=1,prob=N)
      #sample stage SEIR
      stage <- sample(SEIR,size=1,prob=Y[,patch])
      Y[stage,patch] = Y[stage,patch]-1
    } else { #migration
      #sample patch of migrant
      ipatch <- sample(p,size=1,prob=N)
      stage <- sample(SEIR,size=1,prob=Y[,ipatch])
      jpatch <- sample(p,size=1,prob=C[ipatch,])
      
      Y[stage,c(ipatch,jpatch)] <- Y[stage,c(ipatch,jpatch)] + c(-1,1)
    }
  }
    
  output$SEIR <- abind(output$SEIR,Y,along=3)
}

  
    


library(plotrix)




output <- mSEIR()
Time <- output$time
Y <- output$SEIR
p <- dim(Y)[2]

Y <- aperm(Y,perm=c(3,1,2))
par(mfrow=c(p+1,1))

cols <- rainbow(4)
for (i in 1:p){
  stackpoly(Time,Y[,,i],stack=T,col = cols)
  if (i==1){
    pp <- par('usr')
    legend((pp[2]-pp[1])/2,0.9*pp[4],legend=c('S','E','I','R'),fill=cols)
  }
}

# plot(Time,colSums(output$SEIR['I',,]))
plot(Time,Y['I',,])
# lines(Time,output$SEIr['I',,])


