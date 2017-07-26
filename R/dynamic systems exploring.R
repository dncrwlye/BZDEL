
r=2
k=1000
nt_matrix <- matrix(0,50,100)

for(l in 1:20)
{
nt=1
for (i in 1:100)
{
  delta_i <- r * nt* (1- nt/k)
  nt = nt + abs(rnorm(1, mean = 1, sd = .5))*delta_i
  nt_matrix[l,i] <- nt
}
}

plot(nt_matrix[1,], type = 'l')
for(z in 2:20)
{
lines(nt_matrix[z,], type = 'l')
}



ds
nt_ave_matrix <-matrix(0,1,25)
for (r in seq(.5,3, by = .001))
{
r=r
print(r)

k=1000
nt=1
nt_matrix <- matrix(0,1,500)
for (i in 1:500)
{
  delta_i <- r * nt* (1- nt/k)
  nt = nt + delta_i
  nt_matrix[1,i] <- nt
}
plot(nt_matrix[1,])
nt_ave_matrix[1,r]<-mean(nt_matrix[1,])
}

plot(nt_ave_matrix[1,])

###########################################################################

nt_matrix <-matrix(0,1,5000)
k=1000
r=.0205
N_o = 100
x<- seq(1,500, by=1)
for (i in x)
{
  print(i)
  nt <- (N_o * k) / ((k-N_o)*(exp(-r*x[i]))+N_o)
  
  nt_matrix[1,i] <- nt
}

plot(nt_matrix[1,1:500], type='l')



###############################################################################
library(deSolve)

population=function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    # rate of change
    dN=r*nt*(1-nt/k);
       list(c(dN))  
    }) 
}

times=seq(0,10,by=.01);    #i assume this is running it for 50 years
state=c(nt=1)
parameters =c(r=1, k=10)

sim=ode(y=state,times=times,func=population,parms=parameters);

#graphing the SIR model with demographics over 50 years

plot(sim)

####################
#random walk

N <-100
ZZ <-seq(1,1.3,by=.1)
x_storage <- matrix(0,length(ZZ),N+1)

for (M in length(ZZ))
{
  print(M)x
  t<-ZZ[M]
  print(t)
  #print(t)
  x<-3
  #t=100
  seq <- seq(1,N, by = t)
  for(i in seq)
  {
    r <- rbinom(1,1,.4)
    #print(r)
    if (r == 1)
    {
      x<-x + t*rnorm(1,0,t)
    }
    if (r == 0)
    {
      x<- x- t*rnorm(1,0,.1) 
    }
    x_storage[M,i] <- x
  }
  x_storage[M,101]<-mean(x_storage[M,1:100])
}
plot(x_storage[1,], type='l')
for(i in 2:length(ZZ))
{lines(x_storage[i,])}



hist(x_storage[,101])
mean(x_storage[])

####################



x<-matrix(0,1,1000)
for(i in 1:1000)
{
  
  x[,i] <-sum(rbinom(1000, 1, .5))
  
}

hist(x[,])



#####################################
#wolves and rabbits

LotVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*alpha - beta*y*x
    dy = y*-gamma + delta*x*y
    return(list(c(dx, dy)))
  })
}

Pars <- c(alpha = 20, beta = 2, gamma = 20, delta = 2)
State <- c(x = 11, y = 9)
Time <- seq(0, 1000, by = 10)

out <- ode(func = LotVmod, y = State, parms = Pars, times = Time)

plot(out)











