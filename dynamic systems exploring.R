
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


