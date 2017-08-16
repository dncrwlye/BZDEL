a=5
b=1
Omega=5           #system size parameter
C=1/Omega
m=.1              #migration needed to prevent fixation in the event of extinction


xx=1
X=xx

Time=0
Tmax=5

tt=0
Events=c('Birth','Death')

### Gillespie algorithm for simulating logistic growth one-step process

while(tt<Tmax){

  p1=a*xx+m
  p2=b*xx+C*xx^2

  S=p1+p2

  Dt=rexp(1,S)
  Dx=sample(Events,size = 1,prob = c(p1,p2))

  tt=tt+Dt
  Time<-c(Time,tt)
  X <- c(X,xx)

  if(Dx=='Birth'){
    xx=xx+1
  } else {
    xx=xx-1
  }

  Time<-c(Time,tt)
  X <- c(X,xx)

}

###### Stochastic differential equation:
#### dY=(a-b)*Y(1-Y/K)dt+(a+b)*Y(1+Y/K)dW

Dt=0.01
Y <- numeric(Tmax/Dt+1)
Ty <- seq(0,Tmax,by=Dt)
dW <- rnorm(length(Y),sd=sqrt(Dt))
Y[1]=X[1]
for (i in 2:length(Y)){
  Y[i]=Y[i-1]+
    ((a-b)*Y[i-1]-C*Y[i-1]^2)*Dt+
    sqrt((a+b)*Y[i-1]+C*Y[i-1]^2)*dW[i]
}

# plot(Time,X,type='l')
lines(Ty,Y,type='l',col='blue')
lines(Time,X,lwd=2)
lines(c(0,Tmax),c((a-b)/C,(a-b)/C),lwd=2,col='red')

