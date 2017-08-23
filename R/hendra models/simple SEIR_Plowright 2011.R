## simple model of hendra virus dynamics within a camp
## closed population
## dbecker@uga.edu
## 08/01/2017

## clean
rm(list=ls()) 
graphics.off()

## fresh plot window
par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),mfrow=c(1,1),oma=c(0,0,0,0))

## load packages
library(deSolve)
library(png)

## ODE from Plowright et al 2011 (no spatial structure)
hendramod=function(t,y,params){
  
  ## state variables: S, E, I, R
  S=y[1]
  E=y[2]
  I=y[3]
  R=y[4]

  ## equations with parameters as a list
  with(as.list(params),{
    
    dS=(b*(1-((S+E+I+R)/k))) - ((beta*I)*S) - (d*S) + (e*R)
    dE=((beta*I)*S) - (sigma*E) - (d*E)
    dI=(sigma*E) - (gamma*I) - (d*I) 
    dR=(gamma*I) - (d*R) - (e*R)
      
    ## vector of state variables, as list
    dy=c(dS,dE,dI,dR)
    list(dy)})}

## define starting conditions
N0=10000
infecteds=1
prop=0
immune=prop/100*N0
R0=immune
S0=N0-infecteds-immune
I0=infecteds
E0=0
Y0=c(S0,E0,I0,R0)

## set parameters in days
beta=mean(c(2e-5,5e-5))
beta=2e-5
gamma=1/7
sigma=1/6
b=0.4/365
d=1/(10*365)
k=N0/(1-(d/b))

## waning immunity
#e=1/182
e=0 ## lifelong immunity

## parameter vector
pars=c(beta=beta,gamma=gamma,sigma=sigma,b=b,d=d,k=k,e=e)

## time in days
life=(1/d)/365
yrs=life
yrs=5
tmax=365*yrs
time=seq(1,tmax,by=1)

## run model
out=as.data.frame(lsoda(Y0,time,hendramod,pars))  ## store in data frame
names(out)=c("times","S","E","I","R") ## make names easy to call

## population size
out$N=with(out,S+E+I+R)

## derive seroprevalence
out$seroprev=with(out,(I+R)/N)

## plot dynamics
par(mar=c(4.5,4.5,0.5,0.5),oma=c(0,0,0,0))
plot(seroprev~time,col="purple",type="l",lwd=2,data=out,xlab="years after viral introduction",
     ylab="Hendra virus seroprevalence in Pteropids",las=1,xaxt="n")
axis(1,at=seq(0,tmax,l=yrs),labels=seq(1,tmax/365,l=yrs))

plot(I~time,col="purple",type="l",lwd=2,data=out,xlab="years after viral introduction",
     ylab="number of infected Pteropids",las=1,xaxt="n")
axis(1,at=seq(0,tmax,l=yrs),labels=seq(1,tmax/365,l=yrs))
