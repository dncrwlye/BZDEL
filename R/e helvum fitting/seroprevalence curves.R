## e helvum seropositivity to henipavirus
## dbecker@uga.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## read in eido data
eidoset=read.csv("ehelvum seroprev timeseries.csv",header=T)

## seroprevalence per sampling date
library(prevalence)
propCI(x=table(eidoset$Henipavirus)["1"],n=nrow(eidoset),method="wald")

## all prevalence
day_prev=with(eidoset,propCI(table(totaldays,Henipavirus)[,"1"],
                             rowSums(table(totaldays,Henipavirus)),method="wald"))

## get totaldays
day_prev$totaldays=as.numeric(rownames(day_prev))

## fix lower to zero
day_prev$lower=ifelse(day_prev$lower<0,0,day_prev$lower)

## plot
par(mfrow=c(1,1),mar=c(4.5,4.5,0.5,0.5),oma=c(0,0,0,0))
plot(p~totaldays,data=day_prev,las=1,xlab="Julian date from 2007 to 2010",
     ylab="henipavirus seroprevalence in Eidolon helvum, Accra",type="b")

## repeat by week
week_prev=with(eidoset,propCI(table(totalweek,Henipavirus)[,"1"],
                              rowSums(table(totalweek,Henipavirus)),method="wald"))

## get totaldays
week_prev$totalweek=as.numeric(rownames(week_prev))

## fix lower to zero
week_prev$lower=ifelse(week_prev$lower<0,0,week_prev$lower)

## plot
par(mfrow=c(1,1),mar=c(4.5,4.5,0.5,0.5),oma=c(0,0,0,0))
plot(p~totalweek,data=week_prev,las=1,xlab="Julian week from 2007 to 2010",
     ylab="henipavirus seroprevalence in E helvum, Accra",type="b")
