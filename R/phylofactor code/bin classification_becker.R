## phylofactor bin classification script
## Daniel Becker
## daniel.becker3@montana.edu

## clear workspace
rm(list=ls()) 
graphics.off()

## packages
library(readxl) ## alternative if xlsx fails (25 July 2017)
library(plyr)

## load in data
setwd("~/Desktop/BZDEL/Mammalian_Virus_Phylofactorization/Datasets")
pfactor=read_excel("Olival_w_phylobin_final.xlsx")

## fix virus_name_dc
library(dplyr)
library(magrittr)
pfactor_red <- pfactor %>%
  mutate(virus_name_dc = gsub('_[1-9]', ' ', vVirusNameCorrected)) %>%
  ################################################################ the following mutates needs to take their source material from the new variable
  mutate(virus_name_dc = gsub('-[1-9]', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub(' [0-9]', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('_', ' ', virus_name_dc)) %>%
  mutate(virus_name_dc = gsub('-', ' ', virus_name_dc)) %>% 
  mutate(virus_name_dc = gsub('\\s+$', '', virus_name_dc)) %>% 
  mutate(virus_name_dc = tolower(virus_name_dc)) 

## merge in zoonotic virus trait data
setwd("~/Desktop/BZDEL/Mammalian_Virus_Phylofactorization/Olival")
prelease=read.csv("pathogen release data_done.csv",header=T)

## retain just virus_name_dc, pr_vector, pr_slaughter, pr_excrete, pr_source, pr_notes
keeps=c("virus_name_dc","pr_vector","pr_slaughter","pr_excrete","pr_source","pr_notes")
prelease=prelease[keeps]

## merge
data=merge(pfactor_red,prelease,by="virus_name_dc",all.x=T)

## load in human human
hh=read.csv("human to human zoonosis data_DB edit.csv",header=T)
keeps=which(colnames(hh)=="virus_name_dc"):ncol(hh)
hh=hh[keeps]

## rename
names(hh)=c("virus_name_dc","hh_trans","hh_source","hh_notes")

## merge
data=merge(data,hh,by="virus_name_dc",all.x=T)

## just names for bin classification
knames=c(names(hh),names(prelease))
knames=unique(knames)

## add phylofactor bins and zoonotics
keeps=c(knames,"Envelope","IsZoonotic","IsZoonotic.stringent","pfIsZ","pfIsZS","nonZ","pfCBdist")
data=data[keeps]

## convert columns to character
data=data %>%
  mutate_all(as.character)

## make slaughter * excretion
data$pr_se=with(data,as.numeric(pr_slaughter)*as.numeric(pr_excrete))

## convert columns to character
data=data %>%
  mutate_all(as.character)

## convert NA to unknown
data$hh_trans=ifelse(is.na(data$hh_trans),"unknown",data$hh_trans)
data$pr_excrete=ifelse(is.na(data$pr_excrete),"unknown",data$pr_excrete)
data$pr_slaughter=ifelse(is.na(data$pr_slaughter),"unknown",data$pr_slaughter)
data$pr_vector=ifelse(is.na(data$pr_vector),"unknown",data$pr_vector)
data$pr_se=ifelse(is.na(data$pr_se),"unknown",data$pr_se)

## revalue
data$hh_trans=revalue(data$hh_trans,c("no"="no onward transmission",
                        "yes"="onward transmission"))

## classify the zoonoses for IsZoonotic
zdata=data[which(data$IsZoonotic=="TRUE"),]

## sort bins
lvls=as.character(sort(unique(as.character(zdata$pfIsZ))))
zdata$pfIsZ=factor(zdata$pfIsZ,levels=lvls)

## hh_trans
par(mar=c(4.5,5,2.5,0.5),mfrow=c(1,1))
barplot(prop.table(table(zdata$hh_trans,zdata$pfIsZ),margin=2),las=1,legend=T,
        ylab="proportion of zoonotic viruses",
        xlab="phylofactor bin",
        args.legend=list(x="top",bty="n",horiz=T,inset=c(-0.125,-0.125),
                         text.width=c(3,6,0)))

## pr_vector
#zdata$pr_vector=revalue(zdata$pr_vector,c("0"="not vector","1"="vector-borne"))
par(mar=c(4.5,5,2.5,0.5),mfrow=c(1,1))
barplot(prop.table(table(zdata$pr_vector,zdata$pfIsZ),margin=2),las=1,legend=T,
        ylab="proportion of zoonotic viruses",
        xlab="phylofactor bin",
        args.legend=list(x="top",bty="n",horiz=T,inset=c(-0.125,-0.125),
                         text.width=c(3,6,0)))

## tabulate proportions for IsZoonotic
pdata=cbind(t(prop.table(table(zdata$pr_vector,zdata$pfIsZ),margin=2)),
      t(prop.table(table(zdata$pr_excrete,zdata$pfIsZ),margin=2)),
      t(prop.table(table(zdata$pr_slaughter,zdata$pfIsZ),margin=2)),
      t(prop.table(table(zdata$hh_trans,zdata$pfIsZ),margin=2)))
pdata=data.frame(pdata)
pdata=round(pdata,2)

## write to csv
setwd("~/Dropbox (MSU projects)/Phylofactor/bin classifications")
write.csv(pdata,"pdataIsZ.csv")
table(zdata$pfIsZ)

## repeat for zoonotic stringent
sdata=data[which(data$IsZoonotic.stringent=="TRUE"),]

## sort bins
lvls=as.character(sort(unique(as.character(sdata$pfIsZS))))
sdata$pfIsZS=factor(sdata$pfIsZS,levels=lvls)

## tabulate proportions for IsZoonotic
psdata=cbind(t(prop.table(table(sdata$pr_vector,sdata$pfIsZS),margin=2)),
            t(prop.table(table(sdata$pr_excrete,sdata$pfIsZS),margin=2)),
            t(prop.table(table(sdata$pr_slaughter,sdata$pfIsZS),margin=2)),
            t(prop.table(table(sdata$hh_trans,sdata$pfIsZS),margin=2)))
psdata=data.frame(psdata)
psdata=round(psdata,2)
psdata

## write to csv
setwd("~/Dropbox (MSU projects)/Phylofactor/bin classifications")
write.csv(psdata,"pdataIsZS.csv")
table(sdata$pfIsZS)

## repeat for X
cdata=data[which(data$IsZoonotic=="TRUE"),]

## sort bins
cdata$pfCBdist=factor(cdata$pfCBdist,levels=c("Family_Togaviridae","Genus_Flavivirus","Family_Bunyaviridae",
                                              "Family_Picornaviridae","Genus_Orthopoxvirus",))

## tabulate proportions for IsZoonotic
psdata=cbind(t(prop.table(table(sdata$pr_vector,sdata$pfIsZS),margin=2)),
             t(prop.table(table(sdata$pr_excrete,sdata$pfIsZS),margin=2)),
             t(prop.table(table(sdata$pr_slaughter,sdata$pfIsZS),margin=2)),
             t(prop.table(table(sdata$hh_trans,sdata$pfIsZS),margin=2)))
psdata=data.frame(psdata)
psdata=round(psdata,2)

## write to csv
setwd("~/Dropbox (MSU projects)/Phylofactor/bin classifications")
write.csv(psdata,"pdataIsZS.csv")
table(sdata$pfIsZS)

