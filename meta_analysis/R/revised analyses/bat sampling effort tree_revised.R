## code for bat filovirus and HNV sampling effort phylofactor
## daniel.becker3@montana.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(rotl)
library(caper)
library(ape)
library(plyr)

## read in Rdata file
setwd("~/Dropbox (MSU projects)/Spillover postdoc/bat virus meta-analysis")
load("seroprevalence_revised.Rdata")
#setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
#load(file='data/seroprevalence.Rdata')

#load("/Users/danielbecker/Desktop/BZDEL/meta_analysis/data/seroprevalence.Rdata")
data=seroprevalence
data=data.frame(data)
rm(seroprevalence)

## make seroprevalence accessible
data$serop=data$seroprevalence_percentage
data$serop=as.numeric(data$serop)/100

## clean
data$seroprevalence_percentage=NULL

## sample size accessible
data$sample=data$sample_size
data$sample=as.numeric(data$sample)

## make detection
data$detection=ifelse(data$serop>0,1,0)

## convert virus to lower
data$virus=tolower(data$virus)

## standardize virus
table(data$virus)
data$virus=revalue(data$virus,c("bombali virus"="filovirus",
                                "cedar"="henipavirus",
                                "ebola"="filovirus",
                                "hendra"="henipavirus",
                                "nipah"="henipavirus"))
table(data$virus)

## if methodology==isolation
data$detection=ifelse(data$methodology=="isolation",1,data$detection)

## trim data down to just species, detection, and virus
data=data[c("species","virus","detection","sample")]

## fix species names
data$species=revalue(data$species,c("rhinolophus refulgens"="rhinolophus lepidus",
                                    "hipposideros cf caffer"="hipposideros caffer",
                                    "hipposideros cf caffer/ruber"="hipposideros caffer",
                                    "hipposideros cf ruber"="hipposideros ruber",
                                    "natalus lanatus"="natalus mexicanus",
                                    "myonycteris leptodon"="myonycteris torquata",
                                    "tadarida plicata"="chaerephon plicatus",
                                    "pipistrellus cf nanus/nanulus"="neoromicia nanus",
                                    "neoromicia nana"="neoromicia nanus",
                                    "hypsygnathus monstrosus"="hypsignathus monstrosus",
                                    "miniopterus schrebersii"="miniopterus schreibersii",
                                    "lissonycteris angolensis"="myonycteris angolensis",
                                    "pteropus medius"="pteropus giganteus",
                                    "artibeus planirostris"="artibeus jamaicensis"))

## fix mini
data$species=revalue(data$species,c("miniopterus schreibersii"="miniopterus schreibersii schreibersii"))

## record pooled ones
pool=c("bat","bats*",
       "carollia species","chalinobus species","cynopterus species","eonycteris spelaea & rousettus species",
       "epomophorus species","glossophaginae species1","glossophaginae species3",
       "hipposideros species","hypsignathus monstrosus, epomops franqueti, myonycteris torquata",
       "microchiroptere species","micropteropus/nanonycteris",
       "miniopterus species","molossideo species1","mops species",
       "myotis species","nycteris species","nyctophilus species",
       "pipistrellus species",
       "pteropus conspicillatus & dobsonia magna",
       "pteropus species","rhinolophus species","scotorepens species",
       "synconycterus species","tadarida species","unknown species",
       "vespertilionidae species1","vespertilioniformes species",
       "unidentified nycterid bat",
       "unidentified myonycteris/epomophorous bat",
       "unidentified molossid bat",
       "unidentified hipposiderod bat",
       "unidentified epomophorous bat",
       "unidentified chaerephon bat",
       "rhinolophus sp.",
       "nycteris sp.",
       "mops sp.",
       "miniopterus sp.",
       "microchiroptere species",
       "kerivoula sp.",
       "hipposideros sp.",
       "hipposideros ruber/caffer",
       "glossophaginae species1",
       "glossophaginae species3",
       "epomophorus species",
       "eonycteris spelaea & rousettus species",
       "cynopterus species",
       "chalinobus species",
       "chaerephon sp.",
       "carollia species",
       "neoromicia sp.")
pool=sort(unique(pool))

## code to drop these
data$phy_drop=ifelse(data$species%in%pool,"genus pool","keep")

## subset
data=data[-which(data$phy_drop=="genus pool"),]

## fix multiple species
data$species_final=revalue(data$species,c("pteropus alecto, pteropus poliocephalus"="pteropus alecto",
                                          "pteropus alecto, pteropus poliocephalus, pteropus scapulatus"="pteropus alecto",
                                          "pteropus alecto, pteropus scapulatus"="pteropus alecto",
                                          "pteropus conspicillatus, pteropus scapulatus"="pteropus conspicillatus"))
data$species=data$species_final
data$species_final=NULL

## data frame of all species
sdata=data.frame(sort(unique(data$species)))
names(sdata)="species"

## filovirus
fdata=data[which(data$virus=="filovirus"),]
filodata=data.frame(table(fdata$species))
names(filodata)=c("species","filo_samps")

## hnv
hdata=data[which(data$virus=="henipavirus"),]
hnvdata=data.frame(table(hdata$species))
names(hnvdata)=c("species","hnv_samps")

## merge into sdata
sdata=merge(sdata,filodata,by="species",all.x=T)
sdata=merge(sdata,hnvdata,by="species",all.x=T)
rm(hnvdata,filodata)

## total detections
sdata$total_samps=rowSums(sdata[c("filo_samps","hnv_samps")],na.rm=T)

## tabulate total detection
adetect=as.data.frame.matrix(table(data$species,data$detection))
names(adetect)=c("no","yes")
adetect$species=rownames(adetect)
adetect$positive=ifelse(adetect$yes>0,1,0)
adetect=adetect[c("species","positive")]

## tabulate filovirus detection
fdetect=as.data.frame.matrix(table(fdata$species,fdata$detection))
names(fdetect)=c("no","yes")
fdetect$species=rownames(fdetect)
fdetect$filo_positive=ifelse(fdetect$yes>0,1,0)
fdetect=fdetect[c("species","filo_positive")]

## tabulate henipavirus detection
hdetect=as.data.frame.matrix(table(hdata$species,hdata$detection))
names(hdetect)=c("no","yes")
hdetect$species=rownames(hdetect)
hdetect$hnv_positive=ifelse(hdetect$yes>0,1,0)
hdetect=hdetect[c("species","hnv_positive")]

## merge into sdata
sdata=merge(sdata,fdetect,by="species",all.x=T)
sdata=merge(sdata,hdetect,by="species",all.x=T)
sdata=merge(sdata,adetect,by="species",all.x=T)
rm(hdata,fdata,hdetect,fdetect,adetect)

## sdata sort
sdata=sdata[order(sdata$species),]

## obtain host phylogeny, remove spp records
library(rotl)
library(ape)
phy=tnrs_match_names(names=as.character(sdata$species),context_name="Animals",do_approximate_matching=T)

## bind into data frame
sdata=data.frame(sdata,phy)
rm(phy)

## order
sdata=sdata[order(sdata$unique_name,sdata$search_string),]

## remove NA
sdata=sdata[!is.na(sdata$ott_id),]

## rotl the tree
tree=tol_induced_subtree(ott_ids=sdata$ott_id)

## remove ott information from the tips
tree$tip.label=strip_ott_ids(tree$tip.label)

## check binary
is.binary.tree(tree)

## resolve multifurcations
tree=multi2di(tree)
is.binary.tree(tree)

## check ultrametric
is.ultrametric(tree)

## assign branch lengths
tree=compute.brlen(tree,method="Grafen")
is.ultrametric(tree)

## ladderize tree
tree=ladderize(tree)

## plot
par(oma=c(0,0,0,0),mar=c(1,1,1,1),mfrow=c(1,1))
plot(tree,cex=0.5)

## do all our host names match those in the tree?
length(unique(sdata$unique_name))
length(unique(tree$tip.label))

## save previous unique names
sdata$species=sdata$unique_name

## convert names to match
sdata$unique_name=gsub(" ","_",sdata$unique_name)
rm(data)

## which names don't match
sdata$unique_name[!sdata$unique_name%in%tree$tip.label]

## fix sdata and tree
sdata$unique_name=revalue(sdata$unique_name,c("Miniopterus_schreibersii_schreibersii"="Miniopterus_schreibersii"))
tree$tip.label=revalue(tree$tip.label,c("Miniopterus_schreibersii_schreibersii"="Miniopterus_schreibersii"))

## load wilson and reeder bat data
setwd("~/Dropbox (MSU projects)/Spillover postdoc/bat virus meta-analysis/trait data")
bats=read.csv("wilson and reeder_chiroptera.csv",header=TRUE)
#bats = read.csv("/Users/buckcrowley/Dropbox_gmail/Dropbox/bat virus meta-analysis/trait data/wilson and reeder_chiroptera.csv", header =TRUE)

## remove no species
bats=bats[!is.na(bats$Species),]
bats=bats[-which(bats$Species==""),]

## make species
bats$species=paste(bats$Genus,bats$Species,sep=" ")

## remove duplicates
bats=bats[!duplicated(bats$species),]

## sort by species
bats=bats[order(bats$species),]

## fix names to match our phylogeny
bats$species=revalue(bats$species,c("Lissonycteris angolensis"="Myonycteris angolensis",
                                      "Nyctalus plancyi"="Nyctalus velutinus",
                                      "Rhinolophus stheno"="Rhinolophus microglobosus",
                                      "Natalus stramineus"="Natalus stramineus mexicanus",
                                      "Neoromicia somalicus"="Neoromicia malagasyensis"))

## clean up bats data
bats=bats[c("Family","Genus","TypeLocality","species")]

## get locality
Sys.setlocale('LC_ALL','C')
bats$TypeLocality=as.character(bats$TypeLocality)
locs=strsplit(bats$TypeLocality,",")
locs=sapply(locs,function(x) x[1])
bats$loc=locs

## fix values
unique(bats$loc)
bats$region=revalue(bats$loc,c("Indonesia"="Asia",
                                  "Philippines"="Asia",
                                  "Malaysia (N Borneo)"="Asia",
                                  "Papua New Guinea"="Oceania",
                                  "Cameroon"="Africa",
                                  "Borneo"="Asia",
                                  "Sumatra"="Asia",
                                  "India"="Asia",
                                  "Solomon Isls"="Oceania",
                                  "Malaysia"="Asia",
                                  "Madagascar"="Africa",
                                  "Senegal (restricted by K. Andersen"="Africa",
                                  "Phillipines"="Asia",
                                  "Burma"="Asia",
                                  "Angola"="Africa",
                                  "Mozambique"="Africa",
                                  "Gambia"="Africa",
                                  "Sudan"="Africa",
                                  "Ethiopia"="Africa",
                                  "Zanzibar."="Africa",
                                  "South Africa"="Africa",
                                  "Liberia"="Africa",
                                  "Gabon."="Africa",
                                  "Thailand"="Asia",
                                  "Gabon"="Africa",
                                  "Solomon Isls."="Oceania",
                                  "Given by Andersen (1912:790) as \"New Ireland"="Oceania",
                                  "Nigeria"="Africa",
                                  "S\213o Tom\216 and Pr\222ncipe"="Africa",
                                  "Kenya"="Africa",
                                  "N Angola."="Africa",
                                  "Fiji Isls"="Oceania",
                                  "New Caledonia"="Oceania",
                                  "New Guinea"="Oceania",
                                  "Australia"="Australia",
                                  "Seychelles"="Africa",
                                  "Vanuatu"="Oceania",
                               "Papua New Guinea. The type locality was initially given as New Ireland Isl."="Oceania",
                               "Solomon Isl"="Oceania",
                               "Japan"="Asia",
                               "Caroline Isls"="Oceania",
                               "S Burma"="Asia",
                               "Philippines."="Asia",
                               "Comoro Isls"="Africa",
                               "West Pacific"="Asia",
                               "Mascarene Isls"="Africa",
                               "Micronesia"="Oceania",
                               "Australia."="Australia",
                               "Madagascar. Restricted to \"N. and C. Madagascar\" by K. Andersen (1908)."="Africa",
                               "Samoan Isls"="Oceania",
                               "Seychelle Isls"="Africa",
                               "Mariana Isls"="Oceania",
                               "Tonga Isls"="Oceania",
                               "Ualan (= Kosrae; Micronesia)"="Oceania",
                               "New Caledonia (France)."="Oceania",
                               "Tanzania"="Africa",
                               "W Carolines"="Oceania",
                               "Egypt"="Africa",
                               "Uganda"="Africa",
                               "Republic of Congo"="Africa",
                               "Ghana"="Africa",
                               "SE Europe; restricted to Italy by Ellerman et al. (1953:59)."="Europe",
                               "Turkmenistan"="Asia",
                               "Sulawesi"="Asia",
                               "Saudi Arabia"="Middle East",
                               "Japan."="Asia",
                               "Zimbabwe"="Africa",
                               "Italy"="Europe",
                               "France."="Europe",
                               "Taiwan."="Asia",
                               "Guinea"="Africa",
                               "Rwanda"="Africa",
                               "Key-Inseln (= Kai Isls)."="Asia",
                               "Equatorial Guinea"="Africa",
                               "Nepal."="Asia",
                               "Romania"="Europe",
                               "Taiwan"="Asia",
                               "Timor"="Asia",
                               "China"="Asia",
                               "Vietnam"="Asia",
                               "Dem. Rep. Congo"="Africa",
                               "Zambia"="Africa",
                               "Bangladesh"="Asia",
                               "Malaya"="Asia",
                               "Sri Lanka"="Asia",
                               "Pakistan"="Asia",
                               "East Solomon Isls"="Oceania",
                               "Ghana."="Africa",
                               "Sierra Leone"="Africa",
                               "C\231te d\325Ivoire"="Africa",
                               "Eritrea"="Africa",
                               "Singapore"="Asia",
                               "Laos"="Asia",
                               "Sao Tome and Princepe"="Africa",
                               "Iran"="Middle East",
                               "E Madagascar."="Africa",
                               "Ethiopia."="Africa",
                               "Senegal."="Africa",
                               "Oman"="Middle East",
                               "Mauritius."="Africa",
                               "Ecuador"="LA",
                               "Guatemala"="LA",
                               "Costa Rica"="LA",
                               "Panama"="LA",
                               "Brazil"="LA",
                               "Colombia"="LA",
                               "Samoa."="Oceania",
                               "Bismarck Archipelago"="Oceania",
                               "Surinam."="LA",
                               "Trinidad"="LA",
                               "Senegal"="Africa",
                               "\"Guinea\"."="Africa",
                               "Sierra Leone."="Africa",
                               "Somalia"="Africa",
                               "Madagascar."="Africa",
                               "New Zealand"="Oceania",
                               "New Zealand."="Oceania",
                               "Paraguay"="LA",
                               "Guyana",
                               "Cuba"="LA",
                               "Puerto Rico"="LA",
                               "Jamaica."="LA",
                               "Venezuela"="LA",
                               "Guatemala."="LA",
                               "Mexico."="LA",
                               "Mexico"="LA",
                               "Nicaragua"="LA",
                               "Peru"="LA",
                               "Mexico. The type locality was incorrectly changed to Brazil"="LA",
                               "French Guiana"="LA",
                               "Trinidad and Tobago"="LA",
                               "Surinam"="LA",
                               "Bolivia"="LA",
                               "USA"="NA",
                               "Haiti."="LA",
                               "Jamaica"="LA",
                               "Honduras"="LA",
                               "Peru."="LA",
                               "Guadeloupe (Lesser Antilles)"="LA",
                               "Dominica (Lesser Antilles)."="LA",
                               "Not designated in original publication (presumably Jamaica)."="LA",
                               "Costa Rica."="LA",
                               "Brazil; see Carter and Dolan (1978)."="LA",
                               "Honduras (= R\222o Segovia) (McCarthy et al."="LA",
                               "Not designated in original publication (probably Virgin Isls)."="LA",
                               "Surinam (restricted by Thomas"="LA",
                               "Bahamas"="LA",
                               "Dominican Republic"="LA",
                               "Cuba."="LA",
                               "N Sudan"="Africa",
                               "Namibia"="Africa",
                               "Sao Tome and Principe"="Africa",
                               "\"Brazil.\""="LA",
                               "French Guiana."="LA",
                               "Argentina"="LA",
                               "France"="Europe",
                               "Mauritius"="Africa",
                               "\"Central Peru.\""="LA",
                               "Australia. Probably New South Wales"="Australia",
                               "Kazakhstan"="Asia",
                               "Yemen."="Middle East",
                               "Mongolia"="Africa",
                               "Korea"="Asia",
                               "Sweden."="Asia",
                               "Uruguay"="LA",
                               "Haiti"="LA",
                               "Chile"="LA",
                               "Barbados (Lesser Antilles)"="LA",
                               "St. Vincent (Lesser Antilles"="LA",
                               "Liberia."="Africa",
                               "\"India\"."="LA",
                               "Portugal"="Europe",
                               "Northern Italy"="Europe",
                               "Germany"="Europe",
                               "British Papua (= Papua New Guinea)"="Oceania",
                               "Libya"="Africa",
                               "England"="Europe",
                               "Restricted by Kock (1969<i>a</i>) to the Nile Valley between north of Aswan"="Africa",
                               "Austria"="Europe",
                               "Croatia"="Europe",
                               "Spain"="Europe",
                               "Chile."="Europe",
                               "Israel"="Middle East",
                               "Botswana"="Africa",
                               "Sweden. Baag\277e (2001<i>b</i>) indicated that the type locality is probably near Uppsala"="Europe",
                               "Peking (China)."="Asia",
                               "\"Eastern United States\"."="NA",
                               "Russia"="Europe",
                               "Greece"="Europe",
                               "Tajikistan"="Asia",
                               "\"Southern China\"."="Asia",
                               "Tunisia"="Africa",
                               "Germany."="Europe",
                               "California"="NA",
                               "Singapore."="Asia",
                               "Cambodia."="Asia",
                               "NE Angola."="Africa",
                               "\"Eastern coast of South Africa.\""="Africa",
                               "Hungary (Pliocene). See discussion in Hor\207cek et al. (2000)."="Europe",
                               "Martinique (Lesser Antilles)"="LA",
                               "South Australia."="Australia",
                               "Canada"="NA",
                               "\"Southern Brazil.\""="LA",
                               "Armenia"="Asia",
                               "Denmark"="Europe",
                               "Nepal"="Asia",
                               "United States"="NA",
                               "Philippine Isls."="Oceania",
                               "Guyana"="LA",
                               "Cura\215ao"="LA",
                               "Specified as unknown in the original description. Cabrera (1958) restricted the type locality to Lagoa Sanata"="LA",
                               "New Hebrides (= Vanatu)"="Oceania",
                               "\"Am\216rique m\216ridionale.\""="LA",
                               "R\216union Isl (France)."="Africa",
                               "Duke of York Isl"="Oceania",
                               "Madeira Isls"="Europe",
                               "Guinea-Bissau"="Africa",
                               "Moluccas"="Africa",
                               "New Caledonia and Loyalty Isls."="Oceania",
                               "New Caledonia (France)"="Oceania",
                               "Natal"="Africa",
                               "Not definitely identifiable"="Africa"))
unique(bats$region)

## just region
bats=bats[c("species","region")]

## get phylogeny
bphy=tnrs_match_names(names=bats$species,context_name="Animals",do_approximate_matching=F)

## fuck you R, do this 1116/250 times
bs=seq(1,nrow(bats),by=250)

## get sets of names
bphy1=tnrs_match_names(names=bats$species[bs[1]:(bs[2]-1)],
                      context_name="Animals",do_approximate_matching=T)
bphy2=tnrs_match_names(names=bats$species[bs[2]:(bs[3]-1)],
                       context_name="Animals",do_approximate_matching=T)
bphy3=tnrs_match_names(names=bats$species[bs[3]:(bs[4]-1)],
                       context_name="Animals",do_approximate_matching=T)
bphy4=tnrs_match_names(names=bats$species[bs[4]:(bs[5]-1)],
                       context_name="Animals",do_approximate_matching=T)
bphy5=tnrs_match_names(names=bats$species[bs[5]:nrow(bats)],
                       context_name="Animals",do_approximate_matching=T)

## merge together
bphy=rbind.data.frame(bphy1,bphy2,bphy3,bphy4,bphy5)

## bind with bats
bphy=data.frame(bphy,bats)

## order by species
bphy=bphy[order(bphy$species),]

## get tree of all bats from rotl
btree=tol_induced_subtree(ott_ids=bphy$ott_id)

## clean
rm(bphy1,bphy2,bphy3,bphy4,bphy5)
# ## remove bad otts
# bphy=bphy[-which(bphy$ott_id=="687487"),]
# bphy=bphy[-which(bphy$ott_id=="1034143"),]
# bphy=bphy[-which(bphy$ott_id=="756131"),]
# bphy=bphy[-which(bphy$ott_id=="513426"),]
# bphy=bphy[-which(bphy$ott_id=="1055363"),]
# bphy=bphy[-which(bphy$ott_id=="1080286"),]
# bphy=bphy[-which(bphy$ott_id=="342715"),]
# bphy=bphy[-which(bphy$ott_id=="238425"),]

## convert names to match
bphy$unique_name=gsub(" ","_",bphy$unique_name)

## which names in our set don't match Wilson & Reeder?
nrow(bphy)
nrow(sdata)
sdata$unique_name[!sdata$unique_name%in%bphy$unique_name]
length(sdata$unique_name[!sdata$unique_name%in%bphy$unique_name])

## assign surveyed
sdata$surveyed=1

## standardize unique names
sdata$umatch=tolower(sdata$unique_name)
bphy$umatch=tolower(bphy$unique_name)

## match bat species
ourbats=sdata$umatch
table(ourbats%in%bphy$umatch)

## get those regions
oreg=bphy[bphy$umatch%in%ourbats,]
oreg=oreg[c("umatch","region")]

## remove duplicates
oreg=oreg[!duplicated(oreg$umatch),]

## sort by species then merge
oreg=oreg[order(oreg$umatch),]
sdata$region=oreg[,"region"]

## subtract our species from bphy
bphy=bphy[!bphy$umatch%in%ourbats,]

## add appropriate columns to make life easy
bphy$species=bphy$search_string

## add blank columns
st=which(names(sdata)=="species")+1
end=which(names(sdata)=="positive")
blanks=c(names(sdata)[st:end])
#blanks=names(sdata)[st:(length(names(sdata))-1)]
for (i in blanks) {
  bphy[, i]=0
}
bphy$surveyed=0

## convert positive to NA
bphy$positive=NA
bphy$hnv_positive=NA
bphy$filo_positive=NA

## make sure sdata and bphy are same order of columns
bphy=bphy[names(sdata)]

## species lower in both
bphy$species=tolower(bphy$species)
sdata$species=tolower(sdata$species)

## join phy and bphy
batphy=rbind.data.frame(sdata,bphy)

## add underdash for species
batphy$tree_species=batphy$unique_name

## clean
batphy$umatch=NULL

## yes no for filo and hnv
batphy$filo_surv=ifelse(batphy$surveyed==0,NA,ifelse(is.na(batphy$filo_samps)==T,0,1))
batphy$hnv_surv=ifelse(batphy$surveyed==0,NA,ifelse(is.na(batphy$hnv_samps)==T,0,1))
batphy$both_surv=ifelse(batphy$surveyed==0,NA,ifelse(batphy$filo_surv==1 & batphy$hnv_surv==1,1,0))

## export batphy
#setwd("~/Dropbox (MSU projects)/Spillover postdoc/bat virus meta-analysis/bat phylogeny")
setwd("~/Desktop/BZDEL/meta_analysis/R/revised analyses")
write.csv(batphy,"batphy_for_rotl.csv")

## as Rdata
save(batphy,file="batphy_revised.Rdata")

## get tree with otts
bat_tree=tol_induced_subtree(ott_ids=batphy$ott_id)

## remove ott information from the tips
bat_tree$tip.label=strip_ott_ids(bat_tree$tip.label)

## assign branch lengths
bat_tree=compute.brlen(bat_tree,method="Grafen")

## sort batphy into same order as bat_tree
bat_tree=makeLabel(bat_tree)
batphy=batphy[match(bat_tree$tip.label,batphy$tree_species),]
batphy=batphy[!is.na(batphy$tree_species),]

## give row names
rownames(batphy)=batphy$tree_species

## use caper to join (takes a while)
library(caper)
cdata=comparative.data(phy=bat_tree,data=batphy,names.col=tree_species,vcv=T,na.omit=F,warn.dropped=T)

