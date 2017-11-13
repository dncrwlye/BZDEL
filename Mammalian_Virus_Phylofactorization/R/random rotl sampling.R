## get all Chordata
phy=tnrs_match_names("Chordata")
tree=tol_subtree(ott_id=phy$ott_id[1])

## sample tips
set.seed(1)
tips=sample(tree$tip.label,50)

## strip ott id
tips=strip_ott_ids(tips)

## remove underdash
tips=strsplit(tips,"_")
tips=sapply(tips,function(x) paste(x[1],x[2]))

## search again
phy=tnrs_match_names(names=tips,context_name="Animals",do_approximate_matching=F)

## remove NA
phy=phy[!is.na(phy$ott_id),]

## get tree
tree=tol_induced_subtree(ott_ids=phy$ott_id)

## strip ott id
tree$tip.label=strip_ott_ids(tree$tip.label)

## plot nice a pretty
par(oma=c(0,0,0,0),mar=c(0.5,0.5,0.5,0.5))
plot(tree,cex=0.75)
