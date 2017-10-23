## phylofactor figure script
library(rotl)
library(ape)

## make vector of species
specs=c("Tursiops truncatus",
        "Carcharodon carcharias",
        "Ambystoma maculatum",
        "Phocoena phocoena",
        "Homo sapiens",
        "Gorilla gorilla",
        "Balaenoptera musculus",
        "Macropus giganteus",
        "Varanus varius",
        "Thunnus maccoyii",
        "Anguilla rostrata",
        "Salmo trutta",
        "Hyla cinerea",
        "Zalophus californianus",
        "Halichoerus grypus",
        "Eunectes murinus")

## get phylo info
phylo=tnrs_match_names(names=specs,context_name="Animals",do_approximate_matching=T)

## tree
tree=tol_induced_subtree(ott_ids=phylo$ott_id)

## plot
par(oma=c(0,0,0,0),mar=c(0,0,0,0),mfrow=c(1,1))
plot(tree)

## remove ott information from the tips
tree$tip.label=strip_ott_ids(tree$tip.label)
plot(tree)
length(tree$tip.label)

## resolve multifurcations
tree=multi2di(tree)

## ladderize tree?
tree=ladderize(tree)
plot(tree,type="unrooted",label.offset=0.5,cex=0.75)
