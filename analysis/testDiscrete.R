library(ape)
setwd("~/Documents/MPI/ClimateAndLanguage/PHOIBLE_Replication/analysis/")

source("fitPagel.r")

runDiscrete = function(tree,x,y, iQ, dQ){
  
  # fit discrete
  fit.ape = NA
  try(fit.ape<-fitPagel(tree,x,y,iQ=iQ,dQ=dQ,method='ace'))
  return(list(list(fit.ape,table(x,y),qx)))
}


# single-rate model
iQ<-matrix(c(0,1,2,0,3,0,0,2,4,0,0,1,0,4,3,0),4,4,byrow=TRUE)

# differet rates for complex -> simple in humid and dry conditions
dQ =matrix(c(0,1,2,0,3,0,0,5,4,0,0,1,0,4,3,0),4,4,byrow=TRUE) 


tree = read.nexus("../data/grollemund_et_al2015-d-place.NEXUS")
tree.t = read.csv("../data/grollemund_et_al2015-d-place_Taxa.csv", stringsAsFactors = F)

tree$tip.label = tree.t[match(tree$tip.label, tree.t$taxon),]$glottocode

p = read.csv("../data/phoibleTonesAndHumidity.csv", stringsAsFactors = F)
load('../../Phylogenetics/GlottologTree/ANU_Data.rDat')
g = read.csv("../data/glottolog-languoid.csv/languoid.csv", stringsAsFactors = F)
anu$glottocode = g[match(anu$iso, g$hid),]$id
anu = anu[!is.na(anu$glottocode),]
anu = anu[!anu$glottocode %in% p$Glottocode,]
names(anu)[names(anu) %in% c("Number.of.tones","glottocode")] = c("Tones","Glottocode")
comb = rbind(anu[,c("Glottocode","Tones",'specH.mean')],
      p[,c("Glottocode","Tones",'specH.mean')])

tree = drop.tip(tree, tree$tip.label[!tree$tip.label %in% comb$Glottocode])

tones = comb[match(tree$tip.label, comb$Glottocode),]$Tones
names(tones) = comb[match(tree$tip.label, comb$Glottocode),]$Glottocode
humidity = comb[match(tree$tip.label, comb$Glottocode),]$specH.mean
names(humidity) = comb[match(tree$tip.label, comb$Glottocode),]$Glottocode

tones.bin = c("simple","complex")[1+as.numeric(tones>=3)]
names(tones.bin) = names(tones)

# Estimate quantiles from whole distribution
quantilesToRun = seq(from=0.05,to=1,by=0.05)
quantiles = quantile(comb$specH.mean, quantilesToRun)

res = list()

for(qx in quantiles){
    # set dry/humid by quantile
    humidity.bin = c('humid','dry')[1+as.numeric(humidity < qx)]
    names(humidity.bin) = names(humidity)
    
    # run discrete
    resx = runDiscrete(tree,tones.bin, humidity.bin,
                      iQ,dQ)
    res[[length(res)+1]] = resx
    
}

par(mfrow=c(1,2))
DrawRateGraph2(resx[[1]][[1]], independent = T)
title("Independent", line = -1)
DrawRateGraph2(resx[[1]][[1]])
title("Dependent", line = -1)
