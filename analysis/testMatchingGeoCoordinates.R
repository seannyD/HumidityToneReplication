library(fields)
library(ggplot2)
setwd("~/Documents/MPI/ClimateAndLanguage/PHOIBLE_Replication/analysis/")
anu = read.csv("~/Desktop/Stuff/Everett/ANU_numTones_SpecificHumidity_GlottoFams_utf8.csv", stringsAsFactors = F)

glotto = read.csv("../data/glottolog-languoid.csv/languoid.csv", stringsAsFactors = F)

anu$glotto = glotto[match(anu$iso,glotto$hid),]$id

anu$glottoLat = glotto[match(anu$glotto,glotto$id),]$latitude
anu$glottoLong = glotto[match(anu$glotto,glotto$id),]$longitude


anu = anu[complete.cases(anu[,c("Latitude",'glottoLat')]),]

dist = rdist.earth(cbind(anu$Longitude, anu$Latitude), cbind(anu$glottoLong, anu$glottoLat), miles=F)
dist = diag(dist)
names(dist) = anu$glotto


ggplot(data.frame(dist = dist), aes(dist)) + geom_histogram() + scale_x_log10()

# number of points that differ by 100km
sum(dist>1000)

# 95% of comparisons are within:
quantile(dist,0.95)

err = anu[anu$glotto %in% names(dist[dist>1000]),]

# proportion of potentially problematic languages
nrow(err) / nrow(anu)

#List of mismatches
err[,c("glottoLong",'Longitude')]

# draw the mismatched languages
# ANU language points are in red
map(interior = F)
for(i in 1:nrow(err)){
  dx = anu[i,]
  points(x = dx$Longitude, dx$Latitude, pch=1, col=2)
  points(x = dx$glottoLong, dx$glottoLat, pch=1, col=3)
  lines(c(dx$Longitude, dx$glottoLong), y = c(dx$Latitude, dx$glottoLat),
         col = 2)
}
