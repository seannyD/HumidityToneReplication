# The data-wrangling part of this script kindly supplied by Marton Soskuthy

library(dplyr)
library(sp)
library(fields)

hz2bark = function(hz){
  (26.81/(1+(1960/hz))) - 0.53
}

setwd("~/Documents/MPI/ClimateAndLanguage/PHOIBLE_Replication/processing/")
vs <- read.csv("../data/BeckerVowelCorpus.csv", fileEncoding="UTF-8", header=T, 
                 quote='"', sep=",",stringsAsFactors=F)
# for now: manually removing erroneous Hungarian entry...
vs <- vs[!(vs$ISO=="hun" & vs$Method=="Wordlist"),]

# keeping only one dialect / language
# - if there is a "Standard"/"standard"/"stadard"/"Stadard" keep that - not implemented yet!
# - otherwise:
#   + if there are multiple recording styles, preferentially choose
#     o isolation > wordlist > carrier sentence > meaningful sentence > running speech

dialects <- unique(vs[,c("ISO","Dialect","Method")])
dialects$Method <- factor(dialects$Method, 
                          levels=c("Isolation","Wordlist","Carrier sentence",
                                   "Meaningful sentence","Running speech")
)
dialects <- dialects[order(dialects$ISO,dialects$Dialect,dialects$Method),]
dialects <- dialects[!duplicated(dialects$ISO),]
dialects.to.keep <- paste(dialects$ISO, dialects$Dialect, dialects$Method)

vs$combined.label <- paste(vs$ISO, vs$Dialect, vs$Method)
vs <- vs[vs$combined.label %in% dialects.to.keep,]
vs <- vs[!duplicated(paste(vs$Language,vs$Quantity)),]

# only keeping unanimously long / short / uniform subsystems
unique(vs$Quantity)
vs$Quantity <- recode(vs$Quantity, 
                      `0`="Uniform",
                      Conbined="Combined",
                      long="Long")
#vs <- vs[vs$Quantity %in% c("Long","Short","Uniform"),]

# moving to long format
vs.long <- reshape(vs, varying=list(vowel=paste0("V",1:14,"OS"),
                                    ps=paste0("V",1:14,"PS"),
                                    f1=paste0("V",1:14,"F1"),
                                    f2=paste0("V",1:14,"F2"),
                                    f3=paste0("V",1:14,"F3"),
                                    f4=paste0("V",1:14,"F4")),
                   v.names=c("vowel","ps","f1","f2","f3","f4"),
                   timevar="vnum",direction="long")

# reorder according to language / dialect / subsystem (called quantity here)
vs.long <- vs.long[order(vs.long$ISO, vs.long$Dialect, vs.long$Quantity),]

# filtering out lines with no relevant info
vs.long <- vs.long[!is.na(vs.long$f1) & !vs.long$f1=="" & !vs.long$f2=="" & !vs.long$f1=="XXX",]

# converting to numeric cols

vs.long$f1 <- as.numeric(vs.long$f1)
vs.long$f2 <- as.numeric(vs.long$f2)

vs.long$f1.bark = hz2bark(vs.long$f1)
vs.long$f2.bark = hz2bark(vs.long$f2)

write.csv(vs.long,"../data/BeckerVowelCorpus_Long.csv")

plot(vs.long$f1, vs.long$f2)
plot(vs.long$f1.bark, vs.long$f2.bark)

####
d  = vs.long

g = read.csv("../data/glottolog-languoid.csv/languoid.csv")
g$family = g[match(g$family_pk, g$pk),]$name

d$glotto = g[match(d$ISO,g$hid),]$id

d$ISO[is.na(d$glotto)]

d[d$Language=='Akuntsu',]$glotto ='kano1245'
d[d$ISO=='bnh',]$glotto = 'jama1261'
d[d$ISO=='ori',]$glotto = 'oriy1255'
d[d$ISO=='que',]$glotto = 'cald1236'
d[d$ISO=='roa',]$glotto = 'norm1245'

d = d[!is.na(d$glotto),]

d$glotto = as.character(d$glotto)
d$family = g[match(d$glotto,g$id),]$family
d$family = as.character(d$family)
d$family[is.na(d$family)]  = d$glotto[is.na(d$family)]

d$latitude = g[match(d$glotto,g$id),]$latitude
d$longitude = g[match(d$glotto,g$id),]$longitude


d = d[!is.na(d$latitude),]

####
# Autotyp area

load("~/Documents/MPI/Neandertals_Collab/fossil/autotyp.data/autotyp.backbone.rda")
autotyp.backbone$glottolog_LID.2014 =as.character(autotyp.backbone$glottolog_LID.2014)
autotyp.backbone$area =as.character(autotyp.backbone$area)

# Find actual matches by glottoID
d$autotyp.area = autotyp.backbone[match(d$glotto,autotyp.backbone$glottolog_LID.2014),]$area

# split into found and missing
atFound = as.matrix(autotyp.backbone[,c("longitude","latitude")])
atMissing = as.matrix(d[is.na(d$autotyp.area),c("longitude","latitude")])

# For missing, find area of closest langauge
atDist = rdist.earth(atFound,atMissing)
closestArea = autotyp.backbone$area[apply(atDist,2,function(X){which(X==min(X,na.rm=T),arr.ind=T)[1]})]

d[is.na(d$autotyp.area),]$autotyp.area = closestArea




#####
# Link humidity data
h.mean = read.csv("../data/MeanMonthlyMeans.csv")
h.mean = as.matrix(h.mean[,c(2:ncol(h.mean))])

cols = 190
rows = 94

d$latitude.grid = as.numeric(cut(d$latitude,seq(90,-90,length.out=rows)))
d$longitude.grid = as.numeric(cut(d$longitude,seq(-180,180,length.out=cols)))

d$specH.mean = apply(cbind(d$latitude.grid,d$longitude.grid),1,function(X){h.mean[X[2],X[1]]})



######
# Calculate area of vowels
p.area = data.frame(glotto = unique(d$glotto), 
                    Language = unique(d$Language), 
                    area = NA, numVowels = NA)

p.area$specH.mean = d[match(p.area$glotto,d$glotto),]$specH.mean
p.area$autotyp.area = d[match(p.area$glotto,d$glotto),]$autotyp.area
p.area$family = d[match(p.area$glotto,d$glotto),]$family

pdf("../results/VowelPolygons.pdf")
plot(range(d$f1), range(d$f2), type='n',
     xlab="F1 (Hz)", ylab="F2 (Hz)")
cols = rainbow(length(unique(d$glotto)))
cols.i = 0
for(glotto in unique(d$glotto)){
  # get convex hull of the vowel measures
  dx = d[d$glotto==glotto,]
  ch = chull(dx$f1, dx$f2)
  ch = c(ch,ch[1])
  ch.points = cbind(dx$f1, dx$f2)[ch,]
  polygon(ch.points[,1],ch.points[,2],
          border = cols[cols.i])
  cols.i = cols.i + 1
  # calculate area by
  px = Polygon(ch.points)
  
  p.area$area[p.area$glotto==glotto] = px@area
  p.area$numVowels[p.area$glotto==glotto] = nrow(dx)
}
dev.off()

x = tapply(p.area$area, p.area$family, length)>5
cfams = names(x[x])
cx = rep(NA,length(cfams))
names(cx) = cfams
for(f in cfams){
  #cx[f] = cor(p.area[p.area$family==f,]$area,p.area[p.area$family==f,]$specH.mean)
  x = lm(p.area[p.area$family==f,]$area~
       p.area[p.area$family==f,]$specH.mean +
       p.area[p.area$family==f,]$numVowels
     )
  cx[f] = summary(x)$coefficients[2,3]
}

sort(cx)



fx = "Austronesian"
plot(p.area[p.area$family==fx,]$specH.mean,
     p.area[p.area$family==fx,]$area,
     xlim=range(p.area$specH.mean))
text(p.area[p.area$family==fx,]$specH.mean,
     p.area[p.area$family==fx,]$area,
     p.area[p.area$family==fx,]$glotto)

lx = c("tach1250","amha1245",'soma1255')
lx = rev(c("maor1246","paic1239",'nias1242'))
lnamex = paste0(as.character(p.area[match(lx,p.area$glotto),]$Language),
               " (",
               round(p.area[match(lx,p.area$glotto),]$specH.mean,3),")")

ldx= p.area[p.area$glotto %in% lx,]
library(RColorBrewer)
cols = brewer.pal(3,'Set1')
#cols = paste0(cols,"4D")
pdf("../results/AustronesianVowelExample.pdf",
    width = 6, height = 5)
plot(range(d$f1), range(d$f2), type='n',
     xlab="F1 (Hz)", ylab="F2 (Hz)")

i =1
for(l in lx){

dx = d[d$glotto==l,]
ch = chull(dx$f1, dx$f2)
ch = c(ch,ch[1])
ch.points = cbind(dx$f1, dx$f2)[ch,]
polygon(ch.points[,1],ch.points[,2],
        border = cols[i], lwd=3)
points(dx$f1, dx$f2, col=cols[i],pch=(15:17)[i])
i = i +1
}
legend(825, 2500, lnamex, lty=1, lwd=2, col=cols, pch=15:17, bty='n')
dev.off()

write.csv(p.area, "../data/BeckerVowelCorpus_Area.csv")
