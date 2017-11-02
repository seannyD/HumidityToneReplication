library(ggplot2)
library(fields)
library(psych)
library(lme4)
setwd("~/Documents/MPI/ClimateAndLanguage/PHOIBLE_Replication/processing/")


# ANU

g = read.csv("../../Grollemund/data/glottolog/languages-and-dialects-geo.csv", stringsAsFactors = F)

a = read.csv("~/Desktop/Stuff/Everett/ANU_numTones_SpecificHumidity_GlottoFams_utf8.csv", stringsAsFactors = F, fileEncoding = 'utf-8')


a2 = a[,c("Language","OID_","iso","Latitude",'Longitude','Number.of.tones')]

a2 = a2[a2$iso!='',]
a2$Phonemes = NA
a2$Consonants = NA
a2$Vowels = NA
a2$Source = "ANU"
a2$Glottocode = g[match(a2$iso,g$isocodes),]$glottocode
names(a2) = c("Language","InventoryID","iso","Latitude",'Longitude','Tones','Phonemes','Consonants',"Vowels","Source","Glottocode")


# Phoible

p = read.csv("../data/phoibleTonesAndHumidity.csv")

## Combined

combined = p[p$Glottocode %in% a2$Glottocode,]
combined$ANU.Tones = a2[match(combined$Glottocode, a2$Glottocode),]$Tones

write.csv(combined, "../data/phoibleAndANU_combined.csv")

