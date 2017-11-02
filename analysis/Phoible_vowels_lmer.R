library(lme4)
library(sjPlot)
library(ggplot2)
library(caret)
library(car)

setwd("~/Documents/MPI/ClimateAndLanguage/PHOIBLE_Replication/analysis/")

p = read.csv("../data/phoibleVowelsAndHumidty.csv")

p = p[p$Trump==1,]

pp = preProcess(p[,c('Tones','specH.mean')], method="BoxCox")

p$specH.mean.center = bcPower(p$specH.mean, lambda = pp$bc$specH.mean$lambda)

p$specH.mean.center = scale(p$specH.mean.center)
p$specH.mean.center2 = p$specH.mean.center^2
hist(p$specH.mean)
hist(p$specH.mean.center)

#p$specH.mean.center = scale(p$specH.mean)
p$prop.vowels = p$Vowels/(p$Consonants + p$Vowels)
p$prop.vowels.scaled = scale(p$prop.vowels)

p$inventorySize = p$Consonants + p$Vowels
p$inventorySize = scale(p$inventorySize)

hist(p$prop.vowels.scaled)

m0 = lmer(prop.vowels.scaled ~ 1 + (1+specH.mean.center|Family) + (1+specH.mean.center|autotyp.area), data=p)

m1 = lmer(prop.vowels.scaled ~ specH.mean.center + (1+specH.mean.center|Family) + (1+specH.mean.center|autotyp.area), data=p)

m2 = lmer(prop.vowels.scaled ~ specH.mean.center + specH.mean.center2 + (1+specH.mean.center|Family) + (1+specH.mean.center|autotyp.area), data=p)

anova(m0,m1,m2)

summary(m1)


plot(p$prop.vowels, p$specH.mean)

pdf("../results/PropVowels_SpecH.pdf", width=4, height=4)
ggplot(p, aes(x=prop.vowels, y=specH.mean)) + geom_point() +stat_smooth() + xlab("Proportion of vowels") + ylab("Specific Humidity")
dev.off()


x = sjp.lmer(m1,'fe.slope', vars=c("specH.mean.center"), show.scatter = T, show.ci = T, prnt.plot = F) 

#x$plot.list[[1]]$data$x = x$plot.list[[1]]$data$x*
#    attr(p$specH.mean.center,'scaled:scale') + attr(p$specH.mean.center, 'scaled:center')

x$plot.list[[1]]$data$x = p[complete.cases(p[,c("specH.mean.center",'prop.vowels.scaled','Family','autotyp.area')]),]$specH.mean

x$plot.list[[1]]$data$y = x$plot.list[[1]]$data$y*
  attr(p$prop.vowels.scaled,'scaled:scale') + attr(p$prop.vowels.scaled, 'scaled:center')

pdf("../results/PropVowels_SpecH_Estimates.pdf", width=4, height=4)
x$plot.list[[1]] + xlab("Specific Humidity")+ ylab("Vowel Index")
dev.off()  


########

library(MCMCglmm)

p = p[complete.cases(p[,c("Family",'autotyp.area','prop.vowels.scaled','specH.mean.center')]),]

familyRandomEffectsN = 2
areaRandomEffectsN = 2

prior.m3 <- list(
  R=list(V=1, n=1, fix=1),
  G=list(G1=list(V        = diag(familyRandomEffectsN),  # family intercept+slope
                 n        = familyRandomEffectsN,
                 alpha.mu = rep(0, familyRandomEffectsN),
                 alpha.V  = diag(familyRandomEffectsN)*25^2),
         G2=list(V        = diag(areaRandomEffectsN), # area intercept+slope
                 n        = areaRandomEffectsN,
                 alpha.mu = rep(0, areaRandomEffectsN),
                 alpha.V  = diag(areaRandomEffectsN)*25^2)))

m3.mcmcglmm <- MCMCglmm(
  prop.vowels.scaled ~ specH.mean.center,
  ~ us(1 + specH.mean.center):Family +
    us(1 + specH.mean.center):autotyp.area,
  data   = p,
  family = "gaussian",
  prior  = prior.m3,
  thin   =     100,
  burnin =   10000,
  nitt   = 1010000,
  verbose = TRUE)

save(m3.mcmcglmm, file="../results/m3_mcmcglmm_vowels.RDat")
# load("../results/m3_mcmcglmm_vowels.RDat")
plot(m3.mcmcglmm, random=F)
sx = summary(m3.mcmcglmm)

sx

fe = m3.mcmcglmm$Sol[,2]
#fe = sample(fe, 10000)

dx = density(fe)
plot(dx, main='', xlab='Parameter estimate')
abline(v=0)


re = m3.mcmcglmm$VCV
re = as.data.frame(re)
re.area = re$`specH.mean.center:specH.mean.center.autotyp.area`
re.area.d = density(re.area)
plot(re.area.d)



# test total size of inventory

m3 = lmer(prop.vowels.scaled ~ specH.mean.center+ inventorySize + (1+specH.mean.center|Family) + (1+specH.mean.center|autotyp.area), data=p)

anova(m1,m3)

summary(m3)
# inventory size doesn't take away the effect