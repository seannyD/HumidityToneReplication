library(lme4)
library(sjPlot)
library(caret)
library(car)

setwd("~/Documents/MPI/ClimateAndLanguage/PHOIBLE_Replication/analysis/")

p = read.csv("../data/phoibleTonesAndHumidity.csv")
p = p[p$Tones!=1,]

# transform with box-cox

pp = preProcess(p[,c('Tones','specH.mean')], method="BoxCox")

p$specH.mean.center = bcPower(p$specH.mean, lambda = pp$bc$specH.mean$lambda)

p$specH.mean.center = scale(p$specH.mean.center)

hist(p$specH.mean.center)

# Test need for random effects

m0 = glmer(Tones~1  + (1|Family) + (1|autotyp.area), 
           data=p, family=poisson, control = glmerControl(optimizer = 'bobyqa'))

m1 = glmer(Tones~1  + (1|Family) + (1+specH.mean.center|autotyp.area), 
           data=p, family=poisson,
           control = glmerControl(optimizer = 'bobyqa'))

m2 = glmer(Tones~1  + (1+specH.mean.center||Family) + (1+specH.mean.center|autotyp.area), 
           data=p, family=poisson,
           control = glmerControl(optimizer = 'bobyqa'))

anova(m0,m1,m2)

m3 = glmer(Tones~1 +specH.mean.center  + (1+specH.mean.center||Family) + (1+specH.mean.center|autotyp.area), 
           data=p, family=poisson,
           control = glmerControl(optimizer = 'bobyqa'))

anova(m2,m3)
summary(m3)
sjp.glmer(m3, 'eff', vars=c("specH.mean.center"), show.scatter = T, show.ci = T)

x = ranef(m3)
x2 = x$autotyp.area$specH.mean.center
names(x2) = rownames(x$autotyp.area)
pdf("../results/PHOIBLE_ranef.pdf", height=5, width=4)
dotplot(sort(x2), xlab='Humidity Random Effect')
dev.off()

pdf("../results/PHOIBLE_Tones_raw.pdf", width=4, height=4)
ggplot(p, aes(x=as.factor(Tones), y = specH.mean)) + geom_boxplot() +
  xlab("Number of tones") + ylab("Specific Humidity")
dev.off()



#######
library(MCMCglmm)


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
            Tones ~ specH.mean.center,
               ~ idh(1 + specH.mean.center):Family +
                 us(1 + specH.mean.center):autotyp.area,
               data   = p,
               family = "poisson",
               prior  = prior.m3,
               thin   =     100,
               burnin =   10000,
               nitt   = 1010000,
               verbose = TRUE)

save(m3.mcmcglmm, file="../results/m3_mcmcglmm.RDat")
# load("../results/m3_mcmcglmm.RDat")
#plot(m3.mcmcglmm, random=F)


sx = summary(m3.mcmcglmm)

sx

fe = m3.mcmcglmm$Sol[,2]
#fe = sample(fe, 10000)

dx = density(fe)
plot(dx, main='', xlab='Parameter estimate')
abline(v=0)


re = m3.mcmcglmm$VCV
re = as.data.frame(re)
re.area = sample(re$`specH.mean.center:specH.mean.center.autotyp.area`,10000)
re.area.d = density(re.area)
plot(re.area.d)

