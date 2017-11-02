library(lme4)
library(ggplot2)
setwd("~/Documents/MPI/ClimateAndLanguage/PHOIBLE_Replication/analysis/")

d = read.csv("../data/ASJP_and_SlonimskaRoberts_VowelProportions.csv", stringsAsFactors = F)



ggplot(d[!is.na(d$asjp.VProp),], aes(x =asjp.VProp, y = S.VProp)) +
         xlab("ASJP Vowel Proportion") +
         ylab("IDS+ Vowel Proportion") +
         geom_point() +
         geom_smooth()

cor.test(d$asjp.VProp, d$S.VProp)

ggplot(d[!is.na(d$asjp.VProp),], aes(x =specH.mean, y = asjp.VProp)) +
  xlab("Humidity") +
  ylab("ASJP Vowel Proportion") +
  geom_point() +
  geom_smooth()

ggplot(d, aes(x =specH.mean, y = S.VProp)) +
  xlab("Humidity") +
  ylab("IDS+ Vowel Proportion") +
  geom_point() +
  geom_smooth()



par(mfrow = c(1,2))
hist(d$asjp.VProp, main="ASJP")
hist(d$S.VProp,main="IDS+")
par(mfrow=c(1,1))

d$asjp.VProp.scale = scale(d$asjp.VProp)
d$S.VProp.scale = scale(d$S.VProp)


m0 = lmer(asjp.VProp.scale~ 1 +
            (1  | family) +
            (1  | autotyp.area),
          data = d[!is.na(d$specH.mean),])

m1 = lmer(asjp.VProp.scale~ 1 +
            specH.mean + 
            (1  | family) +
            (1  | autotyp.area),
          data = d[!is.na(d$specH.mean),])

anova(m0,m1)

m0 = lmer(S.VProp.scale~ 1 +
          (1 + specH.mean | family) +
          (1 + specH.mean | autotyp.area),
            data = d,
          control = lmerControl(optimizer = 'bobyqa',
                                optCtrl = list(maxfun=50000)))

m1 = lmer(S.VProp.scale~ 1 +
            specH.mean + 
            (1 + specH.mean | family) +
            (1 + specH.mean | autotyp.area),
          data = d,
          control = lmerControl(optimizer = 'bobyqa',
                                optCtrl = list(maxfun=50000)))

anova(m0,m1)
