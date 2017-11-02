
setwd("~/Documents/MPI/ClimateAndLanguage/CHILDES/")

d = read.delim("Data/toneCounts.tab",sep='\t',header=T)

d = d[d$corpus!='PaidoCantonese',]

d$date2 = as.Date(d$date,format="%d-%b-%Y")

d$month = as.numeric(format(d$date2, "%m"))
d$summer = d$month %in% c(7:9)
d$season = cut(d$month, c(0,6,9,13))

d = d[!is.na(d$date2),]

d = d[d$language=="yue , eng",]

h = read.delim("Data/SpecificHumidity_HongKong.tsv", skip=1)
head(h)

h$month = ceiling(h$months.since.1960.01.01 %%12)
meanHumidity = tapply(h$unitless, h$month, mean)


d[,c("T1p","T2p",'T3p','T4p','T5p','T6p')] = t(apply(d[,c("T1","T2",'T3','T4','T5','T6')],1,prop.table))

d$prop.contour = d$T1p + d$T2p + d$T4p + d$T5p

d$humidity = meanHumidity[d$month]

library(ggplot2)

ggplot(d, aes(x=as.factor(month), y = T1p)) + geom_boxplot()


ggplot(d, aes(x=humidity, y = prop.contour)) + geom_point() + geom_smooth()

library(lme4)

m0 = lmer(log(prop.contour) ~ 1 + (1|corpus), data=d)
m1 = lmer(log(prop.contour) ~ humidity + (1|corpus), data=d)
anova(m0,m1)

summary(m1)

plot(1:12,ylim=c(0,0.5), xlab="Month", ylab="Proportion of use")
i = 1
for(x in c("T1p","T2p",'T3p','T4p','T5p','T6p')){
	points(tapply(d[,x],d$month,mean, na.rm=T),col=i, type='l')
	i = i +1
	
}
legend(1,0.5,legend=paste("Tone",1:6),col=1:6,pch=1, ncol=2)


tones = matrix(1:12,ncol=12)
for(i in c("T1","T2",'T3','T4','T5','T6')){
	tones = rbind(tones,tapply(d[,i],d$month,sum,na.rm=T))
}
tones = tones[2:nrow(tones),]
tones2 = prop.table(tones,margin=1)
plot(1:12,ylim=c(0,0.2))
for(i in 1:6){
	lines(tones2[i,],col=i)
}
legend(1,0.2,legend=1:6,col=1:6,pch=1)


