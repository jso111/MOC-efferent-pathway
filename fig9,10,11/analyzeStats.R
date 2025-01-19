#title: Awake mouse analysis of Ari Cobo-Cuan's data
#author: John Oghalai
#date: 9/22/2023
#install.packages("tidyverse")

#setwd("C:/Users/jsogh/Dropbox (Personal)/Cobo-Cuan, Ariadna")

source('compareTuningCurves.R')
library(gsignal)
library(ggplot2)
library(gghighlight)
library(ggsignif)
library(gridExtra)
library(ggpubr)
library(plotrix)
library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)

#####
# Functions
#####
checkTtests<-function(dT,p) {
  # go through each freq and do t-tests on every combination in case the overall p-value is significant
  # and you want to look freq by freq
  xSteps=unique(dT$select2)
  count1=0
  count2=0
  for (x in xSteps) {
    dTsub<-filter(dT,select2==x)
    if (colSums(!is.na(dTsub))[3]>3) {
      tFit1<-pairwise.t.test(dTsub$select3, dTsub$select1, paired=p, p.adj = "bonf")
      if (any(tFit1$p.value<0.05,na.rm = TRUE)){
        print(sprintf('x axis with a significant p-value: %3.1f',x))
        print(tFit1)
        count1=count1+1
      }
    }
  } 
  print(sprintf('Number of x axis steps with a significant p-value: %2.0f',count1))
}

#####
# Load data and clean it
#####
v1 <- tibble(read.csv('v1.txt'))
v2 <- tibble(read.csv('v2.txt'))
v3 <- tibble(read.csv('v3.txt'))
v4 <- tibble(read.csv('v4.txt'))
vBin <- tibble(read.csv('vBin.txt'))
v1AA <- tibble(read.csv('v1AA.txt'))
v2AA <- tibble(read.csv('v2AA.txt'))
v3AA <- tibble(read.csv('v3AA.txt'))

conditions = c("Awake","Anesth","Dead")
protocols = c("Multitone","Singletone","TM")
genotypes = c('WT','Alpha9KO')
pupilSizes=c('Small','Medium','Large')
c=1
p=1

v1$experiment<-factor(v1$experiment)
v1$condition<-factor(v1$condition,levels=conditions)
v1$protocol<-factor(v1$protocol,levels=protocols)
v1$genotype<-factor(v1$genotype,levels=genotypes)
v1$pupilSize<-factor(v1$pupilSize,levels=pupilSizes)
v1$freq<-v1$freq/1000
v1$phi<-v1$phi*2*pi
v2$experiment<-factor(v2$experiment)
v2$condition<-factor(v2$condition,levels=conditions)
v2$protocol<-factor(v2$protocol,levels=protocols)
v2$genotype<-factor(v2$genotype,levels=genotypes)
v2$pupilSize<-factor(v2$pupilSize,levels=pupilSizes)
v2$freq<-v2$freq/1000
v3$experiment<-factor(v3$experiment)
v3$condition<-factor(v3$condition,levels=conditions)
v3$protocol<-factor(v3$protocol,levels=protocols)
v3$genotype<-factor(v3$genotype,levels=genotypes)
v3$pupilSize<-factor(v3$pupilSize,levels=pupilSizes)
v4$experiment<-factor(v4$experiment)
v4$condition<-factor(v4$condition,levels=conditions)
v4$protocol<-factor(v4$protocol,levels=protocols)
v4$genotype<-factor(v4$genotype,levels=genotypes)
v4$pupilSize<-factor(v4$pupilSize,levels=pupilSizes)

vBin$experiment<-factor(vBin$experiment)
vBin$genotype<-factor(vBin$genotype,levels=genotypes)

v1AA$experiment<-factor(v1AA$experiment)
v1AA$condition<-factor(v1AA$condition,levels=conditions)
v1AA$protocol<-factor(v1AA$protocol,levels=protocols)
v1AA$genotype<-factor(v1AA$genotype,levels=genotypes)
v1AA$pupilSize<-factor(v1AA$pupilSize,levels=pupilSizes)
v1AA$freq<-v1AA$freq/1000
v1AA$phi<-v1AA$phi*2*pi
v2AA$experiment<-factor(v2AA$experiment)
v2AA$condition<-factor(v2AA$condition,levels=conditions)
v2AA$protocol<-factor(v2AA$protocol,levels=protocols)
v2AA$genotype<-factor(v2AA$genotype,levels=genotypes)
v2AA$pupilSize<-factor(v2AA$pupilSize,levels=pupilSizes)
v2AA$freq<-v2AA$freq/1000
v3AA$experiment<-factor(v3AA$experiment)
v3AA$condition<-factor(v3AA$condition,levels=conditions)
v3AA$protocol<-factor(v3AA$protocol,levels=protocols)
v3AA$genotype<-factor(v3AA$genotype,levels=genotypes)
v3AA$pupilSize<-factor(v3AA$pupilSize,levels=pupilSizes)

###
# create table that fits the organization of the Linear Mixed-Effects Example.R function 
# and then call it to analyze the data statistically
###
print('')
print('')
print('')
degree=3

# WT and alpha9 gain vs pupilsize
dTab<-rename(v2,id=experiment, cohort=pupilSize,gain=gain1)
dTab<-filter(dTab,condition==conditions[1])
dTab<-filter(dTab,protocol==protocols[1])
dTable<-filter(dTab,genotype==genotypes[1])
dTable<-subset(dTable, select = -c(condition,protocol,genotype,gain2))
dTable<-filter(dTable,freq>8)
filename='AriWT gain vs pupilsize.pdf'
fit<-lmer(gain ~ cohort * stats::poly(freq,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
#print(fit_anova) 
print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
if (last(fit_anova$`Pr(>F)`)<0.05) {
  dT=select(dTable,select=c(cohort,freq,gain))
  p=TRUE
  checkTtests(dT,p) 
}

dTable<-filter(dTab,genotype==genotypes[2])
dTable<-subset(dTable, select = -c(condition,protocol,genotype,gain2))
dTable<-filter(dTable,freq>8)
filename='AriAlpha9 gain vs pupilsize.pdf'
fit<-lmer(gain ~ cohort * stats::poly(freq,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
#print(fit_anova) 
print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
if (last(fit_anova$`Pr(>F)`)<0.05) {
  dT=select(dTable,select=c(cohort,freq,gain))
  p=TRUE
  checkTtests(dT,p) 
}

# WT and alpha9 BF and Q10dB vs pupilsize
dTab<-rename(v3,id=experiment, cohort=pupilSize)
dTab<-filter(dTab,condition==conditions[1])
dTab<-filter(dTab,protocol==protocols[1])
dTable<-filter(dTab,genotype==genotypes[1])
dTable<-subset(dTable, select = -c(condition,protocol,genotype))
dTable<-filter(dTable,level>15)
filename='AriWT BF vs pupilsize.pdf'
fit<-lmer(bf ~ cohort * stats::poly(level,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
#print(fit_anova) 
print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
if (last(fit_anova$`Pr(>F)`)<0.05) {
  dT=select(dTable,select=c(cohort,level,bf))
  p=TRUE
  checkTtests(dT,p) 
}

filename='AriWT Q10dB vs pupilsize.pdf'
fit<-lmer(q ~ cohort * stats::poly(level,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
#print(fit_anova) 
print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
if (last(fit_anova$`Pr(>F)`)<0.05) {
  dT=select(dTable,select=c(cohort,level,q))
  p=TRUE
  checkTtests(dT,p) 
}

dTable<-filter(dTab,genotype==genotypes[2])
dTable<-subset(dTable, select = -c(condition,protocol,genotype))
dTable<-filter(dTable,level>15)
filename='AriAlpha9 BF vs pupilsize.pdf'
fit<-lmer(bf ~ cohort * stats::poly(level,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
#print(fit_anova) 
print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
if (last(fit_anova$`Pr(>F)`)<0.05) {
  dT=select(dTable,select=c(cohort,level,bf))
  p=TRUE
  checkTtests(dT,p) 
}

filename='AriAlpha9 Q10dB vs pupilsize.pdf'
fit<-lmer(q ~ cohort * stats::poly(level,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
#print(fit_anova) 
print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`))) 
if (last(fit_anova$`Pr(>F)`)<0.05) {
  dT=select(dTable,select=c(cohort,level,q))
  p=TRUE
  checkTtests(dT,p) 
}


# WT and alpha9 gain vs AA
dTab<-rename(v2AA,id=experiment, cohort=condition,gain=gain1)
dTable<-filter(dTab,genotype==genotypes[1])
dTable<-subset(dTable, select = -c(protocol,pupilSize,genotype,gain2))
dTable<-filter(dTable,freq>8)
filename='AriWT gain vs AA'
fit<-lmer(gain ~ cohort * stats::poly(freq,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
#print(fit_anova) 
print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`))) 
if (last(fit_anova$`Pr(>F)`)<0.05) {
  dT=select(dTable,select=c(cohort,freq,gain))
  p=TRUE
  checkTtests(dT,p) 
}

dTable<-filter(dTab,genotype==genotypes[2])
dTable<-subset(dTable, select = -c(protocol,pupilSize,genotype,gain2))
dTable<-filter(dTable,freq>8)
filename='AriAlpha9 gain vs AA.pdf'
fit<-lmer(gain ~ cohort * stats::poly(freq,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
#print(fit_anova) 
print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
if (last(fit_anova$`Pr(>F)`)<0.05) {
  dT=select(dTable,select=c(cohort,freq,gain))
  p=TRUE
  checkTtests(dT,p) 
}


# WT and alpha9 BF and Q10dB vs AA
dTab<-rename(v3AA,id=experiment, cohort=condition)
dTable<-filter(dTab,genotype==genotypes[1])
dTable<-subset(dTable, select = -c(pupilSize,protocol,genotype))
dTable<-filter(dTable,level<80)
filename='AriWT BF vs AA.pdf'
fit<-lmer(bf ~ cohort * stats::poly(level,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
#print(fit_anova) 
print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
if (last(fit_anova$`Pr(>F)`)<0.05) {
  dT=select(dTable,select=c(cohort,level,bf))
  p=TRUE
  checkTtests(dT,p) 
}

filename='AriWT Q10dB vs AA.pdf'
fit<-lmer(q ~ cohort * stats::poly(level,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
#print(fit_anova) 
print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
if (last(fit_anova$`Pr(>F)`)<0.05) {
  dT=select(dTable,select=c(cohort,level,q))
  p=TRUE
  checkTtests(dT,p) 
}

dTable<-filter(dTab,genotype==genotypes[2])
dTable<-subset(dTable, select = -c(pupilSize,protocol,genotype))
dTable<-filter(dTable,level<80)
filename='AriAlpha9 BF vs AA.pdf'
fit<-lmer(bf ~ cohort * stats::poly(level,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
#print(fit_anova) 
print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
if (last(fit_anova$`Pr(>F)`)<0.05) {
  dT=select(dTable,select=c(cohort,level,bf))
  p=TRUE
  checkTtests(dT,p) 
}

filename='AriAlpha9 Q10dB vs AA.pdf'
fit<-lmer(q ~ cohort * stats::poly(level,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
#print(fit_anova) 
print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
if (last(fit_anova$`Pr(>F)`)<0.05) {
  dT=select(dTable,select=c(cohort,level,q))
  p=TRUE
  checkTtests(dT,p) 
}


# analyze binned pupil data
degree=3
dTab<-rename(vBin,id=experiment)
dTable<-dTab
filename='Ari_bin gain.pdf'
fit<-lmer(gain ~ genotype * stats::poly(pupilBin,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
print(fit_anova) 
#print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
filename='Ari_bin CF.pdf'
fit<-lmer(cf ~ genotype * stats::poly(pupilBin,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
print(fit_anova) 
#print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
filename='Ari_bin Q10dB.pdf'
fit<-lmer(q ~ genotype * stats::poly(pupilBin,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
print(fit_anova) 
#print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
filename='Ari_bin magBF.pdf'
fit<-lmer(mag ~ genotype * stats::poly(pupilBin,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
print(fit_anova) 
#print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
filename='Ari_bin maghalfBF.pdf'
fit<-lmer(maglow ~ genotype * stats::poly(pupilBin,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
print(fit_anova) 
#print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`))) 
filename='Ari_bin phaseCF.pdf'
fit<-lmer(phase ~ genotype * stats::poly(pupilBin,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
print(fit_anova) 
#print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))


# WT and alpha9 tuning curves vs pupilsize
dTab<-rename(v1,id=experiment, cohort=pupilSize,phase=phi)
dTab<-filter(dTab,condition==conditions[1])
dTab<-filter(dTab,protocol==protocols[1])
dTable<-filter(dTab,genotype==genotypes[1])
dTable<-subset(dTable, select = -c(condition,protocol,genotype))
filename='AriWT TC vs pupilsize.pdf'
compareTuningCurves(dTable,filename,degree)

dTable<-filter(dTab,genotype==genotypes[2])
dTable<-subset(dTable, select = -c(condition,protocol,genotype))
dTable<-filter(dTable,freq<13)
filename='AriAlpha9 TC vs pupilsize.pdf'
compareTuningCurves(dTable,filename,degree)

# WT and alpha9 tuning curves awake vs anesthetized
dTab<-rename(v1AA,id=experiment, cohort=condition,phase=phi)
dTable<-filter(dTab,genotype==genotypes[1])
dTable<-subset(dTable, select = -c(pupilSize,protocol,genotype))
filename='AriWT TC vs AA.pdf'
compareTuningCurves(dTable,filename,degree)

dTable<-filter(dTab,genotype==genotypes[2])
dTable<-subset(dTable, select = -c(pupilSize,protocol,genotype))
dTable<-filter(dTable,freq<13)
filename='AriAlpha9 TC vs AA.pdf'
compareTuningCurves(dTable,filename,degree)




