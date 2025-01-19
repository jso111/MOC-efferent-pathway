#title:  analysis of Yuki's double-cross alpha9/vglut3 data
#author: John Oghalai
#date: 9/29/2023
#install.packages("tidyverse")

#setwd("C:\Users\jsogh\Dropbox (Personal)\Awake mouse paper\stats for yukis data")

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
#    print(x)
    dTsub<-filter(dT,select2==x)
    wt<-dTsub$select3[dTsub$select1==genotypes[1]]
    vGlut<-dTsub$select3[dTsub$select1==genotypes[3]]
    double<-dTsub$select3[dTsub$select1==genotypes[4]]
    if (sum(!is.nan(vGlut))>3&sum(!is.nan(double))>3){
      tFit1<-t.test(vGlut,double,paired = FALSE,na.rm=TRUE)
      if (tFit1$p.value<0.05){
        print(sprintf('x axis with a significant p-value vGLut vs double: %3.1f',x))
        print(tFit1)
        if (sum(!is.nan(wt))>3) {
          tFit2<-t.test(wt,double,paired = FALSE,na.rm=TRUE)
          print('control:wt vs double')
          print(tFit2)
        }
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

genotypes = c('WT','Alpha9KO','VGLUT3KO','VGLUT3KOAlpha9KO')

v1$experiment<-factor(v1$experiment)
v1$genotype<-factor(v1$genotype,levels=genotypes)
v1$freq<-v1$freq/1000
v1$phi<-v1$phi*2*pi
v2$experiment<-factor(v2$experiment)
v2$genotype<-factor(v2$genotype,levels=genotypes)
v2$freq<-v2$freq/1000
v3$experiment<-factor(v3$experiment)
v3$genotype<-factor(v3$genotype,levels=genotypes)
v4$experiment<-factor(v4$experiment)

###
# create table that fits the organization of the Linear Mixed-Effects Example.R function 
# and then call it to analyze the data statistically
###
print('')
print('')
print('')
degree=3

# gain vs genotype
dTable<-rename(v2,id=experiment, cohort=genotype)
dTable<-filter(dTable,freq<13)
dTable<-filter(dTable,freq>8)
filename='Yuki gain vs genotype.pdf'
fit<-lmer(gain ~ cohort * stats::poly(freq,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
#print(fit_anova) 
print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
if (last(fit_anova$`Pr(>F)`)<0.05) {
  dT=select(dTable,select=c(cohort,freq,gain))
  p=FALSE
  checkTtests(dT,p) 
}

# BF and Q10dB vs genotype
dTable<-rename(v3,id=experiment, cohort=genotype)
filename='Yuki BF vs genotype.pdf'
fit<-lmer(bf ~ cohort * stats::poly(level,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
#print(fit_anova) 
print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
if (last(fit_anova$`Pr(>F)`)<0.05) {
  dT=select(dTable,select=c(cohort,level,bf))
  p=FALSE
  checkTtests(dT,p) 
}

filename='Yuki Q10db vs genotype.pdf'
fit<-lmer(q ~ cohort * stats::poly(level,degree) + (1|id), data=dTable)
fit_anova<-anova(fit)
#print(fit_anova) 
print(sprintf('p-value: %5.5f',last(fit_anova$`Pr(>F)`)))
if (last(fit_anova$`Pr(>F)`)<0.05) {
  dT=select(dTable,select=c(cohort,level,q))
  p=FALSE
  checkTtests(dT,p) 
}

# calculate anova for figure parts D,E,F
filename='Yuki maxGain vs genotype.pdf'
dTab_controls<-filter(v4,(genotype==genotypes[1])|(genotype==genotypes[2])|(genotype==genotypes[4]))
fit<-aov(maxGain~genotype, data=v4)
#print(summary(fit)) 
x<-summary(fit)
print(sprintf('p-value: %5.5f',first(x[[1]]$`Pr(>F)`)))
if (first(x[[1]]$`Pr(>F)`)<0.05) {
  vGlut<-v4$maxGain[v4$genotype==genotypes[3]]
  double<-v4$maxGain[v4$genotype==genotypes[4]]
  
  tFit1<-t.test(vGlut,double,paired = FALSE,na.rm=TRUE)
  print(tFit1)
  tFit2<-aov(maxGain~genotype, data=dTab_controls)
  print(summary(tFit2))
}

filename='Yuki sensCF vs genotype.pdf'
fit<-aov(sensCF~genotype, data=v4)
#print(summary(fit)) 
x<-summary(fit)
print(sprintf('p-value: %5.5f',first(x[[1]]$`Pr(>F)`)))
if (first(x[[1]]$`Pr(>F)`)<0.05) {
  vGlut<-v4$sensCF[v4$genotype==genotypes[3]]
  double<-v4$sensCF[v4$genotype==genotypes[4]]
  
  tFit1<-t.test(vGlut,double,paired = FALSE,na.rm=TRUE)
  print(tFit1)
  tFit2<-aov(sensCF~genotype, data=dTab_controls)
  print(summary(tFit2))
}

filename='Yuki sens5k vs genotype.pdf'
fit<-aov(sens5k~genotype, data=v4)
#print(summary(fit)) 
x<-summary(fit)
print(sprintf('p-value: %5.5f',first(x[[1]]$`Pr(>F)`)))
if (first(x[[1]]$`Pr(>F)`)<0.05) {
  vGlut<-v4$sens5k[v4$genotype==genotypes[3]]
  double<-v4$sens5k[v4$genotype==genotypes[4]]
  
  tFit1<-t.test(vGlut,double,paired = FALSE,na.rm=TRUE)
  print(tFit1)
  tFit2<-aov(sens5k~genotype, data=dTab_controls)
  print(summary(tFit2))
}





# tuning curves vs genotype
degree=3
dTab<-rename(v1,id=experiment, cohort=genotype,phase=phi)
dTab<-filter(dTab,freq<10)
dTab<-filter(dTab,freq>4)
dTab<-filter(dTab,level>5)
dTable<-filter(dTab,(cohort==genotypes[3])|(cohort==genotypes[4]))
filename='YukiTC Vglut3 vs doubles.pdf'
compareTuningCurves(dTable,filename,degree)

dTable<-filter(dTab,(cohort==genotypes[1])|(cohort==genotypes[2])|(cohort==genotypes[4]))
filename='YukiTC wt vs alpha9 vs doubles.pdf'
compareTuningCurves(dTable,filename,degree)

# dTable<-filter(dTab,(cohort==genotypes[1])|(cohort==genotypes[4]))
# filename='YukiTC wt vs doubles.pdf'
# compareTuningCurves(dTable,filename,degree)
# 
# dTable<-filter(dTab,(cohort==genotypes[2])|(cohort==genotypes[4]))
# filename='YukiTC alpha9 vs doubles.pdf'
# compareTuningCurves(dTable,filename,degree)



