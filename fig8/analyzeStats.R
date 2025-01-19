#title: Awake mouse analysis of Michele Pei's data
#author: John Oghalai
#date: 9/22/2023
#install.packages("tidyverse")

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
# Load data and clean it
#####
v1 <- tibble(read.csv('v1.txt'))
v2 <- tibble(read.csv('v2.txt'))
v3 <- tibble(read.csv('v3.txt'))
v4 <- tibble(read.csv('v4.txt'))

conditions = c("Awake")
protocols = c("Singletone")
genotypes = c('WT','Alpha9KO')
pupilSizes=c('Small','Medium','Large')
c=1
p=1

# v1$experiment<-factor(v1$experiment)
v1$condition<-factor(v1$condition,levels=conditions)
v1$protocol<-factor(v1$protocol,levels=protocols)
v1$genotype<-factor(v1$genotype,levels=genotypes)
v1$pupilSize<-factor(v1$pupilSize,levels=pupilSizes)
v1$freq<-v1$freq/1000
v2$experiment<-factor(v2$experiment)
v2$condition<-factor(v2$condition,levels=conditions)
v2$protocol<-factor(v2$protocol,levels=protocols)
v2$genotype<-factor(v2$genotype,levels=genotypes)
v2$pupilSize<-factor(v2$pupilSize,levels=pupilSizes)
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


###
# create table that fits the organization of the Linear Mixed-Effects Example.R function
###
degree=3
dTab<-rename(v1,id=experiment, cohort=pupilSize,phase=phi)
dTable<-filter(dTab,genotype==genotypes[1])
dTable<-subset(dTable, select = -c(condition,protocol,genotype))
filename='PeiWT.pdf'
compareTuningCurves(dTable,filename,degree)

dTable<-filter(dTab,genotype==genotypes[2])
dTable<-subset(dTable, select = -c(condition,protocol,genotype))
filename='PeiAlpha9.pdf'
compareTuningCurves(dTable,filename,degree)



  