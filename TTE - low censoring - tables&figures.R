##########################################################################
#                                                                        #
# Script for TTE outcome simulations (low censoring)  - TABLES AND PLOTS #
#                                                                        #
##########################################################################

library(tidyverse)
library(dplyr)
library(xtable)
library(ggpattern)
library(ggpubr)
theme_set(theme_classic())


#LOW CENSORING
load("results.low.over.RData")
load("results.low.well.RData")
load("results.low.under.RData")

c1<-results.low.under
c2<-results.low.well
c3<-results.low.over

data1<-rbind(c1,c2,c3)

x<-data1[c(1,5,9,13,17,21),1:6]
x<-as.vector(t(x))

y<-c(rep("hazard underestimated", 12), rep("hazard well estimated",12) , rep("hazard overestimated",12))
y<-as.factor(y)
yf<-factor(y,levels=c("hazard underestimated", "hazard well estimated","hazard overestimated" ))

z<-rep(c("E","SS"),18)
w<-rep(c(rep(0.9,2),rep(0.7,2),rep(0.5,2)),2)
v<-rep(c(rep("interim",6), rep("overall",6)),3)

data<-data.frame(power=x, assumedH=yf,strategy=z, trueRR=w, est=v)

power<-ggplot(data,aes(x=trueRR,y=power, linetype = strategy))+geom_point(aes(colour=strategy))+
  geom_line(aes(colour=strategy))+
  facet_grid(est ~assumedH)+scale_x_continuous(breaks=c(0.5,0.7,0.9))+
  xlab("true HR")+labs(colour = "IA strategy", linetype = "IA strategy")+
  theme(panel.grid.major.y = element_line(color = "lightgrey",size = 0.25,linetype = 2))+
  theme(text = element_text(size = 15))

ggsave(power,file="power.tte.low.png", height = 5, width = 8.889,dpi=400)

#SAMPLE SIZE 

x<-data1[c(4,8,12,16,20,24),1:6]
x<-as.vector(t(x))

y<-c(rep("hazard underestimated", 12), rep("hazard well estimated",12) , rep("hazard overestimated",12))
y<-as.factor(y)
yf<-factor(y,levels=c("hazard underestimated", "hazard well estimated","hazard overestimated" ))

z<-rep(c("E","SS"),18)
w<-rep(c(rep(0.9,2),rep(0.7,2),rep(0.5,2)),2)
v<-rep(c(rep("interim",6), rep("overall",6)),3)

data<-data.frame(ss=x, assumedH=yf,strategy=z, trueHR=w, est=v)

ss<-ggplot(data,aes(x=trueHR,y=ss,fill=strategy))+
  geom_bar_pattern(aes(pattern=strategy),stat='identity',colour='black',position="dodge")+
  facet_grid(est ~assumedH)+
  scale_x_continuous(breaks=c(0.5,0.7,0.9))+
  xlab("true HR")+
  ylab("expected sample size")+
  labs(fill = "IA strategy")+
  theme(panel.grid.major.y = element_line(color = "lightgrey",size = 0.3, linetype = 2),
        panel.grid.minor.y = element_line(color = "lightgrey",linetype = 2))+
  theme(text = element_text(size = 15))+
  scale_pattern_manual(name="IA strategy",values = c("none", "stripe"))

ggsave(ss,file="ss.tte.low.png", height = 5, width = 8.889,dpi=400)

figure <- ggarrange(power, ss,
                    labels = c("(a)", "(b)"),
                    ncol = 1, nrow = 2)

ggsave(figure,file="tte.plots.low.png", height = 10, width = 8.889,dpi=400)

data1<-data1[c(5,6,3,4,1,2,7,8)]
xtable(data1)

