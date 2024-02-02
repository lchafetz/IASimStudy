#####################################################
#                                                   #
# Script for tables and figures for binary outcomes #
#                                                   #
#####################################################

library(tidyverse)
library(dplyr)
library(xtable)
library(ggpattern)
library(ggpubr)
theme_set(theme_classic())

rm(list=ls())
load("bin_0.1_E.Rdata")
load("bin_0.1_SS.Rdata")
load("bin_0.2_E.Rdata")
load("bin_0.2_SS.Rdata")
load("bin_0.3_E.Rdata")
load("bin_0.3_SS.Rdata")

results_0.1 <- NULL
for (i in 1:4) {
  column_E <- results_0.1_E[, i]
  column_SS <- results_0.1_SS[, i]
  results_0.1 <- cbind(results_0.1, column_E, column_SS)
  results_0.1<-as.data.frame(results_0.1)
}

results_0.2<- NULL
for (i in 1:4) {
  column_E <- results_0.2_E[, i]
  column_SS <- results_0.2_SS[, i]
  results_0.2 <- cbind(results_0.2, column_E, column_SS)
  results_0.2<-as.data.frame(results_0.2)
}

results_0.3<- NULL
for (i in 1:4) {
  column_E <- results_0.3_E[, i]
  column_SS <- results_0.3_SS[, i]
  results_0.3 <- cbind(results_0.3, column_E, column_SS)
  results_0.3<-as.data.frame(results_0.3)
}

all_results<-rbind(results_0.1, results_0.2, results_0.3) %>% as.data.frame(.)
colnames(all_results)=c("E1","E2","E3","E4","S1","S2","S3","S4")
rownames(all_results)<-c("P_1","N_1","D_1","P_2","N_2","D_2",
                         "P1","N1","D1","P2","N2","D2",
                         "P.1","N.1","D.1","P.2","N.2","D.2")


# Probability of stopping plot
x<-all_results[c(7,10,13,16,1,4),c(1:6)]
x<-as.vector(t(x))
y<-c(rep("risk well estimated", 12), rep("risk overestimated",12) , rep("risk underestimated",12))
y<-as.factor(y)
yf<-factor(y,levels=c("risk underestimated", "risk well estimated","risk overestimated" ))
z<-rep(c("E","SS"),18)
w<-rep(c(rep(0.9,2),rep(0.7,2),rep(0.5,2)),2)
v<-rep(c(rep("interim",6), rep("overall",6)),3)
data<-data.frame(power=x, assumedH=yf,strategy=z, trueRR=w, est=v)

power.bin<-ggplot(data,aes(x=trueRR,y=power, color = strategy, linetype = strategy, shape = strategy))+geom_point()+
  geom_line()+facet_grid(est ~assumedH)+scale_x_continuous(breaks=c(0.5,0.7,0.9))+
  xlab("true RR")+labs(colour = "IA strategy", linetype = "IA strategy", shape = "IA strategy")+theme(panel.grid.major.y = element_line(color = "lightgrey",
                                                size = 0.25,linetype = 2))+theme(text = element_text(size = 15))
power.bin
ggsave(power.bin,file="power.bin.png", height = 5, width = 8.889,dpi=400)

# Sample size plot
x<-all_results[c(2,5,8,11,14,17),c(1:6)]
x<-as.vector(t(x))
y<-c(rep("risk underestimated", 12), rep("risk well estimated",12) , rep("risk overestimated",12))
y<-as.factor(y)
yf<-factor(y,levels=c("risk underestimated", "risk well estimated","risk overestimated" ))
z<-rep(c("E","SS"),18)
w<-rep(c(rep(0.9,2),rep(0.7,2),rep(0.5,2)),2)
v<-rep(c(rep("interim",6), rep("overall",6)),3)
data<-data.frame(ss=x, assumedH=yf,strategy=z, trueHR=w, est=v)

ss.bin<-ggplot(data,aes(x=trueHR,y=ss, fill=strategy))+geom_bar_pattern(aes(pattern=strategy),stat='identity',colour='black',position="dodge")+
  facet_grid(est ~assumedH)+scale_x_continuous(breaks=c(0.5,0.7,0.9))+
  xlab("true RR")+ylab("expected sample size")+
  labs(fill = "IA strategy")+
  theme(panel.grid.major.y = element_line(color = "lightgrey",size = 0.3,linetype = 2),
        panel.grid.minor.y = element_line(color = "lightgrey",size = 0.25,linetype = 2))+
  theme(text = element_text(size = 15))+
  scale_pattern_manual(name="IA strategy",values = c("none", "stripe"))

ss.bin
ggsave(ss.bin,file="ss.bin.png", height = 5, width = 8.889,dpi=400)


figure <- ggarrange(power.bin, ss.bin,
                    labels = c("(a)", "(b)"),
                    ncol = 1, nrow = 2)

ggsave(figure,file="bin.plots.png", height = 10, width = 8.889,dpi=400)

all_results<-all_results[c(5,6,3,4,1,2,7,8)]
all_results
xtable(all_results)

