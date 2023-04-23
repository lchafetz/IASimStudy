###################################################
# 
# Script for TTE outcome simulations
# 
####################################################

library(rstanarm)
library(survival)
library(dplyr)

#set parameter values for log(HRs)
beta0<--1  
beta1_90<-log(0.9)
beta1_70<-log(0.7)
beta1_50<-log(0.5)
beta1_1<-log(1)

shape=1.5 #shape parameter

#Reparameterization
sigma0=1/(exp(beta0/shape))
sigma1_90=1/(exp((beta0+beta1_90)/shape))
sigma1_70=1/(exp((beta0+beta1_70)/shape))
sigma1_50=1/(exp((beta0+beta1_50)/shape))
sigma1_1=1/(exp((beta0+beta1_1)/shape))

sigma1_t<-c(sigma1_90,sigma1_70,sigma1_50,sigma1_1)

ta<-3 # accrual period
tf<-1.7 # follow-up period
cutoff<-ta+tf # administrative censoring time
factor<-1.3 # 1.3 for high censoring, change to 0.15 for low censoring

sup<-0.98 # efficacy stopping threshold
n=1018 # sample size for assumed beta_0 of -1 (change sample size depending on censoring factor and assumed beta_0 as defined in Table 4 of paper)
D=196 # number of expected events 
halfwayE<-D/2 # halfway point for strategy E
halfwaySS<-n/2 # halfway point for strategy SS
nsim<-1000 # number of replications 

T<-c(rep(c(0,1),n/2)) #Arm indicator
ID<-1:n

#Simulate entry times 
entry_times<-replicate(nsim,runif(n,0,ta),simplify=FALSE)
for(i in 1:nsim){
  entry_times[[i]]<-sort(entry_times[[i]]) # sort in chronological order 
}

entry_c<-vector("list")
for( i in 1:nsim){
  entry_c[[i]]<-entry_times[[i]][c(TRUE,FALSE)] # entry times for control arm
}

entry_t<-vector("list")
for( i in 1:nsim){
  entry_t[[i]]<-entry_times[[i]][c(FALSE,TRUE)] # entry times for treatment arm 
}

entry_times<-vector("list")
for(i in 1:nsim){
  entry_times[[i]]<-c(rbind(entry_c[[i]],entry_t[[i]]))
}

#Survival times control group
t_c<-replicate(nsim,rweibull(n/2,shape=shape,scale=sigma0),simplify = FALSE)

#Censoring times control group
cens_c<-replicate(nsim,rexp(n/2,factor),simplify=FALSE)

#Observed survival times in control group 
yt_c<-vector("list")
for( i in 1:nsim){
  yt_c[[i]]<-entry_c[[i]]+t_c[[i]]
}

#Observed censoring times in control group 
ycens_c<-vector("list")
for( i in 1:nsim){
  ycens_c[[i]]<-entry_c[[i]]+cens_c[[i]]
}

#Survival times treatment group
t_t<-vector("list",4)
for(p in 1:4){
  t_t[[p]]<-replicate(nsim, rweibull(n/2,shape=shape,scale=sigma1_t[p]),simplify = FALSE)
}

#Censoring times treatment group
cens_t<-vector("list",4)
for(p in 1:4){
  cens_t[[p]]<-replicate(nsim, rexp(n/2,factor),simplify = FALSE)
}

#Observed survival times in treatment group
yt_t<-vector("list",4)
for( p in 1:4){
  yt_t[[p]]<-vector("list",1000)
}

for(p in 1:4){
  for(i in 1:nsim){
    yt_t[[p]][[i]]<-entry_t[[i]]+t_t[[p]][[i]]
  }
}

#Merge control and treatment survival times
Survtimes<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
    Survtimes[[p]][[i]]<-c(rbind(t_c[[i]],t_t[[p]][[i]]))
  }
}

#Observed censoring times in treatment group
ycens_t<-vector("list",4)
for(p in 1:4){
  ycens_t[[p]]<-vector("list",1000)
}

for(p in 1:4){
  for(i in 1:nsim){
    ycens_t[[p]][[i]]<-entry_t[[i]]+cens_t[[p]][[i]]
  }
}

#Trial end date (administrative censoring)
censoring_time<-rep(cutoff,n/2)

#Observation times for control group = min(observed event time, observed censoring time, trial end date)
Obstimes_C<-vector("list",1000)
for(i in 1:nsim){
  Obstimes_C[[i]]<-pmin(yt_c[[i]],ycens_c[[i]], censoring_time)
}

#Observation times for treatment groups
Obstimes_T<-vector("list",4)

for(p in 1:4){
  for( i in 1:nsim)
    Obstimes_T[[p]][[i]]<-pmin(yt_t[[p]][[i]],ycens_t[[p]][[i]],censoring_time)
}

#Merge control and treatment observed times
Obstimes<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
    Obstimes[[p]][[i]]<-c(rbind(Obstimes_C[[i]],Obstimes_T[[p]][[i]]))
  }
}


#Create binary censoring indicator variable for the simulated datasets
Cens_c<-vector("list",4)

for(i in 1:nsim){
  Cens_c[[i]]<- ifelse(ycens_c[[i]]<yt_c[[i]] | censoring_time<yt_c[[i]] , 0, 1)
}

Cens_t<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
    Cens_t[[p]][[i]]<-ifelse(ycens_t[[p]][[i]]<yt_t[[p]][[i]] | censoring_time<yt_t[[p]][[i]], 0, 1)
  }
}

cens<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
    cens[[p]][[i]]<-c(rbind(Cens_c[[i]],Cens_t[[p]][[i]]))
  }
}

#Create datasets for analysis
dataset<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
    dataset[[p]][[i]]<-data.frame(trt=T,entry=entry_times[[i]], Surv_times=Survtimes[[p]][[i]], Obs_times=Obstimes[[p]][[i]],event=cens[[p]][[i]])
    dataset[[p]][[i]]<-arrange(dataset[[p]][[i]],Obs_times)
    dataset[[p]][[i]][,"cum_trt_pts"]<-cumsum(dataset[[p]][[i]]$trt)
    dataset[[p]][[i]][,"cum_ctrl_pts"]<-cumsum(1-dataset[[p]][[i]]$trt)
    dataset[[p]][[i]][,"total_cum_pts"]<-dataset[[p]][[i]]$cum_trt_pts+dataset[[p]][[i]]$cum_ctrl_pts
    dataset[[p]][[i]][,"trt_events"]<-ifelse(dataset[[p]][[i]]$trt==1 & dataset[[p]][[i]]$event==1,1,0)
    dataset[[p]][[i]][,"ctrl_events"]<-ifelse(dataset[[p]][[i]]$trt==0 & dataset[[p]][[i]]$event==1,1,0)
    dataset[[p]][[i]][,"cum_trt_events"]<-cumsum(dataset[[p]][[i]]$trt_events)
    dataset[[p]][[i]][,"cum_ctrl_events"]<-cumsum(dataset[[p]][[i]]$ctrl_events)
    dataset[[p]][[i]][,"cum_events"]<-dataset[[p]][[i]]$cum_trt_events+dataset[[p]][[i]]$cum_ctrl_events
    dataset[[p]][[i]]$id<-ID
  }}

#Compute total number of events in each dataset
max_events<-vector("list",4) 

for(p in 1:4){
  for(i in 1:nsim){
    max_events[[p]][[i]]<-max(dataset[[p]][[i]]$cum_events)
  }} 

#Create datasets for halfway point (strategy E)
data.halfE<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    data.halfE[[p]][[i]]<-transform(dataset[[p]][[i]],event=ifelse(cum_events>halfwayE,0,event))
  }}

for(p in 1:4){
  for(i in 1:nsim){
    data.halfE[[p]][[i]]<-transform(data.halfE[[p]][[i]],total_cum_pts=ifelse(cum_events>halfwayE,0,total_cum_pts))
  }}

obs_times.halfE<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    ifelse(max_events[[p]][[i]]<halfwayE, 
           obs_times.halfE[[p]][[i]]<-max(data.halfE[[p]][[i]][which(data.halfE[[p]][[i]]$cum_events==max_events[[p]][[i]]),"Obs_times"]), 
           obs_times.halfE[[p]][[i]]<-min(data.halfE[[p]][[i]][which(data.halfE[[p]][[i]]$cum_events==halfwayE),"Obs_times"]))
  }}


for(p in 1:4){
  for( i in 1:nsim){
    data.halfE[[p]][[i]]$Surv_times_half<-data.halfE[[p]][[i]]$Surv_times;
    data.halfE[[p]][[i]][data.halfE[[p]][[i]]$cum_events>halfwayE,]$Surv_times_half<-
      obs_times.halfE[[p]][[i]]-data.halfE[[p]][[i]][data.halfE[[p]][[i]]$cum_events>halfwayE,]$entry;
  }}

for(p in 1:4){
  for(i in 1:nsim){
    data.halfE[[p]][[i]]<-data.halfE[[p]][[i]][data.halfE[[p]][[i]]$Surv_times_half>=0,]
  }
}

#Time of IA
obs_times.halfE<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
    ifelse(max_events[[p]][[i]]<halfwayE, 
           obs_times.halfE[[p]][[i]]<-max(data.halfE[[p]][[i]][which(data.halfE[[p]][[i]]$cum_events==max_events[[p]][[i]]),"Obs_times"]), 
           obs_times.halfE[[p]][[i]]<-min(data.halfE[[p]][[i]][which(data.halfE[[p]][[i]]$cum_events==halfwayE),"Obs_times"]))
  }}


data.halfSS<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    data.halfSS[[p]][[i]]<-transform(dataset[[p]][[i]],event=ifelse(total_cum_pts>halfwaySS,0,event))
  }}

obs_times.halfSS<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    obs_times.halfSS[[p]][[i]]<-min(data.halfSS[[p]][[i]][which(data.halfSS[[p]][[i]]$total_cum_pts==halfwaySS),"Obs_times"])
  }}


for(p in 1:4){
  for( i in 1:nsim){
    data.halfSS[[p]][[i]]$Surv_times_half<-data.halfSS[[p]][[i]]$Surv_times;
    data.halfSS[[p]][[i]][data.halfSS[[p]][[i]]$total_cum_pts>halfwaySS,]$Surv_times_half<-
      obs_times.halfSS[[p]][[i]]-data.halfSS[[p]][[i]][data.halfSS[[p]][[i]]$total_cum_pts>halfwaySS,]$entry;
  }
}

for(p in 1:4){
  for(i in 1:nsim){
    data.halfSS[[p]][[i]]<-data.halfSS[[p]][[i]][data.halfSS[[p]][[i]]$Surv_times_half>=0,]
  }
}

for(p in 1:4){
  for(i in 1:nsim){
    data.halfSS[[p]][[i]]<-transform(data.halfSS[[p]][[i]],cum_events=ifelse(total_cum_pts>halfwaySS,0,cum_events))
  }}



#Compute total number of events at halfway point under strategy SS
max_halfSS<-vector("list",4) 
for(p in 1:4){
  for(i in 1:nsim){
    max_halfSS[[p]][[i]]<-max(data.halfSS[[p]][[i]]$cum_events)
  }} 

#Samples sizes at halfway point (Strategy E)
SS.half.E<-vector("list",4) 
for(p in 1:4){
  for(i in 1:nsim){
    SS.half.E[[p]][[i]]<-max(data.halfE[[p]][[i]]$total_cum_pts)
  }
}

#stan_surv model for halfway dataset (Strategy E)
mod.half.E<-vector("list",4) 

for(p in 1:4){
  for(i in 1:nsim){
    mod.half.E[[p]][[i]]<-stan_surv(Surv(Surv_times_half, event) ~ trt, data = data.halfE[[p]][[i]],
                                    basehaz = "weibull", chains = 2, seed = 20220103,iter=1000)
  }}

#stan_surv model for halfway dataset (Strategy SS)
mod.half.SS<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    mod.half.SS[[p]][[i]]<-stan_surv(Surv(Surv_times_half, event) ~ trt, data = data.halfSS[[p]][[i]], 
                                     basehaz = "weibull", chains = 2, seed = 20220103,iter=1000)
  }}


#Posterior probability halfway (strategy E)
param_est.E<-vector("list",4)
post_prob.E<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    param_est.E[[p]][[i]]<-as.data.frame(mod.half.E[[p]][[i]]$stanfit)
    post_prob.E[[p]][i]<-mean(param_est.E[[p]][[i]]$trt<0)
  }
}

#Probability of stopping halfway (strategy E)
prob_stop.E<-c()
for(p in 1:4){
  prob_stop.E[p]<-mean(unlist(post_prob.E[[p]])>sup)
}

#Posterior probability halfway (strategy SS)
param_est.SS<-vector("list",4)
post_prob.SS<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    param_est.SS[[p]][[i]]<-as.data.frame(mod.half.SS[[p]][[i]]$stanfit)
    post_prob.SS[[p]][i]<-mean(param_est.SS[[p]][[i]]$trt<0)
  }
}

#Probability of stopping halfway (Strategy SS)
prob_stop.SS<-c()
for(p in 1:4){
  prob_stop.SS[p]<-mean(unlist(post_prob.SS[[p]])>sup)
}

#Dataset for final analysis (Strategy E)
data.final.E<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
    data.final.E[[p]][[i]]<-transform(dataset[[p]][[i]],event=ifelse(cum_events>D,0,event))
  }}
for(p in 1:4){
  for(i in 1:nsim){
    data.final.E[[p]][[i]]<-transform(data.final.E[[p]][[i]],total_cum_pts=ifelse(cum_events>D,0,total_cum_pts))
  }}

obs_times.final.E<-vector("list",4) # time of final analysis (Strategy E)

for(p in 1:4){
  for(i in 1:nsim){
    ifelse(max_events[[p]][[i]]<D, obs_times.final.E[[p]][[i]]<-max(data.final.E[[p]][[i]][which(data.final.E[[p]][[i]]$cum_events==max_events[[p]][[i]]),"Obs_times"]), obs_times.final.E[[p]][[i]]<-min(data.final.E[[p]][[i]][which(data.final.E[[p]][[i]]$cum_events==D),"Obs_times"]))
  }}

for(p in 1:4){
  for( i in 1:nsim){
    data.final.E[[p]][[i]]$Surv_times_final<-data.final.E[[p]][[i]]$Surv_times;
    if(max_events[[p]][[i]]>D) {data.final.E[[p]][[i]][data.final.E[[p]][[i]]$cum_events>D,]$Surv_times_final<-obs_times.final.E[[p]][[i]]-data.final.E[[p]][[i]][data.final.E[[p]][[i]]$cum_events>D,]$entry}
  }
}

for(p in 1:4){
  for(i in 1:nsim){
    data.final.E[[p]][[i]]<-data.final.E[[p]][[i]][data.final.E[[p]][[i]]$Surv_times_final>=0,]
  }
}


#Dataset for final analysis (Strategy SS)
data.final.SS<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    data.final.SS[[p]][[i]]<-transform(dataset[[p]][[i]],event=ifelse(total_cum_pts>n,0,event))
  }}

obs_times.final.SS<-vector("list",4) # time of final analysis (Strategy SS)

for(p in 1:4){
  for(i in 1:nsim){
    obs_times.final.SS[[p]][[i]]<-min(data.final.SS[[p]][[i]][which(data.final.SS[[p]][[i]]$total_cum_pts==n),"Obs_times"])
  }}

for(p in 1:4){
  for( i in 1:nsim){
    data.final.SS[[p]][[i]]$Surv_times_final<-data.final.SS[[p]][[i]]$Surv_times;
    data.final.SS[[p]][[i]][data.final.SS[[p]][[i]]$total_cum_pts>n,]$Surv_times_final<-
      obs_times.final.SS[[p]][[i]]-data.final.SS[[p]][[i]][data.final.SS[[p]][[i]]$total_cum_pts>n,]$entry
  }
}

for(p in 1:4){
  for(i in 1:nsim){
    data.final.SS[[p]][[i]]<-data.final.SS[[p]][[i]][data.final.SS[[p]][[i]]$Surv_times_final>=0,]
  }
}

for(p in 1:4){
  for(i in 1:nsim){
    data.final.SS[[p]][[i]]<-transform(data.final.SS[[p]][[i]],cum_events=ifelse(total_cum_pts>n,0,cum_events))
  }}


#Total number of participants at final analysis (Strategy SS)
max_SS.final<-vector("list",4) 
for(p in 1:4){
  for(i in 1:nsim){
    max_SS.final[[p]][[i]]<-max(data.final.SS[[p]][[i]]$cum_events)
  }} 

#Total number of participants at final analysis (Strategy SS)
SS.final.E<-vector("list",4) 
for(p in 1:4){
  for(i in 1:nsim){
    SS.final.E[[p]][[i]]<-max(data.final.E[[p]][[i]]$total_cum_pts)
  }
}

#Extract info at IA to update priors (Strategy E)
summary.mod.half.E<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
    summary.mod.half.E[[p]][[i]]<-as.data.frame(mod.half.E[[p]][[i]]$stan_summary)
  }
}

#Stan_surv model for final analysis (Strategy E)
mod.half.E.final<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
    mod.half.E.final[[p]][[i]]<-stan_surv(Surv(Surv_times_final, event) ~ trt, data = data.final.E[[p]][[i]], 
                                          basehaz = "weibull", chains = 2, seed = 20220103,
                                          prior=normal(summary.mod.half.E[[p]][[i]]$mean[2],summary.mod.half.E[[p]][[i]]$sd[2]),iter=1000)
  }}


#Extract info at IA to update priors (Strategy SS)
summary.mod.half.SS<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
    summary.mod.half.SS[[p]][[i]]<-as.data.frame(mod.half.SS[[p]][[i]]$stan_summary)
  }
}


#Stan_surv model for final analysis (Strategy SS)
mod.half.SS.final<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
    mod.half.SS.final[[p]][[i]]<-stan_surv(Surv(Surv_times_final, event) ~ trt, data = data.final.SS[[p]][[i]], 
                                           basehaz = "weibull", chains = 2, seed = 20220103,
                                           prior=normal(summary.mod.half.SS[[p]][[i]]$mean[2],summary.mod.half.SS[[p]][[i]]$sd[2]),iter=1000)
  }}


#Posterior probability of success (Strategy E)
param_est.final.E<-vector("list",4)
post_prob.final.E<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    if(post_prob.E[[p]][i]<sup){
      param_est.final.E[[p]][[i]]<-as.data.frame(mod.half.E.final[[p]][[i]]$stanfit);
      post_prob.final.E[[p]][i]<-mean(param_est.final.E[[p]][[i]]$trt<0)}
    else{post_prob.final.E[[p]][i]<-NA}
  }
}

#Probability of stopping at final analysis (Strategy E)
prob_stop.final.E<-vector("list",4)
for(p in 1:4){
  prob_stop.final.E[[p]]<-mean(na.omit(unlist(post_prob.final.E[[p]])>sup))
}

#Probability of stopping overall (Strategy E)
prob_stop_overall.E<-c()
for(p in 1:p){
  prob_stop_overall.E[p]<-prob_stop.E[[p]]+(1-prob_stop.E[[p]])*prob_stop.final.E[[p]]
}

#Posterior probability of success (Strategy SS)
param_est.final.SS<-vector("list",4)
post_prob.final.SS<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    if(post_prob.SS[[p]][i]<sup){
      param_est.final.SS[[p]][[i]]<-as.data.frame(mod.half.SS.final[[p]][[i]]$stanfit);
      post_prob.final.SS[[p]][i]<-mean(param_est.final.SS[[p]][[i]]$trt<0)}
    else{post_prob.final.SS[[p]][i]<-NA}
  }
}

#Probability of stopping at final analysis (Strategy SS)
prob_stop.final.SS<-vector("list",4)
for(p in 1:4){
  prob_stop.final.SS[[p]]<-mean(na.omit(unlist(post_prob.final.SS[[p]])>sup))
}

#Probability of stopping overall (Strategy SS)
prob_stop_overall.SS<-c()
for(p in 1:p){
  prob_stop_overall.SS[p]<-prob_stop.SS[[p]]+(1-prob_stop.SS[[p]])*prob_stop.final.SS[[p]]
}

############ RESULTS ################

# Define helper function to calculate mean of lists
mean_of_list <- function(lst) {
  mean(unlist(lst))
}

# Calculate values for each row
#Row 1: Probability of stopping halfway (P_1)
row1 <- data.frame(E = round(prob_stop.E, 2), SS = round(prob_stop.SS, 2))
#Row 2: Expected time of IA (T_1)
row2 <- data.frame(E = round(sapply(obs_times.halfE, mean_of_list),1), SS = round(sapply(obs_times.halfSS, mean_of_list),1))
#Row 3: Number of expected events halfway (D_1)
row3 <- data.frame(E = rep(halfwayE, 4), SS = ceiling(sapply(max_halfSS, mean_of_list)))
#Row 4: Expected sample size halfway (N_1)
row4 <- data.frame(E = ceiling(sapply(SS.half.E, mean_of_list)), SS = rep(halfwaySS, 4))
#Row 5: Probability of stopping overall (P_2)
row5 <- data.frame(E = round(prob_stop_overall.E,2), SS = round(prob_stop_overall.SS,2))
#Row 6: Expected time of final analysis (T_2)
row6 <- data.frame(E = round(mapply(function(p, t1, t2) p*t1 + (1-p)*t2, prob_stop.E, sapply(obs_times.halfE, mean_of_list), sapply(obs_times.final.E, mean_of_list)),2),
                   SS = round(mapply(function(p, t1, t2) p*t1 + (1-p)*t2, prob_stop.SS, sapply(obs_times.halfSS, mean_of_list), sapply(obs_times.final.SS, mean_of_list)),2))
#Row 7: Number of expected events overall (D_2)
row7 <- data.frame(E = ceiling(mapply(function(p) p*halfwayE + (1-p)*D, prob_stop.E)),
                   SS = ceiling(mapply(function(p, t1, t2) p*t1 + (1-p)*t2, prob_stop.SS, sapply(max_halfSS, mean_of_list), sapply(max_SS.final, mean_of_list))))
#Row 8: Expected sample size overall (N_2)
row8 <- data.frame(E = ceiling(mapply(function(p, t1, t2) p*t1 + (1-p)*t2, prob_stop.E, sapply(SS.half.E, mean_of_list), sapply(SS.final.E, mean_of_list))),
                   SS = ceiling(mapply(function(p) p*halfwaySS + (1-p)*n, prob_stop.SS)))


r1<-c(unlist(row1$E),unlist(row1$SS))[c(1,5,2,6,3,7,4,8)]
r2<-c(unlist(row2$E),unlist(row2$SS))[c(1,5,2,6,3,7,4,8)]
r3<-c(unlist(row3$E),unlist(row3$SS))[c(1,5,2,6,3,7,4,8)]
r4<-c(unlist(row4$E),unlist(row4$SS))[c(1,5,2,6,3,7,4,8)]
r5<-c(unlist(row5$E),unlist(row5$SS))[c(1,5,2,6,3,7,4,8)]
r6<-c(unlist(row6$E),unlist(row6$SS))[c(1,5,2,6,3,7,4,8)]
r7<-c(unlist(row7$E),unlist(row7$SS))[c(1,5,2,6,3,7,4,8)]
r8<-c(unlist(row8$E),unlist(row8$SS))[c(1,5,2,6,3,7,4,8)]


results<-data.frame(rbind(r1,r2,r3,r4,r5,r6,r7,r8))
rownames(results)<-c("P_1","T_1","D_1","N_1","P_2","T_2","D_2","N_2")
results





