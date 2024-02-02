###################################################
#                                                 #
#       Script for TTE outcome simulations        #
#   High censoring, control hazard well estimated #
###################################################

rm(list=ls())

#Load required packages 
library(rstanarm)
library(survival)
library(dplyr)


set.seed(48290103)

#Set parameter values for log (HRs)
beta0<--1 
beta1 <- log(c(0.9, 0.7, 0.5, 1))
shape=1.5 #shape parameter of Weibull distribution

#Reparameterization
sigma0=1/(exp(beta0/shape))
sigma1 <- sapply(beta0 + beta1, function(x) 1 / exp(x / shape))

ta<-3 #accrual period
tf<-1.7 #follow-up period
cutoff<-ta+tf  #administrative censoring time


phi<-1.3 #1.3 for high dropout rate; Change to 0.15 for low dropout rate
sup<-0.98 #efficacy stopping threshold

n=1018 #sample size for assumed beta_0 of -1 (change sample size depending on phi and assumed beta_0)
D=196 #number of expected events 
halfwayE<-ceiling(D/2) #halfway point for strategy E
halfwaySS<-ceiling(n/2) #halfway point for strategy SS
nsim<-1000 #number of replications 

entry_times <- replicate(nsim, runif(n, 0, ta), simplify = FALSE)
entry_times <- lapply(entry_times, function(x) sort(x))

entry_c <- lapply(entry_times, function(x) x[c(TRUE, FALSE)]) #control arm 
entry_t <- lapply(entry_times, function(x) x[c(FALSE, TRUE)]) #treatment arm

surv_c <- replicate(nsim, rweibull(n / 2, shape = shape, scale = sigma0), simplify = FALSE) #control arm event times
cens_c <- replicate(nsim, rexp(n / 2, phi), simplify = FALSE) #control arm censoring times 
ysurv_c <- lapply(1:nsim, function(i) entry_c[[i]] + surv_c[[i]]) #observed event times
ycens_c <- lapply(1:nsim, function(i) entry_c[[i]] + cens_c[[i]]) #observed censoring times

surv_t <- lapply(sigma1, function(s) replicate(nsim, rweibull(n / 2, shape = shape, scale = s), simplify = FALSE)) #treatment arm survival times
cens_t <- lapply(rep(n/2, 4), function(x) replicate(nsim, rexp(x, phi), simplify = FALSE)) #treatment arm censoring times
ysurv_t <- lapply(1:4, function(p) lapply(1:nsim, function(i) entry_t[[i]] + surv_t[[p]][[i]])) #observed event times
ycens_t <- lapply(1:4, function(p) lapply(1:nsim, function(i) entry_t[[i]] + cens_t[[p]][[i]])) #observed censoring times

censoring_time<-rep(cutoff,n/2)

Obstimes_C <- lapply(1:nsim, function(i) pmin(ysurv_c[[i]], ycens_c[[i]], censoring_time)) #observed times in control arm
Obstimes_T <- lapply(1:4, function(p) lapply(1:nsim, function(i) pmin(ysurv_t[[p]][[i]], ycens_t[[p]][[i]], censoring_time))) #observed times in treatment arm
Obstimes <- lapply(1:4, function(p) lapply(1:nsim, function(i) c(rbind(Obstimes_C[[i]], Obstimes_T[[p]][[i]])))) 

Cens_c <- lapply(1:nsim, function(i) ifelse(ycens_c[[i]] < ysurv_c[[i]] | censoring_time < ysurv_c[[i]], 0, 1)) #censoring indicator variable 
Cens_t <- lapply(1:4, function(p) lapply(1:nsim, function(i) ifelse(ycens_t[[p]][[i]] < ysurv_t[[p]][[i]] | censoring_time < ysurv_t[[p]][[i]], 0, 1)))
Cens <- lapply(1:4, function(p) lapply(1:nsim, function(i) c(rbind(Cens_c[[i]], Cens_t[[p]][[i]]))))

arm<-c(rep(c(0,1),n/2)) #arm indicator
ID<-1:n

#Create datasets for analysis
dataset<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
    dataset[[p]][[i]]<-data.frame(trt=arm,entry=entry_times[[i]], Obs_times=Obstimes[[p]][[i]],event=Cens[[p]][[i]])
    dataset[[p]][[i]]<-arrange(dataset[[p]][[i]],Obs_times)
    dataset[[p]][[i]][,"cum_trt_pts"]<-cumsum(dataset[[p]][[i]]$trt) #cumulative number of patients in trt arm
    dataset[[p]][[i]][,"cum_ctrl_pts"]<-cumsum(1-dataset[[p]][[i]]$trt) #cumulative number of patients in control arm
    dataset[[p]][[i]][,"total_cum_pts"]<-dataset[[p]][[i]]$cum_trt_pts+dataset[[p]][[i]]$cum_ctrl_pts #total cumulative number of patients
    dataset[[p]][[i]][,"trt_events"]<-ifelse(dataset[[p]][[i]]$trt==1 & dataset[[p]][[i]]$event==1,1,0) #event is in treatment arm
    dataset[[p]][[i]][,"ctrl_events"]<-ifelse(dataset[[p]][[i]]$trt==0 & dataset[[p]][[i]]$event==1,1,0) #event is in control arm
    dataset[[p]][[i]][,"cum_trt_events"]<-cumsum(dataset[[p]][[i]]$trt_events) #cumulative number of events in treatment arm
    dataset[[p]][[i]][,"cum_ctrl_events"]<-cumsum(dataset[[p]][[i]]$ctrl_events) #cumulative number of events in control arm
    dataset[[p]][[i]][,"cum_events"]<-dataset[[p]][[i]]$cum_trt_events+dataset[[p]][[i]]$cum_ctrl_events #total cumulative number of events
    dataset[[p]][[i]]$id<-ID
  }}

#Compute total number of events in each dataset
max_events <- lapply(dataset, function(dat) sapply(dat, function(d) max(d$cum_events)))
#compute total number of events in each dataset

#Create datasets for halfway point IA (strategy E)
halfway_row<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
  halfway_row[[p]][[i]]<-which(dataset[[p]][[i]]$cum_events==halfwayE)[1] # row where interim analysis occurs
  halfway_row[[p]][[i]]<-ifelse(is.na(halfway_row[[p]][[i]]), n, halfway_row[[p]][[i]])
  }}

data.halfE<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    data.halfE[[p]][[i]]<-transform(dataset[[p]][[i]],event=ifelse(id>halfway_row[[p]][[i]],0,event),
                                    total_cum_pts = ifelse(id> halfway_row[[p]][[i]], 0, total_cum_pts))
  }}

obs_times.halfE<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    ifelse(max_events[[p]][[i]]<halfwayE, 
           obs_times.halfE[[p]][[i]]<-max(data.halfE[[p]][[i]][which(data.halfE[[p]][[i]]$cum_events==max_events[[p]][[i]]),"Obs_times"]), 
           obs_times.halfE[[p]][[i]]<-min(data.halfE[[p]][[i]][which(data.halfE[[p]][[i]]$cum_events==halfwayE),"Obs_times"]))
    # if total number of events is smaller than halfway point, there is no interim analysis, go to final analysis directly
  }}


for(p in 1:4){
  for( i in 1:nsim){
    data.halfE[[p]][[i]]$Surv_times_half<-data.halfE[[p]][[i]]$Obs_times-data.halfE[[p]][[i]]$entry;
    data.halfE[[p]][[i]][data.halfE[[p]][[i]]$id>halfway_row[[p]][[i]],]$Surv_times_half<-
      obs_times.halfE[[p]][[i]]-data.halfE[[p]][[i]][data.halfE[[p]][[i]]$id>halfway_row[[p]][[i]],]$entry;
    data.halfE[[p]][[i]]<-data.halfE[[p]][[i]][data.halfE[[p]][[i]]$Surv_times_half>=0, ]
  }}

#Create datasets for halfway point IA (strategy SS)

data.halfSS<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
    data.halfSS[[p]][[i]]<-transform(dataset[[p]][[i]], event=ifelse(total_cum_pts>halfwaySS,0,event),
                                                        cum_events = ifelse(total_cum_pts > halfwaySS, 0, cum_events))
  }}

obs_times.halfSS<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    obs_times.halfSS[[p]][[i]]<-min(data.halfSS[[p]][[i]][which(data.halfSS[[p]][[i]]$total_cum_pts==halfwaySS),"Obs_times"])
  }}


for(p in 1:4){
  for( i in 1:nsim){
    data.halfSS[[p]][[i]]$Surv_times_half<-data.halfSS[[p]][[i]]$Obs_times-data.halfSS[[p]][[i]]$entry;
    data.halfSS[[p]][[i]][data.halfSS[[p]][[i]]$total_cum_pts>halfwaySS,]$Surv_times_half<-
      obs_times.halfSS[[p]][[i]]-data.halfSS[[p]][[i]][data.halfSS[[p]][[i]]$total_cum_pts>halfwaySS,]$entry;
    data.halfSS[[p]][[i]]<-data.halfSS[[p]][[i]][data.halfSS[[p]][[i]]$Surv_times_half>=0,]
  }
}


#Number of events at halfway point (Strategy E)
events_halfSS<-vector("list",4) 
for(p in 1:4){
  for(i in 1:nsim){
    events_halfSS[[p]][[i]]<-max(data.halfSS[[p]][[i]]$cum_events)
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
                                    basehaz = "weibull", chains = 2, seed = 482910,iter=2000)
  }}

save(mod.half.E, file="E.half.high.well.RData")


#stan_surv model for halfway dataset (Strategy SS)
mod.half.SS<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    mod.half.SS[[p]][[i]]<-stan_surv(Surv(Surv_times_half, event) ~ trt, data = data.halfSS[[p]][[i]], 
                                     basehaz = "weibull", chains = 2, seed = 482910,iter=2000)
  }}

save(mod.half.SS, file="SS.half.high.well.RData")
# Function to calculate posterior probabilities and stopping probabilities

calc_prob_stop <- function(mod, halfway, sup) {
  obj<-list()
  param_est <- lapply(mod, function(x) lapply(x, function(y) as.data.frame(y$stanfit))) #parameter estimates
  obj$post_prob <- lapply(param_est, function(x) sapply(x, function(y) mean(y$trt < 0))) #posterior probabilities of alternative hypothesis
  obj$prob_stop <- sapply(obj$post_prob, function(x) mean(unlist(x) > sup)) #probability of stopping
  return(obj)
}

# Probability of stopping halfway (strategy E)
prob_E <- calc_prob_stop(mod.half.E, halfwayE, sup)

# Probability of stopping halfway (strategy SS)
prob_SS <- calc_prob_stop(mod.half.SS, halfwaySS, sup)

#Dataset for final analysis (Strategy E AND Strategy SS)
data.final<-dataset

obs_times.final<-vector("list",4) # time of final analysis (Strategy E)

for(p in 1:4){
  for(i in 1:nsim){
    obs_times.final[[p]][[i]]<-max(data.final[[p]][[i]][,"Obs_times"])
  }}

for(p in 1:4){
  for( i in 1:nsim){
    data.final[[p]][[i]]$Surv_times_final<-data.final[[p]][[i]]$Obs_times-data.final[[p]][[i]]$entry
  }
}

#Total number of events at final analysis (Strategy SS)
events_finalSS<-vector("list",4) 
for(p in 1:4){
  for(i in 1:nsim){
    events_finalSS[[p]][[i]]<-max(data.final[[p]][[i]]$cum_events)
  }} 

#Total number of participants at final analysis (Strategy E)
SS.final.E<-vector("list",4) 
for(p in 1:4){
  for(i in 1:nsim){
    SS.final.E[[p]][[i]]<-max(data.final[[p]][[i]]$total_cum_pts)
  }
}

#Extract info at IA to update priors (Strategy E)
summary.mod.half.E<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
    summary.mod.half.E[[p]][[i]]<-as.data.frame(mod.half.E[[p]][[i]]$stan_summary)
  }
}

mod.final.E<-vector("list",4)

#Stan_surv model for final analysis (Strategy E)
for(p in 1:4){
  for(i in 1:nsim){
    mod.final.E[[p]][[i]]<-stan_surv(Surv(Surv_times_final, event) ~ trt, data = data.final[[p]][[i]], 
                                          basehaz = "weibull", chains = 2, seed = 2392810,
                                          prior=normal(summary.mod.half.E[[p]][[i]]$mean[2],summary.mod.half.E[[p]][[i]]$sd[2]),iter=2000)
  }}

save(mod.final.E, file="E.final.high.well.RData")

#Extract info at IA to update priors (Strategy SS)
summary.mod.half.SS<-vector("list",4)
for(p in 1:4){
  for(i in 1:nsim){
    summary.mod.half.SS[[p]][[i]]<-as.data.frame(mod.half.SS[[p]][[i]]$stan_summary)
  }
}

mod.final.SS<-vector("list",4)

#Stan_surv model for final analysis (Strategy SS)
for(p in 1:4){
  for(i in 1:nsim){
    mod.final.SS[[p]][[i]]<-stan_surv(Surv(Surv_times_final, event) ~ trt, data = data.final[[p]][[i]], 
                                           basehaz = "weibull", chains = 2, seed = 2392810,
                                           prior=normal(summary.mod.half.SS[[p]][[i]]$mean[2],summary.mod.half.SS[[p]][[i]]$sd[2]),iter=2000)
  }}

save(mod.final.SS, file="SS.final.high.well.RData")


#Posterior probability of success (Strategy E)
param_est.final.E<-vector("list",4)
post_prob.final.E<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    if(prob_E$post_prob[[p]][i]<sup){
      param_est.final.E[[p]][[i]]<-as.data.frame(mod.final.E[[p]][[i]]$stanfit);
      post_prob.final.E[[p]][i]<-mean(param_est.final.E[[p]][[i]]$trt<0)}
    else{post_prob.final.E[[p]][i]<-NA}
  }
}

#Probability of stopping at final analysis (Strategy E)
prob_stop.final.E<-c()
for(p in 1:4){
  prob_stop.final.E[p]<-mean(na.omit(unlist(post_prob.final.E[[p]])>sup))
}

#Probability of stopping overall (Strategy E)
prob_stop_overall.E<-c()
for(p in 1:p){
  prob_stop_overall.E[p]<-prob_E$prob_stop[p]+(1-prob_E$prob_stop[p])*prob_stop.final.E[p]
}

prob_stop_overall.E

#Posterior probability of success (Strategy SS)
param_est.final.SS<-vector("list",4)
post_prob.final.SS<-vector("list",4)

for(p in 1:4){
  for(i in 1:nsim){
    if(prob_SS$post_prob[[p]][i]<sup){
      param_est.final.SS[[p]][[i]]<-as.data.frame(mod.final.SS[[p]][[i]]$stanfit);
      post_prob.final.SS[[p]][i]<-mean(param_est.final.SS[[p]][[i]]$trt<0)}
    else{post_prob.final.SS[[p]][i]<-NA}
  }
}

#Probability of stopping at final analysis (Strategy SS)
prob_stop.final.SS<-c()
for(p in 1:4){
  prob_stop.final.SS[p]<-mean(na.omit(unlist(post_prob.final.SS[[p]])>sup))
}

#Probability of stopping overall (Strategy SS)
prob_stop_overall.SS<-c()
for(p in 1:p){
  prob_stop_overall.SS[p]<-prob_SS$prob_stop[p]+(1-prob_SS$prob_stop[p])*prob_stop.final.SS[p]
}


############ RESULTS ################

# Define helper function to calculate mean of lists
mean_of_list <- function(lst) {
  mean(unlist(lst))
}

# Calculate values for each row
#Row 1: Probability of stopping halfway (P_1)
row1 <- data.frame(E = round(prob_E$prob_stop, 2), SS = round(prob_SS$prob_stop, 2))
#Row 2: Expected time of IA (T_1)
row2 <- data.frame(E = round(sapply(obs_times.halfE, mean_of_list),1), SS = round(sapply(obs_times.halfSS, mean_of_list),1))
#Row 3: Number of expected events halfway (D_1)
row3 <- data.frame(E = rep(halfwayE, 4), SS = round(sapply(events_halfSS, mean_of_list)))
#Row 4: Expected sample size halfway (N_1)
row4 <- data.frame(E = round(sapply(SS.half.E, mean_of_list)), SS = rep(halfwaySS, 4))
#Row 5: Probability of stopping overall (P_2)
row5 <- data.frame(E = round(prob_stop_overall.E,2), SS = round(prob_stop_overall.SS,2))
#Row 6: Expected time of final analysis (T_2)
row6 <- data.frame(E = round(mapply(function(p, t1, t2) p*t1 + (1-p)*t2, prob_E$prob_stop, sapply(obs_times.halfE, mean_of_list), sapply(obs_times.final, mean_of_list)),1),
                   SS = round(mapply(function(p, t1, t2) p*t1 + (1-p)*t2, prob_SS$prob_stop, sapply(obs_times.halfSS, mean_of_list), sapply(obs_times.final, mean_of_list)),1))
#Row 7: Number of expected events overall (D_2)
row7 <- data.frame(E = round(mapply(function(p) p*halfwayE + (1-p)*D, prob_E$prob_stop)),
                   SS = round(mapply(function(p, t1, t2) p*t1 + (1-p)*t2, prob_SS$prob_stop, sapply(events_halfSS, mean_of_list), sapply(events_finalSS, mean_of_list))))
#Row 8: Expected sample size overall (N_2)
row8 <- data.frame(E = round(mapply(function(p, t1, t2) p*t1 + (1-p)*t2, prob_E$prob_stop, sapply(SS.half.E, mean_of_list), sapply(SS.final.E, mean_of_list))),
                   SS = round(mapply(function(p) p*halfwaySS + (1-p)*n, prob_SS$prob_stop)))


r1<-c(unlist(row1$E),unlist(row1$SS))[c(1,5,2,6,3,7,4,8)]
r2<-c(unlist(row2$E),unlist(row2$SS))[c(1,5,2,6,3,7,4,8)]
r3<-c(unlist(row3$E),unlist(row3$SS))[c(1,5,2,6,3,7,4,8)]
r4<-c(unlist(row4$E),unlist(row4$SS))[c(1,5,2,6,3,7,4,8)]
r5<-c(unlist(row5$E),unlist(row5$SS))[c(1,5,2,6,3,7,4,8)]
r6<-c(unlist(row6$E),unlist(row6$SS))[c(1,5,2,6,3,7,4,8)]
r7<-c(unlist(row7$E),unlist(row7$SS))[c(1,5,2,6,3,7,4,8)]
r8<-c(unlist(row8$E),unlist(row8$SS))[c(1,5,2,6,3,7,4,8)]


results.high.well<-data.frame(rbind(r1,r2,r3,r4,r5,r6,r7,r8))
rownames(results.high.well)<-c("P1","T1","D1","N1","P2","T2","D2","N2")

save(results.high.well,file="results.high.well.RData")
results.high.well
