######################################
#                                    #                                    
# Script for Pfizer trial simulation #  
#                                    #                                    
######################################

library(dplyr)
rm(list = ls())

n<-37000 # Total number of patients in trial

# Solve for true attack rate
lambda_0<-round(solve(n/2*(104/365)+n/2*(104/365)*(1-0.95),94),3) 
lambda_0

# Compute total number of expected events given assumed attack rate of 1.3% and VE of 60%
total_events<- round(n/2*0.013+n/2*0.013*(1-0.6))
total_events

# Proportion of EXPECTED events accrued for first IA (strategy E)
prop<-94/total_events
prop

# Number of patients accrued for first IA (had strategy SS been used)
n_IA<-round(prop*n)
n_IA

# Timing (in days) of IA had strategy SS been used 
days<-round(solve((n/2)*0.013/365+(n/2)*0.013/365*(1-0.60),94))
days


# Patient-years that would have elapsed at IA for strategy SS
round((days/365)*prop*n)


##### STRATEGY E#####
set.seed(287182229)

Nsim<-1000 # Number of simulated trials

pc<-lambda_0 # control arm incidence rate i.e. lambda_0
VE<-0.6 # vaccine efficacy (CHANGE TO 0.60, 0.70, 0.80 or 0.90 DEPENDING ON SCENARIO (or 0 for FPR))
pt<-(1-VE)*pc # vaccine arm incidence rate
n<-37000 # total number of patients in trial
D<-94 # number of events at IA for strategy E

T<-list(); ID<-list()

for(i in 1:Nsim){
  T[[i]]<-rep(c(0,1),n/2)
  ID[[i]]<-c(1:n)
}

EVENT_C<-replicate(Nsim,rbinom(n/2,1,pc), simplify = FALSE)

EVENT_T<-replicate(Nsim,rbinom(n/2,1,pt), simplify = FALSE)

EVENT<-vector("list",Nsim)

for(i in 1:Nsim){
  EVENT[[i]]<-c(rbind(EVENT_C[[i]],EVENT_T[[i]]))
}

dataset<-vector("list",Nsim)

for(i in 1:Nsim){
  dataset[[i]]<-data.frame(id=ID[[i]],trt=T[[i]],event=EVENT[[i]])
  dataset[[i]][,"cum_trt_pts"]<-cumsum(dataset[[i]]$trt) #cumulative number of patients in trt arm
  dataset[[i]][,"cum_ctrl_pts"]<-cumsum(1-dataset[[i]]$trt) #cumulative number of patients in control arm
  dataset[[i]][,"total_cum_pts"]<-dataset[[i]]$cum_trt_pts+dataset[[i]]$cum_ctrl_pts #total cumulative number of patients
  dataset[[i]][,"trt_events"]<-ifelse(dataset[[i]]$trt==1 & dataset[[i]]$event==1,1,0) #event is in treatment arm
  dataset[[i]][,"ctrl_events"]<-ifelse(dataset[[i]]$trt==0 & dataset[[i]]$event==1,1,0) #event is in control arm
  dataset[[i]][,"cum_trt_events"]<-cumsum(dataset[[i]]$trt_events) #cumulative number of events in treatment arm
  dataset[[i]][,"cum_ctrl_events"]<-cumsum(dataset[[i]]$ctrl_events) #cumulative number of events in control arm
  dataset[[i]][,"cum_events"]<-cumsum(dataset[[i]]$event) #total cumulative number of events
}


nc<-vector("list",Nsim); dc<-vector("list",Nsim)
post_a_c<-vector("list",Nsim); post_b_c<-vector("list",Nsim)


for(i in 1:Nsim){
  nc[[i]]<-min(dataset[[i]][which(dataset[[i]]$cum_events==D),"cum_ctrl_pts"])  #find number of control patients when total events=94
  dc[[i]]<-min(dataset[[i]][which(dataset[[i]]$cum_events==D),"cum_ctrl_events"]) #find number of events from control patients when total events=94
  
  #Posterior parameters: 
  post_a_c[[i]]<-1+dc[[i]]
  post_b_c[[i]]<-1+nc[[i]]-dc[[i]]
}

nt<-vector("list",Nsim); dt<-vector("list",Nsim)
post_a_t<-vector("list",Nsim); post_b_t<-vector("list",Nsim)

for(i in 1:Nsim){
  nt[[i]]<-min(dataset[[i]][which(dataset[[i]]$cum_events==D),"cum_trt_pts"])
  dt[[i]]<-min(dataset[[i]][which(dataset[[i]]$cum_events==D), "cum_trt_events"])
  #Posterior parameters:
  post_a_t[[i]]<-1+dt[[i]]
  post_b_t[[i]]<-1+nt[[i]]-dt[[i]]
}

# relative risk
RRsamps<-vector("list",Nsim)

# simulate 2000 RRs at every IA to obtain posterior probability that RR<1
for(i in 1:Nsim){
  RRsamps[[i]]<-rbeta(2000,post_a_t[[i]],post_b_t[[i]])/
    rbeta(2000,post_a_c[[i]],post_b_c[[i]])
}


# posterior prob that RR<1
postprob_E<-vector("list",Nsim)
for(i in 1:Nsim){
  postprob_E[[i]]<-mean(RRsamps[[i]]<0.70) #post prob of success for each simulation
}

# probability of stopping at IA
prob_stop_E<-mean(postprob_E>0.995)


# number of days that would have elapsed
solve(n/2*lambda_0/365+n/2*lambda_0/365*(1-VE),94)


##### STRATEGY SS ##### 
nc<-vector("list",Nsim); dc<-vector("list",Nsim)
post_a_c<-vector("list",Nsim); post_b_c<-vector("list",Nsim)

for(i in 1:Nsim){
  nc[[i]]<-n_IA/2  
  dc[[i]]<-min(dataset[[i]][which(dataset[[i]]$total_cum_pts==n_IA),"cum_ctrl_events"]) #find number of events from control patients when total patients=halfway
  #Posterior parameters: 
  post_a_c[[i]]<-1+dc[[i]]
  post_b_c[[i]]<-1+nc[[i]]-dc[[i]]
}

nt<-vector("list",Nsim); dt<-vector("list",Nsim)
post_a_t<-vector("list",Nsim); post_b_t<-vector("list",Nsim)

for(i in 1:Nsim){
  nt[[i]]<-n_IA/2
  dt[[i]]<-min(dataset[[i]][which(dataset[[i]]$total_cum_pts==n_IA), "cum_trt_events"])
  #Posterior parameters:
  post_a_t[[i]]<-1+dt[[i]]
  post_b_t[[i]]<-1+nt[[i]]-dt[[i]]
}


# relative risk
RRsamps<-vector("list",Nsim)

# simulate 2000 RRs at every IA to obtain posterior probability that RR<1
for(i in 1:Nsim){
  RRsamps[[i]]<-rbeta(2000,post_a_t[[i]],post_b_t[[i]])/
    rbeta(2000,post_a_c[[i]],post_b_c[[i]])
}


# posterior prob that RR<1
postprob_SS<-vector("list",Nsim)

for(i in 1:Nsim){
  postprob_SS[[i]]<-mean(RRsamps[[i]]<0.70) #post prob of success for each simulation
}

# probability of stopping at IA
prob_stop_SS<-mean(postprob_SS>0.995)
prob_stop_SS



