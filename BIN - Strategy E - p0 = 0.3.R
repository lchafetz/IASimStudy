####################################################################
# 
# Script for binary outcome simulations for Strategy E (event-based)
# assumed p_0=0.30
####################################################################

rm(list=ls())
set.seed(382012)

sup=0.98 # efficacy stopping threshold
n<-574 # sample size for assumed p_0=0.30)
D<-148 # number of expected events for assumed p_0=0.30
halfwayE<-ceiling(D/2) # halfway point for strategy E
Nsim<-1000 # number of replications
pc<-0.2 # true control arm rate 
pt<-pc*c(0.9,0.7,0.5,1) # true treatment arm rates 

T<-list(); ID<-list()

for(i in 1:Nsim){
  T[[i]]<-rep(c(0,1),n/2)
  ID[[i]]<-c(1:n)
}

EVENT_C<-replicate(Nsim,rbinom(n/2,1,pc), simplify = FALSE)

EVENT_T<-vector("list",4)

for(j in 1:4){
  EVENT_T[[j]]<-replicate(Nsim,rbinom(n/2,1,pt[j]), simplify = FALSE)
}

EVENT<-vector("list",4)

for(j in 1:4){
  for(i in 1:Nsim){
    EVENT[[j]][[i]]<-c(rbind(EVENT_C[[i]],EVENT_T[[j]][[i]]))
  }}

dataset<-vector("list",4)
for(j in 1:4){
  for(i in 1:Nsim){
    dataset[[j]][[i]]<-data.frame(id=ID[[i]],trt=T[[i]],event=EVENT[[j]][[i]])
    dataset[[j]][[i]][,"cum_trt_pts"]<-cumsum(dataset[[j]][[i]]$trt)
    dataset[[j]][[i]][,"cum_ctrl_pts"]<-cumsum(1-dataset[[j]][[i]]$trt)
    dataset[[j]][[i]][,"total_cum_pts"]<-dataset[[j]][[i]]$cum_trt_pts+dataset[[j]][[i]]$cum_ctrl_pts
    dataset[[j]][[i]][,"trt_events"]<-ifelse(dataset[[j]][[i]]$trt==1 & dataset[[j]][[i]]$event==1,1,0)
    dataset[[j]][[i]][,"ctrl_events"]<-ifelse(dataset[[j]][[i]]$trt==0 & dataset[[j]][[i]]$event==1,1,0)
    dataset[[j]][[i]][,"cum_trt_events"]<-cumsum(dataset[[j]][[i]]$trt_events)
    dataset[[j]][[i]][,"cum_ctrl_events"]<-cumsum(dataset[[j]][[i]]$ctrl_events)
    dataset[[j]][[i]][,"cum_events"]<-cumsum(dataset[[j]][[i]]$event)
  }}

#compute total number of events in each dataset
max_events<-vector("list",4) 
for(j in 1:4){
  for(i in 1:Nsim){
    max_events[[j]][[i]]<-max(dataset[[j]][[i]]$cum_events)
  }} 

nc<-vector("list",4); dc<-vector("list",4)
post_a_c<-vector("list",4); post_b_c<-vector("list",4)

for(p in 1:4){
  for(i in 1:Nsim){
    nc[[p]][i]<-ifelse(max_events[[p]][[i]]>halfwayE, min(dataset[[p]][[i]][which(dataset[[p]][[i]]$cum_events==halfwayE),"cum_ctrl_pts"]),
                       max(dataset[[p]][[i]][which(dataset[[p]][[i]]$cum_events==max_events[[p]][[i]]),"cum_ctrl_pts"]))
    # if the expected number of events halfway is greater than the total number of events in the trial, 
    # then simply take the total number of events (at end of trial)
    dc[[p]][i]<-ifelse(max_events[[p]][[i]]>halfwayE, min(dataset[[p]][[i]][which(dataset[[p]][[i]]$cum_events==halfwayE),"cum_ctrl_events"]),
                       max(dataset[[p]][[i]][which(dataset[[p]][[i]]$cum_events==max_events[[p]][[i]]),"cum_ctrl_events"]))
    
    #Posterior parameters at IA (control group): 
    post_a_c[[p]][i]<-1+dc[[p]][i]
    post_b_c[[p]][i]<-1+nc[[p]][i]-dc[[p]][i]
  }}

nt<-vector("list",4); dt<-vector("list",4)
post_a_t<-vector("list",4); post_b_t<-vector("list",4)

for(p in 1:4){
  for(i in 1:Nsim){
    nt[[p]][i]<-ifelse(max_events[[p]][[i]]>halfwayE, min(dataset[[p]][[i]][which(dataset[[p]][[i]]$cum_events==halfwayE),"cum_trt_pts"]),
                       max(dataset[[p]][[i]][which(dataset[[p]][[i]]$cum_events==max_events[[p]][[i]]),"cum_trt_pts"]))
    dt[[p]][i]<-ifelse(max_events[[p]][[i]]>halfwayE, min(dataset[[p]][[i]][which(dataset[[p]][[i]]$cum_events==halfwayE), "cum_trt_events"]),
                       max(dataset[[p]][[i]][which(dataset[[p]][[i]]$cum_events==max_events[[p]][[i]]),"cum_trt_events"]))
    
    #Posterior parameters at IA (treatment group):
    post_a_t[[p]][i]<-1+dt[[p]][i]
    post_b_t[[p]][i]<-1+nt[[p]][i]-dt[[p]][i]
  }}

#Expected sample size at IA
avg_sample_size<-c()
for(p in 1:4){
  avg_sample_size[p]<-round(mean(unlist(nc[[p]]))+mean(unlist(nt[[p]])))
}

#Probability of stopping at IA
RRsamps<-vector("list",4)

#Obtain 2000 samples of p_t/p_c to obtain posterior probability that RR<1
for(p in 1:4){
  for(i in 1:Nsim){
    RRsamps[[p]][[i]]<-rbeta(2000,post_a_t[[p]][i],post_b_t[[p]][i])/
      rbeta(2000,post_a_c[[p]][i],post_b_c[[p]][i])
  }}

#Posterior probability that RR<1
postprob<-vector("list",4)

for(p in 1:4){
  for(i in 1:Nsim){
    postprob[[p]][[i]]<-mean(RRsamps[[p]][[i]]<1) 
  }}

#Probability of stopping at IA
probstop<-c()

for(p in 1:4){
  probstop[p]<-mean(postprob[[p]]>sup)
}


#Probability of stopping overall 
nc_final<-vector("list",4); dc_final<-vector("list",4)
post_a_c_final<-vector("list",4); post_b_c_final<-vector("list",4)
post_mean_c_final<-vector("list",4)

for(p in 1:4){
  for(i in 1:Nsim){
    if (postprob[[p]][[i]]<sup){
      nc_final[[p]][[i]]<-n/2;
      
      dc_final[[p]][[i]]<-min(dataset[[p]][[i]][which(dataset[[p]][[i]]$cum_ctrl_pts==n/2),"cum_ctrl_events"]);
      
      #Posterior parameters at final analysis (control group): 
      post_a_c_final[[p]][i]<-post_a_c[[p]][[i]]+dc_final[[p]][[i]];
      post_b_c_final[[p]][i]<-post_b_c[[p]][[i]]+nc_final[[p]][[i]]-dc_final[[p]][[i]];
      post_mean_c_final[[p]][i]<-post_a_c_final[[p]][[i]]/(post_a_c_final[[p]][[i]]+post_b_c_final[[p]][[i]])
    }
    else { nc_final[[p]][i]<-NA;
           dc_final[[p]][i]<-NA;
           post_a_c_final[[p]][i]<-NA;
           post_b_c_final[[p]][i]<-NA;
           post_mean_c_final[[p]][i]<-NA
    }}}

nt_final<-vector("list",4); dt_final<-vector("list",4)
post_a_t_final<-vector("list",4); post_b_t_final<-vector("list",4)
post_mean_t_final<-vector("list",4)


for(p in 1:4){
  for(i in 1:Nsim){
    if (postprob[[p]][[i]]<sup){
      nt_final[[p]][[i]]<-n/2;
      dt_final[[p]][[i]]<-min(dataset[[p]][[i]][which(dataset[[p]][[i]]$cum_ctrl_pts==n/2),"cum_trt_events"]);
      
      #Posterior parameters at final analysis (treatment group): 
      post_a_t_final[[p]][i]<-post_a_t[[p]][[i]]+dt_final[[p]][[i]];
      post_b_t_final[[p]][i]<-post_b_t[[p]][[i]]+nt_final[[p]][[i]]-dt_final[[p]][[i]];
      post_mean_t_final[[p]][i]<-post_a_t_final[[p]][[i]]/(post_a_t_final[[p]][[i]]+post_b_t_final[[p]][[i]])
    }
    else { nt_final[[p]][i]<-NA;
           dt_final[[p]][i]<-NA;
           post_a_t_final[[p]][i]<-NA;
           post_b_t_final[[p]][i]<-NA;
           post_mean_t_final[[p]][i]<-NA }}}

#Simulate 2000 p_t/p_c at every IA to obtain posterior probability that RR<1
RRsamps_final<-vector("list",4)

for(p in 1:4){
  for(i in 1:Nsim){
    if(!is.na(post_a_c_final[[p]][[i]])){
      RRsamps_final[[p]][[i]]<-rbeta(2000,post_a_t_final[[p]][[i]],
                                     post_b_t_final[[p]][[i]])/rbeta(2000,post_a_c_final[[p]][[i]],post_b_c_final[[p]][[i]])
    }
    else RRsamps_final[[p]][[i]]<-rep(NA,2000)
  }}

#Expected number of events at final IA
avg_d<-c()
for(p in 1:4){
  avg_d[p]<-round(mean(na.omit(unlist(dc_final[[p]])))+mean(na.omit(unlist(dt_final[[p]]))))
}

#Expected number of events overall
avg_d_overall<-probstop*halfwayE+(1-probstop)*avg_d

#Posterior probability that RR<1 at final analysis
postprob_final<-vector("list",4)

for(p in 1:4){
  for(i in 1:Nsim){
    postprob_final[[p]][i]<-mean(RRsamps_final[[p]][[i]]<1)
  }}


#Probability of stopping at final analysis
probstop_final<-c()
for(p in 1:4){
  probstop_final[p]<-mean(na.omit(unlist(postprob_final[[p]])>sup))
}

probstop_final<-probstop+(1-probstop)*probstop_final


#Total expected sample size at final analysis
avg_sample_size_final<-c()
for(p in 1:4){
  avg_sample_size_final[p]<-mean(na.omit(unlist(nc_final[[p]])))+mean(na.omit(unlist(nt_final[[p]])))
}

#Expected sample size overall
avg_n_overall<-probstop*avg_sample_size+(1-probstop)*avg_sample_size_final

# Results 
r1<-round(probstop,digits=2) # probability of stopping at IA (P_1)
r2<-avg_sample_size # expected sample size at IA (N_1)
r3<-rep(halfwayE,4) # expected number of events halfway  (D_1)
r4<-round(probstop_final, digits=2) # probability of stopping overall (P_2)
r5<-round(avg_n_overall) # expected sample size overall (N_2)
r6<-round(avg_d_overall) # expected number of events overall (D_2)


results_0.3_E<-data.frame(rbind(r1,r2,r3,r4,r5,r6))
rownames(results_0.3_E)<-c("P_1","N_1","D_1","P_2","N_2","D_2")
results_0.3_E

save(results_0.3_E, file="bin_0.3_E.RData")