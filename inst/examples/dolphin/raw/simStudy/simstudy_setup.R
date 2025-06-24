## Simulation study setup ##

source("JSsimfit_fun.R") #load function to simulate from the RPT model

# Simulation design --------------------

# structure of the capture occasions

T_simul <-  c(10,20,30,40) #different considered scenarios for each T

year_start <- vector("list",length = length(T_simul))
for(i in 1:length(T_simul)) year_start[[i]] <- c(1,1+1:i*T_simul[1])

# time lags between subsequent occasions

set.seed(123); t_lag_max <- rep(c(rgeom(T_simul[1]-1,0.05)+1, 240),length(T_simul))
t_lag_max <- t_lag_max[-T_simul[length(T_simul)]][-T_simul[length(T_simul)]]

t_lag_max_monthly <- t_lag_max/30 

t_lag_monthly<- list(t_lag_max_monthly[1:(T_simul[1]-1)],
                     t_lag_max_monthly[1:(T_simul[2]-1)],
                     t_lag_max_monthly[1:(T_simul[3]-1)],
                     t_lag_max_monthly[1:(T_simul[4]-1)])

# recruitment probabilities

set.seed(123); rho_r_max <- c(0.4,rep(0.0025,T_simul[length(T_simul)]-1)) 
rho_r_max[year_start[[length(T_simul)]][-c(1,length(year_start[[length(T_simul)]]))]] <- 0.02 #after a long period without sampling the population, rho_r increases proportional to lags between consequent occasions

set.seed(123); rho_p_max <-  c(0.4,rep(0.005,T_simul[length(T_simul)]-1)) 
rho_p_max[year_start[[length(T_simul)]][-c(1,length(year_start[[length(T_simul)]]))]] <- 0.04 #after a long period without sampling the population, rho_p increases proportional to lags between consequent occasions

rho_r<- list(rho_r_max[1:T_simul[1]],
             rho_r_max[1:T_simul[2]],
             rho_r_max[1:T_simul[3]],
             rho_r_max[1:T_simul[4]])
rho_p<- list(rho_p_max[1:T_simul[1]],
             rho_p_max[1:T_simul[2]],
             rho_p_max[1:T_simul[3]],
             rho_p_max[1:T_simul[4]])
rho_t<- list(rep(0.02,T_simul[1]),
             rep(0.02,T_simul[2]),
             rep(0.02,T_simul[3]),
             rep(0.02,T_simul[4]))

rho<- lapply(1:length(rho_r), function(i) rbind(rho_r[[i]],
                                                rho_p[[i]],
                                                rho_t[[i]]))

# inflation parameters

psi_r<- c(1-prod(1-rho_r[[1]]),
          1-prod(1-rho_r[[2]]),
          1-prod(1-rho_r[[3]]),
          1-prod(1-rho_r[[4]]))
psi_p<- c(1-prod(1-rho_p[[1]]),
          1-prod(1-rho_p[[2]]),
          1-prod(1-rho_p[[3]]),
          1-prod(1-rho_p[[4]]))
psi_t<- c(1-prod(1-rho_t[[1]]),
          1-prod(1-rho_t[[2]]),
          1-prod(1-rho_t[[3]]),
          1-prod(1-rho_t[[4]]))

# survival probabilities

phi <- c(0.01,0.997) #monthly survival

# detection probabilities

mu <- 0 

set.seed(123); tau_max <- rnorm(T_simul[length(T_simul)],0,0.5)
tau_max[T_simul[1]] <- -sum(tau_max[1:(T_simul[1]-1)])
tau_max[T_simul[2]] <- -sum(tau_max[1:(T_simul[2]-1)])
tau_max[T_simul[3]] <- -sum(tau_max[1:(T_simul[3]-1)])
tau_max[T_simul[4]] <- -sum(tau_max[1:(T_simul[4]-1)])

tau<- list(tau_max[1:T_simul[1]],
           tau_max[1:T_simul[2]],
           tau_max[1:T_simul[3]],
           tau_max[1:T_simul[4]]) #tau's sum to 0

p <- list(exp(mu+tau[[1]])/(1+exp(mu+tau[[1]])),
          exp(mu+tau[[2]])/(1+exp(mu+tau[[2]])),
          exp(mu+tau[[3]])/(1+exp(mu+tau[[3]])),
          exp(mu+tau[[4]])/(1+exp(mu+tau[[4]])))

delta <-  0.7

# mixture weights

w <-  c(0.2, 0.45, 0.35) 

# super-population abundance

M_sim <- 500 #we fix the size of the pseudo-population

Nsup_expect<- c(ceiling(M_sim*(w[1]*psi_r[1]+w[2]*psi_p[1]+w[3]*psi_t[1])),
                ceiling(M_sim*(w[1]*psi_r[2]+w[2]*psi_p[2]+w[3]*psi_t[2])),
                ceiling(M_sim*(w[1]*psi_r[3]+w[2]*psi_p[3]+w[3]*psi_t[3])),
                ceiling(M_sim*(w[1]*psi_r[4]+w[2]*psi_p[4]+w[3]*psi_t[4])))

## Expected abundance by group

# barplot(rbind(M_sim*w[1]*psi_r[1:4],
#               M_sim*w[2]*psi_p[1:4],
#               M_sim*w[3]*psi_t[1:4]),
#         ylim=c(0,200),ylab="",main="Expected abundance",xlab="Scenario",names.arg = 1:4,beside=T,col=1:3)
# legend("topleft",c("Resident","Part-time","Transient"),horiz = T,bty="n",pch=15,col=1:3)

# Simulating datasets --------------------------

set.seed(123); seeds <- sample(x=1000:9999,size=9000,replace=F)

K <- 50 #number of datasets for each scenario
sim <- vector("list",K) #it will contain a list of K datasets for a particular scenario

#N.B. variable "scenario" is defined in files "simstudy_JStype.R" and "simstudy_RPT.R"

for(j in 1:K){
  
  set.seed(seeds[j])
  sim[[j]]<- sim.RPT(M.simul = M_sim,
                     T.simul = T_simul[scenario],
                     time.lag = t_lag_monthly[[scenario]],
                     year.start = year_start[[scenario]],
                     rho = rho[[scenario]],
                     phi = phi,
                     mu = mu,
                     tau = tau[[scenario]],
                     delta = delta,
                     w = w,
                     nzeros=500)
  
}

assign(paste0("sim_rpt_",scenario),sim) #we assign this name to the content of sim 
rm(sim)