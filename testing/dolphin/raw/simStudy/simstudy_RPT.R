## Script to fit RPT model on the K simulated datasets ##

scenario <- 1  #changing this value, we change the scenario

source("simStudy/simstudy_setup.R") #load the simulation experiment setup

# Model fitting (parallel computing is required) -------------

require(doParallel)
require(doRNG)

s.init = Sys.time()

set.seed(123) #set a generation code for reproducibility

num_cores <- 5 #SELECT THE NUMBER OF CORES TO USE
cl <- parallel::makeCluster(num_cores, outfile="") 
doParallel::registerDoParallel(cl)

out <-  foreach(j = 1:K, .packages=c("R2jags","rjags"), .export=c(paste0("sim_rpt_",scenario)), .inorder=T) %dorng%{
  
  
  print(paste0("Fitting model to dataset #",j, " from scenario ",scenario))
  #gc()
  temp<- JS.RPT.fit.jags(CR.data.matrix = get(paste0("sim_rpt_",scenario))[[j]]$y,
                         t_lag = get(paste0("sim_rpt_",scenario))[[j]]$time.lag,
                         year_start = get(paste0("sim_rpt_",scenario))[[j]]$year.start,
                         nc = 2,
                         sample = 2e4,
                         burnin = 5e3, 
                         thin = 2,
                         pars.to.save=c("w","p","phi","rho","Nsuper","N.y","clust","loglik_i","mu","delta"))
  
  temp
  
} 
stopCluster(cl)

print(Sys.time()-s.init)

assign(paste0("fit_RPT_",scenario),out) #we assign this name to the content of out
rm(out)

#fold_name = paste0("simStudy/rpt_",scenario,".RData")
#save.image(fold_name) #to save the output