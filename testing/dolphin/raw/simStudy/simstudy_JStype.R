## Script to fit JS-type models on the K simulated datasets ##

scenario <- 1 #changing this value, we change the scenario
modNumber <- 1  #changing this value, we change the model to fit

source("simStudy/simstudy_setup.R") #load the simulation experiment setup

# Model fitting (parallel computing is required) -------------

competingModels<- data.frame(shortname=paste0("mod",1:10),
                             BUGSname=c(rep("phi,p_t.txt",1),
                                        rep("phi,p_tg.txt",3),
                                        rep("phi_g,p_t.txt",3),
                                        rep("phi_g,p_tg.txt",3)),
                             numComp=c(1,rep(2:4,3)))

require(doParallel)
require(doRNG)

s.init = Sys.time()

set.seed(123) #set a generation code for reproducibility

num_cores <- 5 #SELECT THE NUMBER OF CORES TO USE
cl <- parallel::makeCluster(num_cores, outfile="") 
doParallel::registerDoParallel(cl)

out <- foreach(j = 1:K, .packages=c("R2jags","rjags"), .export=c(paste0("sim_rpt_",scenario)), .inorder=T) %dorng%{
  
  
  print(paste0("Fitting model to dataset #",j, " from scenario ",scenario))
  gc()
  temp<- JStype.fit.jags(CR.data.matrix = get(paste0("sim_rpt_",scenario))[[j]]$y,
                         t_lag = get(paste0("sim_rpt_",scenario))[[j]]$time.lag,
                         year_start = get(paste0("sim_rpt_",scenario))[[j]]$year.start,
                         G = competingModels$numComp[modNumber],
                         bugs_model=competingModels$BUGSname[modNumber],
                         nc = 2,
                         sample = 2e4,
                         burnin = 5e3, 
                         thin = 2)
  
  temp
  
} 
stopCluster(cl)

print(Sys.time()-s.init)

assign(paste0("fit_mod",modNumber,"_",scenario),out) #we assign this name to the content of out
rm(out)

#fold_name = paste0("simStudy/rpt_mod",modNumber,"_",scenario,".RData")
#save.image(fold_name) #to save the output