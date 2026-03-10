## Fitting JS-type models on real data ##

# Loading files (a) containing data and (b) functions for model fitting ----

source("realData/dolphins_uneven_3years.R") #(a)
source("JSsimfit_fun.R") #(b)

# JS-type competing models -------------------------------------------------

competingModels <- data.frame(shortname=paste0("mod",1:10),
                              BUGSname=c(rep("phi,p_t.txt",1),
                                         rep("phi,p_tg.txt",3),
                                         rep("phi_g,p_t.txt",3),
                                         rep("phi_g,p_tg.txt",3)),
                              numComp=c(1,rep(2:4,3)))

# Model fitting ------------------------------------------------------------

modNumber <- 1  #change this value to change the model to fit

### THE FOLLOWING PART OF THE CODE MAKES USE OF INTENSIVE COMPUTATIONS,
### FOR WHICH IT IS SUGGESTED TO USE HIGH PERFORMANCE COMPUTING.

set.seed(123)
assign(as.character(competingModels$shortname[modNumber]),
       JStype.fit.jags(CR.data.matrix = datJS.3years.aug,
                       t_lag=time_lags_3years,
                       year_start=year_start_3years,
                       G=competingModels$numComp[modNumber],
                       bugs_model=competingModels$BUGSname[modNumber],
                       nc = 2,
                       sample = 2e4,
                       burnin = 5e3, 
                       thin = 2))

#save.image(paste0("dolph_mod",modNumber,".RData")) #to save the output