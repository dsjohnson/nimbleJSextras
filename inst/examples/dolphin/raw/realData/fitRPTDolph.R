## Fitting RPT model on real data ##

# Loading files (a) containing data and (b) functions for model fitting ----

source("realData/dolphins_uneven_3years.R") #(a)
source("JSsimfit_fun.R") #(b)

# Model fitting ------------------------------------------------------------

### THE FOLLOWING PART OF THE CODE MAKES USE OF INTENSIVE COMPUTATIONS,
### FOR WHICH IT IS SUGGESTED TO USE HIGH PERFORMANCE COMPUTING.

set.seed(123)
res.RPT.fit <- JS.RPT.fit.jags(CR.data.matrix = datJS.3years.aug,
                               t_lag=time_lags_3years,
                               year_start=year_start_3years,
                               nc = 2,
                               sample = 2e4,
                               burnin = 5e3, 
                               thin = 2,
                               pars.to.save = c("w","p","phi","rho","Nsuper","N.y","clust","loglik_i","z","mu","delta"))


#save.image("dolph_RPT.RData") #to save the output