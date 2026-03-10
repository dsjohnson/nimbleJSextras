## Loading data ##

# Packages loading --------------------------------------------------------

require(tidyverse)
require(magrittr)
require(lubridate)

# Load .csv file containing the real data ---------------------------------

datJS.3years <- read.csv("realData/bottlenoseDolph.csv",
                         check.names = "F") #names of the columns are dates

T.3years <- ncol(datJS.3years)

# Computing time lags -----------------------------------------------------

cap_occ <- datJS.3years %>% 
  gather(key="Date", value="Value") %>% 
  arrange(Date) %>%  
  distinct(Date) %>% 
  mutate(TDiffDays = as.numeric(difftime(Date, lag(Date, default = Date[1]), units = 'days'))) #capture occasion

time_lags_3years <- cap_occ$TDiffDays[-1]/365  #T-1 time lags btw capture occasions (survival in yearly scale) 

year_start_3years <- c(1,16,51,T.3years+1) #useful to keep trace of the first capture of each year

# Building the augmented dataset for fitting with PX-DA approach ----------

nzeros <- 500 
datJS.3years.aug <- rbind(as.matrix(datJS.3years),
                          matrix(rep(0,nzeros*T.3years),ncol=T.3years)) #augmented data matrix
M <- nrow(datJS.3years.aug) #number of individuals in the augmented dataset