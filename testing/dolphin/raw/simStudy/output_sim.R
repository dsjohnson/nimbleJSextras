## DISCLAIMER: The .RData output analysed in the current script are available at
## https://drive.google.com/drive/folders/1OMiyAP8iupf6FD_vOhSL_aI9woP9HsD7?usp=sharing .
## Following the link above, please download the folder "outputSim" 
## and add it in the same folder as the current R script.
##
## However, due to the very large size of the RData files, 
## we suggest to directly load the file "simResults.RData"
## to see the output of the code which follows (large-sized objects have been removed).
## N.B. Please, load the following packages first.

#### Required packages

require(HDInterval)
require(tibble)
require(ggplot2)
require(dplyr)

#### Simulation study

modNames = c("RPT",paste0("mod",1:10)) #competing models names

estErrorN = matrix(rep(NA, 4*length(modNames)*50), ncol=4*length(modNames)) #matrix for boxplots of Figure 2

tab1 = matrix(rep(NA,11*20),nrow=11) #matrix which will store Table 1 in the paper
colnames(tab1) = c(sapply(seq(10,40,10), function(i) paste0(c("MAErel","Cov","CIWrel","WAIC","waicPerc"),"_T",i)))
rownames(tab1) = modNames

tab2 = matrix(rep(NA,4*12),nrow=4) #matrix which will store Table 2 in the paper
colnames(tab2) = c(sapply(1:4, function(i) paste0(c("MAE","Cov","CIW"),i)))
rownames(tab2) = paste0("T=",seq(10,40,10))

for(s in 1:4){ #four scenarios considered T=10,20,30,40
  
  ## RPT model
  load(paste0("simStudy/outputSim/rpt_",s,".RData")) #change the directory accordingly
  
  ### TABLE 2
  
  paramSelect = c("phi[1]","phi[2]","delta","mu")
  trueValParam = c(get(paste0("sim_rpt_",s))[[1]]$phi,
                   get(paste0("sim_rpt_",s))[[1]]$delta,
                   get(paste0("sim_rpt_",s))[[1]]$mu) #true values of the parameters (yet introduced in the main text of the paper)
  
  for(p in 1:length(paramSelect)){
    
    # Chain matrices of the selected parameters 
    # (MCMC iterations on rows, replicas on columns)
    
    chainsMat = sapply(get(paste0("fit_RPT_",s)), function(x) x$chains_mat[,paramSelect[p]])
    
    # Posterior estimate (medians) and 95% CI of the parameter of interest for each replica
    
    postEst = apply(chainsMat,2,median)
    postCI = apply(chainsMat,2,hdi,credMass=0.95)
    
    # MAE 
    
    tab2[s,paste0("MAE",p)] = mean(abs(postEst-trueValParam[p]))
    
    # Cov
    
    tab2[s,paste0("Cov",p)] = mean(apply(postCI, 2, function(x) x[1]<=trueValParam[p] & x[2]>=trueValParam[p]))
    
    # CIW
    
    tab2[s,paste0("CIW",p)] = mean(diff(postCI))
    
  }
  
  ### TABLE 1
  
  trueValNsuper = sapply(get(paste0("sim_rpt_",s)), function(x) x[["Nsuper"]]) #vector of the true value of the parameter of interest for each replica
  waic = matrix(rep(NA,length(modNames)*50),ncol=50) #matrix with WAIC for each model (rows) and each replica (columns)
  
  for(m in 1:length(modNames)){
    
    if(modNames[m]!="RPT"){
      load(paste0("simStudy/outputSim/rpt_",modNames[m],"_",s,".RData")) #change the directory accordingly
    }
    
    modelFit = get(paste0("fit_",modNames[m],"_",s))
    rm(list=paste0("fit_",modNames[m],"_",s))
    gc()
    
    chainsMat = sapply(modelFit, function(x) x$chains_mat[,"Nsuper"])
    
    # Posterior estimate (medians) and 95% CI of the parameter of interest for each replica
    
    postEst = apply(chainsMat,2,median)
    postCI = apply(chainsMat,2,hdi,credMass=0.95)
    
    # relative MAE (Nsuper) - N.B. expected value of Nsuper_t's yet given in the main text of the paper
    
    estErrorN[,4*(m-1)+s] = postEst-trueValNsuper #estimation error (Nsuper)
    
    tab1[m,paste0("MAErel_T",seq(10,40,10)[s])] = mean(abs(estErrorN[,4*(m-1)+s]))/Nsup_expect[s] 
    
    # Cov (Nsuper)
    
    tab1[m,paste0("Cov_T",seq(10,40,10)[s])] = mean(mapply(function(x,y,z) x<=z & y>=z,
                                                           postCI[1,],
                                                           postCI[2,],
                                                           trueValNsuper))
    
    # relative CIW (Nsuper)
    
    tab1[m,paste0("CIWrel_T",seq(10,40,10)[s])] = mean(diff(postCI))/Nsup_expect[s]
    
    # median WAIC
    
    waic[m,] = sapply(modelFit, function(x) x$WAIC$estimates["waic",1])
    tab1[m,paste0("WAIC_T",seq(10,40,10)[s])] = median(waic[m,])
    
    rm(modelFit)
    gc()
    
  }
  
  # Percentual of times best waic
  
  tab1[,paste0("waicPerc_T",seq(10,40,10)[s])] = table(factor(apply(waic,2,which.min),levels=1:length(modNames)))/ncol(waic)*100
  
}

write.csv(tab1,"table1.csv")
write.csv(tab2,"table2.csv")

######## Figure 2 #######

jpeg("Figure2.jpg", width = 1200, height = 700, res = 150)

tibble(Model = rep(c("RPT", paste0("M", 1:10)), each = 50*4), ee = c(estErrorN), TT = rep(rep(c("T=10","T=20","T=30","T=40"), each = 50), 11)) %>% 
  mutate(Model = factor(Model, levels = c("RPT", paste0("M", 1:10))),
         Nexp = case_when(TT == "T=10"~170, TT == "T=20"~209, TT == "T=30"~243, TT == "T=40"~271), # expected Nsuper varying T
         ee = ee/Nexp) %>% # compute MAE rel
  ggplot(aes(Model, ee, fill = factor(TT))) +
  geom_boxplot(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = I("grey")) +
  labs(x = "", y = "Relative estimation error", fill = "") +
  scale_fill_manual(values = c("lightblue", "skyblue3", "blue", "darkblue")) +
  theme_bw() +
  theme(legend.position = "top", text = element_text(size = 18))

dev.off()