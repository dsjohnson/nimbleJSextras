## DISCLAIMER: The .RData output analysed in the current script are available at
## https://drive.google.com/drive/folders/1OMiyAP8iupf6FD_vOhSL_aI9woP9HsD7?usp=sharing .
## Following the link above, please download the folder "outputDolph" 
## and add it in the same folder as the current R script.

#### Required packages

require(HDInterval)
require(tibble)
require(ggplot2)
require(dplyr)

##### Real data

modNames = c("RPT",paste0("mod",1:10)) #competing models names

tab3 = matrix(rep(NA,length(modNames)*4),nrow=length(modNames)) #matrix which will store Table 1 in the paper
colnames(tab3) = c("WAIC","Nsuper_median","Nsuper_lb","Nsuper_ub") 
rownames(tab3) = modNames

for(m in 1:length(modNames)){
  
  load(paste0("realData/outputDolph/dolph_",modNames[m],".RData")) #change the directory accordingly
  
  if(modNames[m]=="RPT"){
    
    tab3[m,1] = get(paste0("res.",modNames[m],".fit"))$WAIC$estimates["waic",1] #WAIC
    tab3[m,2] = median(get(paste0("res.",modNames[m],".fit"))$chains_mat[,"Nsuper"]) #Nsuper posterior median
    tab3[m,3:4] = hdi(get(paste0("res.",modNames[m],".fit"))$chains_mat[,"Nsuper"]) #Nsuper 95% posterior HDI 
    
    rm(list=paste0("res.",modNames[m],".fit"))
    gc()
    
  }else{
    
    tab3[m,1] = get(modNames[m])$WAIC$estimates["waic",1] #WAIC
    tab3[m,2] = median(get(modNames[m])$chains_mat[,"Nsuper"]) #Nsuper posterior median
    tab3[m,3:4] = hdi(get(modNames[m])$chains_mat[,"Nsuper"]) #Nsuper 95% posterior HDI 
    
    rm(list=modNames[m])
    gc()
    
  }
  
}

tab3_ord = tab3[order(tab3[,1]),] #sorted by ascending WAIC (Table 3)

write.csv(tab3_ord,"table3.csv")
