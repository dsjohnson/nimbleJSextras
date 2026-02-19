library(tidyverse)
library(here)

setwd(here("testing","turtle","poisson","raw"))

honu_ch <- read.delim("yearly_total_sightings.txt",sep=" ")[,-1]
honu_ch$X2013 <- 0
honu_ch$X2020 <- 0
yr <- colnames(honu_ch) %>% stringr::str_remove("X") %>% as.numeric()
honu_ch <- honu_ch[,order(yr)]
honu_ch <- honu_ch[,which(yr>= 2002)]
honu_ch <- honu_ch[rowSums(honu_ch)>0,]

write_csv(honu_ch, file="../honu_ch_pois.csv")




