library(tidyverse)
library(here)

setwd(here("testing","turtle","bernoulli","raw"))

honu_ch <- read.delim("yearly_detection.txt",sep=" ")[,-1]
id <- read.delim("yearly_detection.txt",sep=" ")[,1]
honu_ch$X2013 <- 0
honu_ch$X2020 <- 0
yr <- colnames(honu_ch) %>% stringr::str_remove("X") %>% as.numeric()
honu_ch <- honu_ch[,order(yr)]
yr <- yr[order(yr)]
# ret_ch <- honu_ch[,which(yr>= 1986 & yr<=2016)]
ret_ch <- honu_ch[,which(yr<=2016)]
# capt <- 1.0*(rowSums(honu_ch[,which(yr<1986)])>0)
out_idx <- rowSums(ret_ch)>0
ret_ch <- ret_ch[out_idx,]
# capt <- capt[out_idx]
# id <- id[out_idx]

skip_forage <- function(x){
  rrr <- rle(x)
  out <- max(rrr$lengths[rrr$values==1])
  return(out>1)
}
idx <- as.logical(1-apply(ret_ch, 1, skip_forage))
ret_ch <- ret_ch[idx, ]
# capt <- capt[idx]
# id <- id[idx]

# first <- rep(0,nrow(honu_ch))
# for(i in 1:nrow(honu_ch)){
#   f[i] <- min(which(honu_ch[i,]==1))
# }

ret_ch <- 2*ret_ch + 1
for(i in 1:nrow(ret_ch)){
    f <- min(which(ret_ch[i,]==3))
    ret_ch[i,f] <- 2
}

write_csv(ret_ch, file="../honu_ch.csv")




