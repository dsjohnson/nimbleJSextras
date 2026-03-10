## DISCLAIMER: The .RData output analysed in the current script are available at
## https://drive.google.com/drive/folders/1OMiyAP8iupf6FD_vOhSL_aI9woP9HsD7?usp=sharing .
## Following the link above, please download the folder "outputDolph" 
## and add it in the same folder as the current R script.

# Packages
remotes::install_github("jschoeley/tricolore")
packs <- c("tidyverse", "magrittr", "lubridate", "ggtern", "tricolore", "cowplot")
load_packs <- sapply(packs, require, character.only = T)
if(any(load_packs==F)) install.packages(packs[!load_packs])
load_packs <- sapply(packs, require, character.only = T)

Sys.setlocale(locale = "C")

## DISCLAIMER!!! ##
# Note that most of the outputs will only work with rpt based models being the best performing in terms of WAIC.
# Hence final results only rely on the best fitting model
###################

# Output analysis ---------------------------------------------------------
# Create object with names of all output files from the real data analysis
out_file <- list.files(path = "realData/outputDolph/", full.names = T)
out_file <- out_file[c(1:2, 4:11, 3)] # reorder: RPT first, than from M1 to M10

load(out_file[1]) # load j-th model

D <- M - nzeros

######### Plots in Section Motivating Example to describe the data ######### 
dat_summ <- datJS.3years %>% as.data.frame() %>% 
  rownames_to_column(var = "id") %>% 
  pivot_longer(cols = -id, names_to = "Date", values_to = "cr") %>% 
  mutate(Date = lubridate::ymd(Date), Year = lubridate::year(Date)) %>% 
  filter(cr == 1) %>% 
  group_by(id) %>% 
  summarise(Date = first(Date)) %>% 
  group_by(Date) %>% count

# Yearly number of new identifications
dat_summ %>% group_by(lubridate::year(Date)) %>% summarise(Tot = sum(n)) %>% ungroup() %>% mutate(Perc = Tot/sum(Tot))
# Monthly by year
dat_summ %>% group_by(lubridate::year(Date), lubridate::month(Date)) %>% summarise(Tot = sum(n)) %>% ungroup %>% arrange(desc(Tot)) 

### Figure 1a: Cumulative discorvery rate ###
jpeg("CumulIdRate.jpg", width = 1200, height = 900, res = 150)
dat_summ %>% 
  ggplot(aes(Date, cumsum(n), size = n)) +
  geom_point() +
  labs(x = "", y = "Cumulative identification", size = "New individuals") +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = "bottom")
dev.off()

### Figure 1b: Total number of captures ###
dat_summ <- datJS.3years %>% as.data.frame() %>% 
  apply(., 1, function(x) {
    oo = names(x)[x == 1]
    c(as.numeric(as.Date(oo[length(oo)]) - as.Date(oo[1])), max(cumsum(unlist(x))))
  })
jpeg("CRHist.jpg", width = 1200, height = 900, res = 150)
tibble(NDays = dat_summ[1,], cr = dat_summ[2,]) %>% 
  ggplot(aes(x = NDays, y = cr)) +
  geom_point(size = 2) +
  labs(x = "N. of living days", y = "Total number of captures") +
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()


# Model results -----------------------------------------------------------
# Name of the model object
mod_name <- ls()[grepl("res\\.|^mod[0-9]$", ls())]
est_mod <- get(mod_name)
rm(list = ls()[grepl("res\\.|^mod[0-9]$", ls())], pos = ".GlobalEnv")

est_mod$WAIC$waic # WAIC

###### Inference on uncaught individuals
oo <- matrix(NA, nrow = nrow(est_mod$chains_mat), ncol = 3)
Zchains <- est_mod$chains_mat[, grepl("z", colnames(est_mod$chains_mat))]
for(i in 1:nrow(oo)){ # Very slow!!!
  Zmat <- matrix(Zchains[i,], nrow = M, byrow = F) # create matrix with latent indicator Z at each iteration
  id_z <- apply(Zmat, 1, function(x) any(x==1)) # Get the id of the individuals who have been part of the population
  chains_app <- est_mod$chains_mat[i, grep("clust", colnames(est_mod$chains_mat))][id_z][-c(1:D)] # Get the labels of the uncaught
  chains_app <- factor(chains_app, levels = 1:3) # Check classification
  oo[i,] <- table(chains_app)/sum(table(chains_app))
  #print(i)
}
# Posterior allocation of the uncaught individuals to each of the sub-populations
colMeans(oo)
Rfast2::colQuantile(oo, .5)


###### Estimates of N_super, N_t and N_gt
(Nsuper <- c(est_mod$chains_mat[, "Nsuper"] %>% HDInterval::hdi(), median(est_mod$chains_mat[, "Nsuper"])))
# Yearly Nsuper
(Nt <- apply(est_mod$chains_mat[, grepl("^N.y", colnames(est_mod$chains_mat))], 2, function(x) c(HDInterval::hdi(x), median(x))))

########### Figure 3a ########### 
summ_Nt <- est_mod$chains_mat[, grepl("^N.y", colnames(est_mod$chains_mat))] %>% as.data.frame() %>% 
  summarise(across(everything(), function(x) c(HDInterval::hdi(x), median(x)))) %>% 
  mutate(type = c("q1", "q2", "med")) %>% 
  pivot_longer(names_to = "Year", values_to = "N", cols = -type) %>% 
  mutate(TT = rep(2018:2020, 3))

jpeg("Nocc_dolph_rptlogitp.jpg", width = 900, height = 600, res = 150)
ggplot() + 
  geom_segment(data = summ_Nt %>% filter(type != "med") %>% spread(type, N), 
               aes(x = TT, xend = TT, y = q1, yend = q2), 
               linetype = "dashed", linewidth = 1, color = I("grey30")) +
  geom_point(data = summ_Nt %>% filter(type == "med"), aes(x = TT, y = N), size = 4) +
  scale_x_continuous(breaks = 2018:2020, labels = 2018:2020) + 
  labs(x = "year", y = expression(paste(hat(N)[y]))) +
  theme_bw() +
  theme(text = element_text(size = 18))
dev.off() 


# Ngt
niter <- nrow(est_mod$chains_mat)
ind_alive_18 <- ind_alive_19 <- ind_alive_20 <- matrix(rep(NA,niter*M), ncol=M)
for(i in 1:M){
  ind_alive_18[,i] <- rowSums(est_mod$chains_mat[,paste0("z[",i,",",1:15,"]")])>0
  ind_alive_19[,i] <- rowSums(est_mod$chains_mat[,paste0("z[",i,",",16:50,"]")])>0
  ind_alive_20[,i] <- rowSums(est_mod$chains_mat[,paste0("z[",i,",",51:87,"]")])>0
}
clust_chains <- est_mod$chains_mat[,paste0("clust[",1:M,"]")]

# Divided by group and year
Nsup_18_R_chain <- rowSums(ind_alive_18*(clust_chains==1))
Nsup_19_R_chain <- rowSums(ind_alive_19*(clust_chains==1))
Nsup_20_R_chain <- rowSums(ind_alive_20*(clust_chains==1))

Nsup_18_P_chain <- rowSums(ind_alive_18*(clust_chains==2))
Nsup_19_P_chain <- rowSums(ind_alive_19*(clust_chains==2))
Nsup_20_P_chain <- rowSums(ind_alive_20*(clust_chains==2))

Nsup_18_T_chain <- rowSums(ind_alive_18*(clust_chains==3))
Nsup_19_T_chain <- rowSums(ind_alive_19*(clust_chains==3))
Nsup_20_T_chain <- rowSums(ind_alive_20*(clust_chains==3))

lapply(list(Nsup_18_R_chain,Nsup_18_P_chain,Nsup_18_T_chain,
            Nsup_19_R_chain,Nsup_19_P_chain,Nsup_19_T_chain,
            Nsup_20_R_chain,Nsup_20_P_chain,Nsup_20_T_chain), 
       summary)
# Put Ngt iterations into a dataframe
oo <- as.data.frame(cbind(Nsup_18_R_chain,Nsup_18_P_chain,Nsup_18_T_chain,
                          Nsup_19_R_chain,Nsup_19_P_chain,Nsup_19_T_chain,
                          Nsup_20_R_chain,Nsup_20_P_chain,Nsup_20_T_chain))

########### Figure 3b ########### 
jpeg("Ngt_dolph_rptlogitp.jpg", width = 900, height = 600, res = 150)
oo %>% 
  pivot_longer(cols = everything(), values_to = "val", names_to = "nam") %>% 
  mutate(Y = str_sub(nam, 6,7), Type = str_sub(nam, 9, 9), 
         Year = case_when(Y == "18" ~ 2018, Y == "19" ~ 2019, Y == "20" ~ 2020),
         Type = case_when(Type == "R" ~ "Residents", Type == "P" ~ "Part-time", Type == "T" ~ "Transients"),
         Type = factor(Type, levels = c("Residents", "Part-time", "Transients"))) %>% 
  ggplot(aes(x = factor(Year), y = val, fill = Type)) +
  geom_boxplot(position = "dodge") +
  scale_fill_manual(values = c("skyblue2", "pink2", "gold")) +
  labs(x = "year", y = expression(paste(hat(N)[gt])), fill = "") +
  theme_bw() +
  theme(legend.position = "top", text = element_text(size = 18))
dev.off()
rm(ind_alive_18, ind_alive_19, ind_alive_20, Nsup_18_P_chain, Nsup_18_R_chain, Nsup_18_T_chain,
   Nsup_19_P_chain, Nsup_19_R_chain, Nsup_19_T_chain, Nsup_20_P_chain, Nsup_20_R_chain, Nsup_20_T_chain)
gc()

##### Posterior estimates of the mixture's weights ############
G <- max(clust_chains) #number of mixture components
#relative frequency distribution of the group labels for each individual
rel_freq_estimated_labels <- matrix(unlist(apply(clust_chains,2,function(x) table(factor(x,levels=1:G))/length(x))), ncol=G,byrow=T)
colMeans(rel_freq_estimated_labels) # Allocation with the cluster labels
colMeans(est_mod$chains_mat[,paste0("w[",1:G,"]")]) # Mixture's weights


########### Figure 4 ########### 
# Clustering observed individuals
clust_post <- clust_chains[,1:nrow(datJS.3years)] %>% set_colnames(rownames(datJS.3years))
clust_post_prob <- clust_post %>% as.data.frame() %>%
  gather(ID, Value) %>%
  group_by(ID) %>%
  summarise(class_1 = mean(Value == 1), class_2 = mean(Value == 2), class_3 = mean(Value == 3))
clust_post_lab <- tibble(ID = rownames(datJS.3years), Group = clust_post_prob[,-1] %>% apply(1, which.max))
nms_mat <- expand.grid(1:nrow(datJS.3years), 1:ncol(datJS.3years))
nms_z <- paste0("z[", nms_mat$Var1, ",", nms_mat$Var2, "]")

# Get latent variables Z
z_post <- est_mod$chains_mat[, nms_z] %>% colMeans() %>% round() %>%
  matrix(., nrow = nrow(datJS.3years), byrow = F) %>% as.data.frame() %>%
  set_colnames(value = colnames(datJS.3years)) %>%
  mutate(Group = clust_post_lab$Group, ID = rownames(datJS.3years))

# Clustered capture histories after classification
ids <- split(z_post[1:D,1:T.3years], 1:nrow(z_post[1:D,])) %>%
  map(function(x){
    if(any(x) == 1){
      id <- which(x == 1)
      id2 <- c(id[1], id[length(id)])
    }else{
      id2 <- names(x[[1]])
    }
    return(id2)
  })

# Intermediate dataset
dapp <- map2_dfr(split(datJS.3years, 1:nrow(datJS.3years)), ids, function(x,y){
  if(!is.null(y)){
    histobs <- x[y[1]:y[2]]
    living <- as.numeric(as.Date(colnames(datJS.3years)[y[2]]) - as.Date(colnames(datJS.3years)[y[1]]))
    oo <- tibble(hist = sum(as.numeric(histobs)), time = living)
  }else{
    oo <- NULL
  }
  return(oo)
  
}) %>% mutate(Group = clust_post_lab$Group[1:D]) %>% bind_cols(clust_post_prob[1:D,])

colors_and_legend <- Tricolore(as.data.frame(dapp) %>% rename(R = class_1, P = class_2, `T` = class_3),
                               'R', 'P', "T", breaks = 100, hue = .6, chroma = .5)  
dapp$rbg_disc <- colors_and_legend$rgb

# Figure 4: left panel
freqp <- dapp %>% 
  ggplot(aes(x = time, y = hist, color = rbg_disc, shape = factor(Group))) +
  geom_point(size = 1.5, alpha = .6) +
  scale_color_identity() +
  labs(y = "Total number of captures", x = "N. of living days") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20), text = element_text(size = 18),
        strip.text = element_text(size = 16), legend.position = "none")

# Intermediate dataset for the remaining plot
dat <- map2(split(datJS.3years, 1:nrow(datJS.3years)), ids, function(x,y){
  histobs <- x[y[1]:y[2]]
  tibble(hist = cumsum(as.numeric(histobs)), time = colnames(datJS.3years)[y[1]:y[2]])
}) %>%
  map2_dfr(rownames(datJS.3years), function(dat, nm) dat %>% mutate(ID = nm)) %>%
  group_by(ID) %>%
  mutate(hist = hist + runif(1, 0.3, 0.6), time = ymd(time) + sample(c(-1,1),1)) %>%
  left_join(clust_post_prob) %>%
  left_join(clust_post_lab) %>%
  ungroup() %>%
  rename(R = class_1, P = class_2, `T` = class_3)

colors_and_legend <- Tricolore(as.data.frame(dat), 'R', 'P', "T", breaks = 100, hue = .6, chroma = .5)  
dat$rbg_disc <- colors_and_legend$rgb
# Figure 4: central panel
p <- dat %>%
  mutate(ID = factor(ID, ordered = F, levels = clust_post_prob %>% 
                       arrange(class_1) %>% pull(ID)),
         Group = ifelse(Group == 1, "Residents", ifelse(Group == 2, "Part-time", "Transients")),
         Group = factor(Group, levels = c("Residents", "Part-time", "Transients"))) %>%
  ggplot(aes(x = time, y = hist, group = ID, color = rbg_disc)) +
  geom_line(linewidth = .1) +
  geom_point(size = .5) +
  scale_color_identity() +
  facet_wrap(~Group) +
  labs(y = "Cumulative capture frequencies", x = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), text = element_text(size = 18), strip.text = element_text(size = 16))

# Create Figure 4 with legend
jpeg("ClassPlot_CR.jpg", width = 1200, height = 400, res = 100)
cowplot::plot_grid(freqp, p, ggplotGrob(colors_and_legend$key), rel_widths = c(.3, .55, .15), nrow = 1)
dev.off()


# Appendix ----------------------------------------------------------------
########### Figure 8 (Appendix) ########### 
jpeg("Posterior_pt_rptlogitp.jpg", width = 1200, height = 700, res = 150)
est_mod$chains_mat[, grepl("^p\\[", colnames(est_mod$chains_mat))] %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  mutate(Iter = 1:nrow(.)) %>% 
  pivot_longer(cols = -Iter, values_to = "value", names_to = "pars") %>% 
  mutate(t = parse_number(pars)) %>% 
  group_by(t) %>% 
  summarise(Q1 = HDInterval::hdi(value)[1], M = median(value), Q2 = HDInterval::hdi(value)[2]) %>% ungroup() %>% 
  mutate(dates = lubridate::ymd(colnames(datJS.3years))) %>% 
  ggplot(aes(x = dates)) +
  geom_segment(aes(xend = dates, y = Q1, yend = Q2), linetype = "dashed", color = I("grey40")) +
  geom_point(aes(y = M), size = 1) +
  labs(x = "", y = expression(paste(hat(p)[NP][","][t]))) +
  scale_x_date(date_labels = "%Y/%b", date_breaks = "3 months") +
  theme_bw() +
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust = .9, vjust = .9))
dev.off()

########### Figure 9 abc (Appendix) ########### 
pick_rho <- 1 # 1,2,3
nms <- c("R", "P", "T")
jpeg(paste0("Posterior_rhot_rptlogitp", pick_rho, ".jpg"), width = 1200, height = 700, res = 150)
est_mod$chains_mat[, grepl("^rho\\[", colnames(est_mod$chains_mat))] %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  mutate(Iter = 1:nrow(.)) %>% 
  pivot_longer(cols = -Iter, values_to = "value", names_to = "pars") %>% 
  mutate(G = str_sub(pars, 5, 5), t = parse_number(str_sub(pars, 7, length(pars)))) %>% 
  group_by(G, t) %>% 
  summarise(Q1 = HDInterval::hdi(value, credMass = .95)[1], M = median(value),
            Q2 = HDInterval::hdi(value, credMass = .95)[2]) %>% ungroup() %>%
  mutate(dates = rep(lubridate::ymd(colnames(datJS.3years)), 3)) %>% 
  filter(G == pick_rho) %>% 
  ggplot(aes(x = dates)) +
  geom_segment(aes(xend = dates, y = Q1, yend = Q2), linetype = "dashed", color = I("grey40")) +
  geom_point(aes(y = M), size = 1) +
  labs(x = "", y = expression(paste(hat(rho)[R][","][t]))) +
  scale_x_date(date_labels = "%Y/%b", date_breaks = "3 months") +
  scale_y_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust = .9, vjust = .9))
dev.off()


##################### TRACEPLOTS AND DENSITIES IN THE APPENDIX ######################################
jpeg("Traceplot_Nsup.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, grepl("Nsuper", colnames(est_mod$chains_mat))] %>% 
  as.data.frame %>% mutate(Iter = rep(1:(niter/2), times = 2), chain = rep(c(1,2), each = (niter/2)) %>% factor()) %>%
  ggplot(aes(x = Iter, y = ., color = chain), alpha = .5) + geom_line() +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = "iteration", y = expression(paste(hat(N)[super]))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()

jpeg("Density_Nsup.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, grepl("Nsuper", colnames(est_mod$chains_mat))] %>% 
  as.data.frame %>% mutate(Iter = rep(1:(niter/2), times = 2), chain = rep(c(1,2), each = (niter/2)) %>% factor()) %>%
  ggplot(aes(x = ., y = after_stat(density), color = chain)) +
  geom_density(position = "identity", linewidth = 1, aes(linetype = chain)) +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = expression(paste(hat(N)[super]))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()

jpeg("Traceplot_N2018.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, 1] %>%  
  as.data.frame %>% mutate(Iter = rep(1:(niter/2), times = 2), chain = rep(c(1,2), each = (niter/2)) %>% factor()) %>%
  ggplot(aes(x = Iter, y = .,color = chain), alpha = .5) + geom_line() +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = "iteration", y = expression(paste(hat(N)[2018]))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()

jpeg("Density_N2018.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, 1] %>% 
  as.data.frame %>% mutate(Iter = rep(1:(niter/2), times = 2), chain = rep(c(1,2), each = (niter/2)) %>% factor()) %>%
  ggplot(aes(x = ., y = after_stat(density), color = chain)) +
  geom_density(position = "identity", linewidth = 1, aes(linetype = chain)) +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = expression(paste(hat(N)[2018]))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()

jpeg("Traceplot_N2019.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, 2] %>% 
  as.data.frame %>% mutate(Iter = rep(1:(niter/2), times = 2), chain = rep(c(1,2), each = (niter/2)) %>% factor()) %>%
  ggplot(aes(x = Iter, y = ., color = chain), alpha = .5) + geom_line() +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = "iteration", y = expression(paste(hat(N)[2019]))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()

jpeg("Density_N2019.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, 2] %>% 
  as.data.frame %>% mutate(Iter = rep(1:(niter/2), times = 2), chain = rep(c(1,2), each = (niter/2)) %>% factor()) %>%
  ggplot(aes(x = ., y = after_stat(density), color = chain)) +
  geom_density(position = "identity", linewidth = 1, aes(linetype = chain)) +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = expression(paste(hat(N)[2019]))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()

jpeg("Traceplot_N2020.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, 3] %>% 
  as.data.frame %>% mutate(Iter = rep(1:(niter/2), times = 2), chain = rep(c(1,2), each = (niter/2)) %>% factor()) %>%
  ggplot(aes(x = Iter, y = ., color = chain), alpha = .5) + geom_line() +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = "iteration", y = expression(paste(hat(N)[2020]))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()

jpeg("Density_N2020.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, 3] %>% 
  as.data.frame %>% mutate(Iter = rep(1:(niter/2), times = 2), chain = rep(c(1,2), each = (niter/2)) %>% factor()) %>%
  ggplot(aes(x = ., y = after_stat(density), color = chain)) +
  geom_density(position = "identity", linewidth = 1, aes(linetype = chain)) +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = expression(paste(hat(N)[2020]))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()

# DELTA
mean(est_mod$chains_mat[, grep("^delta", colnames(est_mod$chains_mat))])
HDInterval::hdi(est_mod$chains_mat[, grep("^delta", colnames(est_mod$chains_mat))])

jpeg("Traceplot_delta.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, grep("^delta", colnames(est_mod$chains_mat))] %>% as.data.frame %>%
  mutate(Iter = rep(1:(niter/2), times = 2), chain = rep(c(1,2), each = (niter/2)) %>% factor()) %>%
  ggplot(aes(x = Iter, y = 1-., color = chain), alpha = .5) + geom_line() +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = "iteration", y = expression(paste(hat(delta)))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()

jpeg("Density_delta.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, grep("^delta", colnames(est_mod$chains_mat))] %>% as.data.frame %>%
  mutate(Iter = rep(1:(niter/2), times = 2), chain = rep(c(1,2), each = (niter/2)) %>% factor()) %>%
  ggplot(aes(x = 1-., y = after_stat(density), color = chain)) +
  geom_density(position = "identity", linewidth = 1, aes(linetype = chain)) +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = expression(paste(hat(delta)))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()

# MU
mean(est_mod$chains_mat[, grep("^mu", colnames(est_mod$chains_mat))])
HDInterval::hdi(est_mod$chains_mat[, grep("^mu", colnames(est_mod$chains_mat))])

jpeg("Traceplot_mu.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, grep("^mu", colnames(est_mod$chains_mat))] %>% as.data.frame %>%
  mutate(Iter = rep(1:(niter/2), times = 2), chain = rep(c(1,2), each = (niter/2)) %>% factor()) %>%
  ggplot(aes(x = Iter, y = .,color = chain), alpha = .5) + geom_line() +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = "iteration", y = expression(paste(hat(mu)))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()

jpeg("Density_mu.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, grep("^mu", colnames(est_mod$chains_mat))] %>% as.data.frame %>%
  mutate(Iter = rep(1:(niter/2), times = 2), chain = rep(c(1,2), each = (niter/2)) %>% factor()) %>%
  ggplot(aes(x = ., y = after_stat(density), color = chain)) +
  geom_density(position = "identity", linewidth = 1, aes(linetype = chain)) +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = expression(paste(hat(mu)))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()


# PHI
colMeans(est_mod$chains_mat[, grep("^phi", colnames(est_mod$chains_mat))])
HDInterval::hdi(est_mod$chains_mat[, grep("^phi", colnames(est_mod$chains_mat))])

jpeg("Traceplot_phi1.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, grep("^phi", colnames(est_mod$chains_mat))][,1] %>% as.data.frame %>%
  mutate(Iter = rep(1:(niter/2), times = 2), chain = rep(c(1,2), each = (niter/2)) %>% factor()) %>%
  ggplot(aes(x = Iter, y = ., color = chain), alpha = .5) + geom_line() +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = "iteration", y = expression(paste(hat(phi)[T]))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()

jpeg("Density_phi1.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, grep("^phi", colnames(est_mod$chains_mat))][,1] %>% as.data.frame %>%
  mutate(Iter = rep(1:(niter/2), times = 2), chain = rep(c(1,2), each = (niter/2)) %>% factor()) %>%
  ggplot(aes(x = ., y = after_stat(density), color = chain)) +
  geom_density(position = "identity", linewidth = 1, aes(linetype = chain)) +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = expression(paste(hat(phi)[T]))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()

jpeg("Traceplot_phi2.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, grep("^phi", colnames(est_mod$chains_mat))][,2] %>% as.data.frame %>%
  mutate(Iter = rep(1:7500, times = 2), chain = rep(c(1,2), each = 7500) %>% factor()) %>%
  ggplot(aes(x = Iter, y = ., color = chain), alpha = .5) + geom_line() +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = "iteration", y = expression(paste(hat(phi)[NT]))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()

jpeg("Density_phi2.jpg", width = 800, height = 600, res = 100)
est_mod$chains_mat[, grep("^phi", colnames(est_mod$chains_mat))][,2] %>% as.data.frame %>%
  mutate(Iter = rep(1:(niter/2), times = 2), chain = rep(c(1,2), each = (niter/2)) %>% factor()) %>%
  ggplot(aes(x = ., y = after_stat(density), color = chain)) +
  geom_density(position = "identity", linewidth = 1, aes(linetype = chain)) +
  scale_color_manual(values = c("black", "grey")) +
  labs(x = expression(paste(hat(phi)[NT]))) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.position = "top")
dev.off()
