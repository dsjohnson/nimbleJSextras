##
## CODE FOR JOLLY-SEBER CAPTURE-RECAPTURE MODEL
##

require(nimble)
require(nimbleEcology)
require(nimbleJSextras)

js_code <- nimbleCode({
  
  #' ---------------------------------------------------------------------------
  #' Entry, survival, and reproductive state
  #' ---------------------------------------------------------------------------
  
  # Entry probability
  # Matern RBF coefficients
  sigma_beta_xi ~ dexp(1)
  for(i in 1:df_xi){beta_xi[i] ~ dnorm(0, sd=1)}
  # xi
  for(t in 1:K){
    log(xi_raw[t]) <- sigma_beta_xi*inprod(B[t,1:df_xi], beta_xi[1:df_xi])
  }
  xi[1:K] <- xi_raw[1:K]/sum(xi_raw[1:K])
  # xi_tilde
  xi_tilde[1] <- xi[1]
  for(t in 2:(K-1)){ xi_tilde[t] <- xi[t]/sum(xi[t:K]) }
  xi_tilde[K] <- 1
  
  # survival (phi)
  logit_phi ~ dt(mu=0, sigma=1.52, df=6)
  phi <- expit(logit_phi)
  
  # Foraging dwell time (i.e., reproductive return time)
  for(j in 1:m){
    logit_tau[j] ~ dt(mu=0, sigma=1.52, df=6)
    tau[j] <- expit(logit_tau[j])
  }
  ft_pdf[1] <- tau[1]
  for(j in 2:m){
    ft_pdf[j] <- tau[j] * prod(1-tau[1:(j-1)])
  }
  for(j in 1:(10-m)){
    ft_pdf[j+m] <- prod(1-tau[1:m]) * (1-tau[m])^(j-1) * tau[m]
  }
  
  #' ---------------------------------------------------------------------------
  #' Capture probability
  #' ---------------------------------------------------------------------------
  for(t in 1:2){
    # Reduction in effort for uncaptured animals set to 1 for first 2 occasions
    logit_p[t] ~ dt(mu=0, sigma=1.52, df=6)
    p[t] <- surv[t]*expit(logit_p[t])
    p0[t] <- p[t]
  }
  for(t in 3:K){
    logit_p[t] ~ dt(mu=0, sigma=1.52, df=6)
    p[t] <- surv[t]*expit(logit_p[t])
    logit_eta[t-2] ~ dt(mu=0, sigma=1.52, df=6)
    p0[t] <- surv[t]*expit(logit_eta[t-2])*p[t]
  }
  
  #' ---------------------------------------------------------------------------
  #' HMM matrices (see helper/*.R for nimble functions to load matrices)
  #' ---------------------------------------------------------------------------
  
  for(t in 1:K){
    Pmats[1:n_states,1:3,t] <- make_honu_P(m)
  }
  
  # pi
  pi[1:n_states] <- make_honu_pi(xi_tilde[1], p[1], m)
  
  # Gamma
  for(t in 1:(K-1)){
    Gamma[1:n_states, 1:n_states, t] <- make_honu_Gamma(xi_tilde[t+1], phi, p0[t+1], p[t+1], tau[1:m])
  }
  
  #' ---------------------------------------------------------------------------
  #' Unconditional detection probability
  #' ---------------------------------------------------------------------------
  pstar <- pstar_ms(
    init = pi[1:n_states],
    probObs = Pmats[1:n_states, 1:3, 1:K],
    probTrans = Gamma[1:n_states, 1:n_states, 1:(K-1)],
    len=K
  )
  
  #' ---------------------------------------------------------------------------
  #' HMM likelihood for observed individuals
  #' ---------------------------------------------------------------------------
  for(i in 1:nobs_u){
    ux[i, 1:K] ~ dJS_ms(
      init = pi[1:n_states],
      probObs = Pmats[1:n_states, 1:3, 1:K],
      probTrans = Gamma[1:n_states, 1:n_states, 1:(K-1)],
      pstar=pstar, weight=w[i], len = K
    )
  }
  
  #' ---------------------------------------------------------------------------
  #' Model for number of total captured individuals
  #' ---------------------------------------------------------------------------
  lambda ~ dgamma(1.0e-6, 1.0e-6)
  n ~ dpois(lambda*pstar)
  
  #' ---------------------------------------------------------------------------
  #' Posterior predicted superpopulationn abundance
  #' ---------------------------------------------------------------------------
  Nu ~ dpois(lambda*(1-pstar))
  N <- n+Nu
  
})
