##
## CODE FOR JOLLY-SEBER CAPTURE-RECAPTURE MODEL
##

require(nimble)
require(nimbleEcology)
require(nimbleJSextras)

js_code <- nimbleCode({

  #' ---------------------------------------------------------------------------
  #' Entry and Transition matrix
  #' ---------------------------------------------------------------------------

  # xi
  # sig_xi_R ~ dexp(1)
  # sig_xi_P ~ dexp(1)
  # sig_xi_T ~ dexp(1)
  # beta_xi_R[1] <- 0
  # beta_xi_P[1] <- 0
  # beta_xi_T[1] <- 0
  # for(j in 2:n_surveys){
  #   beta_xi_R[j] ~ dnorm(0,sd=sig_xi_R)
  #   beta_xi_P[j] ~ dnorm(0,sd=sig_xi_P)
  #   beta_xi_T[j] ~ dnorm(0,sd=sig_xi_T)
  # }
  # xi_R[1:K] <- exp( Mxi[1:K,1:n_surveys] %*% beta_xi_R[1:n_surveys] )
  # xi_P[1:K] <- exp( Mxi[1:K,1:n_surveys] %*% beta_xi_P[1:n_surveys] )
  # xi_T[1:K] <- exp( Mxi[1:K,1:n_surveys] %*% beta_xi_T[1:n_surveys] )
  # for(t in 1:(K-1)){
  #   xi_tilde_R[t] <- xi_R[t]/sum(xi_R[t:K])
  #   xi_tilde_P[t] <- xi_P[t]/sum(xi_P[t:K])
  #   xi_tilde_T[t] <- xi_T[t]/sum(xi_T[t:K])
  # }
  # xi_tilde_R[K] <- 1
  # xi_tilde_P[K] <- 1
  # xi_tilde_T[K] <- 1

  for(j in 1:2){
    rho_R[j] ~ dunif(0,1)
    rho_P[j] ~ dunif(0,1)
    rho_T[j] ~ dunif(0,1)
  }
  xi_tilde_R[1] <- rho_R[1]/(1 - ( (1-rho_R[1])*(1-rho_R[2])^(K-1) ))
  xi_tilde_P[1] <- rho_P[1]/(1 - ( (1-rho_P[1])*(1-rho_P[2])^(K-1) ))
  xi_tilde_T[1] <- rho_T[1]/(1 - ( (1-rho_T[1])*(1-rho_T[2])^(K-1) ))
  for(t in 2:(K-1)){
    xi_tilde_R[t] <- rho_R[2]/(1-(1-rho_R[2])^(K-t+1) )
    xi_tilde_P[t] <- rho_P[2]/(1-(1-rho_P[2])^(K-t+1) )
    xi_tilde_T[t] <- rho_T[2]/(1-(1-rho_T[2])^(K-t+1) )
  }
  xi_tilde_R[K] <- 1
  xi_tilde_P[K] <- 1
  xi_tilde_T[K] <- 1

  # phi
  phi_RP ~ dunif(0,1)
  delta_phi ~ dunif(0,1)
  phi_T <- phi_RP * delta_phi

  # alpha
  mlogit_alpha[1] <- 0
  mlogit_alpha[2] ~ dnorm(0,sd=1)
  mlogit_alpha[3] ~ dnorm(0,sd=1)
  for(j in 1:3){
    alpha[j] <- exp(mlogit_alpha[j])/sum(exp(mlogit_alpha[1:3]))
  }

  # pi
  # Resident
  pi[1] <- alpha[1]*(1-xi_tilde_R[1])
  pi[2] <- alpha[1]*xi_tilde_R[1]
  # Part-time
  pi[3] <- alpha[2]*(1-xi_tilde_P[1])
  pi[4] <- alpha[2]*xi_tilde_P[1]
  # Transient
  pi[5] <- alpha[3]*(1-xi_tilde_T[1])
  pi[6] <- alpha[3]*xi_tilde_T[1]
  # Exit
  pi[7] <- 0

  for(t in 1:(K-1)){
    # Resident
    Gamma[1,1,t] <- 1-xi_tilde_R[t+1]
    Gamma[1,2,t] <- xi_tilde_R[t+1]
    for(j in 3:7){ Gamma[1,j,t] <- 0}
    Gamma[2,1,t] <- 0
    Gamma[2,2,t] <- phi_RP
    for(j in 3:6){ Gamma[2,j,t] <- 0}
    Gamma[2,7,t] <- 1-phi_RP
    # Part-time resident
    for(j in 1:2){ Gamma[3,j,t] <- 0 }
    Gamma[3,3,t] <- 1-xi_tilde_P[t+1]
    Gamma[3,4,t] <- xi_tilde_P[t+1]
    for(j in 5:7){ Gamma[3,j,t] <- 0 }
    for(j in 1:3){ Gamma[4,j,t] <- 0 }
    Gamma[4,4,t] <- phi_RP
    for( j in 5:6){ Gamma[4,j,t] <- 0 }
    Gamma[4,7,t] <- 1-phi_RP
    # Transient
    for(j in 1:4){ Gamma[5,j,t] <- 0 }
    Gamma[5,5,t] <- 1-xi_tilde_T[t+1]
    Gamma[5,6,t] <- xi_tilde_T[t+1]
    Gamma[5,7,t] <- 0
    for(j in 1:5){ Gamma[6,j,t] <- 0 }
    Gamma[6,6,t] <- phi_T
    Gamma[6,7,t] <- 1-phi_T
    # Exit
    for(j in 1:6){ Gamma[7,j,t] <- 0 }
    Gamma[7,7,t] <- 1
  }

  #' ---------------------------------------------------------------------------
  #' Detection matrix
  #' ---------------------------------------------------------------------------

  # P
  mu_p ~ dnorm(0,sd=1)
  sig_p ~ dexp(1)
  delta_p ~ dunif(0,1)
  for(j in 2:n_surveys){
    eps[j] ~ dnorm(0, sd=sig_p)
    logit(p_RT[j]) <- mu_p + eps[j]
    p_P[j] <- p_RT[j] * delta_p
  }
  eps[1] <- -sum(eps[2:n_surveys])
  logit(p_RT[1]) <- mu_p + eps[1]
  p_P[1] <- p_RT[1] * delta_p

  p_RT_exp[1:K] <- Mp[1:K,1:n_surveys] %*% p_RT[1:n_surveys]
  p_P_exp[1:K] <- Mp[1:K,1:n_surveys] %*% p_P[1:n_surveys]


  for(t in 1:K){
    # Resident
    Pmats[t,1] <- 0
    Pmats[t,2] <- p_RT_exp[t]
    # Part-time resident
    Pmats[t,3] <- 0
    Pmats[t,4] <- p_P_exp[t]
    # Transient
    Pmats[t,5] <- 0
    Pmats[t,6] <- p_RT_exp[t]
    # Exit
    Pmats[t,7] <- 0
  }

  #' ---------------------------------------------------------------------------
  #' Unconditional detection probability
  #' ---------------------------------------------------------------------------
  pstar <- pstar_binom(
    init = pi[1:7],
    prob = Pmats[1:K,1:7],
    size = Rt[1:K],
    probTrans = Gamma[1:7, 1:7, 1:(K-1)],
    len = K
  )

  #' ---------------------------------------------------------------------------
  #' HMM likelihood for observed individuals
  #' ---------------------------------------------------------------------------
  for(i in 1:nobs){
    x[i, 1:K] ~ dJS_binom(
      init = pi[1:7],
      prob = Pmats[1:K,1:7],
      size = Rt[1:K],
      probTrans = Gamma[1:7, 1:7, 1:(K-1)],
      pstar = pstar, len=K
    )
  }

  #' ---------------------------------------------------------------------------
  #' Model for number of detected individuals
  #' ---------------------------------------------------------------------------
  lambda ~ dgamma(1.0e-6, 1.0e-6)
  n ~ dpois(lambda*pstar)

  #' ---------------------------------------------------------------------------
  #' Code below is for the posterior predicted abundance
  #' ---------------------------------------------------------------------------
  nu ~ dpois(lambda*(1-pstar))
  Nsuper <- n+nu

  nu_t[1:7, 1:K] <- sample_undet_binom(nu, pi[1:7], Pmats[1:K,1:7], Rt[1:K],
                                       Gamma[1:7, 1:7, 1:(K-1)])

  for(i in 1:nobs){
    det_state[i,1:K] <- sample_det_binom(x[i,1:K], pi[1:7], Pmats[1:K,1:7],
                                         Rt[1:K], Gamma[1:7, 1:7, 1:(K-1)])
  }
  avail_R[1:nobs,1:K] <- det_state[1:nobs,1:K]==2
  avail_P[1:nobs,1:K] <- det_state[1:nobs,1:K]==4
  avail_T[1:nobs,1:K] <- det_state[1:nobs,1:K]==6

  for(t in 1:K){
    N_R[t] <- sum(avail_R[1:nobs, t]) + nu_t[2,t]
    N_P[t] <- sum(avail_P[1:nobs, t]) + nu_t[4,t]
    N_T[t] <- sum(avail_T[1:nobs, t]) + nu_t[6,t]
  }

})
