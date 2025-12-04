##
## CODE FOR JOLLY-SEBER CAPTURE-RECAPTURE MODEL
##

require(nimble)
require(nimbleEcology)
require(nimbleJSextras)

js_code <- nimbleCode({

  p_dot ~ dunif(0,1)
  rho_dot ~ dunif(0,1)
  for(t in 1:K){rho[t] <- rho_dot}
  phi_dot ~ dunif(0,1)
  psi23_dot ~ dunif(0,1)
  psi32_dot ~ dunif(0,1)

  xi[1] <- rho[1]
  for(t in 2:(K-1)){ xi[t] <- rho[t]/(1-prod(1-rho[t:K])) }
  xi[K] <- 1

  #' ---------------------------------------------------------------------------
  #' Unconditional detection probability
  #' ---------------------------------------------------------------------------
  pstar <- pstar_ms(
    init = pi[1:4],
    probObs = Pmats[1:4, 1:2, 1:K],
    probTrans = Gamma[1:4, 1:4, 1:(K-1)],
    len=K
  )

  #' ---------------------------------------------------------------------------
  #' HMM likelihood for observed individuals
  #' ---------------------------------------------------------------------------
  for(i in 1:nobs){
    x[i, 1:K] ~ dJS_ms(
      init = pi[1:4],
      probObs = Pmats[1:4, 1:2, 1:K],
      probTrans = Gamma[1:4, 1:4, 1:(K-1)],
      pstar=pstar, len = K
    )
  }
  n ~ dpois(lambda*pstar)

  #' ---------------------------------------------------------------------------
  #' Initial entry probability
  #' ---------------------------------------------------------------------------
  pi[1] <- 1-xi[1]
  pi[2] <- xi[1]
  pi[3] <- 0
  pi[4] <- 0

  #' ---------------------------------------------------------------------------
  #' Detection probabilities and matrix
  #' ---------------------------------------------------------------------------
  for(t in 1:K){
    Pmats[1,1,t] <- 1
    Pmats[2,1,t] <- 1-p_dot
    Pmats[3,1,t] <- 1
    Pmats[4,1,t] <- 1

    Pmats[1,2,t] <- 0
    Pmats[2,2,t] <- p_dot
    Pmats[3,2,t] <- 0
    Pmats[4,2,t] <- 0
  }

  #' ---------------------------------------------------------------------------
  #' Transition probability matrix
  #' ---------------------------------------------------------------------------
  for(t in 1:(K-1)){
    Gamma[1,1,t] <- 1-xi[t+1]
    Gamma[1,2,t] <- xi[t+1]
    Gamma[1,3,t] <- 0
    Gamma[1,4,t] <- 0

    Gamma[2,1,t] <- 0
    Gamma[2,2,t] <- phi_dot * (1-psi23_dot)
    Gamma[2,3,t] <- phi_dot * psi23_dot
    Gamma[2,4,t] <- 1 - phi_dot

    Gamma[3,1,t] <- 0
    Gamma[3,2,t] <- phi_dot * psi32_dot
    Gamma[3,3,t] <- phi_dot * (1-psi32_dot)
    Gamma[3,4,t] <- 1 - phi_dot

    Gamma[4,1,t] <- 0
    Gamma[4,2,t] <- 0
    Gamma[4,3,t] <- 0
    Gamma[4,4,t] <- 1
  }

  #' lambda prior
  lambda ~ dgamma(1.0e-6, 1.0e-6)

  #' ---------------------------------------------------------------------------
  #' Code below is for the posterior predicted abundance
  #' ---------------------------------------------------------------------------

  nu ~ dpois(lambda*(1-pstar))
  Nsuper <- n+nu

  nu_t[1:4, 1:K] <- sample_undet_ms(n=nu, init=pi[1:4], probObs=Pmats[1:4,1:2,1:K],
                                    probTrans=Gamma[1:4, 1:4, 1:(K-1)])

  for(i in 1:nobs){
    det_state[i,1:K] <- sample_det_ms(x=x[i,1:K], init=pi[1:4],
                                      probObs=Pmats[1:4, 1:2, 1:K],
                                      probTrans=Gamma[1:4, 1:4, 1:(K-1)])
  }
  nest_obs[1:nobs,1:K] <- det_state[1:nobs,1:K]==2
  internest_obs[1:nobs,1:K] <- det_state[1:nobs,1:K]==3

  for(t in 1:K){
    n_nest_obs[t] <- sum(nest_obs[1:nobs, t])
    n_internest_obs[t] <- sum(internest_obs[1:nobs, t])
    N_nest[t] <- n_nest_obs[t] + nu_t[2,t]
    N_internest[t] <- n_internest_obs[t] + nu_t[3,t]
    N[t] <- N_nest[t] + N_internest[t]
  }
})
