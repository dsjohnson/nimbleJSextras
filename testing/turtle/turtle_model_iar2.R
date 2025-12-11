##
## CODE FOR JOLLY-SEBER CAPTURE-RECAPTURE MODEL
##

require(nimble)
require(nimbleEcology)
require(nimbleJSextras)

js_code <- nimbleCode({

  p_dot ~ dunif(0,1)

  sigma_phi ~ dexp(1)
  tau_phi <- 1/sigma_phi
  # logit_phi[1:(K-1)] ~ dcar_normal(adj[1:L], weights[1:L], num[1:(K-1)], tau_phi, 2)
  for(t in 1:(K-1)){
    logit_phi[t] ~ dnorm(mean=0, sd=sigma_phi)
    phi[t] <- expit(logit_phi[t])
    }


  sigma_rho ~ dexp(1)
  tau_rho <- 1/sigma_rho
  # logit_rho[1:K] ~ dcar_normal(adj_rho[1:L_rho], weights_rho[1:L_rho], num_rho[1:K], tau_rho, 2)
  for(t in 1:K){
    logit_rho[t] ~ dnorm(mean=0, sd=sigma_rho)
    rho[t] <- expit(logit_rho[t])
    }


  sigma_psi23 ~ dexp(1)
  tau_psi23 <- 1/sigma_psi23
  # logit_psi23_dot[1:(K-1)] ~ dcar_normal(adj[1:L], weights[1:L], num[1:(K-1)], tau_psi23, 2)
  for(t in 1:(K-1)){
    logit_psi23[t] ~ dnorm(mean=0, sd=sigma_psi23)
    psi23[t] <- expit(logit_psi23[t])
    }


  sigma_psi32 ~ dexp(1)
  tau_psi32 <- 1/sigma_psi32
  # logit_psi32_dot[1:(K-1)] ~ dcar_normal(adj[1:L], weights[1:L], num[1:(K-1)], tau=tau_psi32, 2)
  for(t in 1:(K-1)){
    logit_psi32[t] ~ dnorm(mean=0, sd=sigma_psi32)
    psi32[t] <- expit(logit_psi32[t])
    }

  # xi[1] <- rho[1]
  for(t in 1:K{ xi[t] <- rho[t]/(1-prod(1-rho[t:K])) }
  # xi[K] <- 1

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
    Gamma[2,2,t] <- phi[t] * (1-psi23[t])
    Gamma[2,3,t] <- phi[t] * psi23[t]
    Gamma[2,4,t] <- 1 - phi[t]

    Gamma[3,1,t] <- 0
    Gamma[3,2,t] <- phi[t] * psi[t]
    Gamma[3,3,t] <- phi[t] * (1-psi32[t])
    Gamma[3,4,t] <- 1 - phi[t]

    Gamma[4,1,t] <- 0
    Gamma[4,2,t] <- 0
    Gamma[4,3,t] <- 0
    Gamma[4,4,t] <- 1
  }


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

  #' lambda prior
  lambda ~ dgamma(1.0e-6, 1.0e-6)
  n ~ dpois(lambda*pstar)

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
