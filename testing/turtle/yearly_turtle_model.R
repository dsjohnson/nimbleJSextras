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

  # # rho
  # for(t in 1:K){
  #   for(j in 1:3){
  #     rho[t,j] ~ dbeta(1/K, 2 - t/K)
  #   }
  # }

  for(t in 1:K){
    xi_raw[t] ~ dgamma(1/K, 1)
  }
  xi_norm <- sum(xi_raw[1:K])
  for(t in 1:(K-1)){
    xi_tilde[t] <- xi_raw[t]/sum(xi_raw[t:K])
    xi[t] <- xi_raw[t]/xi_norm
  }
  xi_tilde[K] <- 1
  xi[K] <- xi_raw[K]/xi_norm

  # phi
  sig_eps ~ dexp(1)
  mu_phi ~ dnorm(0,sd=1.5)
  for(t in 1:(K-1)){
    eps[t] ~ dnorm(0,sd=sig_eps)
    logit(phi[t]) <- mu_phi + eps[t]
  }


  # Foraging dwell time
  theta ~ dexp(10)
  zeta[1:nd] <- load_dt_pois(nd, theta, shift=1)

  # pi
  pi[1:n_states] <- load_dt_pi(xi_tilde[1],nd)

  # Gamma
  for(t in 1:(K-1)){
    Gamma[1:n_states, 1:n_states, t] <- load_dt_Gamma(xi_tilde[t+1], phi[t], zeta[1:nd])
  }

  #' ---------------------------------------------------------------------------
  #' Detection matrix
  #' ---------------------------------------------------------------------------

  # P
  # phi
  sig_tau ~ dexp(1)
  mu_p ~ dnorm(0,sd=1.5)
  for(t in 1:K){
    tau[t] ~ dnorm(0,sd=sig_tau)
    p[t] <- (1-(t==skip_surv))*expit(mu_p + tau[t])
  }

  # t = 27 is the 2020 season where no captures occurred due to COVID
  for(t in 1:K){
    Pmats[1:n_states,1:2,t] <- load_dt_Pmats(p[t],nd)
  }


  #' ---------------------------------------------------------------------------
  #' Unconditional detection probability
  #' ---------------------------------------------------------------------------
  pstar <- pstar_ms(
    init = pi[1:n_states],
    probObs = Pmats[1:n_states, 1:2, 1:K],
    probTrans = Gamma[1:n_states, 1:n_states, 1:(K-1)],
    len=K
  )

  #' ---------------------------------------------------------------------------
  #' HMM likelihood for observed individuals
  #' ---------------------------------------------------------------------------
  for(i in 1:nobs_u){
    ux[i, 1:K] ~ dJS_ms(
      init = pi[1:n_states],
      probObs = Pmats[1:n_states, 1:2, 1:K],
      probTrans = Gamma[1:n_states, 1:n_states, 1:(K-1)],
      pstar=pstar, weight=w[i], len = K
    )
    # loglik_i[i] <- dJS_ms(
    #   ux[i, 1:K],
    #   init = pi[1:9],
    #   probObs = Pmats[1:9, 1:2, 1:K],
    #   probTrans = Gamma[1:9, 1:9, 1:(K-1)],
    #   pstar=pstar, weight=w[i], len = K, log=1
    # )
  }
  # loglik <- sum(loglik_i[1:nobs_u])

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
  # nu_t[1:n_states, 1:K] <- sample_undet_ms(n=nu, init=pi[1:n_states], probObs=Pmats[1:n_states,1:2,1:K], probTrans=Gamma[1:n_states, 1:n_states, 1:(K-1)])
  # for(i in 1:nobs){
  #   det_state[i,1:K] <- sample_det_ms(x=x[i,1:K], init=pi[1:n_states],
  #                                     probObs=Pmats[1:n_states, 1:2, 1:K],
  #                                     probTrans=Gamma[1:n_states, 1:n_states, 1:(K-1)])
  # }
  # nest_obs[1:nobs,1:K] <- det_state[1:nobs,1:K]==2
  # avail_obs[1:nobs,1:K] <- 1-((det_state[1:nobs,1:K]==1) + (det_state[1:nobs,1:K]==nd))

  # for(t in 1:K){
  #   N_nest[t] <- sum(nest_obs[1:nobs, t]) + nu_t[2,t]
  #   N[t] <- nu - (nu_t[1,t]+nu_t[nd,t]) + sum(avail_obs[1:nobs, t])
  # }

})
