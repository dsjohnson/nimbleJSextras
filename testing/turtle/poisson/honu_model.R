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


  sig_beta ~ dexp(1)
  beta_int ~ dnorm(0,sd=1.5)
  for(j in 1:m){beta[j]~dnorm(0,sd=sig_beta)}
  for(t in 1:K){ logit(rho[t]) <- mu_rho + inprod(B[t,1:m], beta[1:m])  }
  xi_raw[1] <- rho[1]
  for(t in 2:K){
    xi_raw[t] <- prod(1-rho[1:(t-1)]) * rho[t]
  }
  # for(t in 1:K){
  #   xi_raw[t] ~ dgamma(1/K, 1)
  # }
  xi_norm <- sum(xi_raw[1:K])
  for(t in 1:(K-1)){
    xi_tilde[t] <- xi_raw[t]/sum(xi_raw[t:K])
    xi[t] <- xi_raw[t]/xi_norm
  }
  xi_tilde[K] <- 1
  xi[K] <- xi_raw[K]/xi_norm

  # phi
  sig_delta ~ dexp(1)
  mu_phi ~ dnorm(0,sd=1.5)
  for(j in 1:m){delta[j]~dnorm(0,sd=sig_delta)}
  for(t in 1:(K-1)){ logit(phi[t]) <- mu_phi + inprod(B[t,1:m], delta[1:m])  }

  # Foraging dwell time (theta)
  sig_gamma ~ dexp(1)
  mu_theta ~ dnorm(0,sd=1.5)
  for(j in 1:m){gamma[j]~dnorm(0,sd=sig_beta)}
  for(t in 1:(K-1)){ log(theta[t]) <- mu_theta + inprod(B[t,1:m], gamma[1:m])  }

  # theta ~ dexp(10)
  # zeta[1:nd] <- make_dt_pois(nd, theta, shift=1)

  # pi
  pi[1:n_states] <- make_honu_pi(xi_tilde[1],nd)

  # Gamma
  for(t in 1:(K-1)){
    zeta[t,1:nd] <- make_dt_pois(nd, theta[t], shift=1)
    Gamma[1:n_states, 1:n_states, t] <- make_honu_Gamma(xi_tilde[t+1], phi[t], zeta[t,1:nd])
  }

  #' ---------------------------------------------------------------------------
  #' Detection matrix
  #' ---------------------------------------------------------------------------

  sig_tau1 ~ dexp(1)
  sig_tau2 ~ dexp(1)
  mu_r1 ~ dnorm(0,sd=1.5)
  mu_r2 ~ dnorm(0,sd=1.5)
  for(t in 1:K){
    tau1[t] ~ dnorm(0,sd=sig_tau1)
    tau2[t] ~ dnorm(0,sd=sig_tau2)
    r[t,1] <- surv[t,1]*exp(mu_r1 + tau1[t])
    r[t,2] <- surv[t,2]*exp(mu_r2 + tau1[t] + tau2[t])
  }

  # t = 27 is the 2020 season where no captures occurred due to COVID
  Rmat[1:K,1:n_states] <- make_honu_Rmat(r[1:K,1:2],nd)

  #' ---------------------------------------------------------------------------
  #' Unconditional detection probability
  #' ---------------------------------------------------------------------------
  pstar <- pstar_pois(
    init = pi[1:n_states],
    rate = Rmat[1:K, 1:n_states],
    probTrans = Gamma[1:n_states, 1:n_states, 1:(K-1)],
    len=K
  )

  #' ---------------------------------------------------------------------------
  #' HMM likelihood for observed individuals
  #' ---------------------------------------------------------------------------
  for(i in 1:nobs_u){
    ux[i, 1:K] ~ dJS_pois(
      init = pi[1:n_states],
      rate = Rmat[1:K, 1:n_states],
      probTrans = Gamma[1:n_states, 1:n_states, 1:(K-1)],
      pstar=pstar, weight=w[i], len = K
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
  # nu_t[1:n_states, 1:K] <- sample_undet_ms(n=nu, init=pi[1:n_states], probObs=Pmats[1:n_states,1:2,1:K], probTrans=Gamma[1:n_states, 1:n_states, 1:(K-1)])
  # for(i in 1:nobs){
  #   det_state[i,1:K] <- sample_det_ms(x=x[i,1:K], init=pi[1:n_states],
  #                                     probObs=Pmats[1:n_states, 1:2, 1:K],
  #                                     probTrans=Gamma[1:n_states, 1:n_states, 1:(K-1)])
  # }
  # nest_obs[1:nobs,1:K] <- det_state[1:nobs,1:K]==2
  # avail_obs[1:nobs,1:K] <- 1-((det_state[1:nobs,1:K]==1) + (det_state[1:nobs,1:K]==n_states))
  #
  # for(t in 1:K){
  #   N_nest[t] <- sum(nest_obs[1:nobs, t]) + nu_t[2,t]
  #   N_avail[t] <- nu - (nu_t[1,t]+nu_t[nd,t]) + sum(avail_obs[1:nobs, t])
  # }

})
