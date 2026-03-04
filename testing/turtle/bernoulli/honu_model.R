##
## CODE FOR JOLLY-SEBER CAPTURE-RECAPTURE MODEL
##

require(nimble)
require(nimbleEcology)
require(nimbleJSextras)

js_code <- nimbleCode({

  #' ---------------------------------------------------------------------------
  #' Entry and survival
  #' ---------------------------------------------------------------------------

  # alpha
  # alpha ~ dbeta(6,48)

  # xi and xi_tilde
  beta_rho_int ~ dnorm(0,sd=1.5)
  sig_rho ~ dexp(1)
  for(j in 1:m){ beta_rho[j]~dnorm(0,sd=sig_rho) }
  # for(t in 1:K){ eps_rho[t] ~ dnorm(0,sd=sig_rho) }
  for(t in 1:K){
    logit(rho[t]) <- beta_rho_int + inprod(B[t,1:m], beta_rho[1:m])
    }
  xi_raw[1] <- rho[1]
  for(t in 2:K){
    xi_raw[t] <- prod(1-rho[1:(t-1)]) * rho[t]
  }
  xi_norm <- sum(xi_raw[1:K])
  for(t in 1:(K-1)){
    xi_tilde[t] <- xi_raw[t]/sum(xi_raw[t:K])
    xi[t] <- xi_raw[t]/xi_norm
  }
  xi_tilde[K] <- 1
  xi[K] <- xi_raw[K]/xi_norm

  # phi
  sig_phi ~ dexp(1)
  beta_phi_int ~ dnorm(0,sd=1.5)
  # for(j in 1:m){ beta_phi[j]~dnorm(0,sd=sig_phi) }
  beta_phi ~ dnorm(0,sd=sig_phi)
  for(t in 1:(K-1)){
    logit(phi[t]) <- beta_phi_int + beta_phi*t #inprod(B[t,1:m], beta_phi[1:m])
    }

  # Foraging dwell time (theta)
  sig_theta ~ dexp(1)
  beta_theta_int ~ dnorm(0,sd=1.5)
  # for(j in 1:m) { beta_theta[j]~dnorm(0,sd=sig_theta) }
  beta_theta ~ dnorm(0,sd=sig_theta)
  for(t in 1:(K-1)){
    log(theta[t]) <- beta_theta_int + beta_theta*t #inprod(B[t,1:m], beta_theta[1:m])
    }

  #' ---------------------------------------------------------------------------
  #' Detection
  #' ---------------------------------------------------------------------------

  sig_p ~ dexp(1)
  beta_p[1] ~ dnorm(0,sd=1.5)
  # for(j in 1:mp){ beta_p[j] ~ dnorm(0,sd=sig_p) }
  beta_p[2] ~ dnorm(0,sd=sig_p)
  for(t in 1:K){
    # tau[t] ~ dnorm(0,sd=sig_tau)
    alpha[t] ~ dbeta(5,1)
    p[t,2] <- surv[t]*expit(inprod(X[t,1:2],beta_p[1:2]))
    p[t,1] <- surv[t]*alpha[t]*p[t,2]
  }

  #' ---------------------------------------------------------------------------
  #' HMM matrices
  #' ---------------------------------------------------------------------------

  for(t in 1:K){
    Pmats[1:n_states,1:3,t] <- make_honu_P(nd)
  }

  # pi
  pi[1:n_states] <- make_honu_pi(xi_tilde[1], p[1,1:2], alpha=0, nd)

  # Gamma
  for(t in 1:(K-1)){
    zeta[t,1:nd] <- make_dt_pois(nd, theta[t], shift=1)
    Gamma[1:n_states, 1:n_states, t] <- make_honu_Gamma(xi_tilde[t+1], phi[t], p[t+1,1:2], alpha=0, zeta[t,1:nd])
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
