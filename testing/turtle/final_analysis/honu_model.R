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

  # xi and xi_tilde
  beta_rho_int ~ dnorm(0,sd=1.5)
  sig_rho ~ dexp(10)
  for(j in 1:k_rho){ beta_rho[j]~dnorm(0,sd=sig_rho) }
  for(t in 1:K){
    # eps_rho[t] ~ dnorm(0,sd=sig_rho)
    logit(rho[t]) <- beta_rho_int + inprod(B_rho[t,1:k_rho], beta_rho[1:k_rho])
    # logit(rho[t]) <- beta_rho_int + eps_rho[t]
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
  # sig_phi ~ dexp(1)
  # beta_phi_int ~ dnorm(0,sd=1.5)
  # for(j in 1:k_phi){ beta_phi[j]~dnorm(0,sd=sig_phi) }
  # # beta_phi ~ dnorm(0,sd=sig_phi)
  # for(t in 1:(K-1)){
  #   logit(phi[t]) <- beta_phi_int + inprod(B_phi[t,1:k_phi], beta_phi[1:k_phi])
  #   }
  phi ~ dbeta(1,1)

  # Foraging dwell time (mu)
  sig_mu ~ dexp(10)
  beta_mu_int ~ dnorm(0,sd=10)
  for(j in 1:k_mu) { beta_mu[j]~dnorm(0,sd=sig_mu) }
  for(t in 1:(K-1)){
    log(mu[t]) <- beta_mu_int + inprod(B_mu[t,1:k_mu], beta_mu[1:k_mu])
    # eps_mu[t] ~ dnorm(0,sd=sig_mu)
    # log(mu[t]) <-  beta_mu_int + eps_mu[t]
    }

  #' ---------------------------------------------------------------------------
  #' Detection
  #' ---------------------------------------------------------------------------

  sig_p ~ dexp(10)
  beta_p_int ~ dnorm(0,sd=1.5)
  for(t in 1:K){
    eps_p[t] ~ dnorm(0,sd=sig_p)
    # alpha[t] ~ dbeta(5,1)
    alpha[t] ~ dbeta(1,1)
    p[t,2] <- surv[t]*expit(beta_p_int + eps_p[t])
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
    zeta[t,1:nd] <- make_dt_pois(nd, mu[t], shift=1)
    Gamma[1:n_states, 1:n_states, t] <- make_honu_Gamma(xi_tilde[t+1], phi, p[t+1,1:2], alpha=0, zeta[t,1:nd])
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
