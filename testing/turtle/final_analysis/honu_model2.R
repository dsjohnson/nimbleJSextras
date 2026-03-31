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
  # sig_rho ~ dexp(1)
  for(t in 1:K){
    eps_rho[t] ~ dt(mu=0, sigma=1.52, df=6)
    logit(rho[t]) <- eps_rho[t]
    # rho[t] ~ dbeta(1/K,2-t/K)
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
  eps_phi ~ dt(mu=0, sigma=1.52, df=6)
  logit(phi) <- eps_phi

  # Foraging dwell time (tau)
  for(j in 1:m){
    eps_tau[j] ~ dt(mu=0, sigma=1.52, df=6)
    logit(tau[j]) <- eps_tau[j]
  }
  ft_pdf[1] <- tau[1]
  for(j in 2:m){
    ft_pdf[j] <- tau[j] * prod(1-tau[1:(j-1)])
  }
  for(j in 1:(10-m)){
    ft_pdf[j+m] <- prod(1-tau[1:m]) * (1-tau[m])^(j-1) * tau[m]
  }

  #' ---------------------------------------------------------------------------
  #' Detection
  #' ---------------------------------------------------------------------------
  p[1] <- 0
  eps_a[1] ~ dt(mu=0, sigma=1.52, df=6)
  logit(alpha[1]) <- eps_a[1]
  for(t in 2:K){
    eps_p[t-1] ~ dt(mu=0, sigma=1.52, df=6)
    p[t] <- surv[t]*expit(eps_p[t-1])
    eps_a[t] ~ dt(mu=0, sigma=1.52, df=6)
    logit(alpha[t]) <- eps_a[t]
  }

  #' ---------------------------------------------------------------------------
  #' HMM matrices
  #' ---------------------------------------------------------------------------

  for(t in 1:K){
    Pmats[1:n_states,1:3,t] <- make_honu_P(m)
  }

  # pi
  pi[1:n_states] <- make_honu_pi(xi_tilde[1], alpha[1], m)

  # Gamma
  for(t in 1:(K-1)){
    Gamma[1:n_states, 1:n_states, t] <- make_honu_Gamma(xi_tilde[t+1], phi, alpha[t+1], p[t+1], tau[1:m])
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
