##
## CODE FOR JOLLY-SEBER CAPTURE-RECAPTURE MODEL (BINOMIAL / SIZE = 1)
##

require(nimble)
require(nimbleEcology)
require(nimbleJSextras)

js_code <- nimbleCode({
  
  for(t in 1:K){
    logit_p[t] ~ dnorm(mu_p, sd=sig_p)
    p[t] <- expit(logit_p[t])
  }
  mu_p ~ dnorm(0, sd=1.5)
  sig_p ~ dexp(1)
  
  for(t in 1:(K-1)){
    logit_phi[t] ~ dnorm(mu_phi, sd=sig_phi)
    phi[t] <- expit(logit_phi[t])
  }
  mu_phi ~ dnorm(0, sd=1.5)
  sig_phi ~ dexp(1)
  
  d[1] <- 1
  xi[1] <- 1
  for(t in 2:K){
    log_f[t-1] ~ dnorm(mu_f, sd=sig_f)
    f[t-1] <- exp(log_f[t-1])
    xi[t] <- d[t-1] * f[t-1]
    d[t] <- d[t-1] * (phi[t-1] + f[t-1])
  }
  mu_f ~ dnorm(0, sd=1.5)
  sig_f ~ dexp(1)
  
  for(t in 1:(K-1)){
    xi_tilde[t] <- xi[t] / sum(xi[t:K])
  }
  xi_tilde[K] <- 1
  
  #' ---------------------------------------------------------------------------
  #' Initial entry probability
  #' ---------------------------------------------------------------------------
  pi[1] <- 1 - xi_tilde[1]
  pi[2] <- xi_tilde[1]
  pi[3] <- 0
  
  #' ---------------------------------------------------------------------------
  #' Detection probabilities matrix (K x J active states)
  #' ---------------------------------------------------------------------------
  # for(t in 1:K){
  #   p_mat[t, 1] <- p[t]  # State 2 = Available/Active state
  # }
  
  #' ---------------------------------------------------------------------------
  #' Transition probability matrix
  #' ---------------------------------------------------------------------------
  for(t in 1:(K-1)){
    Gamma[1:n_states, 1:n_states, t] <- load_Gamma(xi_tilde[t+1], phi[t])
  }
  
  #' lambda prior
  lambda ~ dgamma(1.0e-6, 1.0e-6)
  
  #' ---------------------------------------------------------------------------
  #' Unconditional detection probability (Binomial)
  #' ---------------------------------------------------------------------------
  pstar <- pstar_binom(
    init = pi[1:n_states],
    prob = matrix(p[1:K], K, 1),
    size = size_vec[1:K],
    probTrans = Gamma[1:n_states, 1:n_states, 1:(K-1)],
    len = K
  )
  
  #' ---------------------------------------------------------------------------
  #' HMM likelihood for observed individuals (Binomial)
  #' ---------------------------------------------------------------------------
  for(i in 1:nobs){
    x[i, 1:K] ~ dJS_binom(
      init = pi[1:n_states],
      prob = matrix(p[1:K], K, 1),
      size = size_vec[1:K],
      probTrans = Gamma[1:n_states, 1:n_states, 1:(K-1)],
      pstar = pstar, weight = 1, len = K
    )
  }
  n ~ dpois(lambda * pstar)
  
  #' ---------------------------------------------------------------------------
  #' Posterior predicted abundance (Binomial)
  #' ---------------------------------------------------------------------------
  nu ~ dpois(lambda * (1 - pstar))
  Nsuper <- n + nu
  
  nu_t[1:n_states, 1:K] ~ dJSUndet_binom(
    nu = nu,
    init = pi[1:n_states],
    prob = matrix(p[1:K], K, 1),
    size = size_vec[1:K],
    probTrans = Gamma[1:n_states, 1:n_states, 1:(K-1)]
  )
  
  # Sample observed individual latent histories
  for(i in 1:nobs) {
    z[i, 1:K] ~ dJSCond_binom(
      obs = x[i, 1:K],
      init = pi[1:n_states],
      prob =matrix(p[1:K], K, 1),
      size = size_vec[1:K],
      probTrans = Gamma[1:n_states, 1:n_states, 1:(K-1)]
    )
  }
  
  # Calculate occasion-specific abundance (State 2 = Available)
  for(t in 1:K) {
    N[t] <- nu_t[2, t] + sum(z[1:nobs, t] == 2) 
  }
})