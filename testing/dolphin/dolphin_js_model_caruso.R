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

  # alpha
  for(j in 1:3){alpha_raw[j] ~ dgamma(1,1)}
  for(j in 1:3){alpha[j] <- alpha_raw[j]/sum(alpha_raw[1:3])}

  # rho
  for(t in 1:K){
    for(j in 1:3){
      rho[t,j] ~ dbeta(1/K, 2 - t/K)
    }
  }

  # phi
  phi_RP ~ dbeta(1,1)
  phi_T ~ T(dbeta(1,20), 0, phi_RP)
  phi[1] <- phi_RP
  phi[2] <- phi_RP
  phi[3] <- phi_T

  # pi
  # Resident
  pi[1] <- alpha[1]*(1-rho[1,1])
  pi[2] <- alpha[1]*rho[1,1]
  pi[3] <- 0
  # Part-time
  pi[4] <- alpha[2]*(1-rho[1,2])
  pi[5] <- alpha[2]*rho[1,2]
  pi[6] <- 0
  # Transient
  pi[7] <- alpha[3]*(1-rho[1,3])
  pi[8] <- alpha[3]*rho[1,3]
  pi[9] <- 0

  for(t in 1:(K-1)){
    Gamma[1:9, 1:9, t] <- load_mix_Gamma(rho[t+1,1:3], phi[1:3]^h[t])
  }

  #' ---------------------------------------------------------------------------
  #' Detection matrix
  #' ---------------------------------------------------------------------------

  # P
  mu_p_raw ~ dnorm(0,0.1)
  delta_p ~ dbeta(1,1)
  for(t in 1:K){
    tau_raw[t] ~ dnorm(0,4)
  }
  mean_tau <- mean(tau_raw[1:K])
  tau[1:K] <- tau_raw[1:K]-mean_tau
  mu_p <- mu_p_raw + mean_tau
  for(t in 1:K){
    logit(p[t,1]) <- mu_p + tau[t]
    p[t,2] <- p[t,1] * delta_p
    p[t,3] <- p[t,1]
    Pmats[1:9,1:2,t] <- load_mix_P(p[t,1:3])
  }


  #' ---------------------------------------------------------------------------
  #' Unconditional detection probability
  #' ---------------------------------------------------------------------------
  pstar <- pstar_ms(
    init = pi[1:9],
    probObs = Pmats[1:9, 1:2, 1:K],
    probTrans = Gamma[1:9, 1:9, 1:(K-1)],
    len=K
  )

  #' ---------------------------------------------------------------------------
  #' HMM likelihood for observed individuals
  #' ---------------------------------------------------------------------------
  for(i in 1:nobs_u){
    ux[i, 1:K] ~ dJS_ms(
      init = pi[1:9],
      probObs = Pmats[1:9, 1:2, 1:K],
      probTrans = Gamma[1:9, 1:9, 1:(K-1)],
      pstar=pstar, weight=w[i], len = K
    )
    n2ll_id[i] <- dJS_ms_n2ll(ux[i,1:K], nobs, init = pi[1:9],
                           probObs = Pmats[1:9, 1:2, 1:K],
                           probTrans = Gamma[1:9, 1:9, 1:(K-1)],
                           pstar, lambda, weight=w[i])
  }
  n2ll <- sum(n2ll_id[1:nobs_u])

  #' ---------------------------------------------------------------------------
  #' Model for number of detected individuals
  #' ---------------------------------------------------------------------------
  lambda ~ dgamma(1.0e-6, 1.0e-6)
  n ~ dpois(lambda*pstar)

  #' ---------------------------------------------------------------------------
  #' Code below is for the posterior predicted abundance
  #' ---------------------------------------------------------------------------
  nu_raw ~ dpois(lambda*(1-pstar))
  nu_t[1:9, 1:K] <- sample_undet_ms(nu_raw, pi[1:9], Pmats[1:9,1:2,1:K], Gamma[1:9, 1:9, 1:(K-1)])
  nu <- nu_raw - (nu_t[1,K] + nu_t[4,K] + nu_t[7,K])
  Nsuper <- n+nu

  # for(i in 1:nobs){
  #   det_state[i,1:K] <- sample_det_ms(x[i,1:K],pi[1:9],Pmats[1:9,1:2,1:K],Gamma[1:9,1:9,1:(K-1)])
  #   avail_R[i,1:K] <- det_state[i,1:K]==2
  #   avail_P[i,1:K] <- det_state[i,1:K]==5
  #   avail_T[i,1:K] <- det_state[i,1:K]==8
  #   for(j in 1:3){
  #     avail_R_yr[i,j] <- sum(avail_R[i,yidx[j,1]:yidx[j,2]]) > 0
  #     avail_P_yr[i,j] <- sum(avail_P[i,yidx[j,1]:yidx[j,2]]) > 0
  #     avail_T_yr[i,j] <- sum(avail_T[i,yidx[j,1]:yidx[j,2]]) > 0
  #   }
  # }
  #
  # nu_t[1:9, 1:K] <- sample_undet_ms(nu, pi[1:9], Pmats[1:9,1:2,1:K], Gamma[1:9, 1:9, 1:(K-1)])
  # B_R[1] <- nu_t[2,1]
  # B_P[1] <- nu_t[5,1]
  # B_T[1] <- nu_t[8,1]
  # S_R[1] <-0
  # S_P[1] <- 0
  # S_T[1] <- 0
  # for(t in 2:K){
  #   S_R[t] <- nu_t[2,t-1] - (nu_t[3,t]-nu_t[3,t-1])
  #   S_P[t] <- nu_t[5,t-1] - (nu_t[6,t]-nu_t[6,t-1])
  #   S_T[t] <- nu_t[8,t-1] - (nu_t[9,t]-nu_t[9,t-1])
  #   B_R[t] <- nu_t[2,t] - S_R[t]
  #   B_P[t] <- nu_t[5,t] - S_P[t]
  #   B_T[t] <- nu_t[8,t] - S_T[t]
  # }
  # for(j in 1:3){
  #   N_R_yr[j] <- sum(avail_R_yr[1:nobs,j]) + sum(B_R[yidx[j,1]:yidx[j,2]]) + S_R[yidx[j,1]]
  #   N_P_yr[j] <- sum(avail_P_yr[1:nobs,j]) + sum(B_P[yidx[j,1]:yidx[j,2]]) + S_P[yidx[j,1]]
  #   N_T_yr[j] <- sum(avail_T_yr[1:nobs,j]) + sum(B_T[yidx[j,1]:yidx[j,2]]) + S_T[yidx[j,1]]
  # }

})
