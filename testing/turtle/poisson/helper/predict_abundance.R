
predlist <- nimbleList(
  N_nest = double(2),
  N_avail = double(2)
)

predict_abundance <- nimbleFunction(
  run = function(xi=double(2), phi = double(2), theta=double(2), r=double(3),
                 nd=integer(0), nu = double(1), x = double(2)) {
    iter <- dim(xi)[1]
    K <- dim(xi)[2]
    n <- dim(x)[1]
    n_states <- 4+nd

    xi_tilde <- numeric(K)
    Gamma <- nimArray(dim=c(n_states,n_states,K-1))
    pi <- numeric(n_states)
    Rmat <- nimMatrix(K,n_states)
    nu_t <- nimMatrix(nrow=n_states, ncol=K)
    det_state <- nimMatrix(nrow=n, ncol = K)
    nest_obs <- nimMatrix(nrow=n, ncol = K)
    avail_obs <- nimMatrix(nrow=n, ncol = K)
    N_nest <- nimMatrix(nrow=iter, ncol=K)
    N_avail <- nimMatrix(nrow=iter, ncol=K)
    zeta <- numeric(nd)

    for(j in 1:iter){

      ### Load HMM matrices
      for(t in 1:(K-1)){
        xi_tilde[t] <- xi[j,t]/sum(xi[j,t:K])
      }
      xi_tilde[K] <- 1

      Rmat <- make_honu_Rmat(r[j,,],nd)
      for(t in 1:(K-1)){
        zeta <- make_dt_pois(nd, theta[j,t], shift=1)
        Gamma[,,t] <- make_honu_Gamma(xi_tilde[t+1], phi[j,t], zeta)
      }
      pi <- make_honu_pi(xi_tilde[1],nd)

      ### Predict abundance for uncaptured id
      nu_t <- sample_undet_pois(nu[j], pi, Rmat, Gamma)

      ### Predict state for captured id
      for(i in 1:n){
        det_state[i, ] <- sample_det_pois(x[i, ], pi, Rmat, Gamma)
      }
      ### Determine abundnace for nesting and avaialble id
      for(t in 1:K){
        n_nest_obs <- 0
        n_avail_obs <- 0
        for(i in 1:n){
          if(det_state[i, t] == 2) n_nest_obs <- n_nest_obs + 1
          # Check if individual is 'available' (not in state 1 or n_states)
          if(det_state[i, t] > 1 & det_state[i, t] < n_states) {
            n_avail_obs <- n_avail_obs + 1
          }
        }
        # Assign to your result matrices
        N_nest[j, t] <- n_nest_obs + nu_t[2, t]
        N_avail[j, t] <-  n_avail_obs + nu[j] - (nu_t[1, t] + nu_t[n_states, t])
      }

    } # MCMC iteration loop end (j index)

    out <- predlist$new()
    # out$zeta <- zeta
    # out$pi <- pi
    # out$Gamma <- Gamma
    # out$Pmats <- Pmats
    # out$nu_t <- nu_t
    out$N_avail <- N_avail
    out$N_nest <- N_nest
    returnType(predlist())
    return(out)
  }) # function end
