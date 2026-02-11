
predlist <- nimbleList(
  N_nest = double(2),
  N_avail = double(2)
)

# debuglist <- nimbleList(
#   zeta = double(1),
#   pi = double(1),
#   Gamma = double(3),
#   Pmats = double(3),
#   nu_t = double(2)
# )

predict_abundance <- nimbleFunction(
  run = function(xi=double(2), phi = double(2), theta=double(1), p=double(2),
                 nd=integer(0), nu = double(1), x = double(2)) {
    iter <- dim(xi)[1]
    K <- dim(xi)[2]
    n <- dim(x)[1]
    n_states <- 3+nd

    xi_tilde <- numeric(K)
    Gamma <- nimArray(dim=c(n_states,n_states,K-1))
    pi <- numeric(n_states)
    Pmats <- nimArray(dim=c(n_states,2,K))
    nu_t <- nimMatrix(nrow=n_states, ncol=K)
    det_state <- nimMatrix(nrow=n, ncol = K)
    nest_obs <- nimMatrix(nrow=n, ncol = K)
    avail_obs <- nimMatrix(nrow=n, ncol = K)
    N_nest <- nimMatrix(nrow=iter, ncol=K)
    N_avail <- nimMatrix(nrow=iter, ncol=K)
    zeta <- numeric(nd)

    for(j in 1:iter){

      ### Load HMM matrices
      zeta <- make_dt_pois(nd, theta[j], shift=1)
      for(t in 1:(K-1)){
        xi_tilde[t] <- xi[j,t]/sum(xi[j,t:K])
      }
      xi_tilde[K] <- 1
      for(t in 1:(K-1)){
        Gamma[,,t] <- make_dt_Gamma(xi_tilde[t+1], phi[j,t], zeta)
        Pmats[,,t] <- make_dt_Pmats(p[j,t],nd)
      }
      Pmats[,,K] <- make_dt_Pmats(p[j,K],nd)
      pi <- make_dt_pi(xi_tilde[1],nd)

      ### Predict abundance for uncaptured id
      nu_t <- sample_undet_ms(nu[j], pi, Pmats, Gamma)

      ### Predict state for captured id
      for(i in 1:n){
        det_state[i, ] <- sample_det_ms(x[i, ], pi, Pmats, Gamma)
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
